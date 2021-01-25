"""Convert CellDesigner models to SBML-qual with a rather strict semantics.

Copyright (C) 2019-2021 Sylvain.Soliman@inria.fr

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import csv
import os.path
import sys
import xml.etree.ElementTree as etree
from itertools import chain, repeat
from typing import Dict, IO, List, Optional, Tuple, cast  # noqa: F401

from loguru import logger  # type: ignore

import networkx as nx  # type: ignore

from . import version

NS = {
    "sbml": "http://www.sbml.org/sbml/level2/version4",
    "cd": "http://www.sbml.org/2001/ns/celldesigner",
    "sbml3": "http://www.sbml.org/sbml/level3/version1/core",
    "layout": "http://www.sbml.org/sbml/level3/version1/layout/version1",
    "qual": "http://www.sbml.org/sbml/level3/version1/qual/version1",
    "mathml": "http://www.w3.org/1998/Math/MathML",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "dc": "http://purl.org/dc/elements/1.1/",
    "dcterms": "http://purl.org/dc/terms/",
    "vCard": "http://www.w3.org/2001/vcard-rdf/3.0#",
    "bqbiol": "http://biomodels.net/biology-qualifiers/",
    "bqmodel": "http://biomodels.net/model-qualifiers/",
    "xhtml": "http://www.w3.org/1999/xhtml",
}


# TODO use attrs, see https://www.attrs.org/en/latest/why.html#namedtuples
Transition = collections.namedtuple(
    "Transition", ["type", "reactants", "modifiers", "notes", "annotations"]
)


def read_celldesigner(filename: IO):
    """Parse the given file."""
    root = etree.parse(filename).getroot()
    tag = root.tag
    if tag != "{" + NS["sbml"] + "}sbml":
        raise ValueError("Currently limited to SBML Level 2 Version 4")
    model = root.find("sbml:model", NS)
    if model:
        display = model.find("./sbml:annotation/cd:extension/cd:modelDisplay", NS)
    else:
        raise ValueError("Could not find SBML model element")
    if display is None:
        raise ValueError("Could not find CellDesigner modelDisplay element")
    return (
        get_transitions(model, species_info(model)),
        display.get("sizeX"),
        display.get("sizeY"),
    )


def species_info(model):
    """Create a map from species' ids to their attributes."""
    nameconv = {}
    compartments = {}
    # Find all CellDesigner species used later
    for species in chain(
        model.findall(
            "./sbml:annotation/cd:extension/"
            + "cd:listOfComplexSpeciesAliases/"
            + "cd:complexSpeciesAlias",
            NS,
        ),
        model.findall(
            "./sbml:annotation/cd:extension/"
            + "cd:listOfSpeciesAliases/"
            + "cd:speciesAlias",
            NS,
        ),
    ):
        bound = species.find(".//cd:bounds", NS)
        in_complex = species.get("complexSpeciesAlias")
        if bound is None or in_complex is not None:
            continue
        ref_species = species.get("species")
        logger.debug("parsing ref_species: {ref}", ref=ref_species)
        sbml = model.find(
            './sbml:listOfSpecies/sbml:species[@id="' + ref_species + '"]', NS
        )
        if sbml is None:
            continue
        annot = sbml.find("./sbml:annotation", NS)
        if annot is None:
            continue
        classtype = get_text(annot.find(".//cd:class", NS), "PROTEIN")
        if classtype == "DEGRADED":
            continue
        if classtype == "PROTEIN":
            is_receptor = find_protein_type(annot, model) == "RECEPTOR"
        else:
            is_receptor = False
        mods = get_mods(annot.find(".//cd:listOfModifications", NS))
        activity = annot.find(".//cd:structuralState", NS)
        if activity is not None:
            activity = activity.get("structuralState")
        else:
            activity = "inactive"
        name = make_name_precise(sbml.get("name"), classtype, mods)
        comp_id = species.get("compartmentAlias")
        if comp_id in compartments:
            compartment = compartments[comp_id]
        else:
            compartment = find_compartment(comp_id, model)
            compartments[comp_id] = compartment
        species_id = species.get("id")
        nameconv[species_id] = {
            "activity": activity,
            "x": bound.get("x"),
            "y": bound.get("y"),
            "h": bound.get("h"),
            "w": bound.get("w"),
            "transitions": [],
            "name": name,
            "function": name,
            "ref_species": ref_species,
            "type": classtype,
            "modifications": mods,
            "receptor": is_receptor,
            "annotations": annot.find(".//rdf:RDF", NS),
            "compartment": compartment,
        }
        # also store in nameconv the reverse mapping from SBML species to CD
        # species using the corresponding reference protein
        prot_ref = "__" + sbml.get("name")
        if prot_ref in nameconv:
            nameconv[prot_ref].append(species.get("id"))
        else:
            nameconv[prot_ref] = [species.get("id")]
    add_subcomponents_only(nameconv, model)
    return nameconv


def find_protein_type(annotation, model):
    """Look for the cd:protein type for an annotation's reference protein."""
    ref = get_text(annotation.find(".//cd:proteinReference", NS))
    if ref:
        protein = model.find('.//cd:protein[@id="' + ref + '"]', NS)
        if protein is not None:
            return protein.get("type")
    return "GENERIC"


def find_compartment(comp_id, model):
    """Look for the name of the SBML compartment associated to a CD one."""
    if comp_id is None:
        return "default_compartment"
    sbml_id = model.find('.//cd:compartmentAlias[@id="' + comp_id + '"]', NS).get(
        "compartment"
    )
    return model.find('.//sbml:compartment[@id="' + sbml_id + '"]', NS).get("name")


def make_name_precise(name, ctype, mods):
    """Append molecule type and modifications to its cleaned-up name."""
    to_map = {"&": "", "|": "", "!": "", "underscore": ""}
    to_remove = {"sub", "endsub"}
    newname = "_".join(
        map(
            lambda s: to_map[s] if s in to_map else s,
            filter(lambda t: t not in to_remove, name.split("_")),
        )
    ).replace("__", "_")
    if ctype == "PROTEIN":
        basis = [newname]
    else:
        basis = [newname, ctype.lower()]
    return "_".join(basis + mods)


def add_subcomponents_only(nameconv, model: etree.Element):
    """Add annotations to the parent complex.

    For unused CD species (only subcomponents of complexes)
    """
    for species in model.findall(
        "./sbml:annotation/cd:extension/"
        + "cd:listOfIncludedSpecies/cd:species/"
        + "cd:notes/xhtml:html/xhtml:body/"
        + "rdf:RDF/../../../..",
        NS,
    ):
        add_rdf(
            nameconv,
            reference=decomplexify(species.get("id"), model, field="species"),
            new_rdf=species.find(".//rdf:RDF", NS),
        )


def add_rdf(nameconv, reference: str, new_rdf: Optional[etree.Element]):
    """Add the new_rdf element to nameconv[reference]['annotations']."""
    # reference might still be a complex, then we ignore it (is this correct?)
    if new_rdf is None or reference not in nameconv:
        return
    if nameconv[reference]["annotations"] is not None:
        rdfs = new_rdf.find("./rdf:Description", NS)
        if rdfs is None:
            return
        nameconv[reference]["annotations"].find("./rdf:Description", NS).extend(rdfs[:])
    else:
        nameconv[reference]["annotations"] = new_rdf


def get_transitions(model: etree.Element, info):
    """Find all transitions."""
    for trans in model.findall("./sbml:listOfReactions/sbml:reaction", NS):
        logger.debug("parsing reaction: {tid}", tid=trans.get("id"))
        annot = trans.find("./sbml:annotation/cd:extension", NS)
        if annot is None:
            continue
        rtype = get_text(annot.find("./cd:reactionType", NS))
        reacs = [
            decomplexify(reac.get("alias"), model)
            for reac in chain(
                annot.findall("./cd:baseReactants/cd:baseReactant", NS),
                annot.findall("./cd:listOfReactantLinks/cd:reactantLink", NS),
            )
        ]
        prods = [
            decomplexify(prod.get("alias"), model)
            for prod in chain(
                annot.findall("./cd:baseProducts/cd:baseProduct", NS),
                annot.findall("./cd:listOfProductLinks/cd:productLink", NS),
            )
        ]
        mods = [
            (mod.get("type"), decomplexify(mod.get("aliases"), model))
            for mod in annot.findall("./cd:listOfModification/cd:modification", NS)
        ]
        notes = trans.find("./sbml:notes//xhtml:body", NS)
        rdf = trans.find("./sbml:annotation/rdf:RDF", NS)
        # remove degraded
        reacs = list(filter(lambda x: x in info, reacs))
        prods = list(filter(lambda x: x in info, prods))
        # for each product of a reaction, add this reaction as a transition
        # affecting that species
        for species in prods:
            info[species]["transitions"].append(
                Transition(rtype, reacs, mods, notes, rdf)
            )
    return info


def decomplexify(species: str, model: etree.Element, field: str = "id"):
    """Return external complex if there is one.

    or species unchanged otherwise.
    """
    cmplx = model.find(
        "./sbml:annotation/cd:extension/"
        + "cd:listOfSpeciesAliases/"
        + 'cd:speciesAlias[@{field}="{species}"]'.format(field=field, species=species),
        NS,
    )
    if cmplx is None:
        return species
    return cmplx.get("complexSpeciesAlias", species)


def get_text(cd_class: Optional[etree.Element], default=None):
    """Get the text of an XML field if it exists or return a default."""
    if cd_class is not None:
        return cd_class.text
    return default


def get_mods(cd_modifications: etree.Element) -> List[str]:
    """Celldesigner:listOfModifications to list of mods."""
    if not cd_modifications:
        return []
    return [mod.get("state") for mod in cd_modifications.findall("cd:modification", NS)]


def write_qual(
    filename: str, info, width: str, height: str, remove: int = 0, sif: bool = False
):
    # pylint: disable=too-many-arguments, too-many-locals
    """Write the SBML qual with layout file for our model."""
    for name, space in NS.items():
        etree.register_namespace(name, space)
    root = etree.Element(
        "sbml",
        {
            "level": "3",
            "version": "1",
            "layout:required": "false",
            "xmlns": NS["sbml3"],
            "qual:required": "true",
            "xmlns:layout": NS["layout"],
            "xmlns:qual": NS["qual"],
        },
    )
    model = etree.SubElement(root, "model", id="model_id")
    clist = etree.SubElement(model, "listOfCompartments")
    etree.SubElement(clist, "compartment", constant="true", id="comp1")
    llist = etree.SubElement(model, "layout:listOfLayouts")
    layout = etree.SubElement(llist, "layout:layout", id="layout1")
    etree.SubElement(layout, "layout:dimensions", width=width, height=height)
    qlist = etree.SubElement(model, "qual:listOfQualitativeSpecies")
    tlist = etree.SubElement(model, "qual:listOfTransitions")
    graph = nx.DiGraph()
    fix_all_names(info)
    add_transitions(tlist, info, graph)
    remove_connected_components(tlist, info, graph, remove)
    if sif:
        write_sif(filename, info, graph)
    add_qual_species(layout, qlist, info)
    etree.ElementTree(root).write(filename, encoding="utf-8", xml_declaration=True)


def remove_connected_components(
    tlist: etree.Element, info, graph: nx.DiGraph, remove: int
):
    """Remove connected components of size smaller than remove."""
    # because we did not properly NameSpace all transitions, we cannot use
    # find('./qual:transition[@qual:id=]')
    logger.debug("remove value {S}", S=remove)
    transitions = list(tlist)
    ccs = list(nx.connected_components(graph.to_undirected()))
    # add completely isolated nodes
    nodes = graph.nodes()
    ccs.extend({species} for (species, data) in info.items() if species not in nodes)
    logger.debug("CCs: {ccs}", ccs=ccs)
    if remove < 0:
        remove = len(max(ccs, key=len)) - 1
    logger.debug("remove value {S}", S=remove)
    for component in filter(lambda x: len(x) <= remove, ccs):
        logger.debug(
            "removing connected component {component}", component=list(component)
        )
        for species in component:
            logger.debug("removing species {sp}", sp=species)
            del info[species]
            trans_id = "tr_" + species
            trans = next((t for t in transitions if t.get("qual:id") == trans_id), None)
            if trans:
                logger.debug("removing transition {tr}", tr=trans)
                tlist.remove(trans)


def write_sif(sbml_filename: str, info, graph: nx.DiGraph):
    """Write a SIF file with influences.

    http://www.cbmc.it/fastcent/doc/SifFormat.htm
    """
    with open(sbml_filename[:-4] + "sif", "w", encoding="utf-8", newline="") as f:
        print("# Generated by CaSQ v{version}".format(version=version), file=f)
        with open(
            sbml_filename[:-5] + "_raw.sif", "w", encoding="utf-8", newline=""
        ) as fraw:
            print("# Generated by CaSQ v{version}".format(version=version), file=fraw)
            for source, target, sign in graph.edges.data("sign"):
                print(
                    source.replace(" ", "_"),
                    sign.upper(),
                    target.replace(" ", "_"),
                    file=fraw,
                )
                print(
                    info[source]["name"].replace(" ", "_"),
                    sign.upper(),
                    info[target]["name"].replace(" ", "_"),
                    file=f,
                )


def simplify_model(info):
    """Clean the model w.r.t. some active/inactive species."""
    multispecies = delete_complexes_and_store_multispecies(info)
    # pylint: disable=too-many-nested-blocks
    for key, value in multispecies.items():
        for val in value:
            # check that it does not appear in any other reaction than the
            # activation one
            logger.debug("looking at multispecies: {mul}", mul=val)
            active = get_active(val, info)
            if val not in info:
                # val has been deleted just above
                continue
            if active in value:
                if not info[val]["transitions"]:
                    # we know that active is a str here, since it is in value
                    add_rdf(info, cast(str, active), info[val]["annotations"])
                    logger.debug(
                        "deleting {val} [{active} is active for {key}]",
                        val=val,
                        active=active,
                        key=key,
                    )
                    del info[val]
                elif len(info[active]["transitions"]) == 1:
                    reac = info[active]["transitions"][0]
                    if (
                        reac.type == "TRANSPORT"
                        and reac.modifiers == []
                        and reac.reactants == [val]
                    ):
                        info[active]["transitions"] = info[val]["transitions"]
                        add_rdf(info, cast(str, active), info[val]["annotations"])
                        logger.debug(
                            "merging {val} into {active} in transport for {key}]",
                            val=val,
                            active=active,
                            key=key,
                        )
                        del info[val]


def delete_complexes_and_store_multispecies(info):
    """Delete useless species and store multispecies.

    Useless species are ligand binding to a receptor, or otherwise unused
    proteins that bind to form a complex.
    """
    multispecies = {}  # type: Dict[str, List[str]]
    duplicate_nodes = {}
    replacements = {}
    # we have to create the list since info will change during iteration
    for key, value in list(info.items()):
        if key.startswith("__"):
            del info[key]
            logger.debug("deleting complex: {cplx} value: {val}", cplx=key, val=value)
            if len(value) > 1:
                multispecies[key] = value
        # merge nodes that have the same reference species
        elif (
            value["ref_species"] in duplicate_nodes
            and key in info
            and duplicate_nodes[value["ref_species"]] in info
        ):
            into = duplicate_nodes[value["ref_species"]]
            logger.debug(
                "merging {key} into {into} for {ref} ({name})",
                key=key,
                into=into,
                ref=value["ref_species"],
                name=value["name"],
            )
            # put our annotations with those of into
            add_rdf(info, into, info[key]["annotations"])
            # put our transitions with those of into
            info[into]["transitions"].extend(info[key]["transitions"])
            # replace key in other transitions with into
            replacements[key] = into
            del info[key]
        # delete receptors that only contribute to their complex
        elif value["receptor"] and not value["transitions"]:
            logger.debug("{key} is a RECEPTOR (and an input)", key=key)
            active = get_active(key, info)
            if active and [
                trans
                for trans in info[active]["transitions"]
                if key in trans.reactants and trans.type == "HETERODIMER_ASSOCIATION"
            ]:
                logger.debug(
                    "deleting input RECEPTOR {key} that dimerizes to form {active}",
                    key=key,
                    active=active,
                )
                add_rdf(info, cast(str, active), info[key]["annotations"])
                del info[key]
        elif value["type"] == "COMPLEX":
            for trans in value["transitions"]:
                if trans.type == "HETERODIMER_ASSOCIATION":
                    if len(trans.reactants) != 2:
                        continue
                    [reac1, reac2] = trans.reactants
                    if reac1 not in info or reac2 not in info:
                        continue
                    active1 = get_active(reac1, info)
                    active2 = get_active(reac2, info)
                    if (
                        active1 == key
                        and active2 == key
                        and not info[reac1]["receptor"]
                        and not info[reac2]["receptor"]
                    ):
                        logger.debug(
                            "deleting {reac1} and {reac2} for complex {key}",
                            reac1=reac1,
                            reac2=reac2,
                            key=key,
                        )
                        add_rdf(info, key, info[reac1]["annotations"])
                        add_rdf(info, key, info[reac2]["annotations"])
                        info[key]["transitions"].extend(info[reac1]["transitions"])
                        info[key]["transitions"].extend(info[reac2]["transitions"])
                        if info[reac1]["ref_species"] in duplicate_nodes:
                            duplicate_nodes[info[reac1]["ref_species"]] = key
                        if info[reac2]["ref_species"] in duplicate_nodes:
                            duplicate_nodes[info[reac2]["ref_species"]] = key
                        del info[reac1]
                        del info[reac2]
        else:
            duplicate_nodes[value["ref_species"]] = key
    replace_in_transitions(info, replacements)
    return multispecies


def replace_in_transitions(info, replacements):
    """Change transitions in info to reflect replacements."""
    for _species, data in info.items():
        for trans in data["transitions"]:
            for val in trans.reactants.copy():
                if val in replacements:
                    trans.reactants.append(replacements[val])
                    trans.reactants.remove(val)
            for (modtype, mod_list) in trans.modifiers.copy():
                mlist = mod_list.split(",")
                changed = False
                for val in mlist:
                    if val in replacements:
                        mlist.append(replacements[val])
                        mlist.remove(val)
                        changed = True
                if changed:
                    trans.modifiers.append((modtype, ",".join(mlist)))
                    trans.modifiers.remove((modtype, mod_list))


def get_active(val, info):
    """Find who val activates."""
    active = None
    for species, data in info.items():
        if species.startswith("csa") or species.startswith("sa"):
            for trans in data["transitions"]:
                if val in trans.reactants or val in (
                    mod
                    for (_modtype, modifier_list) in trans.modifiers
                    for mod in modifier_list.split(",")
                ):
                    if active is None:
                        active = species
                        logger.debug("{val} activates {active}", val=val, active=active)
                    else:
                        logger.debug(
                            "{val} also activates {active}", val=val, active=species
                        )
                        return False
    return active


def fix_all_names(info):
    """Use more descriptive names."""
    count_names = collections.Counter(value["name"] for value in info.values())
    ambiguous_name = {
        key: count_names[value["name"]] > 1 for key, value in info.items()
    }
    namedict = {}
    for species, data in info.items():
        name = fix_name(data["name"], ambiguous_name[species], data["compartment"])
        activity = data["activity"]
        if ambiguous_name[species]:
            if name in namedict:
                other_id, other_activity = namedict[name]
                if activity == "active":
                    # FIXME what if name_active in namedict?
                    name = name + "_active"
                elif other_activity == "active":
                    info[other_id]["name"] = name + "_active"
                    info[other_id]["function"] = name + "_active"
                else:
                    # FIXME don't know what to do
                    # print(f"active is {activity}, other was {other_activity}")
                    pass
            namedict[name] = (species, activity)
        data["name"] = name
        data["function"] = name


def add_qual_species(layout: etree.Element, qlist: etree.Element, info):
    """Create layout sub-elements and species."""
    llist = etree.SubElement(layout, "layout:listOfAdditionalGraphicalObjects")
    for species, data in info.items():
        glyph = etree.SubElement(
            llist,
            "layout:generalGlyph",
            {"layout:reference": species, "layout:id": species + "_glyph"},
        )
        box = etree.SubElement(glyph, "layout:boundingBox")
        etree.SubElement(
            box, "layout:position", {"layout:x": data["x"], "layout:y": data["y"]}
        )
        etree.SubElement(
            box,
            "layout:dimensions",
            {"layout:height": data["h"], "layout:width": data["w"]},
        )
        if data["transitions"]:
            constant = "false"
        else:
            constant = "true"
        qspecies = etree.SubElement(
            qlist,
            "qual:qualitativeSpecies",
            {
                "qual:maxLevel": "1",
                "qual:compartment": "comp1",
                "qual:name": data["name"],
                "qual:constant": constant,
                "qual:id": species,
            },
        )
        add_annotation(qspecies, data["annotations"])


def fix_name(name: str, ambiguous: bool, compartment: str):
    """Change name for GINSIM compatibility or to remove subscripts.

    Also add compartment to the name in case of ambiguity.
    """
    if ambiguous:
        new_name = name + "_" + compartment
    else:
        new_name = name
    return (
        new_name.replace("_minus_", "-")
        .replace("_plus_", "+")
        .replace("_super", "^")
        .replace("_slash_", "/")
    )


def add_annotation(node: etree.Element, rdf: Optional[etree.Element]):
    """Add a single RDF element as an annotation node."""
    if rdf is not None:
        etree.SubElement(node, "annotation").append(rdf)


def add_transitions(tlist: etree.Element, info, graph: nx.DiGraph):
    """Create transition elements."""
    known = list(info.keys())
    for species, data in info.items():
        if data["transitions"]:
            trans = etree.SubElement(
                tlist, "qual:transition", {"qual:id": "tr_" + species}
            )
            ilist = etree.SubElement(trans, "qual:listOfInputs")
            add_inputs(ilist, data["transitions"], species, known, graph)
            # there might not be any input left after filtering known species
            if not ilist:
                tlist.remove(trans)
                logger.debug(
                    "transition for {species} exists {trans} but has no inputs",
                    trans=data["transitions"],
                    species=species,
                )
                info[species]["transitions"] = []
                add_function_as_rdf(info, species, info[species]["function"])
            else:
                olist = etree.SubElement(trans, "qual:listOfOutputs")
                etree.SubElement(
                    olist,
                    "qual:output",
                    {
                        "qual:qualitativeSpecies": species,
                        "qual:transitionEffect": "assignmentLevel",
                        "qual:id": "tr_{species}_out".format(species=species),
                    },
                )
                flist = etree.SubElement(trans, "qual:listOfFunctionTerms")
                etree.SubElement(flist, "qual:defaultTerm", {"qual:resultLevel": "0"})
                func = etree.SubElement(
                    flist, "qual:functionTerm", {"qual:resultLevel": "1"}
                )
                add_function(func, data["transitions"], known)
                sfunc = mathml_to_ginsim(func.find("./math/*", NS), info)
                info[species]["function"] = sfunc
                add_function_as_rdf(info, species, sfunc)
                add_notes(trans, data["transitions"])
                add_annotations(trans, data["transitions"])
        else:
            add_function_as_rdf(info, species, info[species]["function"])


def add_notes(trans: etree.Element, transitions: List[Transition]):
    """Add all the found notes."""
    notes = etree.SubElement(trans, "notes")
    html = etree.SubElement(notes, "html", xmlns=NS["xhtml"])
    head = etree.SubElement(html, "head")
    etree.SubElement(head, "title")
    body = etree.SubElement(html, "body")
    some_notes = False
    for reaction in transitions:
        if reaction.notes is not None:
            some_notes = True
            reaction.notes.tag = "p"
            for element in reaction.notes.iter():
                for prefix in ("xhtml", "sbml"):
                    prefix_len = len(NS[prefix]) + 2
                    if element.tag.startswith("{" + NS[prefix] + "}"):
                        element.tag = element.tag[prefix_len:]
            body.append(reaction.notes)
    if not some_notes:
        trans.remove(notes)


def add_annotations(trans: etree.Element, transitions: List[Transition]):
    """Add all the found annotations."""
    annotation = etree.SubElement(trans, "annotation")
    rdf = etree.SubElement(annotation, "rdf:RDF")
    for reaction in transitions:
        if reaction.annotations is not None:
            rdf.append(reaction.annotations[0])
    if not rdf:
        trans.remove(annotation)


def add_function(func: etree.Element, transitions: List[Transition], known: List[str]):
    """Add the complete boolean activation function.

    this is an or over all reactions having the target as product.
    For each reaction it can activate if all reactants are present,
    no inhibitor is present, and one of the activators is present.
    """
    math = etree.SubElement(func, "math", xmlns=NS["mathml"])
    # create or node if necessary
    if len(transitions) > 1:
        apply = etree.SubElement(math, "apply")
        etree.SubElement(apply, "or")
    else:
        apply = math
    for reaction in transitions:
        # we assume that only "BOOLEAN_LOGIC_GATE_AND" has multiple modifiers
        # it is also the only modification that has an AND and therefore ends
        # with reactants
        reactants = [reac for reac in reaction.reactants if reac in known]
        reactants.extend(
            [
                mod
                for (modtype, modifier) in reaction.modifiers
                for mod in modifier.split(",")
                if modtype == "BOOLEAN_LOGIC_GATE_AND" and mod in known
            ]
        )
        activators = [
            modifier
            for (modtype, modifier) in reaction.modifiers
            if modtype not in ("INHIBITION", "BOOLEAN_LOGIC_GATE_AND")
            and modifier in known
            and modifier not in reactants
        ]
        inhibitors = [
            modifier
            for (modtype, modifier) in reaction.modifiers
            if modtype == "INHIBITION" and modifier in known
        ]
        # this should only appear when species is of type PHENOTYPE otherwise
        # non-SBGN compliant
        # just swap reactants and inhibitors, there should not be any activator
        if reaction.type == "NEGATIVE_INFLUENCE":
            reactants, inhibitors = inhibitors, reactants
        # create and node if necessary
        if len(reactants) + len(inhibitors) > 1 or (
            activators and (reactants or inhibitors)
        ):
            lapply = etree.SubElement(apply, "apply")
            etree.SubElement(lapply, "and")
        else:
            lapply = apply
        if len(activators) < 2:
            reactants.extend(activators)
        else:
            # create or node if necessary
            inner_apply = etree.SubElement(lapply, "apply")
            etree.SubElement(inner_apply, "or")
            for modifier in activators:
                set_level(inner_apply, modifier, "1")
        for level, modifier in chain(
            zip(repeat("1"), reactants), zip(repeat("0"), inhibitors)
        ):
            set_level(lapply, modifier, level)


def set_level(elt: etree.Element, modifier: str, level: str):
    """Add mathml to element elt such that modifier is equal to level."""
    trigger = etree.SubElement(elt, "apply")
    etree.SubElement(trigger, "eq")
    math_ci = etree.SubElement(trigger, "ci")
    math_ci.text = modifier
    math_cn = etree.SubElement(trigger, "cn", type="integer")
    math_cn.text = level


def negate(sign: str):
    """Change a sign represented as a string."""
    if sign == "negative":
        return "positive"
    return "negative"


def add_inputs(
    ilist: etree.Element,
    transitions: List[Transition],
    species: str,
    known: List[str],
    graph: nx.DiGraph,
):
    """Add all known inputs."""
    index = 0
    modifiers = []  # type: List[Tuple[str, str]]
    graph.add_node(species)
    for reaction in transitions:
        # we use enumerate to get a dummy modtype for reactants
        for modtype, modifier in chain(
            enumerate(reaction.reactants), reaction.modifiers
        ):
            if modtype == "INHIBITION":
                sign = "negative"
            else:
                sign = "positive"
            # this should only appear when species is of type PHENOTYPE
            # otherwise non-SBGN compliant
            if reaction.type == "NEGATIVE_INFLUENCE":
                sign = negate(sign)
            if (modifier, sign) not in modifiers and modifier in known:
                modifiers.append((modifier, sign))
                graph.add_edge(modifier, species, sign=sign)
                etree.SubElement(
                    ilist,
                    "qual:input",
                    {
                        "qual:qualitativeSpecies": modifier,
                        "qual:transitionEffect": "none",
                        "qual:sign": sign,
                        "qual:id": "tr_{species}_in_{index}".format(
                            species=species, index=index
                        ),
                    },
                )
                index += 1


def write_csv(sbml_filename: str, info):
    """Write a csv file with SBML IDs, CD IDs, Names, Formulae, etc."""
    # pylint: disable=invalid-name
    with open(sbml_filename[:-4] + "csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        for species, data in sorted(info.items()):
            writer.writerow(
                [species, data["name"], data["ref_species"], data["function"]]
            )


def mathml_to_ginsim(math: Optional[etree.Element], info) -> str:
    """Convert a MATHML boolean formula into its GinSIM representation."""
    if math is None:
        raise ValueError("Empty math element")
    if math.tag != "apply":
        raise ValueError(etree.tostring(math))
    children = list(math)
    if children[0].tag == "and":
        return "&".join(map(lambda x: mathml_to_ginsim(x, info), children[1:]))
    if children[0].tag == "or":
        return "|".join(map(lambda x: mathml_to_ginsim(x, info), children[1:]))
    if children[0].tag == "eq":
        species = children[1].text
        species = info[species]["name"]
        if species is None:
            species = ""
        if children[2].text == "0":
            return "!" + species
        return species
    raise ValueError(etree.tostring(math))


def add_function_as_rdf(info, species: str, func: str):
    """Add a new RDF element containing the logical function and name."""
    rdf = etree.Element("{{{rdf}}}RDF".format(rdf=NS["rdf"]))
    descr = etree.SubElement(
        rdf,
        "{{{rdf}}}Description".format(rdf=NS["rdf"]),
        attrib={
            "{{{rdf}}}about".format(rdf=NS["rdf"]): "#" + info[species]["ref_species"]
        },
    )
    bqbiol = etree.SubElement(
        descr, "{{{bqbiol}}}isDescribedBy".format(bqbiol=NS["bqbiol"])
    )
    bag = etree.SubElement(bqbiol, "{{{rdf}}}Bag".format(rdf=NS["rdf"]))
    etree.SubElement(
        bag,
        "{{{rdf}}}li".format(rdf=NS["rdf"]),
        attrib={"{{{rdf}}}resource".format(rdf=NS["rdf"]): "urn:casq:function:" + func},
    )
    bqbiol = etree.SubElement(
        descr, "{{{bqbiol}}}isDescribedBy".format(bqbiol=NS["bqbiol"])
    )
    bag = etree.SubElement(bqbiol, "{{{rdf}}}Bag".format(rdf=NS["rdf"]))
    etree.SubElement(
        bag,
        "{{{rdf}}}li".format(rdf=NS["rdf"]),
        attrib={
            "{{{rdf}}}resource".format(rdf=NS["rdf"]): "urn:casq:cdid:"
            + info[species]["ref_species"]
        },
    )
    add_rdf(info, species, rdf)


def main():
    """Run conversion using the CLI given first argument."""
    parser = argparse.ArgumentParser(
        description=" ".join(__doc__.splitlines()[:3]) + " GPLv3"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s v{version}".format(version=version),
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="Display a lot of debug information"
    )
    parser.add_argument(
        "-c",
        "--csv",
        action="store_true",
        help="Store the species information in a separate CSV file",
    )
    parser.add_argument(
        "-s",
        "--sif",
        action="store_true",
        help="Store the influence information in a separate SIF file",
    )
    parser.add_argument(
        "-r",
        "--remove",
        metavar="S",
        type=int,
        default="0",
        help="""Delete connected components in the resulting model if their
        size is smaller than S.
        A negative S leads to keep only the biggest(s) connected component(s)""",
    )
    parser.add_argument(
        "infile",
        type=argparse.FileType("r", encoding="utf-8"),
        nargs="?",
        default=sys.stdin,
        help="CellDesigner File",
    )
    parser.add_argument("outfile", nargs="?", default=sys.stdout, help="SBML-Qual File")
    args = parser.parse_args()

    if not args.debug:
        logger.disable("casq")
    logger.debug("parsing {fname}…", fname=args.infile.name)
    info, width, height = read_celldesigner(args.infile)
    simplify_model(info)
    if args.infile != sys.stdin and args.outfile == sys.stdout:
        args.outfile = os.path.splitext(args.infile.name)[0] + ".sbml"
    write_qual(args.outfile, info, width, height, remove=args.remove, sif=args.sif)
    if args.csv and args.outfile != sys.stdout:
        write_csv(args.outfile, info)


if __name__ == "__main__":  # pragma: no cover
    main()
