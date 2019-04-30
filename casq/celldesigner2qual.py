#!/usr/bin/env python3
"""Convert CellDesigner models to SBML-qual with a rather strict semantics."""

import argparse
import collections
import csv
import os.path
import sys
import xml.etree.ElementTree as etree
from itertools import chain, repeat
from typing import Dict, IO, List, Optional, Tuple, cast

import networkx as nx

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


Transition = collections.namedtuple(
    "Transition", ["type", "reactants", "modifiers", "notes", "annotations"]
)


def read_celldesigner(filename: IO, debug: bool = False):
    """Parse the given file."""
    root = etree.parse(filename).getroot()
    tag = root.tag
    if tag != "{" + NS["sbml"] + "}sbml":
        raise ValueError("Currently limited to SBML Level 2 Version 4")
    model = root.find("sbml:model", NS)
    if model:
        display = model.find("./sbml:annotation/cd:extension/cd:modelDisplay", NS)
    else:
        raise ValueError("Could not find SBML model element")
    if display is None:
        raise ValueError("Could not find CellDesigner modelDisplay element")
    return (
        get_transitions(model, species_info(model, debug), debug),
        display.get("sizeX"),
        display.get("sizeY"),
    )


def species_info(model, debug):
    """Create a map from species' ids to their attributes."""
    nameconv = {}
    # Find all CellDesigner species used later
    for species in chain(
        model.findall(
            "./sbml:annotation/cd:extension/"
            + "cd:listOfComplexSpeciesAliases/"
            + "cd:complexSpeciesAlias[@compartmentAlias]",
            NS,
        ),
        model.findall(
            "./sbml:annotation/cd:extension/"
            + "cd:listOfSpeciesAliases/"
            + "cd:speciesAlias[@compartmentAlias]",
            NS,
        ),
    ):
        # if species.get("complexSpeciesAlias"):
        #     continue
        bound = species.find(".//cd:bounds", NS)
        if bound is None:
            continue
        ref_species = species.get("species")
        if debug:
            print("parsing ref_species:", ref_species, file=sys.stderr)
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
        nameconv[species.get("id")] = {
            "activity": get_text(species.find(".//cd:activity", NS), "inactive"),
            "x": bound.get("x"),
            "y": bound.get("y"),
            "h": bound.get("h"),
            "w": bound.get("w"),
            "transitions": [],
            "name": sbml.get("name"),
            "function": sbml.get("name"),
            "ref_species": ref_species,
            "type": classtype,
            "modifications": get_mods(annot.find(".//cd:listOfModifications", NS)),
            "annotations": annot.find(".//rdf:RDF", NS),
        }
        # also store in nameconv the reverse mapping from SBML species to CD
        # species using the corresponding reference protein
        prot_ref = "__" + sbml.get("name")
        if prot_ref in nameconv:
            nameconv[prot_ref].append(species.get("id"))
        else:
            nameconv[prot_ref] = [species.get("id")]
    add_subcomponents_only(nameconv, model)
    return nameconv


def add_subcomponents_only(nameconv, model: etree.Element):
    """Add annotations to the parent complex.

    For unused CD species (only subcomponents of complexes)
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
    if new_rdf is None:
        return
    if nameconv[reference]["annotations"] is not None:
        rdfs = new_rdf.find("./rdf:Description", NS)
        if rdfs is None:
            return
        nameconv[reference]["annotations"].find("./rdf:Description", NS).extend(rdfs[:])
    else:
        nameconv[reference]["annotations"] = new_rdf


def get_transitions(model: etree.Element, info, debug: bool):
    """Find all transitions."""
    for trans in model.findall("./sbml:listOfReactions/sbml:reaction", NS):
        if debug:
            print("parsing reaction:", trans.get("id"), file=sys.stderr)
        annot = trans.find("./sbml:annotation/cd:extension", NS)
        if annot is None:
            continue
        # currently not used
        rtype = get_text(annot.find("./cd:reactionType", NS))
        reacs = [
            decomplexify(reac.get("alias"), model)
            for reac in annot.findall("./cd:baseReactants/cd:baseReactant", NS)
        ]
        prods = [
            decomplexify(prod.get("alias"), model)
            for prod in annot.findall("./cd:baseProducts/cd:baseProduct", NS)
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
    filename: str,
    info,
    width: str,
    height: str,
    ginsim_names: bool,
    remove: int = 0,
    debug: bool = False,
):
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
    graph = nx.Graph()
    add_transitions(tlist, info, graph)
    remove_connected_components(tlist, info, graph, remove, debug)
    add_qual_species(layout, qlist, info, ginsim_names)
    etree.ElementTree(root).write(filename, encoding="utf-8", xml_declaration=True)


def remove_connected_components(
    tlist: etree.Element, info, graph: nx.Graph, remove: int, debug: bool
):
    """Remove connected components of size smaller than remove."""
    # because we did not properly NameSpace all transitions, we cannot use
    # find('./qual:transition[@qual:id=]')
    transitions = list(tlist)
    ccs = list(nx.connected_components(graph))
    if remove < 0:
        remove = len(max(ccs, key=len)) - 1
    for cc in filter(lambda x: len(x) <= remove, ccs):
        if debug:
            print("removing connected component", list(cc), file=sys.stderr)
        for species in cc:
            if debug:
                print("removing species", species, file=sys.stderr)
            del info[species]
            trans = list(
                filter(lambda t: t.get("qual:id") == "tr_" + species, transitions)
            )
            if trans:
                if debug:
                    print("removing transition", trans[0], file=sys.stderr)
                tlist.remove(trans[0])


def simplify_model(info, debug: bool = False):
    """Clean the model w.r.t. some active/inactive species."""
    multispecies: Dict[str, List[str]] = {}
    for key, value in list(info.items()):
        if not key.startswith("csa") and not key.startswith("sa"):
            del info[key]
            if debug:
                print("deleting complex:", key, "value:", value, file=sys.stderr)
            if len(value) > 1:
                multispecies[key] = value
                # print('multi', key, value)
    for _key, value in multispecies.items():
        for val in value:
            # check that it does not appear in any other reaction than the
            # activation one
            if debug:
                print("looking at multispecies:", val, file=sys.stderr)
            active = None
            for species, data in info.items():
                for trans in data["transitions"]:
                    if val in trans.reactants or val in [
                        mod[0] for mod in trans.modifiers
                    ]:
                        if active is None:
                            active = species
                        else:
                            active = False
                if active is False:
                    break
            if not info[val]["transitions"] and active in value:
                # we know that active is a str here, since it is in value
                add_rdf(info, cast(str, active), info[val]["annotations"])
                # print('deleting {val} [{active} is active for {key}]'.format(
                #     val=val,
                #     active=active,
                #     key=key,
                # ))
                del info[val]


def add_qual_species(
    layout: etree.Element, qlist: etree.Element, info, ginsim_names: bool
):
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
                "qual:name": fix_name(data["name"], species, ginsim_names),
                "qual:constant": constant,
                "qual:id": species,
            },
        )
        add_annotation(qspecies, data["annotations"])


def fix_name(name: str, species: str, ginsim_names: bool):
    """Change name for GINSIM compatibility or to remove subscripts."""
    if ginsim_names:
        # ginsim bug uses name as id
        return name.replace(" ", "_").replace(",", "").replace("/", "_") + "_" + species
    return name.replace("_sub_", "").replace("_endsub_", "")


def add_annotation(node: etree.Element, rdf: Optional[etree.Element]):
    """Add a single RDF element as an annotation node."""
    if rdf is not None:
        etree.SubElement(node, "annotation").append(rdf)


def add_transitions(tlist: etree.Element, info, graph: nx.Graph):
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
    prefix_len = len(NS["xhtml"]) + 2
    for reaction in transitions:
        if reaction.notes is not None:
            some_notes = True
            reaction.notes.tag = "p"
            for element in reaction.notes.getiterator():
                if element.tag.startswith("{" + NS["xhtml"] + "}"):
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


def add_inputs(
    ilist: etree.Element,
    transitions: List[Transition],
    species: str,
    known: List[str],
    graph: nx.Graph,
):
    """Add all known inputs."""
    index = 0
    modifiers: List[Tuple[str, str]] = []
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
            if (modifier, sign) not in modifiers and modifier in known:
                modifiers.append((modifier, sign))
                graph.add_edge(species, modifier)
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
    with open(sbml_filename + ".csv", "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        for species, data in info.items():
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-g", "--ginsim", action="store_true")
    parser.add_argument("-d", "--debug", action="store_true")
    parser.add_argument(
        "-c",
        "--csv",
        action="store_true",
        help="Store the species information in a separate CSV file",
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

    if args.debug:
        print("parsing", args.infile.name, "…", file=sys.stderr)
    info, width, height = read_celldesigner(args.infile, debug=args.debug)
    simplify_model(info, args.debug)
    if args.infile != sys.stdin and args.outfile == sys.stdout:
        args.outfile = os.path.splitext(args.infile.name)[0] + ".sbml"
    write_qual(
        args.outfile,
        info,
        width,
        height,
        ginsim_names=args.ginsim,
        remove=args.remove,
        debug=args.debug,
    )
    if args.csv and args.outfile != sys.stdout:
        write_csv(args.outfile, info)


if __name__ == "__main__":  # pragma: no cover
    main()
