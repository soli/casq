"""Read CellDesigner models.

Copyright (C) 2019, Sylvain.Soliman@inria.fr

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

import collections
import xml.etree.ElementTree as etree
from itertools import chain
from typing import IO, List, Optional

from loguru import logger  # type: ignore

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


def read_celldesigner(fileobj: IO):
    """Parse the given file."""
    root = etree.parse(fileobj).getroot()
    tag = root.tag
    if tag != "{" + NS["sbml"] + "}sbml":
        raise ValueError("Currently limited to SBML Level 2 Version 4")
    model = root.find("sbml:model", NS)
    if model is not None:
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
        (
            to_map[s] if s in to_map else s
            for s in filter(lambda t: t not in to_remove, name.split("_"))
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
            reference=decomplexify(species.get("id", ""), model, field="species"),
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
        logger.debug(
            "adding annotations from {rdfs} to {reference}",
            rdfs=rdfs[:],
            reference=reference,
        )
        description = nameconv[reference]["annotations"].find("./rdf:Description", NS)
        for element in rdfs[:]:
            if element not in description:
                description.append(element)
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
            decomplexify(reac.get("alias", ""), model)
            for reac in chain(
                annot.findall("./cd:baseReactants/cd:baseReactant", NS),
                annot.findall("./cd:listOfReactantLinks/cd:reactantLink", NS),
            )
        ]
        prods = [
            decomplexify(prod.get("alias", ""), model)
            for prod in chain(
                annot.findall("./cd:baseProducts/cd:baseProduct", NS),
                annot.findall("./cd:listOfProductLinks/cd:productLink", NS),
            )
        ]
        mods = [
            (mod.get("type"), decomplexify(mod.get("aliases", ""), model))
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
        + f'cd:speciesAlias[@{field}="{species}"]',
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
    if cd_modifications is None:
        return []
    return [
        mod.get("state", "") for mod in cd_modifications.findall("cd:modification", NS)
    ]
