"""Read SBGN-ML models.

Copyright (C) 2025, issa.kerima-khalil@utoulouse.fr

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

import unicodedata
import xml.etree.ElementTree as etree
from typing import IO

from loguru import logger  # type: ignore

from .readCD import NS, Transition, add_rdf, make_name_precise


def read_sbgnml(fileobj: IO):
    """Parse the given SBGN-ML file, return info and canvas size."""
    root = etree.parse(fileobj).getroot()
    if root.tag != "{" + NS["sbgn"] + "}sbgn":
        raise ValueError("Expected sbgn root element")

    map_info = root.find(".//sbgn:map", namespaces=NS)
    if map_info is None:
        raise ValueError("Could not find sbgn:map element")

    sizeX = map_info.get("width", "1000")
    sizeY = map_info.get("height", "1000")

    nameconv = species_info_sbgn(map_info)
    add_subcomponents_only_sbgn(nameconv, map_info)
    info = get_transitions_sbgn(map_info, nameconv)

    return info, sizeX, sizeY


def species_info_sbgn(map_element):
    """Create a map from species' ids to their attributes."""
    nameconv = {}
    compartments = {}

    for glyph in map_element.findall("sbgn:glyph", namespaces=NS):
        if glyph.get("class") == "compartment":
            cid = glyph.get("id")
            compartments[cid] = {
                "bbox": glyph.find("sbgn:bbox", namespaces=NS),
                "name": glyph.find("sbgn:label", namespaces=NS).get("text")
                if glyph.find("sbgn:label", namespaces=NS) is not None
                else cid,
            }

    for glyph in map_element.findall("sbgn:glyph", namespaces=NS):
        cls = glyph.get("class")
        if cls not in (
            "macromolecule",
            "simple chemical",
            "unspecified entity",
            "complex",
            "nucleic acid feature",
            "phenotype",
            "macromolecule multimer",
        ):
            logger.warning(
                "Skipping glyph with unsupported class = '{}' : id={}",
                cls,
                glyph.get("id"),
            )
            continue

        species_id = glyph.get("id")
        label = glyph.find("sbgn:label", namespaces=NS)
        name = label.get("text") if label is not None else species_id

        bbox = glyph.find("sbgn:bbox", namespaces=NS)
        x = float(bbox.get("x"))
        y = float(bbox.get("y"))
        w = float(bbox.get("w"))
        h = float(bbox.get("h"))

        compartment_name = "default_compartment"
        cx, cy = x + w / 2, y + h / 2
        min_area = float("inf")
        for _cid, comp in compartments.items():
            cbbox = comp["bbox"]
            cx0 = float(cbbox.get("x"))
            cy0 = float(cbbox.get("y"))
            cw = float(cbbox.get("w"))
            ch = float(cbbox.get("h"))
            if cx0 <= cx <= cx0 + cw and cy0 <= cy <= cy0 + ch:
                area = cw * ch
                if area < min_area:
                    min_area = area
                    compartment_name = comp["name"].replace(" ", "_")

        classtype = class_to_type(cls)

        # get the state (activity) from the sub-glyphs of a state variable
        activity = "inactive"  # valeur par défaut
        for state_var in glyph.findall(
            "sbgn:glyph[@class='state variable']", namespaces=NS
        ):
            state_elem = state_var.find("sbgn:state", namespaces=NS)
            if state_elem is not None:
                value = state_elem.get("value")
                if value:
                    activity = value  # "active", "inactive" there is also "P", probably phosphorylation
                    break
        mods = []

        name_clean = make_name_precise(greeks_to_name(name), classtype, mods)
        # ref_species = make_name_precise(name, classtype, [])  # sans mods
        # if classtype == "COMPLEX" and ref_species.endswith("_complex"):
        #    ref_species = ref_species[:-8]
        # annotation = glyph.find("sbgn:annotation", namespaces=NS)
        # rdf = annotation.find(".//rdf:RDF", namespaces=NS) if annotation is not None else None
        rdf = glyph.find(".//rdf:RDF", namespaces=NS)
        logger.debug(
            "Adding entity: id={}, type={}, name={}", species_id, classtype, name_clean
        )
        nameconv[species_id] = {
            "activity": activity,
            "x": str(x),
            "y": str(y),
            "h": str(h),
            "w": str(w),
            "transitions": [],
            "name": name_clean,
            "function": name_clean,
            # "ref_species": species_id,
            "ref_species": f"{name_clean}__{compartment_name}__{activity}",
            "type": classtype,
            "modifications": mods,
            "receptor": False,
            "annotations": rdf,
            "compartment": compartment_name,
        }
        # --- Reverse mapping : key "__Nom" → list of IDs corresponding to this name ---
        prot_ref = "__" + name_clean.rsplit("_", 1)[0]

        if prot_ref in nameconv:
            nameconv[prot_ref].append(species_id)
        else:
            nameconv[prot_ref] = [species_id]
    return nameconv


def get_transitions_sbgn(map_element, info):
    """Find all transitions."""
    port_to_reaction = {}
    reactions = {}
    logic_gates = {}
    port_to_gate = {}

    for gate in map_element.findall("sbgn:glyph[@class='and']", namespaces=NS):
        gid = gate.get("id")
        for port in gate.findall("sbgn:port", namespaces=NS):
            pid = port.get("id")
            port_to_gate[pid] = gid
        logic_gates[gid] = {"type": "and", "inputs": [], "output": None}

    for glyph in map_element.findall("sbgn:glyph", namespaces=NS):
        if glyph.get("class") in {
            "process",
            "association",
            "dissociation",
            "omitted process",
            "uncertain process",
        }:
            rid = glyph.get("id")
            logger.debug("Parsing process: id={}", rid)

            for port in glyph.findall("sbgn:port", namespaces=NS):
                port_to_reaction[port.get("id")] = rid

            glyph_class = glyph.get("class")
            if glyph_class == "association":
                label = "HETERODIMER_ASSOCIATION"
            elif glyph_class == "dissociation":
                label = "DISSOCIATION"
            elif glyph_class == "process":
                label_elem = glyph.find("sbgn:label", namespaces=NS)
                label = label_elem.get("text") if label_elem is not None else "PROCESS"
            else:
                label = "UNKNOWN"
            logger.debug(
                "Process {} (class='{}') assigned label '{}'", rid, glyph_class, label
            )
            # Get notes (HTML or XHTML text) and RDF annotations
            rdf = None
            notes_elem = glyph.find("libsbgn:notes", namespaces=NS)
            # notes_elem = glyph.find("sbgn:notes", namespaces=NS)
            notes = None
            if notes_elem is not None:
                notes = notes_elem.find(".//xhtml:body", namespaces=NS)
            rdf = glyph.find(".//rdf:RDF", namespaces=NS)
            reactions[rid] = {
                "inputs": [],
                "outputs": [],
                "modifiers": [],
                "label": label,
                "notes": notes,
                "annotations": rdf,
            }

    for arc in map_element.findall("sbgn:arc", namespaces=NS):
        source = arc.get("source")
        target = arc.get("target")
        aclass = arc.get("class")
        logger.debug(
            "Parsing arc: class={}, source={}, target={}", aclass, source, target
        )

        source_gate = port_to_gate.get(source, source)
        target_gate = port_to_gate.get(target, target)

        if target_gate in logic_gates:
            logic_gates[target_gate]["inputs"].append(source_gate)

        elif source_gate in logic_gates:
            logic_gates[source_gate]["output"] = target_gate

        elif source in port_to_reaction:
            rid = port_to_reaction[source]
            rxn = reactions[rid]
            if aclass == "production":
                rxn["outputs"].append(target)

        elif target in port_to_reaction or target in reactions:
            if target in port_to_reaction:
                rid = port_to_reaction[target]
            else:
                rid = target
            rxn = reactions[rid]
            if aclass == "consumption":
                rxn["inputs"].append((source, "positive"))
            else:
                modifier_type = classify_arc_modifier(aclass)
                if modifier_type != "UNKNOWN":
                    rxn["modifiers"].append((modifier_type, source))
                else:
                    logger.warning("Unknown modifier arc class: {}", aclass)

        elif target in info:
            modifier_type = classify_arc_modifier(aclass)
            if modifier_type != "UNKNOWN":
                logger.debug(
                    "Creating implicit transition for modifier arc: {} → {}",
                    source,
                    target,
                )
                info[target]["transitions"].append(
                    Transition(
                        type=f"{modifier_type}",
                        reactants=[source],
                        modifiers=[],
                        notes=None,
                        annotations=None,
                    )
                )

    for _lgid, gate in logic_gates.items():
        rid = gate["output"]
        if rid in reactions:
            for mid in gate["inputs"]:
                reactions[rid]["modifiers"].append(("BOOLEAN_LOGIC_GATE_AND", mid))

    # Build a Transition tuple and add it to each product
    for _rid, rxn in reactions.items():
        reactants = [r for r, sign in rxn["inputs"]]
        modifiers = rxn["modifiers"]
        notes = rxn["notes"]
        annotations = rxn["annotations"]
        ttype = rxn["label"]

        for output in rxn["outputs"]:
            if output in info:
                logger.debug(
                    "Assigning transition: reaction={} to product={}", ttype, output
                )
                info[output]["transitions"].append(
                    Transition(
                        type=ttype,
                        reactants=reactants,
                        modifiers=modifiers,
                        notes=notes,
                        annotations=annotations,
                    )
                )

    return info


def find_reaction_for_arc(arc, map_element):
    """Return the reaction to which an arc is connected."""
    source = arc.get("source")
    target = arc.get("target")

    for glyph in map_element.findall("sbgn:glyph[@class='process']", namespaces=NS):
        if glyph.get("id") == source or glyph.get("id") == target:
            return glyph
    return None


# def make_name_precise(name, ctype, mods):
#     if name is None:
#         return ""
#     name=name.replace(" ", "_espace_").replace("/", "_slash_")
#     to_map = {"&": "", "|": "", "!": "", "underscore": ""}
#     to_remove = {"sub", "endsub"}
#     parts = [s for s in name.split("_") if s not in to_remove]
#     newname = "_".join(to_map.get(s, s) for s in parts).replace("__", "_")
#     if ctype == "PROTEIN":
#         basis = [newname]
#     else:
#         basis = [newname, ctype.lower()]
#     return "_".join(basis + mods)

# def add_rdf(nameconv, reference: str, new_rdf: Optional[etree.Element]):
#     if new_rdf is None or reference not in nameconv:
#         return
#     if nameconv[reference]["annotations"] is not None:
#         rdfs = new_rdf.find("./rdf:Description", namespaces=NS)
#         if rdfs is None:
#             return
#         description = nameconv[reference]["annotations"].find("./rdf:Description", namespaces=NS)
#         if description is None:
#             return
#         for element in rdfs:
#             if element not in description:
#                 description.append(element)
#         logger.debug("Merging RDF annotations for {}", reference)
#     else:
#         nameconv[reference]["annotations"] = new_rdf
#         logger.debug("Assigning RDF annotations to {}", reference)


def add_subcomponents_only_sbgn(nameconv, map_element):
    """Add annotations to the parent complex.

    For unused SBGN species (only subcomponents of complexes)
    """
    for glyph in map_element.findall("sbgn:glyph", namespaces=NS):
        parent_ref = glyph.get("parent")
        if parent_ref and parent_ref in nameconv:
            annotation = glyph.find("sbgn:annotation", namespaces=NS)
            rdf = (
                annotation.find(".//rdf:RDF", namespaces=NS)
                if annotation is not None
                else None
            )
            if rdf is not None:
                add_rdf(nameconv, parent_ref, rdf)


def classify_arc_modifier(aclass):
    """Translate SBGN modifications to CS ones."""
    mapping = {
        "stimulation": "STIMULATION",
        "positive influence": "POSITIVE_INFLUENCE",
        "absolute stimulation": "STIMULATION",
        "inhibition": "INHIBITION",
        "negative influence": "NEGATIVE_INFLUENCE",
        "absolute inhibition": "INHIBITION",
        "modulation": "MODULATION",
        "catalysis": "CATALYSIS",
        "necessary stimulation": "NECESSARY_STIMULATION",
        "assignment": "ASSIGNMENT",
        "unknown influence": "UNKNOWN_INFLUENCE",
    }
    return mapping.get(aclass, "UNKNOWN")


def class_to_type(cls: str) -> str:
    """Translate SBGN glyph classes to CS types."""
    cls = cls.lower().replace(" ", "_")
    if cls in (
        # "macromolecule",
        # "simple_chemical",
        "complex",
        "phenotype",
        "unspecified_entity",
        "nucleic_acid_feature",
        # "macromolecule_multimer",
    ):
        return cls.upper()
    return "PROTEIN"


def greeks_to_name(s: str) -> str:
    """Change all greek letters to ascii in s."""
    return "".join(greek_to_name(c) for c in s)


def greek_to_name(c: str) -> str:
    """Change one greek letter to ascii."""
    try:
        greek, size, letter, what, *with_tonos = unicodedata.name(c).split()
    except ValueError:  # not enough values to unpack
        return c
    if (greek, letter) != ("GREEK", "LETTER"):
        return c
    return what.lower() if size == "SMALL" else what.title()
