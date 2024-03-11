"""Simplify CellDesigner models.

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
from typing import cast  # noqa: F401

import networkx as nx  # type: ignore
from loguru import logger  # type: ignore

from .readCD import Transition, add_rdf


def simplify_model(info, upstream, downstream, names_as_ids: bool = False):
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
    fix_all_names(info)
    if names_as_ids:
        use_names_as_ids(info)
    restrict_model(info, upstream, downstream)
    handle_phenotypes(info)


def handle_phenotypes(info):
    """Restructure all reactions targetting a phenotype as a single one."""
    for key, data in list(info.items()):
        if data["type"] != "PHENOTYPE":
            continue
        transitions = data["transitions"]
        modifiers = []
        new_transitions = []
        for t in transitions:
            if len(t.reactants) != 1:
                logger.debug(
                    "ignoring non-unary reaction to phenotype {pheno}",
                    pheno=data["name"],
                )
                new_transitions.append(t)
                continue
            if t.type == "NEGATIVE_INFLUENCE":
                modifiers.append(("INHIBITION", t.reactants[0]))
            else:
                modifiers.append(("CATALYSIS", t.reactants[0]))
        if modifiers:
            new_transitions.append(
                Transition("STATE_TRANSITION", [], modifiers, None, None)
            )
            info[key]["transitions"] = new_transitions


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
                        and not info[reac1]["transitions"]
                        and not info[reac2]["transitions"]
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
                        # info[key]["transitions"].remove(trans)
        else:
            duplicate_nodes[value["ref_species"]] = key
    replace_in_transitions(info, replacements)
    return multispecies


def restrict_model(info, upstream, downstream):
    """Only keep species upstream/downstream of some list of species."""
    name_to_ids = {v["name"]: k for (k, v) in info.items()}
    for name in upstream + downstream:
        if name not in name_to_ids:
            logger.error(name + " was not found, maybe it is ambiguous…")

    if upstream == [] and downstream == []:
        return
    graph = nx.DiGraph()
    for species, data in info.items():
        graph.add_node(species)
        for trans in data["transitions"]:
            for val in trans.reactants:
                graph.add_edge(val, species)
            for _modtype, modlist in trans.modifiers:
                for val in modlist.split(","):
                    graph.add_edge(val, species)
    keep = set()
    for dnname in downstream:
        dn = name_to_ids[dnname]
        keep |= {
            elt for succl in nx.dfs_successors(graph, dn).values() for elt in succl
        } | {dn}
    ggraph = graph.reverse()
    for upname in upstream:
        up = name_to_ids[upname]
        keep |= {
            elt for succl in nx.dfs_successors(ggraph, up).values() for elt in succl
        } | {up}
    for species in list(info.keys()):
        if species not in keep:
            del info[species]


def replace_in_transitions(info, replacements):
    """Change transitions in info to reflect replacements."""
    for _species, data in info.items():
        for trans in data["transitions"]:
            for val in trans.reactants.copy():
                if val in replacements:
                    trans.reactants.append(replacements[val])
                    trans.reactants.remove(val)
            for modtype, mod_list in trans.modifiers.copy():
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
        for old, new in replacements.items():
            data["function"] = data["function"].replace(old, new)


def get_active(val, info):
    """Find who val activates."""
    active = None
    for species, data in info.items():
        if species.startswith("csa") or species.startswith("sa"):
            for trans in data["transitions"]:
                if val in trans.reactants or val in (
                    mod
                    for _modtype, modifier_list in trans.modifiers
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
                if activity == "active" and (name + "_active") not in namedict:
                    name = name + "_active"
                elif other_activity == "active" and (name + "_active") not in namedict:
                    info[other_id]["name"] = name + "_active"
                    info[other_id]["function"] = name + "_active"
                else:
                    newname = name
                    tag = 0
                    while newname in namedict:
                        tag += 1
                        newname = name + "_" + str(tag)
                    name = newname
            namedict[name] = (species, activity)
        data["name"] = name
        data["function"] = name


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


def use_names_as_ids(info):
    """Replace all ids with names."""
    newinfo = {}
    replacements = {}
    for key, data in info.items():
        oname = data["name"]
        name = oname.replace(" ", "_")
        name = "".join(c for c in name if c.isalnum() or c == "_")
        newinfo[name] = data
        newinfo[name]["name"] = name
        replacements[oname] = name
        replacements[key] = name
    info.clear()
    info.update(newinfo)
    replace_in_transitions(info, replacements)
