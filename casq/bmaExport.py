"""Convert CellDesigner models to BMA json.

Copyright (C) 2021 b.hall@ucl.ac.uk

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

import itertools
import json


class booleanFormulaBuilder:
    """Builds a boolean formula."""

    def __init__(self):
        """Init."""
        self.value = "0"
        self.transition = "1"

    def addActivator(self, vid):
        """AddActivator."""
        self.transition = "(min(var({vid}),{current}))".format(
            vid=vid, current=self.transition
        )

    def addInhibitor(self, vid):
        """AddInhibitor."""
        self.transition = "(min(1-var({vid}),{current}))".format(
            vid=vid, current=self.transition
        )

    def addTransition(self):
        """AddTransition."""
        self.transition = "1"

    def addCatalysis(self, vidList):
        """AddCatalysis."""
        base = "0"
        for vid in vidList:
            base = "(max(var({vid}),{base}))".format(vid=vid, base=base)
        self.value = "(min({base},{current}))".format(base=base, current=self.value)

    def finishTransition(self):
        """FinishTransition."""
        self.value = "(max({transition},{current}))".format(
            transition=self.transition, current=self.value
        )
        self.transition = "1"


class multiStateFormulaBuilder:
    """Builds a multistate formula."""

    def __init__(self):
        """Init."""
        self.value = ""

    def addActivator(self, vid):
        """Do nothing."""
        pass

    def addInhibitor(self, vid):
        """Do nothing."""
        pass

    def addCatalysis(self, vidList):
        """Do nothing."""
        pass

    def addTransition(self):
        """Do nothing."""
        pass

    def finishTransition(self):
        """Do nothing."""
        pass


COLOURMAP = {
    0: "BMA_Green",
    1: "BMA_Orange",
    2: "BMA_Purple",
    3: "BMA_Mint",
}


def bma_relationship(source, target, idMap, count, which="Activator"):
    """Return BMA relationship dict."""
    result = {
        "ToVariable": idMap[target],
        "Type": which,
        "FromVariable": idMap[source],
        "Id": next(count),
    }
    return result


def get_relationships(info, idMap, count, granularity, ignoreSelfLoops):
    """Return all BMA relationships."""
    relationships = []
    allFormulae = {}
    for item in info.keys():
        # skip if there are no transitions
        if len(info[item]["transitions"]) == 0:
            continue
        product = item
        if granularity == 1:
            formula = booleanFormulaBuilder()
        else:
            formula = multiStateFormulaBuilder()
        # variables may be missing from the "simplified" model.
        # Test for variable in the ID map before appending
        for transition in info[item]["transitions"]:
            formula.addTransition()
            # reactant
            for reactant in transition[1]:
                if ignoreSelfLoops and reactant == product:
                    continue
                if reactant in idMap:
                    relationships.append(
                        bma_relationship(reactant, product, idMap, count)
                    )
                    formula.addActivator(idMap[reactant])
                else:
                    pass
            # now modifiers
            if len(transition[2]) == 0:
                continue
            modifiers = transition[2]
            # catalysts are a special case
            catalysts = []
            for (impact, m) in modifiers:
                if ignoreSelfLoops and m == product:
                    continue
                if m in idMap:
                    if impact == "CATALYSIS" or impact == "UNKNOWN_CATALYSIS":
                        catalysts.append(idMap[m])
                        relationships.append(bma_relationship(m, product, idMap, count))
                    else:
                        pass
                else:
                    pass
            formula.addCatalysis(catalysts)
            # everything else
            for (impact, m) in modifiers:
                if ignoreSelfLoops and m == product:
                    continue
                if m in idMap:
                    if impact == "UNKNOWN_INHIBITION" or impact == "INHIBITION":
                        relationships.append(
                            bma_relationship(m, product, idMap, count, "Inhibitor")
                        )
                        formula.addInhibitor(idMap[m])
                    elif impact == "CATALYSIS" or impact == "UNKNOWN_CATALYSIS":
                        pass  # deal with this earlier
                    else:
                        # treat all other modifiers as reactants
                        relationships.append(bma_relationship(m, product, idMap, count))
                        formula.addActivator(idMap[m])
                else:
                    pass
            formula.finishTransition()
        allFormulae[item] = formula.value
    return (relationships, allFormulae)


def translateGreek(name):
    """Translate Greek to Latin alphabet."""
    greek_alphabet = "ΑαΒβΓγΔδΕεΖζΗηΘθΙιΚκΛλΜμΝνΞξΟοΠπΡρΣσςΤτΥυΦφΧχΨψΩω"
    latin_alphabet = "AaBbGgDdEeZzHhJjIiKkLlMmNnXxOoPpRrSssTtUuFfQqYyWw"
    greek2latin = str.maketrans(greek_alphabet, latin_alphabet)
    return name.translate(greek2latin)


def depunctuate(name):
    """Replace punctuation by underscores."""
    badChars = " ,-()+:/\\'"
    alternatives = "__________"
    cleanup = str.maketrans(badChars, alternatives)
    return name.translate(cleanup)


def cleanName(name):
    """Remove punctuation and replace Greek letters."""
    noPunctuation = depunctuate(name)
    result = translateGreek(noPunctuation)
    return result


def bma_model_variable(vid, infoVariable, formulaDict, v, granularity):
    """Do something."""
    if v in formulaDict:
        formula = formulaDict[v]
    else:
        # Assume that nodes with no incoming edges are active
        formula = str(granularity)
    result = {
        "Name": cleanName(infoVariable["name"]),
        "Id": vid,
        "RangeFrom": 0,
        "RangeTo": granularity,
        "Formula": formula,
    }
    return result


def bma_layout_variable(vid, infoVariable, fill=None, description=""):
    """Do something else."""
    result = {
        "Id": vid,
        "Name": cleanName(infoVariable["name"]),
        "Type": "Constant",
        "ContainerId": 0,
        "PositionX": float(infoVariable["x"]),
        "PositionY": float(infoVariable["y"]),
        "CellY": 0,
        "CellX": 0,
        "Angle": 0,
        "Description": description,
    }
    if fill is not None:
        result["Fill"] = fill
    return result


def write_bma(filename: str, info, granularity=1, ignoreSelfLoops=False):
    # pylint: disable=too-many-arguments, too-many-locals
    """Write the BMA json with layout file for our model."""
    # granularity must be a non-zero natural
    assert granularity > 0

    # calculate the compartments for colours;
    # four largest compartments are coloured BMA colours, all else default
    compartments = {}
    for k in info.keys():
        location = info[k]["compartment"]
        if location in compartments:
            compartments[location] += 1
        else:
            compartments[location] = 1
    compList = list(compartments.items())
    compList.sort(key=lambda i: -i[1])
    compartmentColour = {}
    for i in range(len(compartments)):
        compartmentColour[compList[i][0]] = COLOURMAP.get(i)

    idGenerator = itertools.count(1)
    idMap = {k: next(idGenerator) for k in info.keys()}

    rm, formula = get_relationships(
        info, idMap, idGenerator, granularity, ignoreSelfLoops
    )

    vm = [
        bma_model_variable(idMap[v], info[v], formula, v, granularity)
        for v in info.keys()
    ]
    vl = [
        bma_layout_variable(
            idMap[v],
            info[v],
            compartmentColour[info[v]["compartment"]],
            info[v]["compartment"],
        )
        for v in info.keys()
    ]

    model = {"Name": "CaSQ-BMA", "Variables": vm, "Relationships": rm}
    layout = {"Variables": vl, "Containers": [], "Description": ""}
    ltl = {"states": [], "operations": []}
    universe = {"Model": model, "Layout": layout, "ltl": ltl}

    json_object = json.dumps(universe, indent=4)
    with open(filename, "w") as outfile:
        outfile.write(json_object)
