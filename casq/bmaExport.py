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

from loguru import logger


class booleanFormulaBuilder:
    """Builds a formula for a boolean network encoded in BMA.

    The formula is a max applied to a series of reactions (transitions),
    so that if at least one reaction is active the variable becomes active
    """

    def __init__(self):
        """Init.

        reactant reflects the species that are required for formation
        modifier reflects the state of the modifiers
        previous is a function for other transitions
        """
        self.modifier = "1"
        self.reactant = "1"
        self.previous = "0"

    def function(self):
        """Return self.previous."""
        return self.previous

    def addActivator(self, vid):
        """Update transition to add an activator.

        A reaction may only take place if all reactants/activators are present.
        """
        self.reactant = "(min(var({vid}),{current}))".format(
            vid=vid, current=self.reactant
        )

    def addInhibitor(self, vid):
        """If any inhibitor is active, the reaction is stopped."""
        self.reactant = "(min(1-var({vid}),{current}))".format(
            vid=vid, current=self.reactant
        )

    def addTransition(self):
        """AddTransition."""
        self.reactant = "1"

    def addCatalysis(self, vidList):
        """All non-reactants, non-inhibitors in casq are treated as catalysts.

        If at least one catalyst is active, the reaction can proceed.
        This is achieved in BMA with a min function
        """
        base = "0"
        for vid in vidList:
            base = "(max(var({vid}),{base}))".format(vid=vid, base=base)
        self.modifier = "(min({base},{current}))".format(
            base=base, current=self.modifier
        )

    def addAnd(self, vidList):
        """All listed elements are required for firing."""
        base = "1"
        for vid in vidList:
            base = "(min(var({vid}),{base}))".format(vid=vid, base=base)
        self.modifier = "(min({base},{current}))".format(
            base=base, current=self.modifier
        )

    def finishTransition(self):
        """Add a single transition formula to the current state.

        The final formula is the min of the transition and the catalyst-modifiers
        The catalyst-modifiers default to 1, the transition defaults to 1
        Resets the transition formula to 1.
        """
        function = "(min({transition},{current}))".format(
            transition=self.reactant, current=self.modifier
        )
        self.previous = "(max({f},{old}))".format(f=function, old=self.previous)
        self.reactant = "1"
        self.modifier = "1"


class multiStateFormulaBuilder:
    """Builds a multistate formula.

    This is more simple as BMA defaults to avg(pos)-avg(neg).
    """

    def __init__(self):
        """Init."""
        self.value = ""

    def function(self):
        """Return self.value."""
        return self.value

    def addActivator(self, vid):
        """Do nothing."""
        pass

    def addInhibitor(self, vid):
        """Do nothing."""
        pass

    def addAnd(self, vid):
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


# hardcoded colour codes so elements still follow BMA colourscheme
# Default pink #ff66cc
COLOURMAP = {0: "#ff66cc", 1: "#33cc00", 2: "#ff9900", 3: "#9966ff", 4: "#00cccc"}


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
        logger.debug(
            item + ", varid = " + str(idMap[item]) + ", name = " + info[item]["name"]
        )
        # skip if there are no transitions
        if len(info[item]["transitions"]) == 0:
            logger.debug(item + "-No transitions")
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
            logger.debug(item + "\tReactants:\t" + str(transition[1]))
            # reactant
            for reactant in transition[1]:
                if ignoreSelfLoops and reactant == product:
                    continue
                if reactant in idMap:
                    if transition.type in (
                        "INHIBITION",
                        "NEGATIVE_INFLUENCE",
                        "UNKNOWN_INHIBITION",
                    ):
                        relationships.append(
                            bma_relationship(
                                reactant, product, idMap, count, "Inhibitor"
                            )
                        )
                        formula.addInhibitor(idMap[reactant])
                    else:
                        relationships.append(
                            bma_relationship(reactant, product, idMap, count)
                        )
                        formula.addActivator(idMap[reactant])
                else:
                    pass
            # now modifiers
            if len(transition[2]) == 0:
                formula.finishTransition()
                continue
            modifiers = transition[2]
            # catalysts are a special case
            catalysts = []
            inhibitors = []
            # List of variables that should be omitted from formulae due to other function
            ignoreList = []
            # everything else
            logger.debug(str(modifiers))
            for (impact, m) in modifiers:
                if ignoreSelfLoops and m == product:
                    continue
                if impact == "BOOLEAN_LOGIC_GATE_AND":
                    logger.debug("Found an AND gate")
                    # indicates that the listed vars will be anded
                    # BAH: should I remove them from the cat code? In this instance its harmless...
                    vidList = [idMap[jtem] for jtem in m.split(",") if jtem in idMap]
                    formula.addAnd(vidList)
                    for jtem in vidList:
                        ignoreList.append(jtem)
                if m in idMap:
                    if impact == "UNKNOWN_INHIBITION" or impact == "INHIBITION":
                        relationships.append(
                            bma_relationship(m, product, idMap, count, "Inhibitor")
                        )

                        formula.addInhibitor(idMap[m])
                        inhibitors.append(idMap[m])
                    else:
                        # treat all other modifiers as catalysts (casq approach)
                        logger.debug(item + "\tFound impact:" + impact)
                        catalysts.append(idMap[m])
                        relationships.append(bma_relationship(m, product, idMap, count))
                else:
                    pass
            logger.debug(item + "\tCatalysts\t" + str(catalysts))
            logger.debug(item + "\tInhibitors\t" + str(inhibitors))
            logger.debug(item + "\tIgnoreList\t" + str(ignoreList))
            # filter catalysts for items to be ignored
            finalCat = [item for item in catalysts if item not in ignoreList]
            if len(finalCat) > 0:
                formula.addCatalysis(finalCat)
            formula.finishTransition()
        allFormulae[item] = formula.function()
    return (relationships, allFormulae)


def translateGreek(name):
    """Translate Greek to Latin alphabet."""
    greek_alphabet = "ΑαΒβΓγΔδΕεΖζΗηΘθΙιΚκΛλΜμΝνΞξΟοΠπΡρΣσςΤτΥυΦφΧχΨψΩω"
    latin_alphabet = "AaBbGgDdEeZzHhJjIiKkLlMmNnXxOoPpRrSssTtUuFfQqYyWw"
    greek2latin = str.maketrans(greek_alphabet, latin_alphabet)
    return name.translate(greek2latin)


def depunctuate(name):
    """Replace punctuation by underscores."""
    badChars = " ,-()+:/\\'[]><"
    alternatives = "______________"
    cleanup = str.maketrans(badChars, alternatives)
    return name.translate(cleanup)


def cleanName(name):
    """Remove punctuation and replace Greek letters."""
    noPunctuation = depunctuate(name)
    result = translateGreek(noPunctuation)
    return result


def bma_model_variable(vid, infoVariable, formulaDict, v, granularity, inputLevel):
    """Return BMA model variable as a dict."""
    if v in formulaDict:
        formula = formulaDict[v]
    else:
        # Assume that nodes with no incoming edges are active
        if inputLevel is None:
            formula = str(granularity)
        else:
            formula = str(inputLevel)
    result = {
        "Name": cleanName(infoVariable["name"]),
        "Id": vid,
        "RangeFrom": 0,
        "RangeTo": granularity,
        "Formula": formula,
    }
    return result


def bma_layout_variable(vid, infoVariable, fill=None, description=""):
    """Return BMA layout variable as a dict."""
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


def write_bma(
    filename: str,
    info,
    granularity=1,
    inputLevel=None,
    ignoreSelfLoops=False,
    colourByCompartment=True,
):
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
        if colourByCompartment:
            compartmentColour[compList[i][0]] = COLOURMAP.get(i)
        else:
            compartmentColour[compList[i][0]] = COLOURMAP.get(0)

    idGenerator = itertools.count(1)
    idMap = {k: next(idGenerator) for k in info.keys()}

    rm, formula = get_relationships(
        info, idMap, idGenerator, granularity, ignoreSelfLoops
    )

    logger.debug(formula)

    vm = [
        bma_model_variable(idMap[v], info[v], formula, v, granularity, inputLevel)
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
