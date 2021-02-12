"""Convert CellDesigner models to BMA json

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

import json

class counter():
    def __init__(self,first):
        self.value = first
    def next(self):
        result = self.value
        self.value += 1
        return(result)

class booleanFormulaBuilder():
    def __init__(self):
        self.value = "0"
        self.transition ="1"
    def addActivator(self,vid):
        self.transition = "(min(var({vid}),{current}))".format(vid = vid, current = self.transition)
    def addInhibitor(self,vid):
        self.transition = "(min(1-var({vid}),{current}))".format(vid = vid, current = self.transition)
    def addTransition(self):
        self.transition = "1"
    def finishTransition(self):
        self.value = "(max({transition},{current}))".format(transition=self.transition, current=self.value)
        self.transition = "1"

class multiStateFormulaBuilder():
    def __init__(self):
        self.value = ""
    def addActivator(self,vid):
        pass
    def addInhibitor(self,vid):
        pass
    def addTransition(self):
        pass
    def finishTransition(self):
        pass

def colourMap(n):
    if   n == 0: return("BMA_Green")
    elif n == 1: return("BMA_Orange")
    elif n == 2: return("BMA_Purple")
    elif n == 3: return("BMA_Mint")
    else: return(None)

def bma_relationship(source,target,idMap,count,which="Activator"):
    result = {
                'ToVariable': idMap[target],
                'Type': which, 
                'FromVariable': idMap[source], 
                'Id': count.next()
                }
    return(result)

def get_relationships(info,idMap,count,granularity):
    relationships = []
    allFormulae = {}
    for item in info.keys():
        #skip if there are no transitions
        if len(info[item]["transitions"])==0: continue
        product = item
        if granularity == 1: 
            formula = booleanFormulaBuilder()
        else:
            formula = multiStateFormulaBuilder()
        #variables may be missing from the "simplified" model. Test for variable in the ID map before appending
        for transition in info[item]["transitions"]:
            formula.addTransition()
            #reactant
            for reactant in transition[1]:
                if reactant in idMap:
                    relationships.append(bma_relationship(reactant,product,idMap,count))
                    formula.addActivator(idMap[reactant])
                else:
                    pass
            #now modifiers
            if len(transition[2]) == 0: continue
            modifiers = transition[2]
            for (impact,m) in modifiers:
                if m in idMap:
                    if impact == "UNKNOWN_INHIBITION" or impact == "INHIBITION" :
                        relationships.append(bma_relationship(m,product,idMap,count,"Inhibitor"))
                        formula.addInhibitor(idMap[m])
                    else:
                        relationships.append(bma_relationship(m,product,idMap,count))
                        formula.addActivator(idMap[m])
                else:
                    pass
            formula.finishTransition()
        allFormulae[item] = formula.value
    return(relationships, allFormulae)

def cleanName(name):
    result= name.replace(' ' , '_').replace(',' , '_').replace('-', '_').replace('(','').replace(')','').replace('+','').replace(':','').replace('/','').replace('\\','')
    return(result)

def bma_model_variable(vid, infoVariable, formulaDict, v, granularity):
    if v in formulaDict:
        formula=formulaDict[v]
    else:
        #Assume that nodes with no incoming edges are active
        formula = str(granularity)
    result = {
        'Name':cleanName(infoVariable["name"]),
        'Id':vid,
        'RangeFrom':0,
        'RangeTo':granularity,
        'Formula': formula
        }
    return(result)

def bma_layout_variable(vid, infoVariable,fill=None):
    result = {
        'Id': vid,
        'Name':cleanName(infoVariable["name"]),
        'Type':"Constant",
        'ContainerId':0,
        'PositionX':float(infoVariable["x"]),
        'PositionY':float(infoVariable["y"]),
        'CellY':0,
        'CellX':0,
        'Angle':0,
        'Description':"",
        }
    if fill != None: result['Fill'] = fill
    return(result)

def write_bma(
    filename: str, info, granularity=1
):
    # pylint: disable=too-many-arguments, too-many-locals
    """Write the BMA json with layout file for our model."""
    #granularity must be a non-zero natural
    assert(granularity>0)

    #calculate the compartments for colours; four largest compartments are coloured BMA colours, all else default
    compartments = {}
    for k in info.keys():
        location = info[k]['compartment']
        if location in compartments:
            compartments[location] += 1
        else:
            compartments[location] = 1
    compList = list(compartments.items())
    compList.sort(key= lambda i: -i[1])
    compartmentColour = {}
    for i in range(len(compartments)):
        compartmentColour[compList[i][0]] = colourMap(i)

    idGenerator = counter(1)
    idMap = dict([(k,idGenerator.next()) for k in info.keys()])

    rm,formula = get_relationships(info, idMap, idGenerator, granularity)

    vm = [bma_model_variable(idMap[v],info[v],formula,v,granularity) for v in info.keys()]
    vl = [bma_layout_variable(idMap[v],info[v],compartmentColour[info[v]['compartment']]) for v in info.keys()]

    model = {"Name":"CaSQ-BMA","Variables":vm,"Relationships":rm}
    layout = {"Variables":vl,"Containers":[],"Description":""}
    ltl = {
        "states": [],
        "operations": []
    }
    universe = {"Model":model,"Layout":layout,'ltl':ltl}

    json_object = json.dumps(universe, indent = 4) 
    with open(filename,"w") as outfile:
        outfile.write(json_object)