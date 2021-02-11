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

def bma_relationship(source,target,idMap,count,which="Activator"):
    result = {
                'ToVariable': idMap[target],
                'Type': which, 
                'FromVariable': idMap[source], 
                'Id': count.next()
                }
    return(result)

def get_relationships(info,idMap,count):
    relationships = []
    formulae = {}
    for item in info.keys():
        #skip if there are no transitions
        if len(info[item]["transitions"])==0: continue
        product = item
        #variables may be missing from the "simplified" model. Test for variable in the ID map before appending
        for transition in info[item]["transitions"]:
            #reactant
            for reactant in transition[1]:
                if reactant in idMap:
                    relationships.append(bma_relationship(reactant,product,idMap,count))
                else:
                    pass
            #now modifiers
            if len(transition[2]) == 0: continue
            modifiers = transition[2]
            for (impact,m) in modifiers:
                if m in idMap:
                    if impact == "UNKNOWN_INHIBITION" or impact == "INHIBITION" :
                        relationships.append(bma_relationship(m,product,idMap,count,"Inhibitor"))
                    else:
                        relationships.append(bma_relationship(m,product,idMap,count))
                else:
                    pass
    return(relationships)

def bma_model_variable(vid, infoVariable):
    result = {
        'Name':infoVariable["name"],
        'Id':vid,
        'RangeFrom':0,
        'RangeTo':1,
        'Formula': ""
        }
    return(result)

def bma_layout_variable(vid, infoVariable):
    result = {
        'Id': vid,
        'Name':infoVariable["name"],
        'Type':"Constant",
        'ContainerId':0,
        'PositionX':float(infoVariable["x"]),
        'PositionY':float(infoVariable["y"]),
        'CellY':0,
        'CellX':0,
        'Angle':0,
        'Description':"",
        
        }
    return(result)

def write_bma(
    filename: str, info
):
    # pylint: disable=too-many-arguments, too-many-locals
    """Write the BMA json with layout file for our model."""
    idGenerator = counter(1)
    idMap = dict([(k,idGenerator.next()) for k in info.keys()])

    vm = [bma_model_variable(idMap[v],info[v]) for v in info.keys()]
    vl = [bma_layout_variable(idMap[v],info[v]) for v in info.keys()]

    rm = get_relationships(info, idMap, idGenerator)

    #add_transitions(tlist, info, graph)
    #add_qual_species(layout, qlist, info)

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