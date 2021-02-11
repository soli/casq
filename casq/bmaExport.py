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
                'ToVariable': idMap[source],
                'Type': which, 
                'FromVariable': idMap[target], 
                'Id': count.next()
                }
    return(result)

def get_relationships(info,idMap,count):
    relationships = []
    #reactants first
    for item in info.keys():
        product = item
        print(item)
        print(info[item])
        print(info[item]["transitions"][0])
        print(info[item]["transitions"][0][1])
        reactants = info[item]["transitions"][0][1]
        for r in reactants:
            relationships.append(bma_relationship(r,product,idMap,count))
    #Now modifiers
    for item in info.keys():
        product = item
        if len(info[item]["transitions"][0][2]) == 0: continue
        modifiers = info[item]["transitions"][0][2]
        print(modifiers)
        for (impact,m) in modifiers:
            if impact == "Inhibition":
                relationships.append(bma_relationship(m,product,idMap,count,"Inhibitor"))
            else:
                relationships.append(bma_relationship(m,product,idMap,count))
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