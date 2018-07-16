#!/usr/bin/env python3
'''
Convert CellDesigner models to SBML-qual with a rather strict semantics
'''
import sys
from itertools import chain, repeat
import xml.etree.ElementTree as etree
import collections


NS = {
    'sbml': 'http://www.sbml.org/sbml/level2/version4',
    'cd': 'http://www.sbml.org/2001/ns/celldesigner',
    'sbml3': 'http://www.sbml.org/sbml/level3/version1/core',
    'layout': 'http://www.sbml.org/sbml/level3/version1/layout/version1',
    'qual': 'http://www.sbml.org/sbml/level3/version1/qual/version1',
    'mathml': 'http://www.w3.org/1998/Math/MathML',
    'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'dc': 'http://purl.org/dc/elements/1.1/',
    'dcterms': 'http://purl.org/dc/terms/',
    'vCard': 'http://www.w3.org/2001/vcard-rdf/3.0#',
    'bqbiol': 'http://biomodels.net/biology-qualifiers/',
    'bqmodel': 'http://biomodels.net/model-qualifiers/',
    'xhtml': 'http://www.w3.org/1999/xhtml'
}


Transition = collections.namedtuple('Transition', [
    'type', 'reactants', 'modifiers', 'notes', 'annotations'])


def read_celldesigner(filename):
    '''main file parsing function'''
    root = etree.parse(filename).getroot()
    tag = root.tag
    if tag != '{' + NS["sbml"] + '}sbml':
        print('Currently limited to SBML Level 2 Version 4')
        exit(1)
    model = root.find('sbml:model', NS)
    display = model.find('./sbml:annotation/cd:extension/cd:modelDisplay', NS)
    return (get_transitions(model, species_info(model)),
            display.get('sizeX'), display.get('sizeY'))


def species_info(model):
    '''create a map from species' ids to their attributes'''
    nameconv = {}
    # Find all CellDesigner species used later
    for species in chain(
            model.findall(
                './sbml:annotation/cd:extension/' +
                'cd:listOfComplexSpeciesAliases/' +
                'cd:complexSpeciesAlias[@compartmentAlias]', NS),
            model.findall(
                './sbml:annotation/cd:extension/' +
                'cd:listOfSpeciesAliases/' +
                'cd:speciesAlias[@compartmentAlias]', NS),
    ):
        bound = species.find('.//cd:bounds', NS)
        ref_species = species.get('species')
        sbml = model.find(
            './sbml:listOfSpecies/sbml:species[@id="' + ref_species + '"]',
            NS)
        annot = sbml.find('./sbml:annotation', NS)
        classtype = get_text(annot.find('.//cd:class', NS), 'PROTEIN')
        if classtype == 'DEGRADED':
            continue
        nameconv[species.get('id')] = {
            'activity': get_text(species.find('.//cd:activity', NS),
                                 'inactive'),
            'x': bound.get('x'),
            'y': bound.get('y'),
            'h': bound.get('h'),
            'w': bound.get('w'),
            'transitions': [],
            'name': sbml.get('name'),
            'type': classtype,
            'modifications': get_mods(annot.find('.//cd:listOfModifications',
                                                 NS)),
            'annotations': annot.find('.//rdf:RDF', NS),
        }
        # also store in nameconv the reverse mapping from SBML species to CD
        # species using the corresponding reference protein
        prot_ref = '__' + sbml.get("name")
        if prot_ref in nameconv:
            nameconv[prot_ref].append(species.get('id'))
        else:
            nameconv[prot_ref] = [species.get('id')]
    add_subcomponents_only(nameconv, model)
    return nameconv


def add_subcomponents_only(nameconv, model):
    '''For unused CD species (only subcomponents of complexes), add their
    annotations to the parent complex'''
    for species in model.findall('./sbml:annotation/cd:extension/' +
                                 'cd:listOfIncludedSpecies/cd:species/' +
                                 'cd:notes/xhtml:html/xhtml:body/' +
                                 'rdf:RDF/../../../..', NS):
        add_rdf(
            nameconv,
            reference=decomplexify(species.get('id'), model, field='species'),
            new_rdf=species.find('.//rdf:RDF', NS))


def add_rdf(nameconv, reference, new_rdf):
    '''Adds the new_rdf element to nameconv[reference]['annotations']'''
    if new_rdf is None:
        return
    if nameconv[reference]['annotations'] is not None:
        nameconv[reference]['annotations'].find(
            './rdf:Description', NS).extend(
                new_rdf.find('./rdf:Description', NS)[:])
    else:
        nameconv[reference]['annotations'] = new_rdf


def get_transitions(model, info):
    '''find all transitions'''
    for trans in model.findall('./sbml:listOfReactions/sbml:reaction', NS):
        annot = trans.find('./sbml:annotation/cd:extension', NS)
        rtype = annot.find('./cd:reactionType', NS).text
        reacs = [decomplexify(reac.get('alias'), model) for reac in
                 annot.findall('./cd:baseReactants/cd:baseReactant', NS)]
        prods = [decomplexify(prod.get('alias'), model) for prod in
                 annot.findall('./cd:baseProducts/cd:baseProduct', NS)]
        mods = [(mod.get('type'),
                 decomplexify(mod.get('aliases'), model)) for mod in
                annot.findall('./cd:listOfModification/cd:modification', NS)]
        notes = trans.find('./sbml:notes//xhtml:body', NS)
        rdf = trans.find('./sbml:annotation/rdf:RDF', NS)
        # remove degraded
        reacs = filter(lambda x: x in info, reacs)
        prods = filter(lambda x: x in info, prods)
        # for each product of a reaction, add this reaction as a transition
        # affecting that species
        for species in prods:
            info[species]['transitions'].append(
                Transition(rtype, reacs, mods, notes, rdf)
            )
    return info


def decomplexify(species, model, field='id'):
    '''return external complex if there is one, or species unchanged
    otherwise'''
    cmplx = model.find('./sbml:annotation/cd:extension/' +
                       'cd:listOfSpeciesAliases/' +
                       'cd:speciesAlias[@{field}="{species}"]'.format(
                           field=field,
                           species=species),
                       NS)
    if cmplx is None:
        return species
    return cmplx.get('complexSpeciesAlias', species)


def get_text(cd_class, default=None):
    '''get the text of an XML field if it exists or return a default'''
    if cd_class is not None:
        return cd_class.text
    return default


def get_mods(cd_modifications):
    '''celldesigner:listOfModifications to list of mods'''
    if not cd_modifications:
        return []
    return [mod.get('state') for mod in
            cd_modifications.findall('cd:modification', NS)]


def write_qual(filename, info, width, height):
    '''write the SBML qual with layout file for our model'''
    for name, space in NS.items():
        etree.register_namespace(name, space)
    root = etree.Element('sbml', {
        'level': '3', 'version': '1', 'layout:required': 'false',
        'xmlns': NS['sbml3'], 'qual:required': 'true',
        'xmlns:layout': NS['layout'], 'xmlns:qual': NS['qual'],
    })
    model = etree.SubElement(root, 'model', id="model_id")
    clist = etree.SubElement(model, 'listOfCompartments')
    etree.SubElement(clist, 'compartment', constant="true", id="comp1")
    llist = etree.SubElement(model, 'layout:listOfLayouts')
    layout = etree.SubElement(llist, 'layout:layout', id="layout1")
    etree.SubElement(layout, 'layout:dimensions', width=width, height=height)
    qlist = etree.SubElement(model, 'qual:listOfQualitativeSpecies')
    add_qual_species(layout, qlist, info)
    tlist = etree.SubElement(model, 'qual:listOfTransitions')
    add_transitions(tlist, info)
    etree.ElementTree(root).write(filename, "UTF-8", xml_declaration=True)


def simplify_model(info):
    '''Cleaning the model w.r.t. some active/inactive species'''
    multispecies = {}
    for key, value in list(info.items()):
        if not key.startswith('csa') and not key.startswith('sa'):
            del info[key]
            if len(value) > 1:
                multispecies[key] = value
                # print('multi', key, value)
    for key, value in multispecies.items():
        for val in value:
            # check that it does not appear in any other reaction than the
            # activation one
            active = None
            for species, data in info.items():
                for trans in data['transitions']:
                    if val in trans.reactants or \
                       val in [mod[0] for mod in trans.modifiers]:
                        if active is None:
                            active = species
                        else:
                            active = False
                if active is False:
                    break
            if not info[val]['transitions'] and active in value:
                add_rdf(info, active, info[val]['annotations'])
                # print('deleting {val} [{active} is active for {key}]'.format(
                #     val=val,
                #     active=active,
                #     key=key,
                # ))
                del info[val]


def add_qual_species(layout, qlist, info):
    '''create layout sub-elements and species'''
    llist = etree.SubElement(layout, 'layout:listOfAdditionalGraphicalObjects')
    for species, data in info.items():
        glyph = etree.SubElement(
            llist, 'layout:generalGlyph',
            {'layout:reference': species, 'layout:id': species + '_glyph'})
        box = etree.SubElement(glyph, 'layout:boundingBox')
        etree.SubElement(
            box, 'layout:position',
            {'layout:x': data['x'], 'layout:y': data['y']})
        etree.SubElement(
            box, 'layout:dimensions',
            {'layout:height': data['h'], 'layout:width': data['w']})
        if data['transitions']:
            constant = "false"
        else:
            constant = "true"
        qspecies = etree.SubElement(
            qlist,
            'qual:qualitativeSpecies',
            {
                'qual:maxLevel': "1",
                'qual:compartment': "comp1",
                'qual:name': fix_name(data['name']),
                'qual:constant': constant,
                'qual:id': species,
            })
        add_annotation(qspecies, data['annotations'])


def fix_name(name):
    '''remove subscripts'''
    return name.replace('_sub_', '').replace('_endsub_', '')


def add_annotation(node, rdf):
    '''add a single RDF element as an annotation node'''
    if rdf is not None:
        etree.SubElement(node, 'annotation').append(rdf)


def add_transitions(tlist, info):
    '''create transition elements'''
    known = list(info.keys())
    for species, data in info.items():
        if data['transitions']:
            trans = etree.SubElement(tlist, 'qual:transition', {
                'qual:id': 'tr_' + species
            })
            ilist = etree.SubElement(trans, 'qual:listOfInputs')
            add_inputs(ilist, data['transitions'], species, known)
            # there might not be any input left after filtering known species
            if not ilist:
                tlist.remove(trans)
            else:
                olist = etree.SubElement(trans, 'qual:listOfOutputs')
                etree.SubElement(olist, 'qual:output', {
                    'qual:qualitativeSpecies': species,
                    'qual:transitionEffect': 'assignmentLevel',
                    'qual:id': 'tr_{species}_out'.format(species=species)
                })
                flist = etree.SubElement(trans, 'qual:listOfFunctionTerms')
                etree.SubElement(flist, 'qual:defaultTerm', {
                    'qual:resultLevel': '0'
                })
                func = etree.SubElement(flist, 'qual:functionTerm', {
                    'qual:resultLevel': '1'
                })
                add_function(func, data['transitions'], known)
                add_notes(trans, data['transitions'])
                add_annotations(trans, data['transitions'])


def add_notes(trans, transitions):
    '''add all the found notes'''
    notes = etree.SubElement(trans, 'notes')
    html = etree.SubElement(notes, 'html', xmlns=NS['xhtml'])
    head = etree.SubElement(html, 'head')
    etree.SubElement(head, 'title')
    body = etree.SubElement(html, 'body')
    some_notes = False
    prefix_len = len(NS['xhtml']) + 2
    for reaction in transitions:
        if reaction.notes is not None:
            some_notes = True
            reaction.notes.tag = 'p'
            for element in reaction.notes.getiterator():
                if element.tag.startswith('{' + NS['xhtml'] + '}'):
                    element.tag = element.tag[prefix_len:]
            body.append(reaction.notes)
    if not some_notes:
        trans.remove(notes)


def add_annotations(trans, transitions):
    '''add all the found annotations'''
    annotation = etree.SubElement(trans, 'annotation')
    rdf = etree.SubElement(annotation, 'rdf:RDF')
    for reaction in transitions:
        if reaction.annotations is not None:
            rdf.append(reaction.annotations[0])
    if not rdf:
        trans.remove(annotation)


def add_function(func, transitions, known):
    '''add the complete boolean activation function

    this is an or over all reactions having the target as product.
    For each reaction it can activate if all reactants are present,
    no inhibitor is present, and one of the activators is present'''
    math = etree.SubElement(func, 'math', xmlns=NS['mathml'])
    # create or node if necessary
    if len(transitions) > 1:
        apply = etree.SubElement(math, 'apply')
        etree.SubElement(apply, 'or')
    else:
        apply = math
    for reaction in transitions:
        # we assume that only "BOOLEAN_LOGIC_GATE_AND" has multiple modifiers
        # it is also the only modification that has an AND and therefore ends
        # with reactants
        reactants = [reac for reac in reaction.reactants if reac in known]
        reactants.extend(
            [mod for (modtype, modifier) in reaction.modifiers
             for mod in modifier.split(',') if
             modtype == 'BOOLEAN_LOGIC_GATE_AND' and mod in known])
        activators = [modifier for (modtype, modifier) in reaction.modifiers
                      if modtype not in ('INHIBITION', 'BOOLEAN_LOGIC_GATE_AND')
                      and modifier in known]
        inhibitors = [modifier for (modtype, modifier) in reaction.modifiers
                      if modtype == 'INHIBITION' and modifier in known]
        # create and node if necessary
        if len(reactants) + len(inhibitors) > 1 or (
                activators and (reactants or inhibitors)):
            lapply = etree.SubElement(apply, 'apply')
            etree.SubElement(lapply, 'and')
        else:
            lapply = apply
        if len(activators) < 2:
            reactants.extend(activators)
        else:
            # create or node if necessary
            inner_apply = etree.SubElement(lapply, 'apply')
            etree.SubElement(inner_apply, 'or')
            for modifier in activators:
                set_level(inner_apply, modifier, '1')
        for level, modifier in chain(zip(repeat('1'), reactants),
                                     zip(repeat('0'), inhibitors)):
            set_level(lapply, modifier, level)


def set_level(elt, modifier, level):
    '''add mathml to element elt such that modifier is equal to level'''
    trigger = etree.SubElement(elt, 'apply')
    etree.SubElement(trigger, 'eq')
    math_ci = etree.SubElement(trigger, 'ci')
    math_ci.text = modifier
    math_cn = etree.SubElement(trigger, 'cn', type='integer')
    math_cn.text = level


def add_inputs(ilist, transitions, species, known):
    '''add all known inputs'''
    index = 0
    modifiers = []
    for reaction in transitions:
        # we use enumerate to get a dummy modtype for reactants
        for modtype, modifier in chain(enumerate(reaction.reactants),
                                       reaction.modifiers):
            if modtype == 'INHIBITION':
                sign = 'negative'
            else:
                sign = 'positive'
            if (modifier, sign) not in modifiers and modifier in known:
                modifiers.append((modifier, sign))
                etree.SubElement(ilist, 'qual:input', {
                    'qual:qualitativeSpecies': modifier,
                    'qual:transitionEffect': 'none',
                    'qual:sign': sign,
                    'qual:id': 'tr_{species}_in_{index}'.format(
                        species=species,
                        index=index,
                    ),
                })
                index += 1


def main():
    '''run conversion using the CLI given first argument'''
    if len(sys.argv) != 3:
        print('Usage: [python3] ' + sys.argv[0] +
              ' <celldesignerinfile.xml> <sbmlqualoutfile.xml>')
        exit(1)
    celldesignerfile = sys.argv[1]
    print('parsing {celldesignerfile}…'.format(
        celldesignerfile=celldesignerfile))
    info, width, height = read_celldesigner(celldesignerfile)
    simplify_model(info)
    write_qual(sys.argv[2], info, width, height)


if __name__ == '__main__':  # pragma: no cover
    main()
