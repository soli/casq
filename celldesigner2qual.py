#!/usr/bin/env python3
'''
Convert CellDesigner models to SBML-qual with a rather strict semantics
'''
import sys
import xml.etree.ElementTree as etree


NS = {
    'sbml': 'http://www.sbml.org/sbml/level2/version4',
    'cd': 'http://www.sbml.org/2001/ns/celldesigner'
}


def read_sbgn(filename):
    '''main file parsing function'''
    root = etree.parse(filename).getroot()
    tag = root.tag
    if tag != f'{{{NS["sbml"]}}}sbml':
        print('Currently limited to SBML Level 2 Version 4')
        exit(1)
    model = root.find('sbml:model', NS)
    return species_conv(model)


def species_conv(model):
    '''create a map from species' ids to their attributes'''
    nameconv = {}
    for species in model.findall('./sbml:listOfSpecies/sbml:species', NS):
        annot = species.find('sbml:annotation', NS)
        cls = get_class(annot.find('.//cd:class', NS))
        mods = get_mods(annot.find('.//cd:listOfModifications', NS))
        nameconv[species.get('id')] = [
            species.get('name'),
            cls,
            mods]
    return nameconv


def get_class(cd_class):
    '''celldesigner:class to class'''
    if cd_class is not None:
        return cd_class.text
    return 'PROTEIN'


def get_mods(cd_modifications):
    '''celldesigner:listOfModifications to list of mods'''
    mods = []
    if cd_modifications:
        for mod in cd_modifications.findall('cd:modification', NS):
            mods.append(mod.get('state'))
    return mods


def main():
    '''run conversion using the CLI given first argument'''
    if len(sys.argv) != 2:
        print(f'Usage: [python3] {sys.argv[0]} <filename.sbgn>')
        exit(1)
    sbgnfile = sys.argv[1]
    print(f'parsing {sbgnfile}…')
    # read_sbgn(sbgnfile)
    print(read_sbgn(sbgnfile))


if __name__ == '__main__':  # pragma: no cover
    main()
