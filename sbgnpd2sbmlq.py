#!/usr/bin/env python3
'''
Convert SBGN-PD models in SBGN-ML to SBML-qual with a rather strict semantics
'''
import sys
import xml.etree.ElementTree as etree


SBGNML = '{http://sbgn.org/libsbgn/0.2}'


def read_sbgn(filename):
    '''main file parsing function'''
    root = etree.parse(filename).getroot()
    tag = root.tag
    if tag != f'{SBGNML}sbgn':
        print('Currently limited to SBGN-ML 0.2')
        exit(1)
    map = root[0]
    if map.get('language') != 'process description':
        print('Currently limited to SBGN PD maps')
        exit(1)
    children = list(map)
    for child in children:
        read_child(child)


def read_child(child):
    '''parses an SBGN-PD glyph or arc'''
    if child.tag == f'{SBGNML}glyph':
        read_glyph(child)
    elif child.tag == f'{SBGNML}arc':
        read_arc(child)
    else:
        print(f'Unknown {child.tag} element in map')
        exit(1)


def read_glyph(glyph):
    cls = glyph.get('class')
    if cls == 'compartment':
        pass
    elif cls == 'source and sink':
        pass
    elif cls == 'process':
        pass # store ports and id
    elif cls in ('macromolecule', 'nucleic acid feature',
                   'simple chemical', 'phenotype', 'complex'):
        print(glyph.find(f'{SBGNML}label').get('text'))
    else:
        print(f'Unknown cls: {cls}')

    if glyph.find(f'{SBGNML}glyph'):
        print('HAS A CHILD ==================')


def read_arc(arc):
    pass # print(arc.get('class'))


def main():
    '''run conversion using the CLI given first argument'''
    if len(sys.argv) != 2:
        print(f'Usage: [python3] {sys.argv[0]} <filename.sbgn>')
        exit(1)
    sbgnfile = sys.argv[1]
    print(f'parsing {sbgnfile}…')
    read_sbgn(sbgnfile)


if __name__ == '__main__':  # pragma: no cover
    main()
