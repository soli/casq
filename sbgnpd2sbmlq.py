#!/usr/bin/env python3
'''
Convert SBGN-PD models in SBGN-ML to SBML-qual with a rather strict semantics
'''
import sys
import xml.etree.ElementTree as etree


def read_sbgn(filename):
    '''main file parsing function'''
    tree = etree.parse(filename)
    tag = tree.getroot().tag
    nspace = tag[1:tag.index('}')]
    return tree, nspace


def main():
    '''run conversion using the CLI given first argument'''
    sbgnfile = sys.argv[0]
    _, nspace = read_sbgn(sbgnfile)
    print(nspace)


if __name__ == '__main__':  # pragma: no cover
    main()
