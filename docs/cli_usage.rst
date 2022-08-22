Command-line usage
==================

Just follow the instructions:

.. command-output:: casq --help

.. ifconfig:: not programoutput_available

   .. code-block::

      $ casq --help
      usage: casq [-h] [-v] [-D] [-c] [-s] [-r S] [-u [UPSTREAM ...]]
                  [-d [DOWNSTREAM ...]] [-b] [-g GRANULARITY] [-i INPUT] [-C]
                  [infile] [outfile]

      Convert CellDesigner models to SBML-qual with a rather strict semantics.
      Copyright (C) 2019-2022 Sylvain.Soliman@inria.fr GPLv3

      positional arguments:
        infile                CellDesigner File
        outfile               SBML-Qual/BMA json File

      optional arguments:
        -h, --help            show this help message and exit
        -v, --version         show program's version number and exit
        -D, --debug           Display a lot of debug information
        -c, --csv             Store the species information in a separate CSV file
        -s, --sif             Store the influence information in a separate SIF file
        -r S, --remove S      Delete connected components in the resulting model if
                              their size is smaller than S. A negative S leads to
                              keep only the biggest(s) connected component(s)
        -u [UPSTREAM ...], --upstream [UPSTREAM ...]
                              Only species upstream of this/these species will be
                              kept
        -d [DOWNSTREAM ...], --downstream [DOWNSTREAM ...]
                              Only species downstream of this/these species will be
                              kept
        -b, --bma             Output to BMA json format
        -g GRANULARITY, --granularity GRANULARITY
                              When exporting to BMA, use this granularity
        -i INPUT, --input INPUT
                              When exporting to BMA, nodes with no input should be
                              set to this value
        -C, --colourConstant  When exporting to BMA, colour all variables pink
                              (defaults to colour by compartment)
