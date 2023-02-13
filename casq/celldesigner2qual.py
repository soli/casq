"""Convert CellDesigner models to SBML-qual with a rather strict semantics.

Copyright (C) 2019, Sylvain.Soliman@inria.fr

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

import argparse
import os.path
import sys
from typing import List

from loguru import logger  # type: ignore

from . import bmaExport, version
from .readCD import read_celldesigner
from .simplify import simplify_model
from .write import write_csv, write_qual


def map_to_model(map_filename: str, model_filename: str, bma=False):
    """Do the full run with defaults arguments."""
    logger.disable("casq")
    with open(map_filename, "r", encoding="utf-8") as f:
        info, width, height = read_celldesigner(f)
    simplify_model(info, [], [])
    if not bma:
        write_qual(model_filename, info, width, height)
    else:
        bmaExport.write_bma(model_filename, info, 1, None, False, True)


def main(argv: List[str] = None):
    """Run conversion using the CLI given first argument."""
    parser = argparse.ArgumentParser(
        description=" ".join(__doc__.splitlines()[:3]) + " GPLv3"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s v{version}",
    )
    parser.add_argument(
        "-D", "--debug", action="store_true", help="Display a lot of debug information"
    )
    parser.add_argument(
        "-c",
        "--csv",
        action="store_true",
        help="Store the species information in a separate CSV (and .bnet) file",
    )
    parser.add_argument(
        "-s",
        "--sif",
        action="store_true",
        help="Store the influence information in a separate SIF file",
    )
    parser.add_argument(
        "-r",
        "--remove",
        metavar="S",
        type=int,
        default="0",
        help="""Delete connected components in the resulting model if their
        size is smaller than S.
        A negative S leads to keep only the biggest(s) connected component(s)""",
    )
    parser.add_argument(
        "-f",
        "--fixed",
        type=argparse.FileType(),
        help="""A CSV file containing input values or knock-ins/knock-outs,
        one per line, with name in the first column and the value in the second.""",
    )
    parser.add_argument(
        "-n",
        "--names",
        action="store_true",
        help="Use the names as IDs in the SBML file",
    )
    if sys.version_info >= (3, 8, 0):
        parser.add_argument(
            "-u",
            "--upstream",
            action="extend",
            nargs="*",
            type=str,
            default=[],
            help="Only species upstream of this/these species will be kept",
        )
        parser.add_argument(
            "-d",
            "--downstream",
            action="extend",
            nargs="*",
            type=str,
            default=[],
            help="Only species downstream of this/these species will be kept",
        )
    else:
        parser.add_argument(
            "-u",
            "--upstream",
            action="append",
            type=str,
            default=[],
            help="Only species upstream of this/these species will be kept",
        )
        parser.add_argument(
            "-d",
            "--downstream",
            action="append",
            type=str,
            default=[],
            help="Only species downstream of this/these species will be kept",
        )
    parser.add_argument(
        "infile",
        type=argparse.FileType("r", encoding="utf-8"),
        nargs="?",
        default=sys.stdin,
        help="CellDesigner File",
    )
    parser.add_argument(
        "-b",
        "--bma",
        action="store_true",
        help="Output to BMA json format",
    )
    parser.add_argument(
        "-g",
        "--granularity",
        type=int,
        default="1",
        help="When exporting to BMA, use this granularity",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=int,
        default=None,
        help="When exporting to BMA, nodes with no input should be set to this value",
    )
    parser.add_argument(
        "-C",
        "--colourConstant",
        action="store_false",
        help="When exporting to BMA, colour all variables pink (defaults to colour by compartment)",
    )
    parser.add_argument(
        "outfile", nargs="?", default=sys.stdout, help="SBML-Qual/BMA json File"
    )
    if argv:
        args = parser.parse_args(argv)
    else:
        args = parser.parse_args()

    if not args.debug:
        logger.disable("casq")
    logger.debug("parsing {fname}…", fname=args.infile.name)
    info, width, height = read_celldesigner(args.infile)
    simplify_model(info, args.upstream, args.downstream, args.names)
    if args.infile != sys.stdin and args.outfile == sys.stdout:
        args.outfile = os.path.splitext(args.infile.name)[0] + ".sbml"
    if args.bma:
        bmaExport.write_bma(
            args.outfile, info, args.granularity, args.input, False, args.colourConstant
        )
    else:
        write_qual(
            args.outfile,
            info,
            width,
            height,
            remove=args.remove,
            sif=args.sif,
            fixed=args.fixed,
        )
    if args.csv and args.outfile != sys.stdout:
        write_csv(args.outfile, info)


if __name__ == "__main__":  # pragma: no cover
    main()
