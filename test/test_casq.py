"""Tests for CaSQ."""

from difflib import ndiff
from filecmp import cmp
from glob import glob
from os import chdir, path, unlink
from unittest.mock import patch

import pytest  # type: ignore

from casq.celldesigner2qual import main, map_to_model
from casq.utils import validate


@pytest.mark.parametrize(
    "infile",
    glob(path.join(str(path.dirname(path.realpath(__file__))), "map_*.xml")),
)
def test_casq_produces_valid_files(tmp_path, infile):
    """Check if the files we produce are valid."""
    outfile = path.join(
        str(tmp_path), path.splitext(path.basename(infile))[0] + ".sbml"
    )
    with patch("sys.argv", ["casq", infile, outfile]):
        main()
    assert validate(outfile) == "OK"

    map_to_model(infile, outfile + "_api")

    assert cmp(outfile, outfile + "_api")


def test_CD_and_SBGNML_similar():
    """Check if the files generated from SBGN and CD are similar."""
    chdir(path.dirname(path.realpath(__file__)))
    with patch("sys.argv", ["casq", "-c", "Tgfb_CD.xml"]):
        main()
    with patch("sys.argv", ["casq", "-c", "Tgfb_sbgnml.xml"]):
        main()

    with open("Tgfb_CD.bnet") as f:
        cdsbml_out = (
            f.read()
            .replace("_empty", "")
            .splitlines(keepends=True)
        )
    with open("Tgfb_sbgnml.bnet") as f:
        sbgnml_out = f.readlines()

    unlink("Tgfb_CD.bnet")
    unlink("Tgfb_CD.csv")
    unlink("Tgfb_CD.sbml")
    unlink("Tgfb_sbgnml.bnet")
    unlink("Tgfb_sbgnml.csv")
    unlink("Tgfb_sbgnml.sbml")

    # import sys
    # sys.stderr.writelines(ndiff(cdsbml_out, sbgnml_out))

    # differences are order of SMAD double-complex (4 lines), missing ubiquitin inhibition from CD file (3 lines)
    assert len([s for s in ndiff(cdsbml_out, sbgnml_out) if not s.startswith(" ")]) == 7
