"""Tests for CaSQ."""

from difflib import ndiff
from filecmp import cmp
from glob import glob
from os import path
from unittest.mock import patch

import pytest  # type: ignore

from casq.celldesigner2qual import main, map_to_model
from casq.utils import validate


@pytest.fixture
def change_test_dir(request, monkeypatch):
    """Change dir to the test directory for the duration of the test."""
    monkeypatch.chdir(request.fspath.dirname)


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


@pytest.mark.parametrize(
    "cdsbml_in",
    glob(path.join(str(path.dirname(path.realpath(__file__))), "*_CD.xml")),
)
def test_CD_and_SBGNML_similar(cdsbml_in, change_test_dir):
    """Check if the files generated from SBGN and CD are similar."""
    sbgnml_in = cdsbml_in[:-6] + "SBGNML.xml"
    with patch("sys.argv", ["casq", "-c", cdsbml_in]):
        main()
    with patch("sys.argv", ["casq", "-c", sbgnml_in]):
        main()

    with open(cdsbml_in[:-3] + "bnet") as f:
        cdsbml_out = (
            f.read()
            .replace("_empty", "")
            .replace("_phosphorylated", "")
            .replace("_ubiquitinated", "")
            .replace("_ion", "")
            .replace("_simple_molecule", "")
            .replace("_rna", "_nucleic_acid_feature")
            .splitlines(keepends=True)
        )
    with open(sbgnml_in[:-3] + "bnet") as f:
        sbgnml_out = f.readlines()

    import sys

    sys.stderr.writelines(ndiff(cdsbml_out, sbgnml_out))

    # differences are order of SMAD double-complex (4 lines), missing ubiquitin inhibition from CD file (3 lines)
    assert len([s for s in ndiff(cdsbml_out, sbgnml_out) if not s.startswith(" ")]) == 7
