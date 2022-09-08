"""Tests for CaSQ."""

from filecmp import cmp
from glob import glob
from os import path
from unittest.mock import patch

from casq.celldesigner2qual import main, map_to_model
from casq.utils import validate

import pytest  # type: ignore


@pytest.mark.parametrize(
    "infile", glob(path.join(str(path.dirname(path.realpath(__file__))), "*.xml")),
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
