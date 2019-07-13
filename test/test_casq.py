"""Tests for CaSQ."""

import glob
import os
import unittest.mock

from casq.celldesigner2qual import main
from casq.utils import validate

import pytest  # type: ignore


@pytest.mark.parametrize(
    "infile",
    glob.glob(os.path.join(str(os.path.dirname(os.path.realpath(__file__))), "*.xml")),
)
def test_casq_produces_valid_files(tmp_path, infile):
    """Check if the files we produce are valid."""
    outfile = os.path.join(
        str(tmp_path), os.path.splitext(os.path.basename(infile))[0] + ".sbml"
    )
    with unittest.mock.patch("sys.argv", ["casq", infile, outfile]):
        main()
    assert validate(outfile) == "OK"
