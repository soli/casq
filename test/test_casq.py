"""Tests for CaSQ."""

import os
import sys

from casq.celldesigner2qual import main
from casq.utils import validate

import pytest  # type: ignore


@pytest.mark.parametrize("infile", ["map_mastcell.xml"])
def test_casq_produces_valid_files(tmp_path, infile):
    """Check if the files we produce are valid."""
    outfile = os.path.join(tmp_path, os.path.splitext(infile)[0] + ".sbml")
    infile = os.path.join(os.path.dirname(os.path.realpath(__file__)), infile)
    store_argv = sys.argv
    sys.argv = ["casq", infile, outfile]
    try:
        main()
    finally:
        sys.argv = store_argv
    assert validate(outfile) == "OK"
