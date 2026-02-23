"""Tests for CaSQ."""

from difflib import ndiff
from filecmp import cmp
from glob import glob
from os import path
from unittest.mock import patch

import pytest  # type: ignore

from casq.celldesigner2qual import main, map_to_model
from casq.utils import validate, validator_unavailable


@pytest.fixture
def change_test_dir(request, monkeypatch):
    """Change dir to the test directory for the duration of the test."""
    monkeypatch.chdir(request.fspath.dirname)


@pytest.mark.skipif(validator_unavailable(), reason="validator unavailable")
@pytest.mark.parametrize(
    "infile",
    glob(path.join(str(path.dirname(path.realpath(__file__))), "map_*.xml")),
)
def test_casq_produces_valid_files_on_maps(tmp_path, infile):
    """Check if the files we produce are valid."""
    outfile = path.join(
        str(tmp_path), path.splitext(path.basename(infile))[0] + ".sbml"
    )
    with patch("sys.argv", ["casq", infile, outfile]):
        main()
    assert validate(outfile) == "OK"

    map_to_model(infile, outfile + "_api")

    assert cmp(outfile, outfile + "_api")


@pytest.mark.skipif(validator_unavailable(), reason="validator unavailable")
@pytest.mark.parametrize(
    "infile",
    glob(path.join(str(path.dirname(path.realpath(__file__))), "R-HSA*.sbgn")),
)
def test_casq_produces_valid_files_on_reactome(tmp_path, infile):
    """Check if the files we produce are valid."""
    outfile = path.join(
        str(tmp_path), path.splitext(path.basename(infile))[0] + ".sbml"
    )
    with patch("sys.argv", ["casq", infile, outfile]):
        main()
    assert validate(outfile) == "OK"


@pytest.mark.parametrize(
    "infile,diffs",
    [
        ("Tgfb", 7),  # order of SMAD cmplx, missing ubiquitin inhibition from CD
        ("E_Prot", 4),  # order for Activity_space_pheno
        ("JNK", 0),
        ("Orf3a", 4),  # order for P65 cmplx and CASP1, two distinct MYD88 in CD
        ("coagulation", 62), # differences due to complex names (_space → _vascular/_host) and order of complexes
        ("Nsp4_Nsp6_inter", 0)
    ],
)
def test_CD_and_SBGNML_similar(infile, diffs, change_test_dir):
    """Check if the files generated from SBGN and CD are similar."""
    cdsbml_in = infile + "_CD.xml"
    sbgnml_in = infile + "_SBGNML.xml"
    with patch("sys.argv", ["casq", "-c", cdsbml_in]):
        main()
    with patch("sys.argv", ["casq", "-c", sbgnml_in]):
        main()

    with open(infile + "_CD.bnet") as f:
        cdsbml_out = (
            f.read()
            .replace("_empty", "")
            .replace("_ion", "")
            .replace("_simple_molecule", "")
            .replace("_rna", "_nucleic_acid_feature")
            .replace(",_space", "")
            .replace("mitochondrial_space_matrix", "mitochondrial_matrix")
            .replace("endoplasmic_space_reticulum", "endoplasmic_reticulum")
            .replace("human_space_host", "human_host")
            .splitlines(keepends=True)
        )
    with open(infile + "_SBGNML.bnet") as f:
        sbgnml_out = (
            f.read()
            .replace("_macromolecule_multimer", "")
            .replace(",_", "_")
            .replace("_unspecified_entity", "_drug")
            .splitlines(keepends=True)
        )

    import sys

    sys.stderr.writelines(ndiff(cdsbml_out, sbgnml_out))

    # differences are order of SMAD double-complex (4 lines), missing ubiquitin inhibition from CD file (3 lines)
    assert (
        len([s for s in ndiff(cdsbml_out, sbgnml_out) if not s.startswith(" ")])
        == diffs
    )


def test_spreadsheet_output():
    """Test our CSV export."""
    filename = path.join(str(path.dirname(path.realpath(__file__))), "Small_Apoptosis")
    with patch("sys.argv", ["casq", "-c", filename + ".xml"]):
        main()

    assert path.isfile(filename + "_Model.csv"), "Model sheet CSV file not generated"
    with open(filename + ".bnet", "r") as f:
        # ignore header
        bnet_lines = f.readlines()[3:]
    with open(filename + "_Species.csv", "r") as f:
        # ignore header
        species_lines = f.readlines()[1:]
    assert len(species_lines) == len(bnet_lines), (
        "Species sheet CSV file and BNET files differ in length"
    )
    non_constant = [line for line in species_lines if ",False," in line]
    with open(filename + "_Transitions.csv", "r") as f:
        # ignore header
        transitions_lines = f.readlines()[1:]
    assert len(transitions_lines) == len(non_constant), (
        "Transitions sheet CSV file and non-constant entries in BNET differ"
    )
    for line in transitions_lines:
        things = line.split(",")
        assert f"{things[0]}, {things[2]}\n" in bnet_lines
