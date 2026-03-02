"""Tests for CaSQ."""

from csv import reader
from difflib import ndiff
from filecmp import cmp
from glob import glob
from os import path
from unittest.mock import patch

import pytest  # type: ignore
from tabularqual import convert_sbml_to_spreadsheet  # ty: ignore[unresolved-import]

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
        ("Tgfb", 4),  # order of SMAD cmplx
        ("E_Prot", 4),  # order for Activity_space_pheno
        ("JNK", 0),
        ("Orf3a", 4),  # order for P65 cmplx and CASP1, TLR3 not receptor in SBGN
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
            .splitlines(keepends=True)
        )
    with open(infile + "_SBGNML.bnet") as f:
        sbgnml_out = (
            f.read()
            .replace("_macromolecule_multimer", "")
            .replace("_Active", "")
            .replace("_active", "")
            .replace("_ion", "")
            .splitlines(keepends=True)
        )

    import sys

    sys.stderr.writelines(ndiff(cdsbml_out, sbgnml_out))

    # differences are order of SMAD double-complex (4 lines), missing ubiquitin inhibition from CD file (3 lines)
    assert (
        len([s for s in ndiff(cdsbml_out, sbgnml_out) if not s.startswith(" ")])
        == diffs
    )


def test_csv_vs_bnet_output():
    """Test our CSV export vs. our .bnet output."""
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
        things = line.replace('"', "").split(",")
        assert f"{things[0]}, {things[1]}\n" in bnet_lines


def test_csv_vs_tabularqual_output():
    """Test our CSV export vs. TabularQual's one."""
    filename = path.join(str(path.dirname(path.realpath(__file__))), "Small_Apoptosis")
    with patch("sys.argv", ["casq", "-c", filename + ".xml"]):
        main()

    convert_sbml_to_spreadsheet(
        filename + ".sbml",
        filename + "_tq",
        output_csv=True,
        use_name=True,
        print_messages=False,
    )

    with open(filename + "_Species.csv", newline="") as f:
        species_lines = list(reader(f))

    with open(filename + "_tq_Species.csv", newline="") as f:
        tq_species_lines = [
            [
                cell.replace(", ", ",").replace("&amp;", "&").replace("&quot;", '"')
                for (i, cell) in enumerate(row)
                if i in (1, 3, 8, 9, 10)
            ]
            for row in reader(f)
        ]

    for i, line in enumerate(species_lines):
        if line[0] in ("M", "Orf3a", "E", "Orf3b", "Orf8a", "N", "Orf6", "Orf9b"):
            # ignore is and isEncodedBy non supported by CaSQ
            continue
        # ignore isDescribedBy Notes1 and Comment
        assert [line[0]] + line[2:-2] == tq_species_lines[i]

    with open(filename + "_Transitions.csv", newline="") as f:
        trans_lines = list(reader(f))
    trans_lines.sort()

    with open(filename + "_tq_Transitions.csv", newline="") as f:
        tq_trans_lines = [
            [cell.replace(" ", "") for (i, cell) in enumerate(row) if i in (2, 4, 6)]
            for row in reader(f)
        ]
    tq_trans_lines.sort()

    for i, line in enumerate(trans_lines):
        if line[0] in ('"BCL2/MCL1/BCL2L1_complex"', '"BAD/BBC3/BCL2L11_complex"'):
            # ignore isEncodedBy non supported by CaSQ
            continue
        # ignore Notes1
        assert line[:-3] + [line[-2]] == tq_trans_lines[i]
