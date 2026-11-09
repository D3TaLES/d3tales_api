import numpy as np
import pytest
from pymatgen.core.periodic_table import Element
from rdkit import Chem

from d3tales_api.Calculators.ocelot_transform import (
    ACParser,
    chiral_stereo_check,
    clean_charges,
    valence_electron,
)


@pytest.mark.unit
@pytest.mark.parametrize("symbol, expected", [
    ("H", 1), ("C", 4), ("N", 5), ("O", 6), ("F", 7),
])
def test_valence_electron_main_group(symbol, expected):
    assert valence_electron(Element(symbol)) == expected


@pytest.mark.unit
def test_chiral_stereo_check_returns_sanitized_molecule():
    mol = Chem.MolFromSmiles("C[C@H](O)F")
    result = chiral_stereo_check(mol)
    assert result is mol
    assert Chem.MolToSmiles(result, isomericSmiles=True)


@pytest.mark.unit
def test_clean_charges_preserves_neutral_molecule():
    mol = Chem.MolFromSmiles("c1ccccc1")
    result = clean_charges(mol)
    assert Chem.MolToSmiles(result) == "c1ccccc1"


@pytest.mark.unit
def test_get_uadu_identifies_unsaturated_atoms():
    ua, du = ACParser.getUADU([4, 2, 1], [3, 2, 0])
    assert ua == [0, 2]
    assert du == [1, 1]


@pytest.mark.unit
def test_get_bonds_returns_only_connected_unsaturated_pairs():
    ac = np.array([
        [0, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [0, 0, 1, 0],
    ])
    assert ACParser.get_bonds([0, 1, 2, 3], ac) == [(0, 1), (1, 2), (2, 3)]


@pytest.mark.unit
def test_get_ua_pairs_returns_empty_tuple_when_no_bonds():
    ac = np.zeros((2, 2), dtype=int)
    assert ACParser.get_UA_pairs([0, 1], ac) == [()]


@pytest.mark.unit
def test_get_ua_pairs_returns_non_overlapping_matching():
    ac = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ])
    pairs = ACParser.get_UA_pairs([0, 1, 2, 3], ac)[0]
    assert {frozenset(pair) for pair in pairs} == {frozenset((0, 1)), frozenset((2, 3))}


@pytest.mark.unit
def test_get_bo_increases_selected_bond_order():
    ac = np.array([[0, 1], [1, 0]], dtype=int)
    bo, final_pairs = ACParser.get_BO(ac, [1, 1], [2, 2], [(0, 1)])
    assert bo.tolist() == [[0, 2], [2, 0]]
    assert final_pairs == ()


@pytest.mark.unit
@pytest.mark.parametrize("atomic_number, valence_electrons, bond_order, expected", [
    (1, 1, 1, 0),
    (5, 3, 3, 0),
    (6, 4, 4, 0),
    (7, 5, 3, 0),
    (8, 6, 2, 0),
    (15, 5, 5, 0),
    (16, 6, 6, 0),
])
def test_get_atomic_charge_neutral_valences(
    atomic_number, valence_electrons, bond_order, expected
):
    assert ACParser.get_atomic_charge(
        atomic_number, valence_electrons, bond_order
    ) == expected


@pytest.mark.unit
def test_valences_not_too_large_accepts_valid_matrix():
    bo = np.array([[0, 1], [1, 0]], dtype=int)
    assert ACParser.valences_not_too_large(bo, [1, 1]) is True


@pytest.mark.unit
def test_valences_not_too_large_rejects_excess_bonds():
    bo = np.array([[0, 2], [2, 0]], dtype=int)
    assert ACParser.valences_not_too_large(bo, [1, 1]) is False


@pytest.mark.unit
def test_parser_initializes_rdkit_molecule_with_expected_atoms():
    parser = ACParser(
        np.array([[0, 1], [1, 0]], dtype=int),
        charge=0,
        atomnumberlist=[6, 8],
        sani=False,
    )
    mol = parser.init_rdmol()
    assert mol.GetNumAtoms() == 2
    assert [atom.GetAtomicNum() for atom in mol.GetAtoms()] == [6, 8]


@pytest.mark.unit
def test_add_bo_to_molecule_creates_single_bond():
    parser = ACParser(
        np.array([[0, 1], [1, 0]], dtype=int),
        charge=0,
        atomnumberlist=[6, 8],
        sani=False,
    )
    mol = parser.init_rdmol()
    result = parser.addBO2mol(
        mol,
        np.array([[0, 1], [1, 0]], dtype=int),
        charged_fragments=False,
        force_single=True,
    )
    assert result.GetNumBonds() == 1
    assert result.GetBondWithIdx(0).GetBondType() == Chem.BondType.SINGLE
