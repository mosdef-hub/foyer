import difflib
import glob
import os
from pkg_resources import resource_filename

import mbuild as mb
from mbuild.examples import Alkane
import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.forcefield import generate_topology
from foyer.forcefield import _check_independent_residues
from foyer.tests.utils import get_fn


FF_DIR = resource_filename('foyer', 'forcefields')
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, '*.xml'))


def test_load_files():
    for ff_file in FORCEFIELDS:
        ff1 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) > 0

        ff2 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) == len(ff2._atomTypes)


def test_duplicate_type_definitions():
    with pytest.raises(ValueError):
        ff4 = Forcefield(name='oplsaa', forcefield_files=FORCEFIELDS)



def test_from_parmed():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2)

    assert sum((1 for at in ethane.atoms if at.type == 'opls_135')) == 2
    assert sum((1 for at in ethane.atoms if at.type == 'opls_140')) == 6
    assert len(ethane.bonds) == 7
    assert all(x.type for x in ethane.bonds)
    assert len(ethane.angles) == 12
    assert all(x.type for x in ethane.angles)
    assert len(ethane.rb_torsions) == 9
    assert all(x.type for x in ethane.dihedrals)

    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    mol2.box_vectors = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2)

    assert ethane.box_vectors == mol2.box_vectors


def test_from_mbuild():
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2)

    assert sum((1 for at in ethane.atoms if at.type == 'opls_135')) == 2
    assert sum((1 for at in ethane.atoms if at.type == 'opls_140')) == 6
    assert len(ethane.bonds) == 7
    assert all(x.type for x in ethane.bonds)
    assert len(ethane.angles) == 12
    assert all(x.type for x in ethane.angles)
    assert len(ethane.rb_torsions) == 9
    assert all(x.type for x in ethane.dihedrals)

def test_write_refs():
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2, references_file='ethane.bib')
    assert os.path.isfile('ethane.bib')

def test_write_refs_multiple():
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(forcefield_files=get_fn('refs-multi.xml'))
    ethane = oplsaa.apply(mol2, references_file='ethane-multi.bib')
    assert os.path.isfile('ethane-multi.bib')
    with open(get_fn('ethane-multi.bib')) as file1:
        with open('ethane-multi.bib') as file2:
            diff = list(difflib.unified_diff(file1.readlines(),
                                             file2.readlines(),
                                             n=0))
    assert not diff

def test_preserve_resname():
    untyped_ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    untyped_resname = untyped_ethane.residues[0].name
    oplsaa = Forcefield(name='oplsaa')
    typed_ethane = oplsaa.apply(untyped_ethane)
    typed_resname = typed_ethane.residues[0].name
    assert typed_resname == untyped_resname

def test_apply_residues():
    from mbuild.examples import Ethane
    ethane = Ethane()
    opls = Forcefield(name='oplsaa')
    typed = opls.apply(ethane, residues='CH3')
    assert len([res for res in typed.residues if res.name == 'CH3']) == 2

def test_from_mbuild_customtype():
    mol2 = mb.load(get_fn('ethane_customtype.pdb'))
    customtype_ff = Forcefield(forcefield_files=get_fn('validate_customtypes.xml'))
    ethane = customtype_ff.apply(mol2)

    assert sum((1 for at in ethane.atoms if at.type == 'C3')) == 2
    assert sum((1 for at in ethane.atoms if at.type == 'Hb')) == 6
    assert len(ethane.bonds) == 7
    assert all(x.type for x in ethane.bonds)
    assert len(ethane.angles) == 12
    assert all(x.type for x in ethane.angles)
    assert len(ethane.rb_torsions) == 9
    assert all(x.type for x in ethane.dihedrals)

def test_improper_dihedral():
    untyped_benzene = pmd.load_file(get_fn('benzene.mol2'), structure=True)
    ff_improper = Forcefield(forcefield_files=get_fn('improper_dihedral.xml'))
    benzene = ff_improper.apply(untyped_benzene, assert_dihedral_params=False)
    assert len(benzene.dihedrals) == 18
    assert len([dih for dih in benzene.dihedrals if dih.improper]) == 6
    assert len([dih for dih in benzene.dihedrals if not dih.improper]) == 12

def test_residue_map():
    ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    ethane *= 2
    oplsaa = Forcefield(name='oplsaa')
    topo, NULL = generate_topology(ethane)
    topo_with = oplsaa.run_atomtyping(topo, use_residue_map=True)
    topo_without = oplsaa.run_atomtyping(topo, use_residue_map=False)
    assert all([a.id for a in topo_with.atoms()][0])
    assert all([a.id for a in topo_without.atoms()][0])
    struct_with = pmd.openmm.load_topology(topo_with, oplsaa.createSystem(topo_with))
    struct_without = pmd.openmm.load_topology(topo_without, oplsaa.createSystem(topo_without))
    for atom_with, atom_without in zip(struct_with.atoms, struct_without.atoms):
        assert atom_with.type == atom_without.type
        b_with = atom_with.bond_partners
        b_without = atom_without.bond_partners
        assert [a0.type for a0 in b_with] == [a1.type for a1 in b_without]
        assert [a0.idx for a0 in b_with] == [a1.idx for a1 in b_without]

def test_independent_residues_molecules():
    """Test to see that _check_independent_residues works for molecules."""
    butane = Alkane(4)
    structure = butane.to_parmed()
    topo, NULL = generate_topology(structure)
    assert _check_independent_residues(topo)
    structure = butane.to_parmed(residues=['RES', 'CH3'])
    topo, NULL = generate_topology(structure)
    assert not _check_independent_residues(topo)

def test_independent_residues_atoms():
    """Test to see that _check_independent_residues works for single aotms."""
    argon = mb.Compound()
    argon.name = 'Ar'
    structure = argon.to_parmed()
    topo, NULL = generate_topology(structure)
    assert _check_independent_residues(topo)

def test_topology_precedence():
    """Test to see if topology precedence is properly adhered to.

    This test uses a force field file where bond, angle, and dihedral
    parameters are present with different counts of `type` definitions.
    It checks that:
        1. The parameters with the higher number of `type` definitions
           are assigned (because they are given the highest precedence)
        2. That if multiple definitions exist with the same number of
           `type` definitions, that the convention from OpenMM is followed
           whereby the definitions that occurs earliest in the XML is
           assigned.
    """
    ethane = mb.load(get_fn('ethane.mol2'))
    ff = Forcefield(forcefield_files=get_fn('ethane-topo-precedence.xml'))
    typed_ethane = ff.apply(ethane)

    assert len([bond for bond in typed_ethane.bonds
                if round(bond.type.req, 2) == 1.15]) == 6
    assert len([bond for bond in typed_ethane.bonds
                if round(bond.type.req, 2) == 1.6]) == 1
    assert len([angle for angle in typed_ethane.angles
                if round(angle.type.theteq, 3) == 120.321]) == 6
    assert len([angle for angle in typed_ethane.angles
                if round(angle.type.theteq, 3) == 97.403]) == 6
    assert len([rb for rb in typed_ethane.rb_torsions
                if round(rb.type.c0, 3) == 0.287]) == 9

@pytest.mark.parametrize("ff_filename,kwargs", [
    ("ethane-angle-typo.xml", {"assert_angle_params": False}),
    ("ethane-dihedral-typo.xml", {"assert_dihedral_params": False})
])
def test_missing_topo_params(ff_filename, kwargs):
    """Test that the user is notified if not all topology parameters are found."""
    ethane = mb.load(get_fn('ethane.mol2'))
    oplsaa_with_typo = Forcefield(forcefield_files=get_fn(ff_filename))
    with pytest.raises(Exception):
        ethane = oplsaa_with_typo.apply(ethane)
    with pytest.warns(UserWarning):
        ethane = oplsaa_with_typo.apply(ethane, **kwargs)

def test_overrides_space():
    ethane = mb.load(get_fn('ethane.mol2'))
    ff = Forcefield(forcefield_files=get_fn('overrides-space.xml'))
    typed_ethane = ff.apply(ethane)
    assert typed_ethane.atoms[0].type == 'CT3'
