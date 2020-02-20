import difflib
import glob
import os
from pkg_resources import resource_filename

from lxml import etree as ET

import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.forcefield import generate_topology
from foyer.forcefield import _check_independent_residues
from foyer.exceptions import FoyerError, ValidationWarning
from foyer.tests.utils import get_fn, register_mock_request
from foyer.utils.io import has_mbuild


FF_DIR = resource_filename('foyer', 'forcefields')
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, 'xml/*.xml'))

RESPONSE_BIB_ETHANE_JA962170 = """@article{Jorgensen_1996,
	doi = {10.1021/ja9621760},
	url = {https://doi.org/10.1021%2Fja9621760},
	year = 1996,
	month = {jan},
	publisher = {American Chemical Society ({ACS})},
	volume = {118},
	number = {45},
	pages = {11225--11236},
	author = {William L. Jorgensen and David S. Maxwell and Julian Tirado-Rives},
	title = {Development and Testing of the {OPLS} All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids},
	journal = {Journal of the American Chemical Society}
}"""

RESPONSE_BIB_ETHANE_JP0484579="""@article{Jorgensen_2004,
	doi = {10.1021/jp0484579},
	url = {https://doi.org/10.1021%2Fjp0484579},
	year = 2004,
	month = {oct},
	publisher = {American Chemical Society ({ACS})},
	volume = {108},
	number = {41},
	pages = {16264--16270},
	author = {William L. Jorgensen and Jakob P. Ulmschneider and Julian Tirado-Rives},
	title = {Free Energies of Hydration from a Generalized Born Model and an All-Atom Force Field},
	journal = {The Journal of Physical Chemistry B}
}"""


def test_load_files():
    for ff_file in FORCEFIELDS:
        ff1 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) > 0

        ff2 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) == len(ff2._atomTypes)


def test_duplicate_type_definitions():
    with pytest.raises(ValueError):
        ff4 = Forcefield(name='oplsaa', forcefield_files=FORCEFIELDS)


def test_missing_type_definitions():
    with pytest.raises(FoyerError):
        FF = Forcefield()
        ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
        FF.apply(ethane)

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

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_from_mbuild():
    import mbuild as mb
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

@pytest.mark.parametrize("mixing_rule", ['lorentz', 'geometric'])
@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_comb_rule(mixing_rule):
    import mbuild as mb
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2, combining_rule=mixing_rule)
    assert ethane.combining_rule == mixing_rule

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_write_refs(requests_mock):
    import mbuild as mb
    register_mock_request(mocker=requests_mock,
                          url='http://api.crossref.org/',
                          path='works/10.1021/ja9621760/transform/application/x-bibtex',
                          headers={'accept': 'application/x-bibtex'},
                          text=RESPONSE_BIB_ETHANE_JA962170)
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(mol2, references_file='ethane.bib')
    assert os.path.isfile('ethane.bib')
    with open(get_fn('ethane.bib')) as file1:
        with open('ethane.bib') as file2:
            diff = list(difflib.unified_diff(file1.readlines(),
                                             file2.readlines(),
                                             n=0))
    assert not diff

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_write_refs_multiple(requests_mock):
    import mbuild as mb
    register_mock_request(mocker=requests_mock,
                          url='http://api.crossref.org/',
                          path='works/10.1021/ja9621760/transform/application/x-bibtex',
                          headers={'accept': 'application/x-bibtex'},
                          text=RESPONSE_BIB_ETHANE_JA962170)
    register_mock_request(mocker=requests_mock,
                          url='http://api.crossref.org/',
                          path='works/10.1021/jp0484579/transform/application/x-bibtex',
                          headers={'accept': 'application/x-bibtex'},
                          text=RESPONSE_BIB_ETHANE_JP0484579)
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

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_write_bad_ref(requests_mock):
    import mbuild as mb
    register_mock_request(mocker=requests_mock,
                          url='http://api.crossref.org/',
                          path='works/10.1021/garbage_bad_44444444jjjj/transform/application/x-bibtex',
                          headers={'accept': 'application/x-bibtex'},
                          status_code=404)
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(forcefield_files=get_fn('refs-bad.xml'))
    with pytest.warns(UserWarning):
        ethane = oplsaa.apply(mol2, references_file='ethane.bib')

def test_preserve_resname():
    untyped_ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    untyped_resname = untyped_ethane.residues[0].name
    oplsaa = Forcefield(name='oplsaa')
    typed_ethane = oplsaa.apply(untyped_ethane)
    typed_resname = typed_ethane.residues[0].name
    assert typed_resname == untyped_resname

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_apply_residues():
    import mbuild.recipes
    propane = mbuild.recipes.Alkane(n=3)

    opls = Forcefield(name='oplsaa')
    typed = opls.apply(propane, residues='CH3')
    assert len([res for res in typed.residues if res.name == 'CH3']) == 2

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_from_mbuild_customtype():
    import mbuild as mb
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

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_urey_bradley():
    import mbuild as mb
    system = mb.Compound()
    first = mb.Particle(name='_CTL2',pos=[-1,0,0])
    second = mb.Particle(name='_CL', pos=[0,0,0])
    third = mb.Particle(name='_OBL', pos=[1,0,0])
    fourth = mb.Particle(name='_OHL', pos = [0,1,0])

    system.add([first, second, third, fourth])

    system.add_bond((first,second))
    system.add_bond((second, third))
    system.add_bond((second, fourth))

    ff = Forcefield(forcefield_files=[get_fn('charmm36_cooh.xml')])
    struc = ff.apply(system, assert_angle_params=False, assert_dihedral_params=False,
            assert_improper_params=False)
    assert len(struc.angles) == 3
    assert len(struc.urey_bradleys) ==2

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_charmm_improper():
    import mbuild as mb
    system = mb.Compound()
    first = mb.Particle(name='_CTL2',pos=[-1,0,0])
    second = mb.Particle(name='_CL', pos=[0,0,0])
    third = mb.Particle(name='_OBL', pos=[1,0,0])
    fourth = mb.Particle(name='_OHL', pos = [0,1,0])

    system.add([first, second, third, fourth])

    system.add_bond((first,second))
    system.add_bond((second, third))
    system.add_bond((second, fourth))

    ff = Forcefield(forcefield_files=[get_fn('charmm36_cooh.xml')])
    struc = ff.apply(system, assert_angle_params=False, assert_dihedral_params=False,
            assert_improper_params=False)
    assert len(struc.impropers) == 1
    assert len(struc.dihedrals) == 0

def test_residue_map():
    ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    ethane *= 2
    oplsaa = Forcefield(name='oplsaa')
    map_with = oplsaa.run_atomtyping(ethane, use_residue_map=True)
    map_without = oplsaa.run_atomtyping(ethane, use_residue_map=False)
    assert all([a['atomtype'] for a in map_with.values()])
    assert all([a['atomtype'] for a in map_without.values()])
    struct_with = ethane
    struct_without = ethane
    oplsaa._apply_typemap(struct_with, map_with)
    oplsaa._apply_typemap(struct_without, map_without)
    for atom_with, atom_without in zip(struct_with.atoms, struct_without.atoms):
        assert atom_with.type == atom_without.type
        b_with = atom_with.bond_partners
        b_without = atom_without.bond_partners
        assert [a0.type for a0 in b_with] == [a1.type for a1 in b_without]
        assert [a0.idx for a0 in b_with] == [a1.idx for a1 in b_without]


@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_independent_residues_molecules():
    """Test to see that _check_independent_residues works for molecules."""
    import mbuild.recipes
    butane = mbuild.recipes.Alkane(4)
    structure = butane.to_parmed()
    assert _check_independent_residues(structure)
    structure = butane.to_parmed(residues=['RES', 'CH3'])
    assert not _check_independent_residues(structure)

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_independent_residues_atoms():
    """Test to see that _check_independent_residues works for single aotms."""
    import mbuild as mb
    argon = mb.Compound()
    argon.name = 'Ar'
    structure = argon.to_parmed()
    assert _check_independent_residues(structure)

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
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
    import mbuild as mb
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

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
@pytest.mark.parametrize("ff_filename,kwargs", [
    ("ethane-angle-typo.xml", {"assert_angle_params": False}),
    ("ethane-dihedral-typo.xml", {"assert_dihedral_params": False})
])
def test_missing_topo_params(ff_filename, kwargs):
    """Test that the user is notified if not all topology parameters are found."""
    import mbuild as mb
    ethane = mb.load(get_fn('ethane.mol2'))
    oplsaa_with_typo = Forcefield(forcefield_files=get_fn(ff_filename))
    with pytest.raises(Exception):
        ethane = oplsaa_with_typo.apply(ethane)
    with pytest.warns(UserWarning):
        ethane = oplsaa_with_typo.apply(ethane, **kwargs)

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_overrides_space():
    import mbuild as mb
    ethane = mb.load(get_fn('ethane.mol2'))
    ff = Forcefield(forcefield_files=get_fn('overrides-space.xml'))
    typed_ethane = ff.apply(ethane)
    assert typed_ethane.atoms[0].type == 'CT3'

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_allow_empty_def():
    import mbuild as mb
    ethane = mb.load(get_fn('ethane.mol2'))
    with pytest.warns(ValidationWarning):
        ff = Forcefield(forcefield_files=get_fn('empty_def.xml'))
    ff.apply(ethane)

@pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
def test_assert_bonds():
    import mbuild as mb
    ff = Forcefield(name='trappe-ua')

    derponium = mb.Compound()
    at1 = mb.Particle(name='H')
    at2 = mb.Particle(name='O')
    at3 = mb.Particle(name='_CH4')

    derponium.add([at1, at2, at3])
    derponium.add_bond((at1, at2))
    derponium.add_bond((at2, at3))

    with pytest.raises(Exception):
        ff.apply(derponium)
    thing = ff.apply(derponium, assert_bond_params=False, assert_angle_params=False)
    assert any(b.type is None for b in thing.bonds)

def test_apply_subfuncs():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    oplsaa = Forcefield(name='oplsaa')

    ethane = oplsaa.apply(mol2)

    typemap = oplsaa.run_atomtyping(mol2, use_residue_map=False)
    oplsaa._apply_typemap(mol2, typemap)
    ethane2 = oplsaa.parametrize_system(mol2)

    # Note: Check ParmEd issue #1067 to see if __eq__ is implemented
    # assert ethane == ethane2
    assert ethane.box == ethane2.box
    assert ethane.positions == ethane2.positions
    for a1, a2 in zip(ethane.atoms, ethane2.atoms):
        assert a1.name == a2.name
        assert a1.idx == a2.idx
        assert a1.atom_type == a2.atom_type

    for b1, b2 in zip(ethane.bonds, ethane2.bonds):
        assert b1.atom1.atom_type == b2.atom1.atom_type
        assert b1.atom2.atom_type == b2.atom2.atom_type
        assert b1.type == b2.type

    for ang1, ang2 in zip(ethane.angles, ethane2.angles):
        assert ang1.type == ang2.type

@pytest.mark.parametrize("filename", ['ethane.mol2', 'benzene.mol2'])
def test_write_xml(filename):
    mol = pmd.load_file(get_fn(filename), structure=True)
    oplsaa = Forcefield(name='oplsaa')
    typed = oplsaa.apply(mol)

    typed.write_foyer(filename='opls-snippet.xml', forcefield=oplsaa, unique=True)
    oplsaa_partial = Forcefield('opls-snippet.xml')
    typed_by_partial = oplsaa_partial.apply(mol)

    for adj in typed.adjusts:
        type1 = adj.atom1.atom_type
        type2 = adj.atom1.atom_type
        sigma_factor_pre = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
        epsilon_factor_pre = adj.type.epsilon / ((type1.epsilon * type2.epsilon) ** 0.5)

    for adj in typed_by_partial.adjusts:
        type1 = adj.atom1.atom_type
        type2 = adj.atom1.atom_type
        sigma_factor_post = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
        epsilon_factor_post = adj.type.epsilon / ((type1.epsilon * type2.epsilon) ** 0.5)

    assert sigma_factor_pre == sigma_factor_post
    assert epsilon_factor_pre == epsilon_factor_post

    # Do it again but with an XML including periodic dihedrals
    mol = pmd.load_file(get_fn(filename), structure=True)
    oplsaa = Forcefield(get_fn('oplsaa-periodic.xml'))
    typed = oplsaa.apply(mol)

    typed.write_foyer(filename='opls-snippet.xml', forcefield=oplsaa, unique=True)
    oplsaa_partial = Forcefield('opls-snippet.xml')
    typed_by_partial = oplsaa_partial.apply(mol)

    for adj in typed.adjusts:
        type1 = adj.atom1.atom_type
        type2 = adj.atom1.atom_type
        sigma_factor_pre = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
        epsilon_factor_pre = adj.type.epsilon / ((type1.epsilon * type2.epsilon) ** 0.5)

    for adj in typed_by_partial.adjusts:
        type1 = adj.atom1.atom_type
        type2 = adj.atom1.atom_type
        sigma_factor_post = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
        epsilon_factor_post = adj.type.epsilon / ((type1.epsilon * type2.epsilon) ** 0.5)

    assert sigma_factor_pre == sigma_factor_post
    assert epsilon_factor_pre == epsilon_factor_post

@pytest.mark.parametrize("filename", ['ethane.mol2', 'benzene.mol2'])
def test_write_xml_multiple_periodictorsions(filename):
    cmpd = pmd.load_file(get_fn(filename), structure=True)
    ff = Forcefield(forcefield_files=get_fn('oplsaa_multiperiodicitytorsion.xml'))
    typed_struc = ff.apply(cmpd, assert_dihedral_params=False)
    typed_struc.write_foyer(filename='multi-periodictorsions.xml', forcefield=ff, unique=True)

    partial_ff = Forcefield(forcefield_files='multi-periodictorsions.xml')
    typed_by_partial = partial_ff.apply(cmpd, assert_dihedral_params=False)

    assert len(typed_struc.bonds) == len(typed_by_partial.bonds)
    assert len(typed_struc.angles) == len(typed_by_partial.angles)
    assert len(typed_struc.dihedrals) == len(typed_by_partial.dihedrals)

    root = ET.parse('multi-periodictorsions.xml')
    periodic_element = root.find('PeriodicTorsionForce')
    assert 'periodicity2' in periodic_element[0].attrib
    assert 'k2' in periodic_element[0].attrib
    assert 'phase2' in periodic_element[0].attrib

@pytest.mark.parametrize("filename", ['ethane.mol2', 'benzene.mol2'])
def test_load_xml(filename):
    mol = pmd.load_file(get_fn(filename), structure=True)
    if filename == 'ethane.mol2':
        ff = Forcefield(get_fn('ethane-multiple.xml'))
    else:
        ff = Forcefield(name='oplsaa')
    typed = ff.apply(mol)
    typed.write_foyer(filename='snippet.xml', forcefield=ff, unique=True)

    generated_ff = Forcefield('snippet.xml')

def test_write_xml_overrides():
    #Test xml_writer new overrides and comments features
    mol = pmd.load_file(get_fn('styrene.mol2'), structure=True)
    oplsaa = Forcefield(name='oplsaa')
    typed = oplsaa.apply(mol, assert_dihedral_params=False)
    typed.write_foyer(filename='opls-styrene.xml', forcefield=oplsaa, unique=True)
    styrene = ET.parse('opls-styrene.xml')
    atom_types = styrene.getroot().find('AtomTypes').findall('Type')
    for item in atom_types:
        attributes = item.attrib
        if attributes['name'] == 'opls_145':
            assert attributes['overrides'] == 'opls_142'
            assert str(item.xpath('comment()')) in {'[<!--Note: original overrides="opls_141,opls_142"-->]',
                                                    '[<!--Note: original overrides="opls_142,opls_141"-->]'}
        elif attributes['name'] == 'opls_146':
            assert attributes['overrides'] == 'opls_144'
            assert str(item.xpath('comment()')) == '[<!--Note: original overrides="opls_144"-->]'

def test_load_metadata():
    lj_ff = Forcefield(get_fn('lj.xml'))
    assert lj_ff.version == '0.4.1'
    assert lj_ff.name == 'LJ'

    lj_ff = Forcefield(forcefield_files=[get_fn('lj.xml'), get_fn('lj2.xml')])
    assert lj_ff.version == ['0.4.1', '4.8.2']
    assert lj_ff.name == ['LJ', 'JL']
