import difflib
import glob
import os

import gmso
import mbuild as mb
import pytest
from pkg_resources import resource_filename

from foyer.exceptions import FoyerError
from foyer.general_forcefield import Forcefield
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn, register_mock_request

FF_DIR = resource_filename("foyer", "forcefields")
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, "xml/*.xml"))

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

RESPONSE_BIB_ETHANE_JP0484579 = """@article{Jorgensen_2004,
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


class TestGeneralForcefield(BaseTest):
    @pytest.fixture(scope="session")
    def oplsaa(self):
        return Forcefield(name="oplsaa", strict=False)

    @pytest.mark.parametrize("ff_file", FORCEFIELDS)
    def test_load_files(self, ff_file):
        ff1 = Forcefield(forcefield_files=ff_file, strict=False)
        assert len(ff1.ff.atom_types) > 0

        ff2 = Forcefield(forcefield_files=ff_file, strict=False)
        assert len(ff1.ff.atom_types) == len(ff2.ff.atom_types)

    """ Relies on https://github.com/mosdef-hub/gmso/pull/526
    def test_duplicate_type_definitions():
        with pytest.raises(ValueError):
            ff4 = Forcefield(name='oplsaa', forcefield_files=FORCEFIELDS, strict=False)
    """

    def test_missing_type_definitions(self):
        with pytest.raises(FoyerError):
            FF = Forcefield()
            ethane = mb.load(get_fn("ethane.mol2"), backend="parmed")
            FF.apply(ethane, assert_improper_params=False)

    def test_unsupported_backend(self):
        with pytest.raises(FoyerError, match=r"Backend not supported"):
            FF = Forcefield(name="oplsaa", backend="void")

    def test_from_gmso(self, oplsaa):
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")
        top = gmso.external.from_mbuild(mol2)
        ethane = oplsaa.apply(top, assert_improper_params=False)

        assert (
            sum((1 for at in ethane.sites if at.atom_type.name == "opls_135"))
            == 2
        )
        assert (
            sum((1 for at in ethane.sites if at.atom_type.name == "opls_140"))
            == 6
        )
        assert len(ethane.bonds) == 7
        assert all(x.bond_type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.angle_type for x in ethane.angles)
        assert len(ethane.dihedrals) == 9
        assert all(x.dihedral_type for x in ethane.dihedrals)

        """
        Skip test for box information until mbuild box overhaul PR is completed

        mol2 = mb.load(get_fn('ethane.mol2'), backend='parmed')
        mol2.box_vectors = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        oplsaa = Forcefield(name='oplsaa', strict=False)
        ethane = oplsaa.apply(mol2, assert_improper_params=False)

        assert ethane.box_vectors == mol2.box_vectors
        """
    def test_remove_untyped(self):
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")
        
        # test removal of untyped each class of connection seperately
        oplsaa_bond = Forcefield(forcefield_files=get_fn('ethane-missing_bond.xml'),
                                 strict=False)
        oplsaa_angle = Forcefield(forcefield_files=get_fn('ethane-missing_angle.xml'),
                                 strict=False)
        oplsaa_dihedral = Forcefield(forcefield_files=get_fn('ethane-missing_dihedral.xml'),
                                 strict=False)

        ethane1 = oplsaa_bond.apply(mol2,
                                    assert_improper_params=False,
                                    assert_bond_params=False,
                                    remove_untyped_connections=True)

        ethane2 = oplsaa_angle.apply(mol2,
                                     assert_improper_params=False,
                                     assert_angle_params=False,
                                     remove_untyped_connections=True)

        ethane3 = oplsaa_dihedral.apply(mol2,
                                        assert_improper_params=False,
                                        assert_dihedral_params=False,
                                        remove_untyped_connections=True)
    
        assert ethane1.n_bonds == 1
        assert ethane2.n_angles == 6
        assert ethane3.n_dihedrals == 0

    def test_from_mbuild(self, oplsaa):
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")
        ethane = oplsaa.apply(mol2, assert_improper_params=False)

        assert (
            sum((1 for at in ethane.sites if at.atom_type.name == "opls_135"))
            == 2
        )
        assert (
            sum((1 for at in ethane.sites if at.atom_type.name == "opls_140"))
            == 6
        )
        assert len(ethane.bonds) == 7
        assert all(x.bond_type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.angle_type for x in ethane.angles)
        assert len(ethane.dihedrals) == 9
        assert all(x.dihedral_type for x in ethane.dihedrals)

    @pytest.mark.parametrize("mixing_rule", ["lorentz", "geometric"])
    def test_comb_rule(self, mixing_rule, oplsaa):
        mol2 = mb.load(get_fn("ethane.mol2"))
        ethane = oplsaa.apply(
            mol2, combining_rule=mixing_rule, assert_improper_params=False
        )
        assert ethane.combining_rule == mixing_rule

    def test_write_refs(self, requests_mock, oplsaa):
        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/ja9621760/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            text=RESPONSE_BIB_ETHANE_JA962170,
        )
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")
        ethane = oplsaa.apply(
            mol2, references_file="ethane.bib", assert_improper_params=False
        )
        assert os.path.isfile("ethane.bib")
        with open(get_fn("ethane.bib")) as file1:
            with open("ethane.bib") as file2:
                diff = list(
                    difflib.unified_diff(
                        file1.readlines(), file2.readlines(), n=0
                    )
                )
        assert not diff

    def test_write_refs_multiple(self, requests_mock):
        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/ja9621760/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            text=RESPONSE_BIB_ETHANE_JA962170,
        )
        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/jp0484579/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            text=RESPONSE_BIB_ETHANE_JP0484579,
        )
        mol2 = mb.load(get_fn("ethane.mol2"))
        oplsaa = Forcefield(
            forcefield_files=get_fn("refs-multi.xml"), strict=False
        )
        ethane = oplsaa.apply(
            mol2,
            references_file="ethane-multi.bib",
            assert_improper_params=False,
        )
        assert os.path.isfile("ethane-multi.bib")
        with open(get_fn("ethane-multi.bib")) as file1:
            with open("ethane-multi.bib") as file2:
                diff = list(
                    difflib.unified_diff(
                        file1.readlines(), file2.readlines(), n=0
                    )
                )
        assert not diff

    def test_write_bad_ref(self, requests_mock):
        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/garbage_bad_44444444jjjj/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            status_code=404,
        )
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")
        oplsaa = Forcefield(
            forcefield_files=get_fn("refs-bad.xml"), strict=False
        )
        with pytest.warns(UserWarning):
            ethane = oplsaa.apply(
                mol2, references_file="ethane.bib", assert_improper_params=False
            )

    """
    These XML files missed the whole nonbonded force section

    def test_from_mbuild_customtype():
        mol2 = mb.load(get_fn('ethane_customtype.pdb'))
        customtype_ff = Forcefield(forcefield_files=get_fn('validate_customtypes.xml'), strict=False)
        ethane = customtype_ff.apply(mol2, assert_improper_params=False)

        assert sum((1 for at in ethane.sites if at.atom_type.name == 'C3')) == 2
        assert sum((1 for at in ethane.sites if at.atom_type.name == 'Hb')) == 6
        assert len(ethane.bonds) == 7
        assert all(x.bond_type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.angle_type for x in ethane.angles)
        assert len(ethane.dihedrals) == 9
        assert all(x.dihedral_type for x in ethane.dihedrals)

    def test_improper_dihedral():
        untyped_benzene = mb.load(get_fn('benzene.mol2'), backend='parmed')
        ff_improper = Forcefield(forcefield_files=get_fn('improper_dihedral.xml'), strict=False)
        benzene = ff_improper.apply(untyped_benzene, assert_dihedral_params=False, assert_improper_params=False)
        assert len(benzene.dihedrals) == 18
        assert len([dih for dih in benzene.dihedrals if dih.improper]) == 6
        assert len([dih for dih in benzene.dihedrals if not dih.improper]) == 12
    """

    def test_urey_bradley(self):
        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = Forcefield(
            forcefield_files=[get_fn("charmm36_cooh.xml")], strict=False
        )
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        assert len(struc.angles) == 3
        assert len(struc.angle_types) == 3  # 1 harmonic, 2 Urey Bradley

    def test_charmm_improper(self):
        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = Forcefield(
            forcefield_files=[get_fn("charmm36_cooh.xml")], strict=False
        )
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
            remove_untyped_connections =False
        )
        assert len(struc.impropers) == 1
        assert len(struc.dihedrals) == 0

    ''' To be implemented -> Lookup connection types with mixed atomtype-atomclass
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
        ethane = mb.load(get_fn('ethane.mol2'), backend='parmed')
        ff = Forcefield(forcefield_files=get_fn('ethane-topo-precedence.xml'), strict=False)
        typed_ethane = ff.apply(ethane, assert_improper_params=False)
        # Need to work on the units of these test
        assert len([bond for bond in typed_ethane.bonds
                    if round(float(bond.bond_type.parameters['r_eq'].value), 3) == 0.115]) == 6
        assert len([bond for bond in typed_ethane.bonds
                    if round(float(bond.bond_type.parameters['r_eq'].value), 2) == 0.16]) == 1
        assert len([angle for angle in typed_ethane.angles
                    if round(float(angle.angle_type.parameters['theta_eq'].value), 3) == 120.321]) == 6
        assert len([angle for angle in typed_ethane.angles
                    if round(float(angle.angle_type.parameters['theta_eq'].value), 3) == 97.403]) == 6
        assert len([rb for rb in typed_ethane.dihedral
                    if round(float(rb.dihedral_type.parameters['c0'].value), 3) == 0.287]) == 9
    '''

    @pytest.mark.parametrize(
        "ff_filename,kwargs",
        [
            ("ethane-angle-typo.xml", {"assert_angle_params": False}),
            ("ethane-dihedral-typo.xml", {"assert_dihedral_params": False}),
        ],
    )
    def test_missing_topo_params(self, ff_filename, kwargs):
        """Test that the user is notified if not all topology parameters are found."""

        ethane = mb.load(get_fn("ethane.mol2"))
        oplsaa_with_typo = Forcefield(
            forcefield_files=get_fn(ff_filename), strict=False
        )
        with pytest.raises(Exception):
            ethane = oplsaa_with_typo.apply(
                ethane, assert_improper_params=False
            )
        with pytest.warns(UserWarning):
            ethane = oplsaa_with_typo.apply(
                ethane, assert_improper_params=False, **kwargs
            )

    def test_assert_bonds(self):
        ff = Forcefield(name="trappe-ua", strict=False)

        derponium = mb.Compound()
        at1 = mb.Particle(name="H")
        at2 = mb.Particle(name="O")
        at3 = mb.Particle(name="_CH4")

        derponium.add([at1, at2, at3])
        derponium.add_bond((at1, at2))
        derponium.add_bond((at2, at3))

        with pytest.raises(Exception):
            ff.apply(derponium, assert_improper_params=False)
        thing = ff.apply(
            derponium,
            assert_bond_params=False,
            assert_angle_params=False,
            assert_improper_params=False,
            remove_untyped_connections=False,
        )
        assert any(b.bond_type is None for b in thing.bonds)

    def test_apply_subfuncs(self, oplsaa):
        mol2 = mb.load(get_fn("ethane.mol2"), backend="parmed")

        ethane = oplsaa.apply(mol2, assert_improper_params=False)

        typemap = oplsaa._run_atomtyping(mol2, use_residue_map=False)

        ethane2 = oplsaa._parametrize(
            mol2, typemap=typemap, assert_improper_params=False
        )

        assert ethane.box == ethane2.box
        assert (ethane.positions == ethane2.positions).all
        for a1, a2 in zip(ethane.sites, ethane2.sites):
            assert a1.name == a2.name
            assert ethane.get_index(a1) == ethane2.get_index(a2)
            assert a1.atom_type == a2.atom_type

        for b1, b2 in zip(ethane.bonds, ethane2.bonds):
            assert (
                b1.connection_members[0].atom_type
                == b2.connection_members[0].atom_type
            )
            assert (
                b1.connection_members[1].atom_type
                == b2.connection_members[1].atom_type
            )
            assert b1.bond_type == b2.bond_type

    def test_non_zero_charge(self, oplsaa):
        compound = mb.load("C1=CC=C2C(=C1)C(C3=CC=CC=C3O2)C(=O)O", smiles=True)
        with pytest.warns(UserWarning):
            oplsaa.apply(
                compound,
                assert_dihedral_params=False,
                assert_improper_params=False,
            )

    """
    @pytest.mark.parametrize("filename", ['ethane.mol2', 'benzene.mol2'])
    def test_write_xml(filename):
        mol = mb.load(get_fn(filename), backend='parmed')
        oplsaa = Forcefield(name='oplsaa', strict=False)
        typed = oplsaa.apply(mol, assert_improper_params=False)

        typed.write_foyer(filename='opls-snippet.xml', forcefield=oplsaa, unique=True)
        oplsaa_partial = Forcefield('opls-snippet.xml', strict=False)
        typed_by_partial = oplsaa_partial.apply(mol, assert_improper_params=False)

        for i in range(len(typed.sites)):
            atype1 = typed.sites[i].atom_type
            atype2 = typed_by_partial.sites[i].atom_type
            assert atype1.expression == atype2.expression
            assert atype1.parameters == atype2.parameters

        for i in range(len(typed.bonds)):
            btype1 = typed.bonds[i].bond_type
            btype2 = typed_by_partial.bonds[i].bond_type
            assert btype1.expression == btype2.expression
            assert btype1.parameters == btype2.parameters

        # Do it again but with an XML including periodic dihedrals
        mol = mb.load(get_fn(filename), backend='parmed')
        oplsaa = Forcefield(get_fn('oplsaa-periodic.xml'), strict=False)
        typed = oplsaa.apply(mol, assert_improper_params=False)

        typed.write_foyer(filename='opls-snippet.xml', forcefield=oplsaa, unique=True)
        oplsaa_partial = Forcefield('opls-snippet.xml', strict=False)
        typed_by_partial = oplsaa_partial.apply(mol, assert_improper_params=False)

        for i in range(len(typed.sites)):
            atype1 = typed.sites[i].atom_type
            atype2 = typed_by_partial.sites[i].atom_type
            assert atype1.expression == atype2.expression
            assert atype1.parameters == atype2.parameters

        for i in range(len(typed.bonds)):
            btype1 = typed.bonds[i].bond_type
            btype2 = typed_by_partial.bonds[i].bond_type
            assert btype1.expression == btype2.expression
            assert btype1.parameters == btype2.parameters

    @pytest.mark.parametrize("filename", ['ethane.mol2', 'benzene.mol2'])
    def test_write_xml_multiple_periodictorsions(filename):
        cmpd = mb.load(get_fn(filename), backend='parmed')
        ff = Forcefield(forcefield_files=get_fn('oplsaa_multiperiodicitytorsion.xml'), strict=False)
        typed_struc = ff.apply(cmpd, assert_dihedral_params=False, assert_improper_params=False)
        typed_struc.write_foyer(filename='multi-periodictorsions.xml', forcefield=ff, unique=True)

        partial_ff = Forcefield(forcefield_files='multi-periodictorsions.xml', strict=False)
        typed_by_partial = partial_ff.apply(cmpd, assert_dihedral_params=False, assert_improper_params=False)

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
        mol = mb.load(get_fn(filename), backend='parmed')
        if filename == 'ethane.mol2':
            ff = Forcefield(get_fn('ethane-multiple.xml'), strict=False)
        else:
            ff = Forcefield(name='oplsaa', strict=False)
        typed = ff.apply(mol, assert_improper_params=False)
        typed.write_foyer(filename='snippet.xml', forcefield=ff, unique=True)

        generated_ff = Forcefield('snippet.xml', strict=False)

    def test_write_xml_overrides():
        #Test xml_writer new overrides and comments features
        mol = mb.load(get_fn('styrene.mol2'), backend='parmed')
        oplsaa = Forcefield(name='oplsaa', strict=False)
        typed = oplsaa.apply(mol, assert_dihedral_params=False, assert_improper_params=False)
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
        lj_ff = Forcefield(get_fn('lj.xml'), strict=False)
        assert lj_ff.version == '0.4.1'
        assert lj_ff.name == 'LJ'

        lj_ff = Forcefield(forcefield_files=[get_fn('lj.xml'), get_fn('lj2.xml')])
        assert lj_ff.version == ['0.4.1', '4.8.2']
        assert lj_ff.name == ['LJ', 'JL']
    """
