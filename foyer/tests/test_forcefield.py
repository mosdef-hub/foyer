import difflib
import glob
import importlib.resources as resources
import os
from typing import List

import parmed as pmd
import pytest
from lxml import etree as ET
from parmed.gromacs.gromacstop import _Defaults

from foyer import Forcefield, forcefields
from foyer.exceptions import (
    FoyerError,
    UnimplementedCombinationRuleError,
    ValidationWarning,
)
from foyer.forcefield import (
    _check_independent_residues,
    _structure_from_residue,
)
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn, register_mock_request
from foyer.utils.io import has_mbuild

FF_DIR = resources.files("foyer").joinpath("forcefields")
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


class TestForcefield(BaseTest):
    @pytest.mark.parametrize("ff_file", FORCEFIELDS)
    def test_load_files(self, ff_file):
        ff1 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) > 0

        ff2 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) == len(ff2._atomTypes)

    def test_duplicate_type_definitions(self):
        with pytest.raises(ValueError):
            Forcefield(name="oplsaa", forcefield_files=FORCEFIELDS)

    def test_missing_type_definitions(self):
        with pytest.raises(FoyerError):
            FF = Forcefield()
            ethane = pmd.load_file(get_fn("ethane.mol2"), structure=True)
            FF.apply(ethane)

    def test_from_parmed(self, oplsaa):
        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        ethane = oplsaa.apply(mol2)

        assert sum((1 for at in ethane.atoms if at.type == "opls_135")) == 2
        assert sum((1 for at in ethane.atoms if at.type == "opls_140")) == 6
        assert len(ethane.bonds) == 7
        assert all(x.type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.type for x in ethane.angles)
        assert len(ethane.rb_torsions) == 9
        assert all(x.type for x in ethane.dihedrals)

        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        mol2.box_vectors = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        ethane = oplsaa.apply(mol2)

        assert ethane.box_vectors == mol2.box_vectors

    def test_structure_meta(self, oplsaa):
        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        ethane = oplsaa.apply(mol2)

        assert isinstance(ethane.defaults, _Defaults)
        assert ethane.defaults.nbfunc == 1
        assert ethane.defaults.comb_rule == 3
        assert ethane.defaults.fudgeLJ == oplsaa.lj14scale
        assert ethane.defaults.fudgeQQ == oplsaa.coulomb14scale
        assert ethane.defaults.gen_pairs == "yes"

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_from_mbuild(self, oplsaa):
        import mbuild as mb

        mol2 = mb.load(get_fn("ethane.mol2"))
        ethane = oplsaa.apply(mol2)

        assert sum((1 for at in ethane.atoms if at.type == "opls_135")) == 2
        assert sum((1 for at in ethane.atoms if at.type == "opls_140")) == 6
        assert len(ethane.bonds) == 7
        assert all(x.type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.type for x in ethane.angles)
        assert len(ethane.rb_torsions) == 9
        assert all(x.type for x in ethane.dihedrals)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_write_refs(self, requests_mock, oplsaa):
        import mbuild as mb

        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/ja9621760/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            text=RESPONSE_BIB_ETHANE_JA962170,
        )
        mol2 = mb.load(get_fn("ethane.mol2"))
        oplsaa.apply(mol2, references_file="ethane.bib")
        assert os.path.isfile("ethane.bib")
        with open(get_fn("ethane.bib")) as file1:
            with open("ethane.bib") as file2:
                diff = list(
                    difflib.unified_diff(file1.readlines(), file2.readlines(), n=0)
                )
        assert not diff

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_write_refs_multiple(self, requests_mock):
        import mbuild as mb

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
        oplsaa = Forcefield(forcefield_files=get_fn("refs-multi.xml"))
        oplsaa.apply(mol2, references_file="ethane-multi.bib")
        assert os.path.isfile("ethane-multi.bib")
        with open(get_fn("ethane-multi.bib")) as file1:
            with open("ethane-multi.bib") as file2:
                diff = list(
                    difflib.unified_diff(file1.readlines(), file2.readlines(), n=0)
                )
        assert not diff

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_write_bad_ref(self, requests_mock):
        import mbuild as mb

        register_mock_request(
            mocker=requests_mock,
            url="http://api.crossref.org/",
            path="works/10.1021/garbage_bad_44444444jjjj/transform/application/x-bibtex",
            headers={"accept": "application/x-bibtex"},
            status_code=404,
        )
        mol2 = mb.load(get_fn("ethane.mol2"))
        oplsaa = Forcefield(forcefield_files=get_fn("refs-bad.xml"))
        with pytest.warns(UserWarning):
            oplsaa.apply(mol2, references_file="ethane.bib")

    def test_preserve_resname(self, oplsaa):
        untyped_ethane = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        untyped_resname = untyped_ethane.residues[0].name
        typed_ethane = oplsaa.apply(untyped_ethane)
        typed_resname = typed_ethane.residues[0].name
        assert typed_resname == untyped_resname

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_apply_residues(self, oplsaa):
        import mbuild.recipes

        propane = mbuild.recipes.Alkane(n=3)

        typed = oplsaa.apply(propane, residues="CH3")
        assert len([res for res in typed.residues if res.name == "CH3"]) == 2

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_from_mbuild_customtype(self):
        import mbuild as mb

        mol2 = mb.load(get_fn("ethane_customtype.pdb"))
        customtype_ff = Forcefield(forcefield_files=get_fn("validate_customtypes.xml"))
        ethane = customtype_ff.apply(mol2)

        assert sum((1 for at in ethane.atoms if at.type == "C3")) == 2
        assert sum((1 for at in ethane.atoms if at.type == "Hb")) == 6
        assert len(ethane.bonds) == 7
        assert all(x.type for x in ethane.bonds)
        assert len(ethane.angles) == 12
        assert all(x.type for x in ethane.angles)
        assert len(ethane.rb_torsions) == 9
        assert all(x.type for x in ethane.dihedrals)

    def test_improper_dihedral(self):
        untyped_benzene = pmd.load_file(get_fn("benzene.mol2"), structure=True)
        ff_improper = Forcefield(forcefield_files=get_fn("improper_dihedral.xml"))
        benzene = ff_improper.apply(untyped_benzene, assert_dihedral_params=False)
        assert len(benzene.dihedrals) == 18
        assert len([dih for dih in benzene.dihedrals if dih.improper]) == 6
        assert len([dih for dih in benzene.dihedrals if not dih.improper]) == 12

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_urey_bradley(self):
        import mbuild as mb

        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = Forcefield(forcefield_files=[get_fn("charmm36_cooh.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        assert len(struc.angles) == 3
        assert len(struc.urey_bradleys) == 2

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_charmm_improper(self):
        import mbuild as mb

        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = Forcefield(forcefield_files=[get_fn("charmm36_cooh.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        assert len(struc.impropers) == 1
        assert len(struc.dihedrals) == 0

    def test_residue_map(self, oplsaa):
        ethane = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        ethane *= 2
        map_with = oplsaa.run_atomtyping(ethane, use_residue_map=True)
        map_without = oplsaa.run_atomtyping(ethane, use_residue_map=False)
        assert all([a["atomtype"] for a in map_with.values()])
        assert all([a["atomtype"] for a in map_without.values()])
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
    def test_independent_residues_molecules(self):
        """Test to see that _check_independent_residues works for molecules."""
        import mbuild.recipes

        butane = mbuild.recipes.Alkane(4)
        structure = butane.to_parmed()
        assert _check_independent_residues(structure)
        structure = butane.to_parmed(residues=["RES", "CH3"])
        assert not _check_independent_residues(structure)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_independent_residues_atoms(self):
        """Test to see that _check_independent_residues works for single aotms."""
        import mbuild as mb

        argon = mb.Compound()
        argon.name = "Ar"
        structure = argon.to_parmed()
        assert _check_independent_residues(structure)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_topology_precedence(self):
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

        ethane = mb.load(get_fn("ethane.mol2"))
        ff = Forcefield(forcefield_files=get_fn("ethane-topo-precedence.xml"))
        typed_ethane = ff.apply(ethane)

        assert (
            len(
                [bond for bond in typed_ethane.bonds if round(bond.type.req, 2) == 1.15]
            )
            == 6
        )
        assert (
            len([bond for bond in typed_ethane.bonds if round(bond.type.req, 2) == 1.6])
            == 1
        )
        assert (
            len(
                [
                    angle
                    for angle in typed_ethane.angles
                    if round(angle.type.theteq, 3) == 120.321
                ]
            )
            == 6
        )
        assert (
            len(
                [
                    angle
                    for angle in typed_ethane.angles
                    if round(angle.type.theteq, 3) == 97.403
                ]
            )
            == 6
        )
        assert (
            len(
                [rb for rb in typed_ethane.rb_torsions if round(rb.type.c0, 3) == 0.287]
            )
            == 9
        )

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    @pytest.mark.parametrize(
        "ff_filename, kwargs",
        [
            ("ethane-angle-typo.xml", {"assert_angle_params": False}),
            ("ethane-dihedral-typo.xml", {"assert_dihedral_params": False}),
        ],
    )
    def test_missing_topo_params(self, ff_filename, kwargs):
        """Test that the user is notified if not all topology parameters are found."""
        import mbuild as mb

        ethane = mb.load(get_fn("ethane.mol2"))
        oplsaa_with_typo = Forcefield(forcefield_files=get_fn(ff_filename))
        with pytest.raises(Exception):
            ethane = oplsaa_with_typo.apply(ethane)
        with pytest.warns(UserWarning):
            ethane = oplsaa_with_typo.apply(ethane, **kwargs)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_overrides_space(self):
        import mbuild as mb

        ethane = mb.load(get_fn("ethane.mol2"))
        ff = Forcefield(forcefield_files=get_fn("overrides-space.xml"))
        typed_ethane = ff.apply(ethane)
        assert typed_ethane.atoms[0].type == "CT3"

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_allow_empty_def(self):
        import mbuild as mb

        ethane = mb.load(get_fn("ethane.mol2"))
        with pytest.warns(ValidationWarning):
            ff = Forcefield(forcefield_files=get_fn("empty_def.xml"))

        with pytest.warns(UserWarning):
            typed = ff.apply(ethane)
        assert typed.defaults is None

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_assert_bonds(self):
        import mbuild as mb

        ff = Forcefield(name="trappe-ua")

        derponium = mb.Compound()
        at1 = mb.Particle(name="H")
        at2 = mb.Particle(name="O")
        at3 = mb.Particle(name="_CH4")

        derponium.add([at1, at2, at3])
        derponium.add_bond((at1, at2))
        derponium.add_bond((at2, at3))

        with pytest.raises(Exception):
            ff.apply(derponium)
        thing = ff.apply(derponium, assert_bond_params=False, assert_angle_params=False)
        assert any(b.type is None for b in thing.bonds)

    def test_apply_subfuncs(self, oplsaa):
        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)

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

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_non_zero_charge(self, oplsaa):
        import mbuild as mb

        compound = mb.load("C1=CC=C2C(=C1)C(C3=CC=CC=C3O2)C(=O)O", smiles=True)
        with pytest.warns(UserWarning):
            oplsaa.apply(compound, assert_dihedral_params=False)

    @pytest.mark.parametrize("filename", ["ethane.mol2", "benzene.mol2"])
    def test_write_xml(self, filename, oplsaa):
        mol = pmd.load_file(get_fn(filename), structure=True)
        typed = oplsaa.apply(mol)

        typed.write_foyer(
            filename="opls-snippet.xml",
            name="oplsaa-snippet",
            version="0.1.0",
            forcefield=oplsaa,
            unique=True,
        )
        oplsaa_partial = Forcefield("opls-snippet.xml")
        assert oplsaa_partial.name == "oplsaa-snippet"
        assert oplsaa_partial.version == "0.1.0"
        assert oplsaa_partial.combining_rule == "geometric"
        typed_by_partial = oplsaa_partial.apply(mol)

        for adj in typed.adjusts:
            type1 = adj.atom1.atom_type
            type2 = adj.atom1.atom_type
            sigma_factor_pre = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
            epsilon_factor_pre = adj.type.epsilon / (
                (type1.epsilon * type2.epsilon) ** 0.5
            )

        for adj in typed_by_partial.adjusts:
            type1 = adj.atom1.atom_type
            type2 = adj.atom1.atom_type
            sigma_factor_post = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
            epsilon_factor_post = adj.type.epsilon / (
                (type1.epsilon * type2.epsilon) ** 0.5
            )

        assert sigma_factor_pre == sigma_factor_post
        assert epsilon_factor_pre == epsilon_factor_post

        # Do it again but with an XML including periodic dihedrals
        mol = pmd.load_file(get_fn(filename), structure=True)
        oplsaa = Forcefield(get_fn("oplsaa-periodic.xml"))
        typed = oplsaa.apply(mol)

        typed.write_foyer(filename="opls-snippet.xml", forcefield=oplsaa, unique=True)
        oplsaa_partial = Forcefield("opls-snippet.xml")
        typed_by_partial = oplsaa_partial.apply(mol)

        for adj in typed.adjusts:
            type1 = adj.atom1.atom_type
            type2 = adj.atom1.atom_type
            sigma_factor_pre = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
            epsilon_factor_pre = adj.type.epsilon / (
                (type1.epsilon * type2.epsilon) ** 0.5
            )

        for adj in typed_by_partial.adjusts:
            type1 = adj.atom1.atom_type
            type2 = adj.atom1.atom_type
            sigma_factor_post = adj.type.sigma / ((type1.sigma + type2.sigma) / 2)
            epsilon_factor_post = adj.type.epsilon / (
                (type1.epsilon * type2.epsilon) ** 0.5
            )

        assert sigma_factor_pre == sigma_factor_post
        assert epsilon_factor_pre == epsilon_factor_post

    @pytest.mark.parametrize("filename", ["ethane.mol2", "benzene.mol2"])
    def test_write_xml_multiple_periodictorsions(self, filename):
        cmpd = pmd.load_file(get_fn(filename), structure=True)
        ff = Forcefield(forcefield_files=get_fn("oplsaa_multiperiodicitytorsion.xml"))
        typed_struc = ff.apply(cmpd, assert_dihedral_params=False)
        typed_struc.write_foyer(
            filename="multi-periodictorsions.xml", forcefield=ff, unique=True
        )

        partial_ff = Forcefield(forcefield_files="multi-periodictorsions.xml")
        typed_by_partial = partial_ff.apply(cmpd, assert_dihedral_params=False)

        assert len(typed_struc.bonds) == len(typed_by_partial.bonds)
        assert len(typed_struc.angles) == len(typed_by_partial.angles)
        assert len(typed_struc.dihedrals) == len(typed_by_partial.dihedrals)

        root = ET.parse("multi-periodictorsions.xml")
        periodic_element = root.find("PeriodicTorsionForce")
        assert "periodicity2" in periodic_element[0].attrib
        assert "k2" in periodic_element[0].attrib
        assert "phase2" in periodic_element[0].attrib

    @pytest.mark.parametrize("filename", ["ethane.mol2", "benzene.mol2"])
    def test_load_xml(self, filename, oplsaa):
        mol = pmd.load_file(get_fn(filename), structure=True)
        if filename == "ethane.mol2":
            ff = Forcefield(get_fn("ethane-multiple.xml"))
        else:
            ff = oplsaa
        typed = ff.apply(mol)
        typed.write_foyer(filename="snippet.xml", forcefield=ff, unique=True)

        Forcefield("snippet.xml")

    def test_write_xml_overrides(self, oplsaa):
        # Test xml_writer new overrides and comments features
        mol = pmd.load_file(get_fn("styrene.mol2"), structure=True)
        typed = oplsaa.apply(mol, assert_dihedral_params=False)
        typed.write_foyer(filename="opls-styrene.xml", forcefield=oplsaa, unique=True)
        styrene = ET.parse("opls-styrene.xml")
        atom_types = styrene.getroot().find("AtomTypes").findall("Type")
        for item in atom_types:
            attributes = item.attrib
            if attributes["name"] == "opls_145":
                assert attributes["overrides"] == "opls_142"
                assert str(item.xpath("comment()")) in {
                    '[<!--Note: original overrides="opls_141,opls_142"-->]',
                    '[<!--Note: original overrides="opls_142,opls_141"-->]',
                }
            elif attributes["name"] == "opls_146":
                assert attributes["overrides"] == "opls_144"
                assert (
                    str(item.xpath("comment()"))
                    == '[<!--Note: original overrides="opls_144"-->]'
                )

    def test_load_metadata(self):
        lj_ff = Forcefield(get_fn("lj.xml"))
        assert lj_ff.version == "0.4.1"
        assert lj_ff.name == "LJ"

        lj_ff = Forcefield(forcefield_files=[get_fn("lj.xml"), get_fn("lj2.xml")])
        assert lj_ff.version == ["0.4.1", "4.8.2"]
        assert lj_ff.name == ["LJ", "JL"]

    def test_load_metadata_single_xml(self):
        from_xml_ff = Forcefield(forcefield_files=get_fn("lj.xml"))
        assert from_xml_ff.version == "0.4.1"
        assert from_xml_ff.name == "LJ"

    def test_load_metadata_list_xml(self):
        from_xml_ff = Forcefield(forcefield_files=[get_fn("lj.xml"), get_fn("lj2.xml")])
        assert isinstance(from_xml_ff.version, List)
        assert isinstance(from_xml_ff.name, List)
        assert all([x in from_xml_ff.version for x in ["0.4.1", "4.8.2"]])
        assert all([x in from_xml_ff.name for x in ["JL", "LJ"]])

        with pytest.raises(FoyerError):
            Forcefield(forcefield_files=[get_fn("lj.xml"), get_fn("lj3.xml")])

    def test_load_metadata_from_internal_forcefield_plugin_loader(self):
        from_xml_ff = forcefields.load_OPLSAA()
        assert from_xml_ff.version == "0.1.0"
        assert from_xml_ff.name == "OPLS-AA"
        assert from_xml_ff.combining_rule == "geometric"

    def test_load_metadata_from_internal_name(self):
        from_xml_ff = Forcefield(name="oplsaa")
        assert from_xml_ff.version == "0.1.0"
        assert from_xml_ff.name == "OPLS-AA"
        assert from_xml_ff.combining_rule == "geometric"

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_no_overlap_residue_atom_overlap(self):
        import mbuild as mb

        mol1 = mb.load("CCC", smiles=True)
        mol2 = mb.load("COC", smiles=True)

        mol1.name = "CCC"
        mol2.name = "COC"

        box = mb.fill_box([mol1, mol2], n_compounds=[2, 2], overlap=0.01, density=700)

        all_substructures = []
        structure = box.to_parmed(residues=["CCC", "COC"])
        for res_id, res in enumerate(structure.residues):
            all_substructures.append(_structure_from_residue(res, structure))

        residue_idx_per_atom = map(lambda x: x.residue.idx, structure.atoms)
        num_unique_residue_indices = len(set([*residue_idx_per_atom]))
        num_residues = len(structure.residues)

        assert num_residues == num_unique_residue_indices

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_unknown_combining_rule(self, oplsaa):
        import mbuild as mb

        oplsaa._combining_rule = "bogus"

        benzene = mb.load("c1ccccc1", smiles=True)

        with pytest.raises(
            UnimplementedCombinationRuleError, match="bogus is not impl"
        ):
            oplsaa.apply(structure=benzene)

        # Reset the combining rule to its original value; re-using the same fixture
        # is probably more performant than re-loading it on each use
        oplsaa._combining_rule = "geometric"

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_combining_rule_in_forcefield_overrides_apply_arg(self, oplsaa):
        """Test that the combining rule specified in a Forcefield object
        overrides the value of the combining_rule argument passed to .apply()"""
        import mbuild as mb

        benzene = mb.load("c1ccccc1", smiles=True)

        out = oplsaa.apply(structure=benzene)

        assert out.combining_rule == "geometric"

        for adj in out.adjusts:
            if adj.atom1.name == "C" and adj.atom2.name == "H":
                break

        found_14_sigma = adj.type.sigma

        expected_14_sigma = (adj.atom1.sigma * adj.atom2.sigma) ** 0.5

        assert abs(found_14_sigma - expected_14_sigma) < 1e-10

    @pytest.mark.parametrize(
        ("ff_name", "expected_combining_rule", "expected_14_sigma"),
        [
            ("benzene_lb.xml", "lorentz", 3.0),
            ("benzene_geo.xml", "geometric", 2.82842712474619),
        ],
    )
    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    def test_combining_rule(self, ff_name, expected_combining_rule, expected_14_sigma):
        import mbuild as mb

        ff = Forcefield(get_fn(ff_name))
        benzene = mb.load("c1ccccc1", smiles=True)

        out = ff.apply(benzene)

        assert out.combining_rule == expected_combining_rule

        for adj in out.adjusts:
            if adj.atom1.name == "C" and adj.atom2.name == "H":
                break

        found_14_sigma = adj.type.sigma

        assert abs(found_14_sigma - expected_14_sigma) < 1e-10
