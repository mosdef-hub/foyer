import lark
import parmed as pmd
import pytest

from foyer.exceptions import FoyerError
from foyer.forcefield import Forcefield
from foyer.smarts import SMARTS
from foyer.smarts_graph import SMARTSGraph
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn
from foyer.topology_graph import TopologyGraph


class TestSMARTS(BaseTest):
    @pytest.fixture(scope="session")
    def rule_match(self, smarts_parser):
        def _rule_match(top, typemap, smart, result):
            rule = SMARTSGraph(
                name="test",
                parser=smarts_parser,
                smarts_string=smart,
                typemap=typemap,
            )
            assert bool(list(rule.find_matches(top, typemap))) is result

        return _rule_match

    @pytest.fixture(scope="session")
    def rule_match_count(self, smarts_parser):
        def _rule_match_count(top, typemap, smart, count):
            rule = SMARTSGraph(
                name="test",
                parser=smarts_parser,
                smarts_string=smart,
                typemap=typemap,
            )
            assert len(list(rule.find_matches(top, typemap))) is count

        return _rule_match_count

    def test_ast(self, smarts_parser):
        ast = smarts_parser.parse("O([H&X1])(H)")
        assert ast.data == "start"
        assert ast.children[0].data == "atom"
        assert ast.children[0].children[0].data == "atom_symbol"
        assert str(ast.children[0].children[0].children[0]) == "O"

    @pytest.mark.parametrize(
        "pattern", ["[#6][#1](C)H", "[O;X2]([C;X4](F)(*)(*))[C;X4]"]
    )
    def test_parse(self, pattern, smarts_parser):
        assert smarts_parser.parse(pattern)

    def test_uniqueness(self, rule_match):
        mol2 = pmd.load_file(get_fn("uniqueness_test.mol2"), structure=True)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in mol2.atoms
        }

        mol2_graph = TopologyGraph.from_parmed(mol2)

        rule_match(mol2_graph, typemap, "[#6]1[#6][#6][#6][#6][#6]1", False)
        rule_match(mol2_graph, typemap, "[#6]1[#6][#6][#6][#6]1", False)
        rule_match(mol2_graph, typemap, "[#6]1[#6][#6][#6]1", True)

    def test_ringness(self, rule_match):
        ring_mol2 = pmd.load_file(get_fn("ring.mol2"), structure=True)
        ring_mol2_graph = TopologyGraph.from_parmed(ring_mol2)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in ring_mol2.atoms
        }

        rule_match(ring_mol2_graph, typemap, "[#6]1[#6][#6][#6][#6][#6]1", True)

        not_ring_mol2 = pmd.load_file(get_fn("not_ring.mol2"), structure=True)
        not_ring_mol2_graph = TopologyGraph.from_parmed(not_ring_mol2)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in not_ring_mol2.atoms
        }

        rule_match(not_ring_mol2_graph, typemap, "[#6]1[#6][#6][#6][#6][#6]1", False)

    def test_fused_ring(self, smarts_parser):
        mol2 = pmd.load_file(get_fn("fused.mol2"), structure=True)
        mol2_graph = TopologyGraph.from_parmed(mol2)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in mol2.atoms
        }

        rule = SMARTSGraph(
            name="test",
            parser=smarts_parser,
            smarts_string="[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2",
            typemap=typemap,
        )

        match_indices = list(rule.find_matches(mol2_graph, typemap))
        assert 3 in match_indices
        assert 4 in match_indices
        assert len(match_indices) == 2

    def test_ring_count(self, smarts_parser):
        # Two rings
        fused = pmd.load_file(get_fn("fused.mol2"), structure=True)
        fused_graph = TopologyGraph.from_parmed(fused)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in fused.atoms
        }
        rule = SMARTSGraph(
            name="test",
            parser=smarts_parser,
            smarts_string="[#6;R2]",
            typemap=typemap,
        )

        match_indices = list(rule.find_matches(fused_graph, typemap))
        for atom_idx in (3, 4):
            assert atom_idx in match_indices
        assert len(match_indices) == 2

        rule = SMARTSGraph(
            name="test",
            parser=smarts_parser,
            smarts_string="[#6;R1]",
            typemap=typemap,
        )
        match_indices = list(rule.find_matches(fused_graph, typemap))
        for atom_idx in (0, 1, 2, 5, 6, 7, 8, 9):
            assert atom_idx in match_indices
        assert len(match_indices) == 8

        # One ring
        ring = pmd.load_file(get_fn("ring.mol2"), structure=True)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in ring.atoms
        }

        ring_graph = TopologyGraph.from_parmed(ring)
        rule = SMARTSGraph(
            name="test",
            parser=smarts_parser,
            smarts_string="[#6;R1]",
            typemap=typemap,
        )
        match_indices = list(rule.find_matches(ring_graph, typemap))
        for atom_idx in range(6):
            assert atom_idx in match_indices
        assert len(match_indices) == 6

    def test_precedence_ast(self, smarts_parser):
        ast1 = smarts_parser.parse("[C,H;O]")
        ast2 = smarts_parser.parse("[O;H,C]")
        assert ast1.children[0].children[0].data == "weak_and_expression"
        assert ast2.children[0].children[0].data == "weak_and_expression"

        assert ast1.children[0].children[0].children[0].data == "or_expression"
        assert ast2.children[0].children[0].children[1].data == "or_expression"

        ast1 = smarts_parser.parse("[C,H&O]")
        ast2 = smarts_parser.parse("[O&H,C]")
        assert ast1.children[0].children[0].data == "or_expression"
        assert ast2.children[0].children[0].data == "or_expression"

        assert ast1.children[0].children[0].children[1].data == "and_expression"
        assert ast2.children[0].children[0].children[0].data == "and_expression"

    def test_precedence(self, rule_match_count):
        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in mol2.atoms
        }

        mol2_graph = TopologyGraph.from_parmed(mol2)

        checks = {
            "[C,O;C]": 2,
            "[C&O;C]": 0,
            "[!C;O,C]": 0,
            "[!C&O,C]": 2,
        }

        for smart, result in checks.items():
            rule_match_count(mol2_graph, typemap, smart, result)

    def test_not_ast(self, smarts_parser):
        checks = {
            "[!C;!H]": "weak_and_expression",
            "[!C&H]": "and_expression",
            "[!C;H]": "weak_and_expression",
            "[!C]": "not_expression",
        }

        for smart, grandchild in checks.items():
            ast = smarts_parser.parse(smart)
            assert ast.children[0].children[0].data == grandchild

        illegal_nots = ["[!CH]", "[!C!H]"]
        for smart in illegal_nots:
            with pytest.raises(lark.UnexpectedInput):
                smarts_parser.parse(smart)

    def test_not(self, rule_match_count):
        mol2 = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in mol2.atoms
        }
        mol2_graph = TopologyGraph.from_parmed(mol2)

        checks = {
            "[!O]": 8,
            "[!#5]": 8,
            "[!C]": 6,
            "[!#6]": 6,
            "[!C&!H]": 0,
            "[!C;!H]": 0,
        }
        for smart, result in checks.items():
            rule_match_count(mol2_graph, typemap, smart, result)

    def test_hexa_coordinated(self):
        ff = Forcefield(forcefield_files=get_fn("pf6.xml"))
        mol2 = pmd.load_file(get_fn("pf6.mol2"), structure=True)

        pf6 = ff.apply(mol2)

        types = [a.type for a in pf6.atoms]
        assert types.count("P") == 1
        assert types.count("F1") == 2
        assert types.count("F2") == 2
        assert types.count("F3") == 2

        assert len(pf6.bonds) == 6
        assert all(bond.type for bond in pf6.bonds)

        assert len(pf6.angles) == 15
        assert all(angle.type for angle in pf6.angles)

    def test_optional_names_bad_syntax(self):
        bad_optional_names = ["_C", "XXX", "C"]
        with pytest.raises(FoyerError):
            SMARTS(optional_names=bad_optional_names)

    def test_optional_names_good_syntax(self):
        good_optional_names = ["_C", "_CH2", "_CH"]
        SMARTS(optional_names=good_optional_names)

    def test_optional_name_parser(self):
        optional_names = ["_C", "_CH2", "_CH"]
        S = SMARTS(optional_names=optional_names)
        ast = S.parse("_CH2_C_CH")
        symbols = [a.children[0] for a in ast.find_data("atom_symbol")]
        for name in optional_names:
            assert name in symbols

    def test_bond_order(self, rule_match_count):
        import mbuild as mb

        smiles_string = "C=CC"  # propene
        cpd = mb.load(smiles_string, smiles=True)
        structure = cpd.to_parmed()
        structure.bonds[0].order = 2.0

        typemap = {
            atom.idx: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for atom in structure.atoms
        }

        mol2_graph = TopologyGraph.from_parmed(structure)

        checks = {
            "[C](=C)(H)(H)": 1,
            "[C](-C)=C": 1,
            "[C](-C)(H)(H)H": 1,
        }

        for smart, result in checks.items():
            rule_match_count(mol2_graph, typemap, smart, result)
