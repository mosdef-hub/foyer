import parmed as pmd
import lark
import pytest

from foyer.exceptions import FoyerError
from foyer.forcefield import generate_topology, Forcefield
from foyer.smarts_graph import SMARTSGraph
from foyer.smarts import SMARTS
from foyer.tests.utils import get_fn
from foyer.tests.base_test import BaseTest


PARSER = SMARTS()
class TestSMARTS(BaseTest):
    
    
    def test_ast(self):
        ast = PARSER.parse('O([H&X1])(H)')
        assert ast.data == "start"
        assert ast.children[0].data == "atom"
        assert ast.children[0].children[0].data == "atom_symbol"
        assert str(ast.children[0].children[0].children[0]) == "O"
    
    
    def test_parse(self):
        smarts = ['[#6][#1](C)H',
                  '[O;X2]([C;X4](F)(*)(*))[C;X4]']
        for pattern in smarts:
            ast = PARSER.parse(pattern)
    
    
    def test_uniqueness(self):
        mol2 = pmd.load_file(get_fn('uniqueness_test.mol2'), structure=True)
        top, _ = generate_topology(mol2)
    
        _rule_match(top, '[#6]1[#6][#6][#6][#6][#6]1', False)
        _rule_match(top, '[#6]1[#6][#6][#6][#6]1', False)
        _rule_match(top, '[#6]1[#6][#6][#6]1', True)
    
    
    def test_ringness(self):
        ring = pmd.load_file(get_fn('ring.mol2'), structure=True)
        top, _ = generate_topology(ring)
        _rule_match(top, '[#6]1[#6][#6][#6][#6][#6]1', True)
    
        not_ring = pmd.load_file(get_fn('not_ring.mol2'), structure=True)
        top, _ = generate_topology(not_ring)
        _rule_match(top, '[#6]1[#6][#6][#6][#6][#6]1', False)
    
    
    def test_fused_ring(self):
        fused = pmd.load_file(get_fn('fused.mol2'), structure=True)
        top, _ = generate_topology(fused)
        rule = SMARTSGraph(name='test', parser=PARSER,
                           smarts_string='[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2')
    
        match_indices = list(rule.find_matches(top))
        assert 3 in match_indices
        assert 4 in match_indices
        assert len(match_indices) == 2
    
    
    def test_ring_count(self):
        # Two rings
        fused = pmd.load_file(get_fn('fused.mol2'), structure=True)
        top, _ = generate_topology(fused)
        rule = SMARTSGraph(name='test', parser=PARSER,
                           smarts_string='[#6;R2]')
    
        match_indices = list(rule.find_matches(top))
        for atom_idx in (3, 4):
            assert atom_idx in match_indices
        assert len(match_indices) == 2
    
        rule = SMARTSGraph(name='test', parser=PARSER,
                           smarts_string='[#6;R1]')
        match_indices = list(rule.find_matches(top))
        for atom_idx in (0, 1, 2, 5, 6, 7, 8, 9):
            assert atom_idx in match_indices
        assert len(match_indices) == 8
    
        # One ring
        ring = pmd.load_file(get_fn('ring.mol2'), structure=True)
        top, _ = generate_topology(ring)
    
        rule = SMARTSGraph(name='test', parser=PARSER,
                           smarts_string='[#6;R1]')
        match_indices = list(rule.find_matches(top))
        for atom_idx in range(6):
            assert atom_idx in match_indices
        assert len(match_indices) == 6
    
    
    def test_precedence_ast(self):
        ast1 = PARSER.parse('[C,H;O]')
        ast2 = PARSER.parse('[O;H,C]')
        assert ast1.children[0].children[0].data == 'weak_and_expression'
        assert ast2.children[0].children[0].data == 'weak_and_expression'
    
        assert ast1.children[0].children[0].children[0].data == 'or_expression'
        assert ast2.children[0].children[0].children[1].data == 'or_expression'
    
        ast1 = PARSER.parse('[C,H&O]')
        ast2 = PARSER.parse('[O&H,C]')
        assert ast1.children[0].children[0].data == 'or_expression'
        assert ast2.children[0].children[0].data == 'or_expression'
    
        assert ast1.children[0].children[0].children[1].data == 'and_expression'
        assert ast2.children[0].children[0].children[0].data == 'and_expression'
    
    
    def test_precedence(self):
        mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
        top, _ = generate_topology(mol2)
    
        checks = {'[C,O;C]': 2,
                  '[C&O;C]': 0,
                  '[!C;O,C]': 0,
                  '[!C&O,C]': 2,
                  }
    
        for smart, result in checks.items():
            _rule_match_count(top, smart, result)
    
    
    def test_not_ast(self):
        checks = {'[!C;!H]': 'weak_and_expression',
                  '[!C&H]': 'and_expression',
                  '[!C;H]': 'weak_and_expression',
                  '[!C]': 'not_expression'}
    
        for smart, grandchild in checks.items():
            ast = PARSER.parse(smart)
            assert ast.children[0].children[0].data == grandchild
    
        illegal_nots = ['[!CH]', '[!C!H]']
        for smart in illegal_nots:
            with pytest.raises(lark.UnexpectedInput):
                PARSER.parse(smart)
    
    
    def test_not(self):
        mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
        top, _ = generate_topology(mol2)
    
        checks = {'[!O]': 8,
                  '[!#5]': 8,
                  '[!C]': 6,
                  '[!#6]': 6,
                  '[!C&!H]': 0,
                  '[!C;!H]': 0,
                  }
        for smart, result in checks.items():
            _rule_match_count(top, smart, result)
    
    
    def test_hexa_coordinated(self):
        ff = Forcefield(forcefield_files=get_fn('pf6.xml'))
        mol2 = pmd.load_file(get_fn('pf6.mol2'), structure=True)
    
        pf6 = ff.apply(mol2)
    
        types = [a.type for a in pf6.atoms]
        assert types.count('P') == 1
        assert types.count('F1') == 2
        assert types.count('F2') == 2
        assert types.count('F3') == 2
    
        assert len(pf6.bonds) == 6
        assert all(bond.type for bond in pf6.bonds)
    
        assert len(pf6.angles) == 15
        assert all(angle.type for angle in pf6.angles)
    
    
    def test_optional_names_bad_syntax(self):
        bad_optional_names = ['_C', 'XXX', 'C']
        with pytest.raises(FoyerError):
            S = SMARTS(optional_names=bad_optional_names)
    
    
    def test_optional_names_good_syntax(self):
        good_optional_names = ['_C', '_CH2', '_CH']
        S = SMARTS(optional_names=good_optional_names)
    
    
    def test_optional_name_parser(self):
        optional_names = ['_C', '_CH2', '_CH']
        S = SMARTS(optional_names=optional_names)
        ast = S.parse('_CH2_C_CH')
        symbols = [a.children[0] for a in ast.find_data('atom_symbol')]
        for name in optional_names:
            assert name in symbols

def _rule_match(top, smart, result):
    rule = SMARTSGraph(name='test', parser=PARSER, smarts_string=smart)
    assert bool(list(rule.find_matches(top))) is result


def _rule_match_count(top, smart, count):
    rule = SMARTSGraph(name='test', parser=PARSER, smarts_string=smart)
    assert len(list(rule.find_matches(top))) is count
