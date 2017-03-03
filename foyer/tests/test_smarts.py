import parmed as pmd
import plyplus
import pytest

from foyer.forcefield import generate_topology, Forcefield
from foyer.rule import Rule
from foyer.smarts import SMARTS
from foyer.tests.utils import get_fn


PARSER = SMARTS()


def _rule_match(atom, smart, result):
    rule = Rule('test', parser=PARSER, smarts_string=smart)
    assert rule.matches(atom) is result


def test_ast():
    ast = PARSER.parse('O([H&X1])(H)')
    assert ast.head == "start"
    assert ast.tail[0].head == "atom"
    assert ast.tail[0].tail[0].head == "atom_symbol"
    assert ast.tail[0].tail[0].head == "atom_symbol"
    assert str(ast.tail[0].tail[0].tail[0]) == "O"


def test_parse():
    smarts = ['[#6][#1](C)H',
              '[O;X2]([C;X4](F)(*)(*))[C;X4]']
    for pattern in smarts:
        ast = PARSER.parse(pattern)


def test_uniqueness():
    mol2 = pmd.load_file(get_fn('uniqueness_test.mol2'), structure=True)
    top, _ = generate_topology(mol2)

    atom1 = next(top.atoms())
    _rule_match(atom1, '[#6]1[#6][#6][#6][#6][#6]1', False)
    _rule_match(atom1, '[#6]1[#6][#6][#6][#6]1', False)
    _rule_match(atom1, '[#6]1[#6][#6][#6]1', True)


def test_ringness():
    ring = pmd.load_file(get_fn('ring.mol2'), structure=True)
    top, _ = generate_topology(ring)
    atom1 = next(top.atoms())
    _rule_match(atom1, '[#6]1[#6][#6][#6][#6][#6]1', True)

    not_ring = pmd.load_file(get_fn('not_ring.mol2'), structure=True)
    top, _ = generate_topology(not_ring)
    atom1 = next(top.atoms())
    _rule_match(atom1, '[#6]1[#6][#6][#6][#6][#6]1', False)


def test_fused_ring():
    fused = pmd.load_file(get_fn('fused.mol2'), structure=True)
    top, _ = generate_topology(fused)
    atoms = list(top.atoms())
    rule = Rule('test', parser=PARSER,
                smarts_string='[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2')

    assert rule.matches(atoms[2]) is False
    for _ in range(10):  # Traversal order during matching is stochastic.
        assert rule.matches(atoms[3]) is True
        assert rule.matches(atoms[4]) is True
    assert rule.matches(atoms[5]) is False


def test_precedence_ast():
    ast1 = PARSER.parse('[C,H;O]')
    ast2 = PARSER.parse('[O;H,C]')
    assert ast1.tail[0].tail[0].head == 'weak_and_expression'
    assert ast2.tail[0].tail[0].head == 'weak_and_expression'

    assert ast1.tail[0].tail[0].tail[0].head == 'or_expression'
    assert ast2.tail[0].tail[0].tail[1].head == 'or_expression'

    ast1 = PARSER.parse('[C,H&O]')
    ast2 = PARSER.parse('[O&H,C]')
    assert ast1.tail[0].tail[0].head == 'or_expression'
    assert ast2.tail[0].tail[0].head == 'or_expression'

    assert ast1.tail[0].tail[0].tail[1].head == 'and_expression'
    assert ast2.tail[0].tail[0].tail[0].head == 'and_expression'


def test_precedence():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    top, _ = generate_topology(mol2)
    atom1 = next(top.atoms())

    checks = {'[C,H;C]': True,
              '[C&H;C]': False,
              '[!C;H,C]': False,
              '[!C&H,C]': True,
              }

    for smart, result in checks.items():
        _rule_match(atom1, smart, result)


def test_not_ast():
    checks = {'[!C;!H]': 'weak_and_expression',
              '[!C&H]': 'and_expression',
              '[!C;H]': 'weak_and_expression',
              '[!C]': 'not_expression'}

    for smart, tail2head in checks.items():
        ast = PARSER.parse(smart)
        assert ast.tail[0].tail[0].head == tail2head

    illegal_nots = ['[!CH]', '[!C!H]']
    for smart in illegal_nots:
        with pytest.raises(plyplus.ParseError):
            PARSER.parse(smart)


def test_not():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    top, _ = generate_topology(mol2)
    atom1 = next(top.atoms())

    checks = {'[!O]': True,
              '[!#5]': True,
              '[!C]': False,
              '[!#6]': False,
              }
    for smart, result in checks.items():
        _rule_match(atom1, smart, result)


def test_hexa_coordinated():
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

