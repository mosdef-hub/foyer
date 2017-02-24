import parmed as pmd

from foyer.forcefield import generate_topology
from foyer.rule import Rule
from foyer.smarts import SMARTS
from foyer.tests.utils import get_fn


PARSER = SMARTS()

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
    top = generate_topology(mol2)

    atom1 = next(top.atoms())
    rule = Rule('test', parser=PARSER, smarts_string='[#6]1[#6][#6][#6][#6][#6]1')
    assert rule.matches(atom1) is False
    rule = Rule('test', parser=PARSER, smarts_string='[#6]1[#6][#6][#6][#6]1')
    assert rule.matches(atom1) is False
    rule = Rule('test', parser=PARSER, smarts_string='[#6]1[#6][#6][#6]1')
    assert rule.matches(atom1) is True


def test_ringness():
    ring = pmd.load_file(get_fn('ring.mol2'), structure=True)
    top = generate_topology(ring)
    atom1 = next(top.atoms())
    rule = Rule('test', parser=PARSER, smarts_string='[#6]1[#6][#6][#6][#6][#6]1')
    assert rule.matches(atom1) is True

    not_ring = pmd.load_file(get_fn('not_ring.mol2'), structure=True)
    top = generate_topology(not_ring)
    atom1 = next(top.atoms())
    rule = Rule('test', parser=PARSER, smarts_string='[#6]1[#6][#6][#6][#6][#6]1')
    assert rule.matches(atom1) is False


def test_fused_ring():
    fused = pmd.load_file(get_fn('fused.mol2'), structure=True)
    top = generate_topology(fused)
    atoms = list(top.atoms())
    rule = Rule('test', parser=PARSER,
                smarts_string='[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2')

    assert rule.matches(atoms[2]) is False
    for _ in range(10):  # Traversal order during matching is stochastic.
        assert rule.matches(atoms[3]) is True
        assert rule.matches(atoms[4]) is True
    assert rule.matches(atoms[5]) is False


def test_not():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    top = generate_topology(mol2)
    atom1 = next(top.atoms())

    rule = Rule('test', parser=PARSER, smarts_string='[!O]')
    assert rule.matches(atom1) is True
    rule = Rule('test', parser=PARSER, smarts_string='[!#5]')
    assert rule.matches(atom1) is True

    rule = Rule('test', parser=PARSER, smarts_string='[!C]')
    assert rule.matches(atom1) is False
    rule = Rule('test', parser=PARSER, smarts_string='[!#6]')
    assert rule.matches(atom1) is False

