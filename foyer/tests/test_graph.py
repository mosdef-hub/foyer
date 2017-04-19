import parmed as pmd

from foyer.forcefield import generate_topology
from foyer.smarts_graph import SMARTSGraph
from foyer.tests.utils import get_fn


TEST_BANK = [
    'O([H&X1])(H)',
    '[O;X2]([C;X4](F)(*)(*))[C;X4]',
    '[#6][#1](C)H',
    '[#6][#6][#6][#6][#6][#6]',
    '[#6]1[#6][#6][#6][#6][#6]1',
    '[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2',
]


def test_init():
    """Initialization test. """
    for smarts in TEST_BANK:
        graph = SMARTSGraph(smarts)
        atoms = graph.ast.select('atom')
        for n, atom in enumerate(atoms):
            assert n in graph.nodes()


def test_lazy_cycle_finding():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    top, _ = generate_topology(mol2)

    rule = SMARTSGraph(smarts_string='[C]')
    list(rule.find_matches(top))
    assert not any([hasattr(a, 'cycles') for a in top.atoms()])

    ring_tokens = ['R1', 'r6']
    for token in ring_tokens:
        rule = SMARTSGraph(smarts_string='[C;{}]'.format(token))
        list(rule.find_matches(top))
        assert all([hasattr(a, 'cycles') for a in top.atoms()])
