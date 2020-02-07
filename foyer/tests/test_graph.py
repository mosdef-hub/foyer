import parmed as pmd

from foyer.smarts_graph import SMARTSGraph, _prepare_atoms
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
        atoms = graph.ast.find_data('atom')
        for n, atom in enumerate(atoms):
            assert n in graph.nodes()


def test_lazy_cycle_finding():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    typemap = {atom.idx: {'whitelist': set(), 'blacklist': set(), 
                            'atomtype': None}     
                            for atom in mol2.atoms}

    rule = SMARTSGraph(smarts_string='[C]', typemap=typemap)
    list(rule.find_matches(mol2, typemap))
    assert not any(['cycles' in typemap[a.idx] for a in mol2.atoms])

    ring_tokens = ['R1', 'r6']
    for token in ring_tokens:
        rule = SMARTSGraph(smarts_string='[C;{}]'.format(token),
                            typemap=typemap)
        list(rule.find_matches(mol2, typemap))
        assert all(['cycles' in typemap[a.idx] for a in mol2.atoms])


def test_cycle_finding_multiple():
    mol2 = pmd.load_file(get_fn('fullerene.pdb'), structure=True)
    typemap = {atom.idx: {'whitelist': set(), 'blacklist': set(), 
                            'atomtype': None}     
                            for atom in mol2.atoms}

    _prepare_atoms(mol2, typemap, compute_cycles=True)
    cycle_lengths = [list(map(len, typemap[atom.idx]['cycles'])) 
                        for atom in mol2.atoms]
    expected = [5, 6, 6]
    assert all(sorted(lengths) == expected for lengths in cycle_lengths)
