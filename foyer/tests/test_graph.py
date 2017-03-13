from foyer.smarts_graph import SMARTSGraph


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
        for atom in atoms:
            assert id(atom) in graph.nodes()
