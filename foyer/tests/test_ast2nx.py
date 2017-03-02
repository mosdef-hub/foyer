from foyer.ast2nx import Ast2Nx
from foyer.smarts import SMARTS

PARSER = SMARTS()
TEST_BANK = [
    'O([H&X1])(H)',  # 0
    '[O;X2]([C;X4](F)(*)(*))[C;X4]',  # 1
    '[#6][#1](C)H',  # 2
    '[#6][#6][#6][#6][#6][#6]',  # 3
    '[#6]1[#6][#6][#6][#6][#6]1',  # 4
    '[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2'  # 5
]


def test_Ast2Nx():
    # nodes test
    for cases in TEST_BANK:
        ast = PARSER.parse(cases)
        atoms = ast.select('atom')
        A2N = Ast2Nx(cases)
        for atom in A2N.NetworkX.node:
            assert atom[0] in atoms


def test_graph_isomorphism():
    # test for graph isomorphism
    assert Ast2Nx(TEST_BANK[0]) == Ast2Nx(TEST_BANK[0])
    assert not Ast2Nx(TEST_BANK[0]) == Ast2Nx(TEST_BANK[1])
    assert Ast2Nx(TEST_BANK[3]) == Ast2Nx(TEST_BANK[3])
    assert not Ast2Nx(TEST_BANK[3]) == Ast2Nx(TEST_BANK[4])
    assert not Ast2Nx(TEST_BANK[3]) == Ast2Nx(TEST_BANK[5])
    

def test_subgraph_isomorphism():
    # test for graph isomorphism
    assert Ast2Nx(TEST_BANK[0]) in Ast2Nx(TEST_BANK[0])
    assert not Ast2Nx(TEST_BANK[0]) in Ast2Nx(TEST_BANK[1])
    assert Ast2Nx(TEST_BANK[3]) in Ast2Nx(TEST_BANK[3])
    assert Ast2Nx(TEST_BANK[3]) in Ast2Nx(TEST_BANK[5])
    assert not Ast2Nx(TEST_BANK[5]) in Ast2Nx(TEST_BANK[3])
    assert not Ast2Nx(TEST_BANK[3]) in Ast2Nx(TEST_BANK[4])
    assert Ast2Nx(TEST_BANK[4]) in Ast2Nx(TEST_BANK[5])
    assert Ast2Nx(TEST_BANK[3]) in Ast2Nx(TEST_BANK[5])
