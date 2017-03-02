from foyer.AST2NetworkX import AST2NetworkX
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


def test_AST2NX():
    # nodes test
    for cases in TEST_BANK:
        ast = PARSER.parse(cases)
        atoms = ast.select('atom')
        A2N = AST2NetworkX(cases)
        for atom in A2N.NetworkX.node:
            assert atom[0] in atoms


def test_Graph_Isomorphism():
    # test for graph isomorphism
    assert AST2NetworkX(TEST_BANK[0]) == AST2NetworkX(TEST_BANK[0])
    assert not AST2NetworkX(TEST_BANK[0]) == AST2NetworkX(TEST_BANK[1])
    assert AST2NetworkX(TEST_BANK[3]) == AST2NetworkX(TEST_BANK[3])
    assert not AST2NetworkX(TEST_BANK[3]) == AST2NetworkX(TEST_BANK[4])
    assert not AST2NetworkX(TEST_BANK[3]) == AST2NetworkX(TEST_BANK[5])
    

def test_Subgraph_Isomorphism():
    # test for graph isomorphism
    assert AST2NetworkX(TEST_BANK[0]) in AST2NetworkX(TEST_BANK[0])
    assert not AST2NetworkX(TEST_BANK[0]) in AST2NetworkX(TEST_BANK[1])
    assert AST2NetworkX(TEST_BANK[3]) in AST2NetworkX(TEST_BANK[3])
    assert AST2NetworkX(TEST_BANK[3]) in AST2NetworkX(TEST_BANK[5])
    assert not AST2NetworkX(TEST_BANK[5]) in AST2NetworkX(TEST_BANK[3])
    assert not AST2NetworkX(TEST_BANK[3]) in AST2NetworkX(TEST_BANK[4])
    assert AST2NetworkX(TEST_BANK[4]) in AST2NetworkX(TEST_BANK[5])
    assert AST2NetworkX(TEST_BANK[3]) in AST2NetworkX(TEST_BANK[5])


test_AST2NX()
test_Graph_Isomorphism()
test_Subgraph_Isomorphism()
