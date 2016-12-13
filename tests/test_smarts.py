import os

import parmed as pmd
from pkg_resources import resource_filename

from foyer.forcefield import generate_topology
from foyer.rule import Rule
from foyer.smarts import parse

resource_dir = resource_filename('foyer', '../opls_validation')

def test_uniqueness():
    mol2 = pmd.load_file(os.path.join(resource_dir, 'uniqueness_test.mol2'),
                         structure=True)
    top = generate_topology(mol2)

    atom1 = next(top.atoms())
    rule = Rule('test', '[#6]1[#6][#6][#6][#6][#6]1')
    result = rule.matches(atom1)
    assert not result

def test_ringness():
    ring = pmd.load_file(os.path.join(resource_dir, 'ring.mol2'),
                         structure=True)
    top = generate_topology(ring)
    atom1 = next(top.atoms())
    rule = Rule('test', '[#6]1[#6][#6][#6][#6][#6]1')
    assert rule.matches(atom1) == True

    not_ring = pmd.load_file(os.path.join(resource_dir, 'not_ring.mol2'),
                         structure=True)
    top = generate_topology(not_ring)
    atom1 = next(top.atoms())
    rule = Rule('test', '[#6]1[#6][#6][#6][#6][#6]1')
    assert rule.matches(atom1) == False

def test_fused_ring():
    fused= pmd.load_file(os.path.join(resource_dir, 'fused.mol2'),
                         structure=True)
    top = generate_topology(fused)
    atoms = list(top.atoms())
    rule = Rule('test', '[#6]12[#6][#6][#6][#6][#6]1[#6][#6][#6][#6]2')

    assert rule.matches(atoms[2]) == False
    for _ in range(10):  # Traversal order during matching is stochastic.
        assert rule.matches(atoms[3]) == True
        assert rule.matches(atoms[4]) == True
    assert rule.matches(atoms[5]) == False

def test_parse():
    smarts = "[O;X2]([C;X4](F)(*)(*))[C;X4]"
    ast = parse(smarts)