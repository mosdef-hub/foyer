import os

import parmed as pmd
from pkg_resources import resource_filename

from foyer.forcefield import generate_topology
from foyer.rule import Rule


def test_uniqueness():
    resource_dir = resource_filename('foyer', '../opls_validation')

    mol2 = pmd.load_file(os.path.join(resource_dir, 'uniqueness_test.mol2'),
                         structure=True)
    top = generate_topology(mol2)

    atom1 = next(top.atoms())
    rule = Rule('test', '[#6]1[#6][#6][#6][#6][#6]1')
    result = rule.matches(atom1)
    assert not result
