import glob
import os

import numpy as np
import parmed as pmd
import sys
from pkg_resources import resource_filename
import pytest

from foyer.forcefield import Forcefield
from foyer.atomtyper import find_atomtypes
from foyer.tests.topology import Topology


class TestOPLS(object):

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    # Please update this file if you implement atom typing for a test case.
    # You can automatically update the files by running the below function
    # `find_correctly_implemented`.
    implemented_tests_path = os.path.join(os.path.dirname(__file__),
                                          'implemented_opls_tests.txt')

    correctly_implemented = [line.split() for line in open(implemented_tests_path)]
    correctly_implemented_mol_names = {x[0] for x in correctly_implemented}
    correctly_implemented_top_files = {x[1] for x in correctly_implemented}

    def find_correctly_implemented(self):
        with open(self.implemented_tests_path, 'a') as fh:
            for name, (top_path, gro_path) in Topology.available().items():
                try:
                    mol_name = self.test_atomtyping(name)
                except:
                    sys.exc_info()[0]
                    continue
                else:
                    print(' *** SUCCESS {} {}\n'.format(mol_name, basename))
                    basename = os.path.basename(top_path)
                    if basename not in self.correctly_implemented_top_files:

                        fh.write('{} {}\n'.format(mol_name, basename))

    @pytest.mark.parametrize('name', correctly_implemented_top_files)
    def test_atomtyping(self, name):

        print("loading {}".format(name))
        structure = Topology.by_name(name)
        print(structure)

        known_opls_types = [atom.type for atom in structure.atoms]

        print("Typing {} ({})...".format(structure.title, name))

        ff = Forcefield.by_name('oplsaa')
        ff.apply(structure, debug=False)

        generated_opls_types = list()
        for i, atom in enumerate(structure.atoms):
            message = ('Found multiple or no OPLS types for atom {} in {} ({}): {}\n'
                       'Should be atomtype: {}'.format(
                i, structure.title, name, atom.type, known_opls_types[i]))
            assert atom.type, message

            generated_opls_types.append(atom.type)
        both = zip(generated_opls_types, known_opls_types)

        n_types = np.array(range(len(generated_opls_types)))
        known_opls_types = np.array(known_opls_types)
        generated_opls_types = np.array(generated_opls_types)

        non_matches = np.array([a != b for a, b in both])
        message = "Found inconsistent OPLS types in {} ({}): {}".format(
            structure.title, name,
            list(zip(n_types[non_matches],
                     generated_opls_types[non_matches],
                     known_opls_types[non_matches])))
        assert not non_matches.any(), message
        return structure.title

    def test_full_parametrization(self):
        top = os.path.join(self.resource_dir, 'benzene.top')
        gro = os.path.join(self.resource_dir, 'benzene.gro')
        ff = os.path.join(self.resource_dir, 'ff.ff/forcefield.itp')
        structure = pmd.load_file(top, xyz=gro)
        parametrized = apply_forcefield(structure, forcefield=ff, debug=False)

        assert sum((1 for at in parametrized.atoms if at.type == 'opls_145')) == 6
        assert sum((1 for at in parametrized.atoms if at.type == 'opls_146')) == 6
        assert len(parametrized.bonds) == 12
        assert all(x.type for x in parametrized.bonds)
        assert len(parametrized.angles) == 18
        assert all(x.type for x in parametrized.angles)
        assert len(parametrized.rb_torsions) == 24
        assert all(x.type for x in parametrized.dihedrals)


if __name__ == "__main__":
    test_class = TestOPLS()

    # mol = 'acetophenone'
    # mol = 'benzyl-alcohol'
    # top_path = test_class.find_topfile_by_mol_name(mol)
    # test_class.test_atomtyping(top_path)

    test_class.find_correctly_implemented()

