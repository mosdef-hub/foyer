import glob
import os

from six import string_types
from pkg_resources import resource_filename
import pytest

from foyer.utils.io import load_top_opls
from foyer.forcefield import prepare_atoms
from foyer.atomtyper import find_atomtypes


class TestOPLS():

    def test_atomtyping(self, only_run=None):
        resource_dir = resource_filename('foyer', '../opls_validation')
        top_files = glob.glob(os.path.join(resource_dir, '*.top'))

        # Please update this file if you implement atom typing for a test case.
        implemented_tests_path = os.path.join(os.path.dirname(__file__), 'implemented_opls_tests.txt')
        correctly_implemented = [line.strip() for line in open(implemented_tests_path)]

        for top in top_files:
            top_name = os.path.split(top)[-1]
            system, known_opls_types, mol_name = load_top_opls(top, only_run)

            if only_run and only_run != mol_name:
                continue
            elif mol_name not in correctly_implemented:
                continue

            print("Typing {} ({})...".format(mol_name, top_name))
            prepare_atoms(system)
            find_atomtypes(list(system.atoms), forcefield='OPLS-AA', debug=False)

            generated_opls_types = list()

            for i, atom in enumerate(system.atoms):
                message = ('Found multiple or no OPLS types for atom {} in {} ({}): {}\n'
                           'Should be atomtype: {}'.format(
                    i, mol_name, top_name, atom.atomtype[0], known_opls_types[i]))
                assert isinstance(atom.atomtype[0], string_types), message
                generated_opls_types.append(atom.atomtype[0])
            both = zip(generated_opls_types, known_opls_types)
            message = "Found inconsistent OPLS types in {} ({}): {}".format(
                mol_name, top_name, list(zip(range(len(generated_opls_types)),
                                        generated_opls_types,
                                        known_opls_types)))

            assert all([a == b for a, b in both]), message
            print("Passed.\n")

    @pytest.mark.skipif(True, reason='Not implemented yet')
    def test_full_parameterization(self):
        #ethane.save('ethane.gro', forcefield='opls-aa')
        pass


if __name__ == "__main__":
    # TestOPLS().test_atomtyping('benzene')
    TestOPLS().test_atomtyping()

    #from mbuild.examples.ethane.ethane import Ethane
    #TestOPLS().test_full_parameterization(Ethane())
