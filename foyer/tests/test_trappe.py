import glob
import itertools as it
import os

import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.tests.utils import atomtype

TRAPPE_UA = Forcefield(name='trappe-ua')

TRAPPE_TESTFILES_DIR = resource_filename('foyer', 'trappe_validation')

class TestTraPPE(object):

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    mol2_files = glob.glob(os.path.join(TRAPPE_TESTFILES_DIR, '*/*.mol2'))

    # Please update this file if you implement atom typing for a test case.
    # You can automatically update the files by running the below function
    # `find_correctly_implemented`.
    implemented_tests_path = os.path.join(os.path.dirname(__file__),
                                          'implemented_trappe_tests.txt')
    with open(implemented_tests_path) as f:
        correctly_implemented = [line.strip() for line in f]

    def find_correctly_implemented(self):
        with open(self.implemented_tests_path, 'a') as fh:
            for mol_path in it.chain(self.mol2_files):
                _, mol_file = os.path.split(mol_path)
                mol_name, ext = os.path.splitext(mol_file)
                try:
                   self.test_atomtyping(mol_name)
                except Exception as e:
                    print(e)
                    continue
                else:
                    if mol_name not in self.correctly_implemented:
                        fh.write('{}\n'.format(mol_name))

    @pytest.mark.parametrize('mol_name', correctly_implemented)
    def test_atomtyping(self, mol_name, testfiles_dir=TRAPPE_TESTFILES_DIR):
        files = glob.glob(os.path.join(testfiles_dir, mol_name, '*'))
        for mol_file in files:
            _, ext = os.path.splitext(mol_file)
            mol2_path = os.path.join(testfiles_dir, mol_name, mol_file)
            structure = pmd.load_file(mol2_path, structure=True)
        atomtype(structure, TRAPPE_UA, non_atomistic=True)

if __name__ == '__main__':
    TestTraPPE().find_correctly_implemented()
