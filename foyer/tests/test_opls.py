import glob
import itertools as it
import os

import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.exceptions import MissingForceError
from foyer.tests.utils import atomtype

OPLSAA = Forcefield(name='oplsaa')

OPLS_TESTFILES_DIR = resource_filename('foyer', 'opls_validation')

class TestOPLS(object):

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    top_files = glob.glob(os.path.join(OPLS_TESTFILES_DIR, '*/*.top'))
    mol2_files = glob.glob(os.path.join(OPLS_TESTFILES_DIR, '*/*.mol2'))

    # Please update this file if you implement atom typing for a test case.
    # You can automatically update the files by running the below function
    # `find_correctly_implemented`.
    implemented_tests_path = os.path.join(os.path.dirname(__file__),
                                          'implemented_opls_tests.txt')
    with open(implemented_tests_path) as f:
        correctly_implemented = [line.strip() for line in f]

    def find_correctly_implemented(self):
        with open(self.implemented_tests_path, 'a') as fh:
            for mol_path in it.chain(self.top_files, self.mol2_files):
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
    def test_atomtyping(self, mol_name, testfiles_dir=OPLS_TESTFILES_DIR):
        files = glob.glob(os.path.join(testfiles_dir, mol_name, '*'))
        for mol_file in files:
            _, ext = os.path.splitext(mol_file)
            if ext == '.top':
                top_filename = '{}.top'.format(mol_name)
                gro_filename = '{}.gro'.format(mol_name)
                top_path = os.path.join(testfiles_dir, mol_name, top_filename)
                gro_path = os.path.join(testfiles_dir, mol_name, gro_filename)
                structure = pmd.load_file(top_path, xyz=gro_path, parametrize=False)
            elif ext == '.mol2':
                mol2_path = os.path.join(testfiles_dir, mol_name, mol_file)
                structure = pmd.load_file(mol2_path, structure=True)
        atomtype(structure, OPLSAA)

    def test_full_parametrization(self):
        top = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.top')
        gro = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.gro')
        structure = pmd.load_file(top, xyz=gro)
        parametrized = OPLSAA.apply(structure)

        assert sum((1 for at in parametrized.atoms if at.type == 'opls_145')) == 6
        assert sum((1 for at in parametrized.atoms if at.type == 'opls_146')) == 6
        assert len(parametrized.bonds) == 12
        assert all(x.type for x in parametrized.bonds)
        assert len(parametrized.angles) == 18
        assert all(x.type for x in parametrized.angles)
        assert len(parametrized.rb_torsions) == 24
        assert all(x.type for x in parametrized.dihedrals)
        assert parametrized.combining_rule == 'geometric'

    def test_get_parameters_atoms(self):
        atom_params = OPLSAA.get_parameters("atoms", "opls_145")
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_get_parameters_bonds(self):
        bond_params = OPLSAA.get_parameters("bonds", ["opls_760", "opls_145"])
        assert bond_params["length"] == 0.146
        assert bond_params["k"] == 334720.0

    def test_get_parameters_angle(self):
        angle_params = OPLSAA.get_parameters(
            "angles", ["opls_166", "opls_772", "opls_167"]
        )
        assert angle_params["theta"] == 2.09439510239
        assert angle_params["k"] == 585.76

    def test_get_parameters_rb_proper(self):
        proper_params = OPLSAA.get_parameters(
            "rb_propers", ["opls_215", "opls_215", "opls_235", "opls_269"]
        )

        assert proper_params["c0"] == 2.28446
        assert proper_params["c1"] == 0.0
        assert proper_params["c2"] == -2.28446
        assert proper_params["c1"] == 0.0
        assert proper_params["c4"] == 0.0
        assert proper_params["c5"] == 0.0

    def test_get_parameters_wildcard(self):
        proper_params = OPLSAA.get_parameters(
            "rb_propers", ["", "opls_235", "opls_544", ""]
        )

        assert proper_params["c0"] == 30.334
        assert proper_params["c1"] == 0.0
        assert proper_params["c2"] == -30.334
        assert proper_params["c1"] == 0.0
        assert proper_params["c4"] == 0.0
        assert proper_params["c5"] == 0.0

    def test_opls_missing_force(self):
        with pytest.raises(MissingForceError):
            OPLSAA.get_parameters('periodic_propers', key=['a', 'b', 'c', 'd'])


if __name__ == '__main__':
    TestOPLS().find_correctly_implemented()
