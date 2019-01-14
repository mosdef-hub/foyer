import glob
import itertools as it
import os
from scipy.constants import N_A, k, calorie
import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.tests.utils import atomtype

from numpy.testing import assert_almost_equal

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
        files = glob.glob(os.path.join(testfiles_dir, mol_name, '*.mol2'))
        for mol_file in files:
            _, ext = os.path.splitext(mol_file)
            mol2_path = os.path.join(testfiles_dir, mol_name, mol_file)
            structure = pmd.load_file(mol2_path, structure=True)
        atomtype(structure, TRAPPE_UA, non_atomistic=True)

@pytest.fixture
def trappe_ref(request):
    # load raw trappe csv file and process slightly
    container = {
        '(pseudo)atom': [],
        'stretch': [],
        'bend': [],
        'torsion': [],
    }
    # split csv file into sections
    csv_file = os.path.join(TRAPPE_TESTFILES_DIR,
                            request.param)
    with open(csv_file, 'r') as inf:
        for line in inf:
            if line.startswith('#'):
                section_clue = line.split(",")[1]
                section = container[section_clue]
                line = next(inf)
                while not line == '\n':
                    # slurp until newline hit
                    section.append(line)
                    try:
                        line = next(inf)
                    except StopIteration:
                        break

    return container


@pytest.fixture
def mol2file(request):
    # load a mol2 file and apply foyer's trappe-ua forcefield
    fname = request.param
    molfile = glob.glob(os.path.join(TRAPPE_TESTFILES_DIR, fname, '*'))[0]
    structure = pmd.load_file(molfile, structure=True)

    return TRAPPE_UA.apply(structure,
                           assert_angle_params=False,
                           assert_dihedral_params=False,
    )


def parmed2trappe(val):
    # convert kcal/mol (parmed units) to K/k_B (trappe units)
    return val * calorie / (k * N_A / 1000.)

def angle_parmed2trappe(val):
    # trappe angle units need doubling
    return 2 * parmed2trappe(val)

def assert_dihedrals_match(foyer_d, reference_params):
    # Check that a dihedral from foyer matches the text in the csv file
    # foyer_d - parmed rb_torsion
    # reference_params - dihedral parameters from trappe csv file
    def conv1(t):
        # convert RB to original Trappe dihedral form:
        # V = c0 + c1(1 + cos phi) + c1(1 - cos 2phi) + c3(1 + cos 3phi)
        R0, R1, R2, R3 = t.c0, t.c1, t.c2, t.c3
        c0 = - R0 + R1 + R2 + R3
        c1 = - R1 - 3/4. * R3
        c2 = - 0.5 * R2
        c3 = - 1/4. * R3

        return c0, c1, c2, c3

    def conv2(t):
        # convert RB to second Trappe dihedral form:
        # V = c0 + c1 cos phi + c2 cos 2phi + c3 cos 3phi + c4 cos 4phi
        R0, R1, R2, R3, R4 = t.c0, t.c1, t.c2, t.c3, t.c4
        print(list(map(lambda x: x * calorie, [R0, R1, R2, R3, R4])))
        c0 = - R0 + 0.5 * R2 + 3/8. * R4
        c1 = - R1 - 3/4. * R3
        c2 = 0.5 * R2 + 0.5 * R4
        c3 = - 1/4. * R3
        c4 = 1/8. * R4

        return [parmed2trappe(v) for v in [c0, c1, c2, c3, c4]]

    reference_params = list(map(float, reference_params))
    # check parameters in *foyer_d* match the raw text from trappe
    # trappe's params are in one of two functional forms.....
    if len(reference_params) == 5:
        # using later form, has an extra column
        foyer_params = conv2(foyer_d.type)
    elif len(reference_params) == 4:
        # using original dihedral form
        foyer_params = conv1(foyer_d.type)

    assert_almost_equal(foyer_params, reference_params, decimal=3)

class TestTrappeReferences(object):
    @staticmethod
    def check_atoms(ref, struc):
        ref_atoms = sorted(ref['(pseudo)atom'],
                           key=lambda x: int(x.split(',')[0]))
        for ref_atom, foyer_atom in zip(ref_atoms, struc.atoms):
            _, _, _, epsilon, sigma, charge = ref_atom.strip().split(',')

            assert_almost_equal(foyer_atom.charge, float(charge))
            assert_almost_equal(foyer_atom.sigma, float(sigma))
            assert_almost_equal(parmed2trappe(foyer_atom.epsilon),
                                float(epsilon),
                                decimal=3)

    @staticmethod
    def check_bonds(ref, struc):
        ref_bonds = ref['stretch']
        # organise foyer's bonds into dict of indices
        bond_dict = {}
        for b in struc.bonds:
            i, j = sorted((b.atom1.idx, b.atom2.idx))
            bond_dict[i+1, j+1] = b

        # iterate reference file bonds
        # lookup foyer bond and compare length
        for ref_bond in ref_bonds:
            _, desc, _, length = ref_bond.strip().split(',')
            i, j = map(int, desc.strip("\"\'").split('-'))
            foyer_bond = bond_dict[i, j]

            assert_almost_equal(foyer_bond.type.req, float(length))
            # trappe is fixed bond length, no point in K comparison

    @staticmethod
    def check_angles(ref, struc):
        ref_angles = ref['bend']
        assert len(ref_angles) == len(struc.angles)

        angle_dict = {}
        for a in struc.angles:
            i, j, k = a.atom1.idx, a.atom2.idx, a.atom3.idx
            if k < i:
                i, j, k = k, j, i
            angle_dict[i+1, j+1, k+1] = a

        for ref_angle in ref_angles:
            _, desc, _, angle, k_angle = ref_angle.strip().split(',')
            i, j, k = map(int, desc.strip("\"\'").split('-'))

            foyer_angle = angle_dict[i, j, k]

            assert_almost_equal(foyer_angle.type.theteq,
                                float(angle),
                                decimal=4)
            assert_almost_equal(angle_parmed2trappe(foyer_angle.type.k),
                                float(k_angle),
                                decimal=3)

    @staticmethod
    def check_dihedrals(ref, struc):
        ref_dihedrals = ref['torsion']

        assert len(ref_dihedrals) == len(struc.rb_torsions)

        dih_dict = {}
        for d in struc.rb_torsions:
            i, j, k, l = d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx

            if l < i:
                i, j, k, l = l, k, j, i
            dih_dict[i+1, j+1, k+1, l+1] = d

        for ref_d in ref_dihedrals:
            _, desc, _, *params = ref_d.strip().split(',')
            i, j, k, l = map(int, desc.strip("\"\'").split('-'))
            if l < i:
                i, j, k, l = l, k, j, i

            foyer_d = dih_dict[i, j, k, l]

            assert_dihedrals_match(foyer_d, params)

    @pytest.mark.parametrize('trappe_ref,mol2file', [
        ('butadiene/trappe_parameters_94.csv', 'butadiene'),
        ('isoprene/trappe_parameters_95.csv', 'isoprene'),
        ('methyl_acrylate/trappe_parameters_96.csv', 'methyl_acrylate'),
    ], indirect=True)
    def test_trappe_reference(self, trappe_ref, mol2file):
        self.check_atoms(trappe_ref, mol2file)

        self.check_bonds(trappe_ref, mol2file)

        self.check_angles(trappe_ref, mol2file)

        self.check_dihedrals(trappe_ref, mol2file)


if __name__ == '__main__':
    TestTraPPE().find_correctly_implemented()
