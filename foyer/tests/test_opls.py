import os
import pytest
from foyer.forcefield import Forcefield
from foyer.tests.topology import Topology


class TestOPLS(object):

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    implemented_tests_path = os.path.join(os.path.dirname(__file__),
                                          'opls_test_cases.txt')

    test_topologies = [line.split() for line in open(implemented_tests_path)]
    test_topologies_mol_names = [x[0] for x in test_topologies]
    test_topologies_names = [x[1] for x in test_topologies]
    test_topologies_expected_results = [x[2]=='True' for x in test_topologies]


    @pytest.mark.parametrize('name,expected_result', zip(test_topologies_names, test_topologies_expected_results))
    def test_atomtyping(self, name, expected_result):

        if expected_result is not True:
            return

        known_structure = Topology.by_name(name, parameterized=True)
        untyped_structure = Topology.by_name(name, parameterized=False)
        typed_structure = Forcefield.by_name('oplsaa').apply(untyped_structure, in_place=False, debug=False)

        atomtyping_errors = find_atomtyping_errors(typed_structure, known_structure)

        assert not atomtyping_errors, atomtyping_errors


def find_atomtyping_errors(structure, known_structure):
    atomtyping_errors = list()

    for i, (atom, known_atom) in enumerate(zip(structure.atoms, known_structure.atoms)):
        if not atom.type:
            atomtyping_errors.append((i, 'Found no OPLS type for atom {} in {}: {}\n'
                                         'Should be atomtype: {}'.format(
                i, structure.title, atom.type, known_atom.type)))

        elif atom.type != known_atom.type:
            atomtyping_errors.append(
                (i, 'Found inconsistent OPLS type in {}: {} != {}'.format(structure.title,
                                                                               atom.type, known_atom.type)))
    return atomtyping_errors
