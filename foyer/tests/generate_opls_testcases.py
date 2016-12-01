import glob
import os

import numpy as np
import parmed as pmd
import sys
from pkg_resources import resource_filename
import pytest

from foyer.forcefield import Forcefield
from foyer.atomtyper import find_atomtypes
from foyer.tests.test_structure import TestStructure
from foyer.tests.test_opls import find_atomtyping_errors
# tmpdir.chdir()

opls_test_case_list_path = os.path.join(os.path.dirname(__file__),
                                      'opls_test_cases.txt')

def create_mol2s():
    for name, (mol2_path, top_path, gro_path) in TestStructure.available().items():
        if not os.path.isfile(mol2_path):
            known_structure = TestStructure.by_name(name, parameterized=True)
            known_structure.save(mol2_path, format='MOL2', overwrite=True)

def update_opls_test_case_list():
    with open(opls_test_case_list_path, 'w') as fh:
        for name, (mol2_path, top_path, gro_path) in TestStructure.available().items():

            try:
                known_structure = TestStructure.by_name(name, parameterized=True)
                untyped_structure = TestStructure.by_name(name, parameterized=False)
                typed_structure = Forcefield.by_name('oplsaa').apply(untyped_structure, in_place=False, debug=False)

                atomtyping_errors = find_atomtyping_errors(typed_structure, known_structure)

                if atomtyping_errors:
                    print('ERROR: {}'.format(atomtyping_errors))
                    expected_result = False
                else:
                    print(' *** SUCCESS {}\n'.format(name))
                    expected_result = True
            except:
                print('ERROR: {}'.format(atomtyping_errors))
                expected_result = False

            basename = os.path.basename(top_path[:-4])

            fh.write('{} {} {}\n'.format(known_structure.title, basename, expected_result))



if __name__ == '__main__':
    create_mol2s()
    update_opls_test_case_list()
