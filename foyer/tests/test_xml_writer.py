import parmed as pmd
import pytest
import os

from pkg_resources import resource_filename
from foyer import Forcefield
from foyer.xml_writer import write_foyer


OPLS_TESTFILES_DIR = resource_filename('foyer', 'opls_validation')

def test_write_xml():
    top = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.top')
    gro = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.gro')
    structure = pmd.load_file(top, xyz=gro)
    forcefield = Forcefield(name='oplsaa')
    param_struc = forcefield.apply(structure)

    param_struc.write_foyer('test.xml', forcefield=forcefield)

def test_load_xml():
    top = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.top')
    gro = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.gro')
    structure = pmd.load_file(top, xyz=gro)
    forcefield = Forcefield(name='oplsaa')
    param_struc = forcefield.apply(structure)

    param_struc.write_foyer('test.xml', forcefield=forcefield)

    generated_ff = Forcefield('test.xml')
