import glob
import itertools as it
import os

import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.tests.utils import atomtype
from foyer.xml_writer import write_foyer

def test_write_xml(filename, ff_file):
   structure = pmd.loadfile(filename)
   forcefield = Forcefield(ff_file)

   structure.write_foyer('test.xml', forcefield=forcefield)

def test_load_xml():
   structure = pmd.loadfile(filename)
   forcefield = Forcefield(ff_file)

   structure.write_foyer('test.xml', forcefield=forcefield)

   generated_ff = Forcefield('text.xml')
