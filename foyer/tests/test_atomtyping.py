import os

import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.exceptions import FoyerError
from foyer.tests.utils import get_fn
from foyer.tests.base_test import BaseTest

OPLS_TESTFILES_DIR = resource_filename('foyer', 'opls_validation')


class TestAtomTyping(BaseTest):
    def test_missing_overrides(self):
        top = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.top')
        gro = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.gro')
        structure = pmd.load_file(top, xyz=gro)
    
        forcefield = Forcefield(get_fn('missing_overrides.xml'))
        with pytest.raises(FoyerError):
            forcefield.apply(structure)
    
    
    def test_missing_definition(self):
        structure = pmd.load_file(get_fn('silly_chemistry.mol2'), structure=True)
    
        forcefield = Forcefield(get_fn('missing_overrides.xml'))
        with pytest.raises(FoyerError):
            forcefield.apply(structure)
