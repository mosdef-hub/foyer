import os

import parmed as pmd
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield
from foyer.exceptions import FoyerError
from foyer.tests.utils import get_fn

OPLS_TESTFILES_DIR = resource_filename('foyer', 'opls_validation')


def test_missing_overrides():
    top = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.top')
    gro = os.path.join(OPLS_TESTFILES_DIR, 'benzene/benzene.gro')
    structure = pmd.load_file(top, xyz=gro)

    forcefield = Forcefield(get_fn('missing_overrides.xml'))
    with pytest.raises(FoyerError):
        forcefield.apply(structure)


def test_missing_definition():
    structure = pmd.load_file(get_fn('silly_chemistry.mol2'), structure=True)

    forcefield = Forcefield(get_fn('missing_overrides.xml'))
    with pytest.raises(FoyerError):
        forcefield.apply(structure)
