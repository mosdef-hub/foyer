import glob
import os
from pkg_resources import resource_filename
import pytest

from foyer import Forcefield


FF_DIR = resource_filename('foyer', 'forcefields')
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, '*.xml'))

def test_load_files():
    ff1 = Forcefield(forcefield_files=FORCEFIELDS)
    assert len(ff1._atomTypes) > 0

    ff2 = Forcefield(forcefield_files=FORCEFIELDS[0])
    assert len(ff1._atomTypes) == len(ff2._atomTypes)

    ff3 = Forcefield(name='oplsaa')
    assert len(ff1._atomTypes) == len(ff3._atomTypes)

def test_duplicate_type_definitions():
    with pytest.raises(ValueError):
        ff4 = Forcefield(name='oplsaa', forcefield_files=FORCEFIELDS)
