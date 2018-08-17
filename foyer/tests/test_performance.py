import mbuild as mb
import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.tests.utils import get_fn


@pytest.mark.timeout(1)
def test_fullerene():
    fullerene = pmd.load_file(get_fn('fullerene.pdb'), structure=True)
    forcefield = Forcefield(get_fn('fullerene.xml'))
    forcefield.apply(fullerene, assert_dihedral_params=False)


@pytest.mark.timeout(15)
def test_surface():
    surface = mb.load(get_fn('silica.mol2'))
    forcefield = Forcefield(get_fn('opls-silica.xml'))
    forcefield.apply(surface)


@pytest.mark.timeout(45)
def test_polymer():
    peg100 = mb.load(get_fn('peg100.mol2'))
    forcefield = Forcefield(name='oplsaa')
    forcefield.apply(peg100)
