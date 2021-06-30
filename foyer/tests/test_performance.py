import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn
from foyer.utils.io import has_mbuild


class TestPerformance(BaseTest):
    @pytest.mark.timeout(1)
    def test_fullerene(self):
        fullerene = pmd.load_file(get_fn("fullerene.pdb"), structure=True)
        forcefield = Forcefield(get_fn("fullerene.xml"))
        forcefield.apply(fullerene, assert_dihedral_params=False)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    @pytest.mark.timeout(15)
    def test_surface(self):
        import mbuild as mb

        surface = mb.load(get_fn("silica.mol2"))
        forcefield = Forcefield(get_fn("opls-silica.xml"))
        forcefield.apply(surface, assert_bond_params=False)

    @pytest.mark.skipif(not has_mbuild, reason="mbuild is not installed")
    @pytest.mark.timeout(45)
    def test_polymer(self, oplsaa):
        import mbuild as mb

        peg100 = mb.load(get_fn("peg100.mol2"))
        oplsaa.apply(peg100)
