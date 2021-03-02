import pytest
import parmed as pmd

from foyer import Forcefield
from foyer.exceptions import FoyerError
from foyer.tests.utils import get_fn

from foyer.tests.base_test import BaseTest


class TestRunAtomTyping(BaseTest):

    @pytest.fixture(scope='session')
    def missing_overrides_ff(self):
        return Forcefield(get_fn('missing_overrides.xml'))

    def test_missing_overrides(self, opls_validation_benzene, missing_overrides_ff):
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(opls_validation_benzene)

    def test_missing_definition(self, missing_overrides_ff):
        structure = pmd.load_file(get_fn('silly_chemistry.mol2'), structure=True)
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(structure)
