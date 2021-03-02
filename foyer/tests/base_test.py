from pathlib import Path
from pkg_resources import resource_filename

import parmed as pmd
import pytest

from foyer.smarts import SMARTS
from foyer import Forcefield

OPLS_TEST_FILE_DIR = Path(resource_filename("foyer", "opls_validation")).resolve()


class BaseTest:
    @pytest.fixture(scope="session")
    def opls_validation_benzene(self):
        top = str(OPLS_TEST_FILE_DIR / "benzene/benzene.top")
        gro = str(OPLS_TEST_FILE_DIR / "benzene/benzene.gro")
        structure = pmd.load_file(top, xyz=gro)
        return structure

    @pytest.fixture(scope='session')
    def smarts_parser(self):
        return SMARTS()

    @pytest.fixture(scope='session')
    def oplsaa(self):
        return Forcefield(name='oplsaa')
