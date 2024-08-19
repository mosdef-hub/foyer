from pathlib import Path

import parmed as pmd
import pytest
from pkg_resources import resource_filename

from foyer import forcefields
from foyer.smarts import SMARTS

OPLS_TEST_FILE_DIR = Path(resource_filename("foyer", "opls_validation")).resolve()


class BaseTest:
    @pytest.fixture(scope="session")
    def opls_validation_benzene(self):
        top = str(OPLS_TEST_FILE_DIR / "benzene/benzene.top")
        gro = str(OPLS_TEST_FILE_DIR / "benzene/benzene.gro")
        structure = pmd.load_file(top, xyz=gro)
        return structure

    @pytest.fixture(scope="session")
    def smarts_parser(self):
        return SMARTS()

    @pytest.fixture(scope="session")
    def oplsaa(self):
        return forcefields.load_OPLSAA()

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()
