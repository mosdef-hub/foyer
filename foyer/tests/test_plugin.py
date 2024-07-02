import pytest

import foyer
from foyer.tests.base_test import BaseTest


class TestPlugin(BaseTest):
    def test_basic_import(self):
        assert "forcefields" in dir(foyer)

    def test_loading_forcefields(self):
        for func in dir(foyer.forcefields):
            if "load_" in func and "__" not in func:
                eval("foyer.forcefields." + func)()

    def test_load_forcefield(self):
        foyer.forcefields.get_forcefield(name="oplsaa")
        foyer.forcefields.get_forcefield(name="trappe-ua")
        with pytest.raises(ValueError):
            foyer.forcefields.get_forcefield("bogus_name")
