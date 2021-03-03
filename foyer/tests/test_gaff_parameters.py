import pytest

from foyer import forcefields
from foyer.exceptions import MissingParametersError


@pytest.mark.skipif(
    condition=not "load_GAFF" in dir(forcefields), reason="GAFF Plugin is not installed"
)
class TestGAFFParameters:
    @pytest.fixture(scope="session")
    def gaff(self):
        return forcefields.load_GAFF()

    def test_missing_group(self, gaff):
        with pytest.raises(ValueError):
            gaff.get_parameters("missing", key=[])

    def test_non_string_keys(self, gaff):
        with pytest.raises(TypeError):
            gaff.get_parameters("atoms", key=1)

    def test_bond_parameters_gaff(self, gaff):
        bond_params = gaff.get_parameters("bonds", ["br", "ca"])
        assert bond_params["length"] == 0.19079000000000002
        assert bond_params["k"] == 219827.35999999996

    def test_bond_params_reversed(self, gaff):
        assert gaff.get_parameters("bonds", ["ca", "br"]) == gaff.get_parameters(
            "bonds", ["ca", "br"]
        )

    def test_missing_bond_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("bonds", ["str1", "str2"])

    def test_angle_parameters(self, gaff):
        angle_params = gaff.get_parameters("angles", ["f", "c1", "f"])
        assert angle_params["theta"] == 3.141592653589793
        assert angle_params["k"] == 487.0176

    def test_missing_angle_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("angles", ["1", "2", "3"])

    def test_periodic_proper_parameters(self, gaff):
        periodic_proper_params = gaff.get_parameters(
            "periodic_propers", ["c3", "c", "sh", "hs"]
        )
        assert periodic_proper_params["periodicity"] == [2.0, 1.0]
        assert periodic_proper_params["k"] == [9.414, 5.4392000000000005]
        assert periodic_proper_params["phase"] == [3.141592653589793, 3.141592653589793]

    def test_periodic_improper_parameters(self, gaff):
        periodic_improper_params = gaff.get_parameters(
            "periodic_impropers", ["c", "", "o", "o"]
        )
        assert periodic_improper_params["periodicity"] == [2.0]
        assert periodic_improper_params["k"] == [4.6024]
        assert periodic_improper_params["phase"] == [3.141592653589793]

    def test_proper_params_missing(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("periodic_impropers", ["a", "b", "c", "d"])
