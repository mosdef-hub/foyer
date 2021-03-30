import pytest
import numpy as np

from foyer import forcefields
from foyer.forcefield import get_available_forcefield_loaders
from foyer.exceptions import MissingParametersError, MissingForceError


@pytest.mark.skipif(
    condition="load_GAFF"
    not in map(lambda func: func.__name__, get_available_forcefield_loaders()),
    reason="GAFF Plugin is not installed",
)
class TestForcefieldParameters:
    @pytest.fixture(scope="session")
    def gaff(self):
        return forcefields.load_GAFF()

    @pytest.fixture(scope="session")
    def opls(self):
        return forcefields.load_OPLSAA()

    def test_gaff_missing_group(self, gaff):
        with pytest.raises(ValueError):
            gaff.get_parameters("missing", key=[])

    def test_gaff_non_string_keys(self, gaff):
        with pytest.raises(TypeError):
            gaff.get_parameters("atoms", key=1)

    def test_gaff_bond_parameters_gaff(self, gaff):
        bond_params = gaff.get_parameters("harmonic_bonds", ["br", "ca"])
        assert np.isclose(bond_params["length"], 0.19079)
        assert np.isclose(bond_params["k"], 219827.36)

    def test_gaff_bond_params_reversed(self, gaff):
        assert gaff.get_parameters("harmonic_bonds", ["ca", "br"]) == gaff.get_parameters(
            "harmonic_bonds", ["ca", "br"]
        )

    def test_gaff_missing_bond_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("harmonic_bonds", ["str1", "str2"])

    def test_gaff_angle_parameters(self, gaff):
        angle_params = gaff.get_parameters("harmonic_angles", ["f", "c1", "f"])
        assert np.allclose(
            [angle_params["theta"], angle_params["k"]], [3.141592653589793, 487.0176]
        )

    def test_gaff_angle_parameters_reversed(self, gaff):
        assert np.allclose(
            list(gaff.get_parameters("harmonic_angles", ["f", "c2", "ha"]).values()),
            list(gaff.get_parameters("harmonic_angles", ["ha", "c2", "f"]).values()),
        )

    def test_gaff_missing_angle_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("harmonic_angles", ["1", "2", "3"])

    def test_gaff_periodic_proper_parameters(self, gaff):
        periodic_proper_params = gaff.get_parameters(
            "periodic_propers", ["c3", "c", "sh", "hs"]
        )
        assert np.allclose(periodic_proper_params["periodicity"], [2.0, 1.0])
        assert np.allclose(periodic_proper_params["k"], [9.414, 5.4392000000000005])
        assert np.allclose(
            periodic_proper_params["phase"], [3.141592653589793, 3.141592653589793]
        )

    def test_gaff_periodic_proper_parameters_reversed(self, gaff):
        assert np.allclose(
            list(
                gaff.get_parameters(
                    "periodic_propers", ["c3", "c", "sh", "hs"]
                ).values()
            ),
            list(
                gaff.get_parameters(
                    "periodic_propers", ["hs", "sh", "c", "c3"]
                ).values()
            ),
        )

    def test_gaff_periodic_improper_parameters(self, gaff):
        periodic_improper_params = gaff.get_parameters(
            "periodic_impropers", ["c", "", "o", "o"]
        )
        assert np.allclose(periodic_improper_params["periodicity"], [2.0])
        assert np.allclose(periodic_improper_params["k"], [4.6024])
        assert np.allclose(periodic_improper_params["phase"], [3.141592653589793])

    def test_gaff_periodic_improper_parameters_reversed(self, gaff):
        assert np.allclose(
            list(
                gaff.get_parameters("periodic_impropers", ["c", "", "o", "o"]).values()
            ),
            list(
                gaff.get_parameters("periodic_impropers", ["c", "o", "", "o"]).values()
            ),
        )

    def test_gaff_proper_params_missing(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("periodic_impropers", ["a", "b", "c", "d"])

    def test_opls_get_parameters_atoms(self, opls):
        atom_params = opls.get_parameters("atoms", "opls_145")
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_atoms_list(self, opls):
        atom_params = opls.get_parameters("atoms", ["opls_145"])
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_atom_class(self, opls):
        atom_params = opls.get_parameters("atoms", "CA", keys_are_atom_classes=True)
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_bonds(self, opls):
        bond_params = opls.get_parameters("harmonic_bonds", ["opls_760", "opls_145"])
        assert bond_params["length"] == 0.146
        assert bond_params["k"] == 334720.0

    def test_opls_get_parameters_bonds_reversed(self, opls):
        assert np.allclose(
            list(opls.get_parameters("harmonic_bonds", ["opls_760", "opls_145"]).values()),
            list(opls.get_parameters("harmonic_bonds", ["opls_145", "opls_760"]).values()),
        )

    def test_opls_get_parameters_bonds_atom_classes_reversed(self, opls):
        assert np.allclose(
            list(opls.get_parameters("harmonic_bonds", ["C_2", "O_2"], True).values()),
            list(opls.get_parameters("harmonic_bonds", ["O_2", "C_2"], True).values()),
        )

    def test_opls_get_parameters_angle(self, opls):
        angle_params = opls.get_parameters(
            "harmonic_angles", ["opls_166", "opls_772", "opls_167"]
        )
        assert np.allclose(
            [angle_params["theta"], angle_params["k"]], [2.0943950239, 585.76]
        )

    def test_opls_get_parameters_angle_reversed(self, opls):
        assert np.allclose(
            list(
                opls.get_parameters(
                    "harmonic_angles", ["opls_166", "opls_772", "opls_167"]
                ).values()
            ),
            list(
                opls.get_parameters(
                    "harmonic_angles", ["opls_167", "opls_772", "opls_166"]
                ).values()
            ),
        )

    def test_opls_get_parameters_angle_atom_classes(self, opls):
        angle_params = opls.get_parameters(
            "harmonic_angles", ["CA", "C_2", "CA"], keys_are_atom_classes=True
        )

        assert np.allclose(
            [angle_params["theta"], angle_params["k"]], [2.09439510239, 711.28]
        )

    def test_opls_get_parameters_angle_atom_classes_reversed(self, opls):
        assert np.allclose(
            list(
                opls.get_parameters(
                    "harmonic_angles", ["CA", "C", "O"], keys_are_atom_classes=True
                ).values()
            ),
            list(
                opls.get_parameters(
                    "harmonic_angles", ["O", "C", "CA"], keys_are_atom_classes=True
                ).values()
            ),
        )

    def test_opls_get_parameters_rb_proper(self, opls):
        proper_params = opls.get_parameters(
            "rb_propers", ["opls_215", "opls_215", "opls_235", "opls_269"]
        )
        assert np.allclose(
            [
                proper_params["c0"],
                proper_params["c1"],
                proper_params["c2"],
                proper_params["c3"],
                proper_params["c4"],
                proper_params["c5"],
            ],
            [2.28446, 0.0, -2.28446, 0.0, 0.0, 0.0],
        )

    def test_get_parameters_rb_proper_reversed(self, opls):
        assert np.allclose(
            list(
                opls.get_parameters(
                    "rb_propers", ["opls_215", "opls_215", "opls_235", "opls_269"]
                ).values()
            ),
            list(
                opls.get_parameters(
                    "rb_propers", ["opls_269", "opls_235", "opls_215", "opls_215"]
                ).values()
            ),
        )

    def test_opls_get_parameters_wildcard(self, opls):
        proper_params = opls.get_parameters(
            "rb_propers", ["", "opls_235", "opls_544", ""]
        )

        assert np.allclose(
            [
                proper_params["c0"],
                proper_params["c1"],
                proper_params["c2"],
                proper_params["c3"],
                proper_params["c4"],
                proper_params["c5"],
            ],
            [30.334, 0.0, -30.334, 0.0, 0.0, 0.0],
        )

    def test_opls_missing_force(self, opls):
        with pytest.raises(MissingForceError):
            opls.get_parameters("periodic_propers", key=["a", "b", "c", "d"])
