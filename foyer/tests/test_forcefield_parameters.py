import numpy as np
import pytest

from foyer import Forcefield, forcefields
from foyer.exceptions import MissingForceError, MissingParametersError
from foyer.forcefield import get_available_forcefield_loaders
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn
from foyer.utils.check_xml_dihedrals_RB_to_OPLS import (
    _test_xml_dihedrals,
    run_xml_test,
    test_xml_dihedral_rb_to_opls,
)


@pytest.mark.skipif(
    condition="load_GAFF"
    not in map(lambda func: func.__name__, get_available_forcefield_loaders()),
    reason="GAFF Plugin is not installed",
)
class TestForcefieldParameters(BaseTest):
    @pytest.fixture(scope="session")
    def gaff(self):
        return forcefields.load_GAFF()

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
        assert gaff.get_parameters(
            "harmonic_bonds", ["ca", "br"]
        ) == gaff.get_parameters("harmonic_bonds", ["ca", "br"])

    def test_gaff_missing_bond_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("harmonic_bonds", ["str1", "str2"])

    def test_gaff_angle_parameters(self, gaff):
        angle_params = gaff.get_parameters("harmonic_angles", ["f", "c1", "f"])
        assert np.allclose(
            [angle_params["theta"], angle_params["k"]],
            [3.141592653589793, 487.0176],
        )

    def test_gaff_angle_parameters_reversed(self, gaff):
        assert np.allclose(
            list(
                gaff.get_parameters(
                    "harmonic_angles", ["f", "c2", "ha"]
                ).values()
            ),
            list(
                gaff.get_parameters(
                    "harmonic_angles", ["ha", "c2", "f"]
                ).values()
            ),
        )

    def test_gaff_missing_angle_parameters(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("harmonic_angles", ["1", "2", "3"])

    def test_gaff_periodic_proper_parameters(self, gaff):
        periodic_proper_params = gaff.get_parameters(
            "periodic_propers", ["c3", "c", "sh", "hs"]
        )
        assert np.allclose(periodic_proper_params["periodicity"], [2.0, 1.0])
        assert np.allclose(
            periodic_proper_params["k"], [9.414, 5.4392000000000005]
        )
        assert np.allclose(
            periodic_proper_params["phase"],
            [3.141592653589793, 3.141592653589793],
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
        assert np.allclose(
            periodic_improper_params["phase"], [3.141592653589793]
        )

    def test_gaff_periodic_improper_parameters_reversed(self, gaff):
        assert np.allclose(
            list(
                gaff.get_parameters(
                    "periodic_impropers", ["c", "", "o", "o"]
                ).values()
            ),
            list(
                gaff.get_parameters(
                    "periodic_impropers", ["c", "o", "", "o"]
                ).values()
            ),
        )

    def test_gaff_proper_params_missing(self, gaff):
        with pytest.raises(MissingParametersError):
            gaff.get_parameters("periodic_impropers", ["a", "b", "c", "d"])

    def test_gaff_scaling_factors(self, gaff):
        assert gaff.lj14scale == 0.5
        assert np.isclose(gaff.coulomb14scale, 0.833333333)

    def test_opls_get_parameters_atoms(self, oplsaa):
        atom_params = oplsaa.get_parameters("atoms", "opls_145")
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_atoms_list(self, oplsaa):
        atom_params = oplsaa.get_parameters("atoms", ["opls_145"])
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_atom_class(self, oplsaa):
        atom_params = oplsaa.get_parameters(
            "atoms", "CA", keys_are_atom_classes=True
        )
        assert atom_params["sigma"] == 0.355
        assert atom_params["epsilon"] == 0.29288

    def test_opls_get_parameters_bonds(self, oplsaa):
        bond_params = oplsaa.get_parameters(
            "harmonic_bonds", ["opls_760", "opls_145"]
        )
        assert bond_params["length"] == 0.146
        assert bond_params["k"] == 334720.0

    def test_opls_get_parameters_bonds_reversed(self, oplsaa):
        assert np.allclose(
            list(
                oplsaa.get_parameters(
                    "harmonic_bonds", ["opls_760", "opls_145"]
                ).values()
            ),
            list(
                oplsaa.get_parameters(
                    "harmonic_bonds", ["opls_145", "opls_760"]
                ).values()
            ),
        )

    def test_opls_get_parameters_bonds_atom_classes_reversed(self, oplsaa):
        assert np.allclose(
            list(
                oplsaa.get_parameters(
                    "harmonic_bonds", ["C_2", "O_2"], True
                ).values()
            ),
            list(
                oplsaa.get_parameters(
                    "harmonic_bonds", ["O_2", "C_2"], True
                ).values()
            ),
        )

    def test_opls_get_parameters_angle(self, oplsaa):
        angle_params = oplsaa.get_parameters(
            "harmonic_angles", ["opls_166", "opls_772", "opls_167"]
        )
        assert np.allclose(
            [angle_params["theta"], angle_params["k"]], [2.0943950239, 585.76]
        )

    def test_opls_get_parameters_angle_reversed(self, oplsaa):
        assert np.allclose(
            list(
                oplsaa.get_parameters(
                    "harmonic_angles", ["opls_166", "opls_772", "opls_167"]
                ).values()
            ),
            list(
                oplsaa.get_parameters(
                    "harmonic_angles", ["opls_167", "opls_772", "opls_166"]
                ).values()
            ),
        )

    def test_opls_get_parameters_angle_atom_classes(self, oplsaa):
        angle_params = oplsaa.get_parameters(
            "harmonic_angles", ["CA", "C_2", "CA"], keys_are_atom_classes=True
        )

        assert np.allclose(
            [angle_params["theta"], angle_params["k"]], [2.09439510239, 711.28]
        )

    def test_opls_get_parameters_angle_atom_classes_reversed(self, oplsaa):
        assert np.allclose(
            list(
                oplsaa.get_parameters(
                    "harmonic_angles",
                    ["CA", "C", "O"],
                    keys_are_atom_classes=True,
                ).values()
            ),
            list(
                oplsaa.get_parameters(
                    "harmonic_angles",
                    ["O", "C", "CA"],
                    keys_are_atom_classes=True,
                ).values()
            ),
        )

    def test_opls_get_parameters_rb_proper(self, oplsaa):
        proper_params = oplsaa.get_parameters(
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

    def test_get_parameters_rb_proper_reversed(self, oplsaa):
        assert np.allclose(
            list(
                oplsaa.get_parameters(
                    "rb_propers",
                    ["opls_215", "opls_215", "opls_235", "opls_269"],
                ).values()
            ),
            list(
                oplsaa.get_parameters(
                    "rb_propers",
                    ["opls_269", "opls_235", "opls_215", "opls_215"],
                ).values()
            ),
        )

    def test_opls_get_parameters_wildcard(self, oplsaa):
        proper_params = oplsaa.get_parameters(
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

    def test_opls_missing_force(self, oplsaa):
        with pytest.raises(MissingForceError):
            oplsaa.get_parameters("periodic_propers", key=["a", "b", "c", "d"])

    def test_opls_scaling_factors(self, oplsaa):
        assert oplsaa.lj14scale == 0.5
        assert oplsaa.coulomb14scale == 0.5

    def test_missing_scaling_factors(self):
        ff = Forcefield(forcefield_files=(get_fn("validate_customtypes.xml")))
        with pytest.raises(AttributeError):
            assert ff.lj14scale
        with pytest.raises(AttributeError):
            assert ff.coulomb14scale

    # check xml files for dihedral conversions from RB to OPLS
    # Note: the non-exact conversions are most likely due to the OPLS equation
    # not using the f0/2 term (i.e., f0 = 0), regardless if the analytical
    # conversion provides a non-zero f0 term.
    def test_oplsaa_xml_dihedral_rb_to_opls(self):
        oplsaa_non_exact_conversion_list = [
            [
                ["CT", "OS", "P", "OS"],
                [1.046, 3.138, 10.0416, -4.184, 0.0, 0.0],
                [20.0832, 0.0, -10.0416, 2.092, -0.0],
                10.0416,
                False,
            ]
        ]

        pass_oplsaa_rb_to_opls = test_xml_dihedral_rb_to_opls(
            "oplsaa",
            "all_oplsaa_rb_to_opls_errors.txt",
            oplsaa_non_exact_conversion_list,
            error_tolerance_rb_to_opls=1e-4,
        )
        assert pass_oplsaa_rb_to_opls is True

    def test_trappeua_xml_dihedral_rb_to_opls(self):
        trappeua_non_exact_conversion_list = [
            [
                ["CH3", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH2", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH2", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH2", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH2", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH2", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH2", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH2", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH2", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH2", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["C", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["C", "CH2", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["C", "CH2", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["C", "CH2", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH", "CH", "CH2"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH2", "CH", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH", "CH", "CH"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH", "CH", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["C", "CH", "CH", "C"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                False,
            ],
            [
                ["CH3", "CH", "O", "H"],
                [2.51338, -5.97885, -0.52315, 5.78421, 0.0, 0.0],
                [3.59118, 3.281385, 0.52315, -2.892105, -0.0],
                1.79559,
                False,
            ],
            [
                ["CH2", "CH", "O", "H"],
                [2.51338, -5.97885, -0.52315, 5.78421, 0.0, 0.0],
                [3.59118, 3.281385, 0.52315, -2.892105, -0.0],
                1.79559,
                False,
            ],
            [
                ["CH", "CH", "O", "H"],
                [2.51338, -5.97885, -0.52315, 5.78421, 0.0, 0.0],
                [3.59118, 3.281385, 0.52315, -2.892105, -0.0],
                1.79559,
                False,
            ],
            [
                ["C", "CH", "O", "H"],
                [2.51338, -5.97885, -0.52315, 5.78421, 0.0, 0.0],
                [3.59118, 3.281385, 0.52315, -2.892105, -0.0],
                1.79559,
                False,
            ],
        ]

        pass_trappeua_rb_to_opls = test_xml_dihedral_rb_to_opls(
            "trappe-ua",
            "all_trappeua_rb_to_opls_errors.txt",
            trappeua_non_exact_conversion_list,
            error_tolerance_rb_to_opls=1e-4,
        )
        assert pass_trappeua_rb_to_opls is True

    def test_xml_dihedrals_function_error_not_float(self):
        xml_file_directory_and_filename = "test_ff"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = "s"

        with pytest.raises(
            TypeError,
            match=f"The error_tolerance_rb_to_opls variable must be a float, "
            f"is type {type(error_tolerance_rb_to_opls)}.",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    # error_tolerance_rb_to_opls range is 1e-1 to  1e-10
    def test_xml_dihedrals_function_error_too_high(self):
        xml_file_directory_and_filename = "test_ff"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = 1.01e-1

        with pytest.raises(
            ValueError,
            match=f"The error_tolerance_rb_to_opls variable is not 1e-10 \<= "
            f"\( entered value is {error_tolerance_rb_to_opls} \) \<= 1e-1.",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    def test_xml_dihedrals_function_error_too_low(self):
        xml_file_directory_and_filename = "test_ff"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = 0.999e-10

        with pytest.raises(
            ValueError,
            match=f"The error_tolerance_rb_to_opls variable is not 1e-10 \<= "
            f"\( entered value is {error_tolerance_rb_to_opls} \) \<= 1e-1.",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    def test_xml_dihedrals_function_ff_xml_not_good_extension(self):
        xml_file_directory_and_filename = "test_ff.out"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = 1e-4

        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct "
            "foyer FF name and not a path to a FF file. "
            "If you are entering a path to a FF file, "
            "please use the forcefield_files variable with the "
            "proper XML extension \(.xml\).",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    def test_xml_dihedrals_function_user_defined_ff_xml_not_exist(self):
        xml_file_directory_and_filename = "test_ff.xml"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = 1e-4

        with pytest.raises(
            ValueError,
            match="Please make sure you are entering the correct foyer FF path, "
            "including the FF file name.xml "
            "If you are using the pre-build FF files in foyer, "
            "only use the string name without any extension.",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    def test_xml_dihedrals_function_foyer_std_ff_xml_not_exist(self):
        xml_file_directory_and_filename = "test_ff"
        output_file_name = "test_filename.txt"
        error_tolerance_rb_to_opls = 1e-4

        with pytest.raises(
            ValueError,
            match="Please make sure you are entering the correct foyer FF name "
            "without the .xml extension.",
        ):
            _test_xml_dihedrals(
                xml_file_directory_and_filename,
                output_file_name,
                error_tolerance_rb_to_opls=error_tolerance_rb_to_opls,
            )

    def test_run_xml_test_input_not_dict(self):
        xml_and_error_file_dict = ["test_ff.xml", "test_filename.txt"]

        with pytest.raises(
            TypeError,
            match=f"The xml_and_error_file_dict variable must be a dict, "
            f"is type {type(xml_and_error_file_dict)}.",
        ):
            run_xml_test(xml_and_error_file_dict)

    def test_run_xml_test_input_key_not_str(self):
        xml_key = 1
        xml_value = "test_filename.txt"
        xml_and_error_file_dict = {xml_key: xml_value}

        with pytest.raises(
            TypeError,
            match=f"The xml_and_error_file_dict keys are not all strings, "
            f"and has the type {type(xml_key)}.",
        ):
            run_xml_test(xml_and_error_file_dict)

    def test_run_xml_test_input_value_not_str(self):
        xml_key = "test_ff.xml"
        xml_value = 1
        xml_and_error_file_dict = {xml_key: xml_value}

        with pytest.raises(
            TypeError,
            match=f"The xml_and_error_file_dict values are not all strings, "
            f"and has the type {type(xml_value)}.",
        ):

            run_xml_test(xml_and_error_file_dict)

    def test_xml_dihedral_rb_to_opls_non_exact_conversion_list_str(self):
        ff_filename = "test_ff.xml"
        error_filename = "test_filename.txt"
        non_exact_conversion_list = "str"

        with pytest.raises(
            TypeError,
            match=f"The non_exact_conversion_list variables is not formated correctly. "
            f"Please see the test_xml_dihedral_rb_to_opls function for the "
            f"proper format infomation.",
        ):
            test_xml_dihedral_rb_to_opls(
                ff_filename, error_filename, non_exact_conversion_list
            )

    def test_xml_dihedral_rb_to_opls_non_exact_conversion_list_wrong_type_str(
        self,
    ):
        ff_filename = "test_ff.xml"
        error_filename = "test_filename.txt"
        non_exact_conversion_list = [
            [
                ["CH3", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
                "False",
            ]
        ]

        with pytest.raises(
            TypeError,
            match=f"The non_exact_conversion_list variables is not formated correctly. "
            f"Please see the test_xml_dihedral_rb_to_opls function for the "
            f"proper format infomation.",
        ):
            test_xml_dihedral_rb_to_opls(
                ff_filename, error_filename, non_exact_conversion_list
            )

    def test_xml_dihedral_rb_to_opls_non_exact_conversion_list_4_length(self):
        ff_filename = "test_ff.xml"
        error_filename = "test_filename.txt"
        non_exact_conversion_list = [
            [
                ["CH3", "CH2", "CH", "CH3"],
                [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
                [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
                2.08734,
            ]
        ]

        with pytest.raises(
            TypeError,
            match=f"The non_exact_conversion_list variables is not formated correctly. "
            f"Please see the test_xml_dihedral_rb_to_opls function for the "
            f"proper format infomation.",
        ):
            test_xml_dihedral_rb_to_opls(
                ff_filename, error_filename, non_exact_conversion_list
            )
