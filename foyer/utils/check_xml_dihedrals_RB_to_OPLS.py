import numpy as np
import os
import xml.etree.ElementTree as ET
from foyer.forcefields import forcefields
from foyer.forcefields.forcefields import get_forcefield
from mbuild.utils.conversion import RB_to_OPLS
from warnings import warn


def _test_xml_dihedrals(xml_file_directory_and_filename, output_file_name,
                     error_tolerance_rb_to_opls=1e-4,
                     ):
    r"""In in file function for the conversion of the Ryckaert-Bellemans (RB) type to OPLS
    type dihedrals (i.e., the mBuild RB_to_OPLS function) for errors and generates ouput
    files with the errors for the single xml file provided, containing conversion errors
    and their info.

    .. math::
    RB_{torsions} &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                  &= c_4*cos(psi)^4 + c_5*cos(5*psi)^5

    .. math::
    Exact_OPLS_torsions &= \frac{f_0}{2} + \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                        &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    .. math::
    Standard_ OPLS_torsions &= \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                           &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    where :math:`psi = t - pi = t - 180 degrees`

    xml_file_directory_and_filename : str
        A string with the RB torsions xml file path, file names, and
        extensions (i.e., '.xml'), or the just the name of the standard foyer force
        without the extension.
    output_file_name : str
         The string which will be the names error output files for the xml_file_list.
    error_tolerance_rb_to_opls : float (range 1e-1 to 1e-10), default=1e-4
        The acceptable absolute tolerance between the RB to OPLS conversion.
        This function converts it to an absolute tolerance.
        The working range is 1e-1 to 1e-10, since the calculated values are rounded
        down to the working range plus 2 decimal places (1e-1 -> 3 decimal places and
        2e-10 -> 11 decimal places).

    Returns
    -------
    Generates text (.txt) file containing the dihedral errors in the mBuild
    RB_to_OPLS function conversion process for the specified force field
    xml file. These files contain the following columns and have the specified meaning:

    dihedral_atoms_or_beads: list of strings [str, ..., str]
        A list of the foyer atoms or beads types in the dihedral, in order.
    rb_constants_c0_c1_c2_c3_c4_c5: list of floats [float, ..., float]
        A list of the RB-torsion values for the dihedral, which includes
        c0, c1, c2, c3, c4, c5.
    opls_constants_f0_f1_f2_f3_f4_rounded_6_decimals : list of floats [float, ..., float]
        A list of the OPLS dihedral values for the dihedral, which includes
        f0, f1, f2, f3, f4.
        NOTE: The standard OPLS dihedral function does not include
        the f0 constant. However, the f0 is included as it is required for the exact analytical conversion.
        Otherwise, for the standard OPLS dihedral, the f0 term must be 0
        for the solution to be analytically correct.
    max_abs_delta_rb_to_opls_calc : float
        The maximum difference between the RB and standard OPLS dihedral functions
        for 100 equally spaced points over the span of 2 Pi.
    exact_rb_to_opls_possible : bool
        This states whether the conversion is possible from RB to OPLS, since
        OPLS does not include the f0 constant. Therefore, if f0/2 is equal to zero
        within the set tolerance (error_tolerance_rb_to_opls), it will be True.
        Otherwise, it will be False.
    """
    if not isinstance(error_tolerance_rb_to_opls, float):
        raise TypeError(
                f"The error_tolerance_rb_to_opls variable must be a float, "
                f"is type {type(error_tolerance_rb_to_opls)}."
            )
    else:
        if not 1e-10 <= error_tolerance_rb_to_opls <= 1e-1:
            print(f"printed  The error_tolerance_rb_to_opls variable is not 1e-10 <= "
                   f"( entered value is {error_tolerance_rb_to_opls} ) <= 1e-1.")
            raise ValueError(
                f"The error_tolerance_rb_to_opls variable is not 1e-10 <= "
                   f"( entered value is {error_tolerance_rb_to_opls} ) <= 1e-1."
            )

    if (
            os.path.splitext(xml_file_directory_and_filename)[1] == ".xml"
    ):
        foyer_packaged_ff = False

    elif (
            os.path.splitext(xml_file_directory_and_filename)[1] == ""
    ):
        foyer_packaged_ff = True

    else:
        print_error_message = (
            r"Please make sure you are entering the correct "
            "foyer FF name and not a path to a FF file. "
            "If you are entering a path to a FF file, "
            "please use the forcefield_files variable with the "
            "proper XML extension (.xml)."
        )
        raise ValueError(print_error_message)

    if foyer_packaged_ff is False:
        ff_full_path_and_filename = xml_file_directory_and_filename
        try:
            ff_xml = ET.parse(ff_full_path_and_filename)
        except:
            print_error_message = (
                "Please make sure you are entering the correct foyer FF path, "
                "including the FF file name.xml "
                "If you are using the pre-build FF files in foyer, "
                "only use the string name without any extension."
            )
            raise ValueError(print_error_message)

    elif foyer_packaged_ff is True:
        ff_full_path_and_filename = (
                forcefields.get_ff_path()[0]
                + "/xml/"
                + xml_file_directory_and_filename
                + ".xml"
        )
        try:
            ff_xml = ET.parse(ff_full_path_and_filename)
        except:
            print_error_message = (
                "Please make sure you are entering the correct foyer FF name "
                "without the .xml extension."
            )
            raise ValueError(print_error_message)

    with open(output_file_name, "w") as data_opls_rb:
        ff_root = ff_xml.getroot()
        rb_torsionForce_proper_root = ff_root.findall('RBTorsionForce/Proper')

        count_no = 0
        for child in rb_torsionForce_proper_root:
            count_no += 1
            class1_iter = child.attrib['class1']
            class2_iter = child.attrib['class2']
            class3_iter = child.attrib['class3']
            class4_iter = child.attrib['class4']

            class_list_iter = [class1_iter,
                               class2_iter,
                               class3_iter,
                               class4_iter]

            c0_iter = float(child.attrib['c0'])
            c1_iter = float(child.attrib['c1'])
            c2_iter = float(child.attrib['c2'])
            c3_iter = float(child.attrib['c3'])
            c4_iter = float(child.attrib['c4'])
            c5_iter = float(child.attrib['c5'])

            cX_list_iter = [c0_iter,
                            c1_iter,
                            c2_iter,
                            c3_iter,
                            c4_iter,
                            c5_iter]

            fX_list_iter = RB_to_OPLS(c0_iter,
                                      c1_iter,
                                      c2_iter,
                                      c3_iter,
                                      c4_iter,
                                      c5_iter,
                                      error_tolerance = error_tolerance_rb_to_opls,
                                      error_if_outside_tolerance=False
                                      )

            f0_iter = fX_list_iter[0]
            f1_iter = fX_list_iter[1]
            f2_iter = fX_list_iter[2]
            f3_iter = fX_list_iter[3]
            f4_iter = fX_list_iter[4]

            fx_list_iter_6_decimals = [np.round(fX_list_iter[0], decimals=6),
                                       np.round(fX_list_iter[1], decimals=6),
                                       np.round(fX_list_iter[2], decimals=6),
                                       np.round(fX_list_iter[3], decimals=6),
                                       np.round(fX_list_iter[4], decimals=6)
                                       ]

            # generate a 100 list for the functions over range of 2 * Pi
            no_points = 100
            radian_range = 2 * np.pi
            rb_torsion_list = []
            rb_to_opls_list = []
            rad_list = []
            abs_differene_opls_rb_list = []
            for rad_fract in range(0, no_points):
                rad_iter = rad_fract * radian_range / no_points
                rad_list.append(rad_iter)

                rad_psi_iter = rad_iter + np.pi

                # RB torsions
                rb_torsion_iter = c0_iter \
                                  + c1_iter * (np.cos(rad_psi_iter)) ** 1 \
                                  + c2_iter * (np.cos(rad_psi_iter)) ** 2 \
                                  + c3_iter * (np.cos(rad_psi_iter)) ** 3 \
                                  + c4_iter * (np.cos(rad_psi_iter)) ** 4 \
                                  + c5_iter * (np.cos(rad_psi_iter)) ** 5
                rb_torsion_list.append(rb_torsion_iter)


                # RB torsions converted to standard OPLS version via alternate function;
                # the f0 is removed here. It is not a part of the standard OPLS dihedral
                # since it has to be zero for the standard OPLS dihedral to be
                # an exact conversion.
                # If f0/2 is added to the rb_to_opls_iter variable, the RB to OPLS conversion
                # should be analytically correct within the machine precision tolerance.
                rb_to_opls_iter = (1 / 2) * (
                                             + f1_iter * (1 + np.cos(rad_iter))
                                             + f2_iter * (1 - np.cos(2 * rad_iter))
                                             + f3_iter * (1 + np.cos(3 * rad_iter))
                                             + f4_iter * (1 - np.cos(4 * rad_iter))
                                             )
                rb_to_opls_list.append(rb_to_opls_iter)

                abs_differene_opls_rb_iter = abs(rb_to_opls_iter - rb_torsion_iter)
                abs_differene_opls_rb_list.append(abs_differene_opls_rb_iter)

            if count_no == 1:
                title_to_output = "{:35s}{:60s}{:60s}{:35s}{:35s}\n" \
                                  "".format("dihedral_atoms_or_beads",
                                            "rb_constants_c0_c1_c2_c3_c4_c5",
                                            "opls_constants_f0_f1_f2_f3_f4_rounded_6_decimals",
                                            "max_abs_delta_rb_to_opls_calc",
                                            "exact_rb_to_opls_possible"
                                            )

                data_opls_rb.write(title_to_output)

            # Is the iteration RB_to_opls conversion possible
            if bool(np.isclose(c5_iter, 0, atol= error_tolerance_rb_to_opls, rtol=0)) is True  \
                    and bool(np.isclose(f0_iter, 0, atol= error_tolerance_rb_to_opls, rtol=0)) is True:
                rb_to_opls_convertable_iter = True

            else:
                rb_to_opls_convertable_iter = False

            max_abs_diff_new_opls_rb = np.round(np.max(abs_differene_opls_rb_list),
                                                decimals=int(np.log10((1/error_tolerance_rb_to_opls))+2)
                                                )
            if max_abs_diff_new_opls_rb > error_tolerance_rb_to_opls:
                warning_to_output = "WARNING: The {} atoms with the {} RB constants not an exact conversion " \
                                    "for the RB_to_OPLS conversion!. " \
                                    "Max diff = {}, f0 = {}.\n" \
                                    "".format(class_list_iter,
                                              cX_list_iter,
                                              max_abs_diff_new_opls_rb,
                                              f0_iter
                                              )
                text_to_output = "{:35s}{:60s}{:60s}{:35s}{:35s}" \
                                 "\n".format(str(class_list_iter),
                                             str(cX_list_iter),
                                             str(fx_list_iter_6_decimals),
                                             str(max_abs_diff_new_opls_rb),
                                             str(rb_to_opls_convertable_iter)
                                             )

                warn(warning_to_output)
                data_opls_rb.write(text_to_output)

    data_opls_rb.close()


def run_xml_test(xml_and_error_file_dict, error_tolerance_rb_to_opls=1e-4):
    r"""In in file function for the conversion of the Ryckaert-Bellemans (RB) type to OPLS
    type dihedrals (i.e., the mBuild RB_to_OPLS function) for errors and generates ouput
    files with the errors for each xml file provided, containing conversion errors
    and their info.

    .. math::
    RB_{torsions} &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                  &= c_4*cos(psi)^4 + c_5*cos(5*psi)^5

    .. math::
    Exact_OPLS_torsions &= \frac{f_0}{2} + \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                        &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    .. math::
    Standard_ OPLS_torsions &= \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                           &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    where :math:`psi = t - pi = t - 180 degrees`


    Parameters
    ----------
    xml_and_error_file_dict : dict, {str, str} and {xml_file: output_file_name, ...}
        A dictionary that contains the key as the xml_file and the value as the output_file_name.
        The xml_file (key) is the RB torsions xml file paths with file names and
        extensions (i.e., '.xml'), or the just the name of the standard foyer force
        without the extension.
        The output_file_name (value) is the names error output files that will be
        generated.
    error_tolerance_rb_to_opls : float (range 1e-1 to 1e-10), default=1e-4
        The acceptable absolute tolerance between the RB to OPLS conversion.
        This function converts it to an absolute tolerance.
        The working range is 1e-1 to 1e-10, since the calculated values are rounded
        down to the working range plus 2 decimal places (1e-1 -> 3 decimal places and
        2e-10 -> 11 decimal places).

    Returns
    -------
    Generates text (.txt) files containing the dihedral errors in the
    RB to OPLS function conversion process. The number of files generated
    is equal to the number of foyer force field xml files provided.
    These files contain the following columns and have the specified meaning.

    dihedral_atoms_or_beads: list of strings [str, ..., str]
        A list of the foyer atoms or beads types in the dihedral, in order.
    rb_constants_c0_c1_c2_c3_c4_c5: list of floats [float, ..., float]
        A list of the RB-torsion values for the dihedral, which includes
        c0, c1, c2, c3, c4, c5.
    opls_constants_f0_f1_f2_f3_f4_rounded_6_decimals : list of floats [float, ..., float]
        A list of the OPLS dihedral values for the dihedral, which includes
        f0, f1, f2, f3, f4.
        NOTE: The standard OPLS dihedral function does not include
        the f0 constant. However, the f0 is included as it is required for the exact analytical conversion.
        Otherwise, for the standard OPLS dihedral, the f0 term must be 0
        for the solution to be analytically correct.
    max_abs_delta_rb_to_opls_calc : float
        The maximum difference between the RB and standard OPLS dihedral functions
        for 100 equally spaced points over the span of 2 Pi.
    exact_rb_to_opls_possible : bool
        This states whether the conversion is possible from RB to OPLS, where
        OPLS does not include the f0 constant. Therefore, if f0/2 is equal to zero
        within the set tolerance (error_tolerance_rb_to_opls), it will be True.
        Otherwise, it will be False.
    """
    if not isinstance(xml_and_error_file_dict, dict):
        raise TypeError(
            f"The xml_and_error_file_dict variable must be a dict, "
            f"is type {type(xml_and_error_file_dict)}."
        )

    xml_file_list = list(xml_and_error_file_dict.keys())
    output_file_name_list = []
    for val_iter in xml_file_list:
        output_file_name_list.append(xml_and_error_file_dict[val_iter])

    for i_iter in range(0, len(xml_file_list)):
        if not isinstance(xml_file_list[i_iter], str):
            raise TypeError(f"The xml_and_error_file_dict keys are not all strings, "
                            f"and has the type {type(xml_file_list[i_iter])}."
                            )

        if not isinstance(output_file_name_list[i_iter], str):
            raise TypeError(f"The xml_and_error_file_dict values are not all strings, "
                            f"and has the type {type(output_file_name_list[i_iter])}."
                            )

        _test_xml_dihedrals(xml_file_list[i_iter],
                            output_file_name_list[i_iter],
                            error_tolerance_rb_to_opls
                            )


def test_xml_dihedral_rb_to_opls(xml_filename,
                                 error_filename,
                                 non_exact_conversion_list,
                                 error_tolerance_rb_to_opls=1e-4):
    """
    Test that all the xml dihedrals either pass the RB to standard OPLS conversion
    via the mBuild RB_to_OPLS function, or they are in the know list that is not
    a direct analytical conversion.


    .. math::
    RB_{torsions} &= c_0 + c_1*cos(psi) + c_2*cos(psi)^2 + c_3*cos(psi)^3 + \\
                  &= c_4*cos(psi)^4 + c_5*cos(5*psi)^5

    .. math::
    Exact_OPLS_torsions &= \frac{f_0}{2} + \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                        &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    .. math::
    Standard_ OPLS_torsions &= \frac{f_1}{2}*(1+cos(t)) + \frac{f_2}{2}*(1-cos(2*t)) + \\
                           &= \frac{f_3}{2}*(1+cos(3*t)) + \frac{f_4}{2}*(1-cos(4*t))

    where :math:`psi = t - pi = t - 180 degrees`

    Parameters
    ----------
    xml_filename : str
        The RB torsions RB torsions xml file paths with file names and
        extensions (i.e., '.xml'), or just the name of the standard foyer force
        without the extension.
    error_filename : str
        The error_filename is the name error output file that will be read and
        compared to the non_exact_conversion_list data.
    non_exact_conversion_list : list
        This is a list of the known dihedrals that do not analyitically pass
        the RB to standard OPLS conversion via the mBuild RB_to_OPLS function.

        This is a list of the following five data in the list, in order (i.e., [(1), (2), (3), (4), (5)]).
        --- (1) dihedral_atoms_or_beads: list of strings [str, ..., str]
        A list of the foyer atoms or beads types in the dihedral, in order.
        --- (2) rb_constants_c0_c1_c2_c3_c4_c5: list of floats [float, ..., float]
        A list of the RB-torsion values for the dihedral, which includes
        c0, c1, c2, c3, c4, c5.
        --- (3) opls_constants_f0_f1_f2_f3_f4_rounded_6_decimals : list of floats [float, ..., float]
        A list of the OPLS dihedral values for the dihedral, which includes
        f0, f1, f2, f3, f4.
        NOTE: The standard OPLS dihedral function does not include
        the f0 constant. However, the f0 is included as it is required for the exact analytical conversion.
        Otherwise, for the standard OPLS dihedral, the f0 term must be 0
        for the solution to be analytically correct.
        --- (4)  max_abs_delta_rb_to_opls_calc : float
        The maximum difference between the RB and standard OPLS dihedral functions
        for 100 equally spaced points over the span of 2 Pi.
        --- (5)  exact_rb_to_opls_possible : bool
        This states whether the conversion is possible from RB to OPLS, where
        OPLS does not include the f0 constant. Therefore, if f0/2 is equal to zero
        within the set tolerance (error_tolerance_rb_to_opls), it will be True.
        Otherwise, it will be False.
        Example list : [
        [
            ['CH3', 'CH2', 'CH', 'CH3'],
            [3.28629, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
            [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
            2.08734,
            False
        ], ....,
        [
            ['CH2', 'CH2', 'CH', 'CH2'],
            [0.0, 7.44211, 1.85995, -14.67569, 0.0, 0.0],
            [-4.17468, 7.129315, -1.85995, 7.337845, -0.0],
            2.08734,
            True
        ]
        ]
    error_tolerance_rb_to_opls : float (range 1e-1 to 1e-10), default=1e-4
        The acceptable absolute tolerance between the RB to OPLS conversion.
        This function converts it to an absolute tolerance.
        The working range is 1e-1 to 1e-10, since the calculated values are rounded
        down to the working range plus 2 decimal places (1e-1 -> 3 decimal places and
        2e-10 -> 11 decimal places).

    Returns
    -------
    passed_test : bool
        True if all the RB to standard OPLS dihedral conversions via the mBuild RB_to_OPLS
        function are analytical correct or are in the existing list known
        non-analytical correct conversions list (i.e., non_exact_conversion_list).
        False if at least one of the RB to standard OPLS dihedral conversions via the mBuild RB_to_OPLS
        function are not analytical correct or are not in the existing list known
        non-analytical correct conversions list (i.e., non_exact_conversion_list).
    """
    passed_test = True

    non_exact_conversion_list_error_txt = f"The non_exact_conversion_list variables is not formated correctly. "\
                                          f"Please see the test_xml_dihedral_rb_to_opls function for the " \
                                          f"proper format infomation."
    if isinstance(non_exact_conversion_list, list):
        len_non_exact_conversion_list = len(non_exact_conversion_list)
        for error_list_i in non_exact_conversion_list:
            if len(error_list_i) == 5:
                if not isinstance( error_list_i[0], list) \
                        or not isinstance( error_list_i[1], list) \
                        or not isinstance( error_list_i[2], list) \
                        or not isinstance( error_list_i[3], float) \
                        or not isinstance( error_list_i[4], bool):
                    raise TypeError(non_exact_conversion_list_error_txt)
            else:
                raise TypeError(non_exact_conversion_list_error_txt)
    else:
        raise TypeError(non_exact_conversion_list_error_txt)

            # put the non_exact_conversion_list in the same string format with the same
    # spacing as the printed/read files.
    correct_file_strings_list = []
    for j_iter in non_exact_conversion_list:
        correct_file_strings_list.append(
            '{:35s}{:60s}{:60s}{:35s}{:35s}\n'.format(str(j_iter[0]),
                                                      str(j_iter[1]),
                                                      str(j_iter[2]),
                                                      str(j_iter[3]),
                                                      str(j_iter[4])
                                                      )
        )

    run_xml_test({xml_filename: error_filename},
                 error_tolerance_rb_to_opls=error_tolerance_rb_to_opls
                 )

    with open(error_filename, "r") as fp:
        out_gomc = fp.readlines()
        for i_iter, line in enumerate(out_gomc):
            split_line = line.split()

            # check the values for the dihedrals
            if i_iter != 0 and len(split_line) > 0 \
                    and line not in correct_file_strings_list:
                    passed_test = False

    return passed_test