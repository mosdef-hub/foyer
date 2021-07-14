Applying a force field
======================

The main method you will use in **foyer** is the
``Forcefield.apply()`` method. There are a few important arguments
you should understand.

The first several are the ``assert_bond_params``
``assert_angle_params``, ``assert_dihedral_params``, and
``assert_improper_params``. These arguments require that the
supplied force field has parameters for every bond, angle, dihedral,
and improper in the system. In most cases, if you get an error,
it means that your force field is missing parameters for one of the
bonds/angles/dihedrals/impropers in the system. This could be
because the parameters are missing or because the atom-typing
(i.e., the SMARTS strings) are incorrect. These arguments are ``True`` by default,
with the exception of ``assert_improper_params``. In all cases, it
is wise to verify that the final files you generate have the expected
number of bonds/angles/dihedrals/impropers for your system.

The other important optional argument is the ``combining_rule`` option,
which is ``"lorentz"`` (Lorentz-Berthelot) by default. The other valid
option is ``"geometric"``, if your force field uses geometric combining rules.

