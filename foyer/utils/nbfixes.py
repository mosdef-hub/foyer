"""Support non-bonded fixes for various interactions."""

from parmed import Structure


def apply_nbfix(struct, atom_type1, atom_type2, sigma, epsilon):
    """Apply a single nbfix to a particular interaction.

    Parameters
    ----------
    struct : parmed.structure.Structure
        The ParmEd structure to which this nbfix will be applied.
    atom_type1 : str
        The name of the first atom type in the nbfix pair
    atom_type2 : str
        The name of the second atom type in the nbfix pair
    sigma : float
        The sigma of the cross-interaction in kcal/mol
    epsilon : float
        The epsilon of the cross-interaction in Angstroms

    Returns
    -------
    struct : parmed.structure.Structure
        The input structure with the nbfix applied.
    """
    struct_copy = struct.copy(cls=Structure, split_dihedrals=True)

    atypes_name = set(a.atom_type.name for a in struct_copy.atoms)
    if atom_type1 not in atypes_name or atom_type2 not in atypes_name:
        raise ValueError(
            "Atom types {} and {} not found " "in structure.".format(
                atom_type1, atom_type2
            )
        )

    # Calculate rmin from sigma because parmed uses it internally
    rmin = sigma * 2 ** (1.0 / 6.0)

    atom_types = list(a.atom_type for a in struct_copy.atoms)
    for atom_type in sorted(atom_types, key=lambda a: a.name):
        if atom_type.name == atom_type1:
            atom_type.add_nbfix(atom_type2, rmin, epsilon)
        if atom_type.name == atom_type2:
            atom_type.add_nbfix(atom_type1, rmin, epsilon)

    return struct_copy
