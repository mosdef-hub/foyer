import os

from intermol.gromacs.gromacs_parser import load_gromacs


def load_top_opls(top_path, mol_name=None):
    """Load a gromacs .top file parameterized with OPLS types. """
    if mol_name is not None:
        with open(top_path) as top_file:
            if mol_name not in top_file.read():
                return None, None, None

    split_path = os.path.split(top_path)
    filename = split_path[-1]
    gro_file = '{}-gas.gro'.format(filename[:-4])
    gro_path = os.path.join(split_path[0], gro_file)

    system = load_gromacs(top_path, gro_path)
    opls_types = [atom.atomtype[0] for atom in system.atoms]
    mol_name = [name for name in system.molecule_types][0]

    return system, opls_types, mol_name
