from __future__ import print_function

import os

import mdtraj as md


def load_top_opls(toppath):
    """Load a gromacs .top file parameterized with OPLS types. """
    split_path = os.path.split(toppath)
    filename = split_path[-1]
    pdb_file = "{}-gas.pdb".format(filename[:-4])
    pdb_path = os.path.join(split_path[0], pdb_file)
    traj = md.load(pdb_path)

    opls_types = list()
    with open(toppath, 'r') as top:
        current = None
        for line in top:
            if ';' in line:
                line = line[:line.index(';')]
            stripped = line.strip()
            if stripped.startswith('*') or len(stripped) == 0:
                continue
            elif stripped.startswith('['):
                current = stripped[1:-1].strip()
            elif current == 'moleculetype':
                mol_name = stripped.split()[0]
            elif current == 'atoms':
                opls_string = stripped.split()[1]
                if 'opls' not in opls_string:
                    # print("Found non 'opls' type in {} ({}): {}.".format(
                    #     mol_name, filename, opls_string))
                    # print("Ignoring file. Need to come up with proper support.\n")
                    return
                #opls_types.append(opls_string.split('_')[1])
                opls_types.append(opls_string)
            elif current == 'bonds':
                atom1 = traj.topology.atom(int(stripped.split()[0]) - 1)
                atom2 = traj.topology.atom(int(stripped.split()[1]) - 1)
                traj.topology.add_bond(atom1, atom2)
    return traj.topology, opls_types, mol_name