import os

import numpy as np
import parmed as pmd
from pkg_resources import resource_filename
from foyer.atomtyper import find_atomtypes
from foyer.forcefield import load

from foyer.smarts import Smarts


if __name__ == '__main__':

    resource_dir = resource_filename('foyer', '../opls_validation')

    top_path = os.path.join(resource_dir, '64-17-5.top')
    gro_path = os.path.join(resource_dir, '64-17-5-gas.gro')

    structure = pmd.gromacs.GromacsTopologyFile(top_path, xyz=gro_path,
                                                parametrize=False)
    structure.title = structure.title.replace(' GAS', '')
    known_opls_types = [atom.type for atom in structure.atoms]

    forcefield_resource_dir = resource_filename('foyer', 'oplsaa')

    forcefield = load(os.path.join(forcefield_resource_dir, 'oplsaa.xml'))

    find_atomtypes(structure.atoms, forcefield, debug=False)

    for atom in structure.atoms:
        print("Atom {} is {}".format(atom, atom.type))