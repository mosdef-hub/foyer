import os

import parmed as pmd
from pkg_resources import resource_filename
from foyer import Forcefield



def test_ethanol():
    resource_dir = resource_filename('foyer', '../opls_validation')

    top_path = os.path.join(resource_dir, '64-17-5.top')
    gro_path = os.path.join(resource_dir, '64-17-5-gas.gro')

    ethanol = pmd.gromacs.GromacsTopologyFile(top_path, xyz=gro_path,
                                              parametrize=False)
    ethanol.title = ethanol.title.replace(' GAS', '')

    known_opls_types = [atom.type for atom in ethanol.atoms]

    oplsaa = Forcefield.by_name('oplsaa')

    typed_ethanol = oplsaa.apply(ethanol)


    for atom in typed_ethanol:
        print("Atom {} is {}".format(atom, atom.type))