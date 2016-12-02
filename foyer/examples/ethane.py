import os
from foyer.forcefield import Forcefield, generate_topology
import parmed as pmd


if __name__ == '__main__':
    mol2_path = os.path.join(os.path.dirname(__file__), 'ethane.mol2')
    ethane = pmd.load_file(mol2_path, structure=True)

    oplsaa = Forcefield.by_name('oplsaa')

    ethane = oplsaa.apply(ethane)

    # topology = generate_topology(ethane)
    # omm_system = oplsaa.createSystem(topology=topology)
    #
    # ethane = pmd.openmm.load_topology(topology=topology, system=omm_system)
    # import ipdb; ipdb.set_trace()

    # print atom types
    print("Atoms:")
    for atom in ethane.atoms:
        print('Atom {} is typed as {}'.format(atom, atom.type))

    # print bonds
    print("Bonds:")
    for bond in ethane.bonds:
        print('{} '.format(bond))

    # print angles
    print("Angles:")
    for angle in ethane.angles:
        print('{} '.format(angle))    # print angles

    # print dihedrals
    print("Dihedrals:")
    for dihedral in ethane.dihedrals:
        print('{} '.format(dihedral))

