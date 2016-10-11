import os
from foyer.forcefield import Forcefield
import parmed as pmd

if __name__ == '__main__':
    # load the methane.mol2 file
    mol2_path = os.path.join(os.path.dirname(__file__), 'ethane.mol2')
    ethane = pmd.load_file(mol2_path, structure=True)

    # apply forcefield
    Forcefield.by_name('oplsaa').apply(ethane)

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