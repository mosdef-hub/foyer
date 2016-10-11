import os
from foyer.forcefield import Forcefield
import parmed as pmd

if __name__ == '__main__':
    # load the methane.mol2 file
    mol2_path = os.path.join(os.path.dirname(__file__), 'methane.mol2')
    methane = pmd.load_file(mol2_path, structure=True)

    # apply forcefield
    Forcefield.by_name('oplsaa').apply(methane)

    # print atom types
    for atom in methane.atoms:
        print('Atom {} is typed as {}'.format(atom, atom.type))
