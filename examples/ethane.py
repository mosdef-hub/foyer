import parmed as pmd
from foyer import Forcefield
from foyer.tests.utils import get_fn

if __name__ == '__main__':
    mol2_path = get_fn('ethane.mol2')
    untyped_ethane = pmd.load_file(mol2_path, structure=True)
    oplsaa = Forcefield(name='oplsaa')
    ethane = oplsaa.apply(untyped_ethane, references_file='ethane.bib')

    print("Atoms:")
    for atom in ethane.atoms:
        print('Atom {} is typed as {}'.format(atom, atom.type))

    print("Bonds:")
    for bond in ethane.bonds:
        print('{} '.format(bond))

    print("Angles:")
    for angle in ethane.angles:
        print('{} '.format(angle))

    print("Dihedrals:")
    for dihedral in ethane.dihedrals:
        print('{} '.format(dihedral))

    # Save to GROMACS
    ethane.save('ethane.gro')
    ethane.save('ethane.top')

    # Within the `Forcefield.apply` method, an intermediate OpenMM system is
    # created. If you wish to use OpenMM, e.g. for use of potential forms not
    # yet supported by ParmEd, you can simply stop the conversion process
    # after the OpenMM System creation by directly invoking the internal
    # method calls:
    from foyer.forcefield import generate_topology
    omm_topology, positions = generate_topology(untyped_ethane)
    omm_system = oplsaa.createSystem(topology=omm_topology)
