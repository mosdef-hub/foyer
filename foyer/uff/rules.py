from foyer.atomtyper import (
    Element, NeighborCount, NeighborsExactly, Whitelist)


# -------------- #
# House of rules #
# -------------- #

@Element('C')
@NeighborCount(4)
@NeighborsExactly('H', 4)
@Whitelist('C_3')
def uff_C_3(atom):
    """ """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist('H_')
def uff_H_(atom):
    """ """
    return True


if __name__ == "__main__":

    from foyer.atomtyper import find_atomtypes
    from foyer.forcefield import prepare_atoms

    from mbuild.examples.methane.methane import Methane

    m = Methane()
    # m = Ethane()

    traj = m.to_trajectory()
    prepare_atoms(traj.top)
    find_atomtypes(traj.top._atoms, forcefield='UFF')

    for atom in traj.top._atoms:
        print("Atom name={}, opls_type={}".format(atom.name, atom.atomtype))
