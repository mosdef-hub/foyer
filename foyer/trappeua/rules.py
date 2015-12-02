from foyer.atomtyper import (
    Element, NeighborCount, NeighborsAtLeast, NeighborsExactly, Whitelist,
    Blacklist, check_atom, InWhitelist)
from foyer.chemical_groups import benzene, dioxolane13


# -------------- #
# House of rules #
# -------------- #

@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 4)
@Whitelist('C_C')
def trappeua_C_C(atom):
    return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 3)
@Whitelist('CH_C')
def trappeua_CH_C(atom):
    return True


@Element('C')
@NeighborCount(2)
@NeighborsExactly('C', 2)
@Whitelist('CH2_C')
def trappeua_CH2_C(atom):
    return True


@Element('C')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist('CH3_C')
def trappeua_CH3(atom):
    return True


@Element('C')
@NeighborCount(2)
@NeighborsExactly('C', 1)
@NeighborsExactly('O', 1)
@Whitelist('CH2_O')
def trappeua_CH2_O(atom):
    return True


@Element('O')
@NeighborCount(2)
@NeighborsExactly('C', 2)
@Whitelist('O_C')
def trappeua_O_C(atom):
    return True


@Element('O')
@NeighborCount(2)
@NeighborsExactly('C', 2)
@Whitelist('O_est')
@Blacklist('O_C')
def trappeua_O_est(atom):
    for neighbor in atom.neighbors:
        if check_atom(neighbor, 'C2_est'):
            return True
    return False

@Element('O')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist('O2_est')
def trappeua_O2_est(atom):
    return True

@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 1)
@NeighborsExactly('O', 2)
@Whitelist('C2_est')
def trappeua_C2_est(atom):
    return True



if __name__ == "__main__":
    import mbuild as mb
    from foyer.atomtyper import find_atomtypes
    from foyer.forcefield import prepare_atoms

    # m = Methane()
    # m = Ethane()
    pass
