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
