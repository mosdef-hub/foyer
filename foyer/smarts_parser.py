
if __name__ == '__main__':

    from topos import Topos
    topos = Topos()

    ethanol = topos.load_topo("64-17-5")

    from forcefield import load

    oplsaa = load('ff/ff.xml')

    from atomtyper import find_atomtypes
    find_atomtypes(ethanol.atoms, oplsaa)

    for atom in ethanol.atoms:
        print("{}: whitelist: {}, blacklist: {}".format(atom, atom.whitelist, atom.blacklist))