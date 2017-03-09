from warnings import warn

from oset import oset as OrderedSet

from foyer.exceptions import FoyerError
from foyer.smarts_graph import SMARTSGraph


def find_atomtypes(topology, forcefield, max_iter=10):
    """Determine atomtypes for all atoms.

    Parameters
    ----------
    topology : simtk.openmm.app.Topology
        The topology that we are trying to atomtype.
    forcefield : foyer.Forcefield
        The forcefield object.
    max_iter : int, optional, default=10
        The maximum number of iterations.

    """

    for atom in topology.atoms():
        atom.whitelist = OrderedSet()
        atom.blacklist = OrderedSet()

    rules = _load_rules(forcefield)
    _iterate_rules(rules, topology, max_iter=max_iter)
    _resolve_atomtypes(topology)


def _load_rules(forcefield):
    """Load atomtyping rules from a forcefield into SMARTSGraphs. """
    rules = dict()
    for rule_name, smarts in forcefield.atomTypeDefinitions.items():
        overrides = forcefield.atomTypeOverrides.get(rule_name)
        rules[rule_name] = SMARTSGraph(smarts_string=smarts,
                                       parser=forcefield.parser,
                                       name=rule_name,
                                       overrides=set(overrides))
    return rules


def _iterate_rules(rules, topology, max_iter):
    """Iteratively run all the rules until the white- and backlists converge.

    Parameters
    ----------
    rules : dict
        A dictionary mapping rule names (typically atomtype names) to
        SMARTSGraphs that evaluate those rules.
    topology : simtk.openmm.app.Topology
        The topology that we are trying to atomtype.
    max_iter : int
        The maximum number of iterations.

    """
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for rule in rules.values():
            matches = rule.find_matches(topology)
            assert False
            if rule.name not in atom.whitelist:
                if rule.matches(atom):
                    atom.whitelist.add(rule.name)
                    atom.blacklist |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")


def _resolve_atomtypes(topology):
    """Determine the final atomtypes from the white- and blacklists. """
    for i, atom in enumerate(topology.atoms()):
        atomtype = [rule_name for rule_name in atom.whitelist - atom.blacklist]

        if len(atomtype) == 1:
            atom.id = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError("Found multiple types for atom {} ({}): {}.".format(
                i, atom.element.name, atomtype))
        else:
            raise FoyerError("Found no types for atom {} ({}).".format(
                i, atom.element.name))
