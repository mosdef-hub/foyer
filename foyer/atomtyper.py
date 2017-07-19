from warnings import warn

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
    rules = _load_rules(forcefield)

    # Only consider rules for elements found in topology
    subrules = dict()
    system_elements = {a.element.symbol for a in topology.atoms()}
    for key,val in rules.items():
        atom = val.node[0]['atom']
        if len(atom.select('atom_symbol')) == 1 and not atom.select('not_expression'):
            try:
                element = atom.select('atom_symbol').strees[0].tail[0]
            except IndexError:
                try:
                    atomic_num = atom.select('atomic_num').strees[0].tail[0]
                    element = pt.Element[int(atomic_num)]
                except IndexError:
                    element = None
        else:
            element = None
        if element is None or element in system_elements:
            subrules[key] = val
    rules = subrules

    _iterate_rules(rules, topology, max_iter=max_iter)
    _resolve_atomtypes(topology)


def _load_rules(forcefield):
    """Load atomtyping rules from a forcefield into SMARTSGraphs. """
    rules = dict()
    for rule_name, smarts in forcefield.atomTypeDefinitions.items():
        overrides = forcefield.atomTypeOverrides.get(rule_name)
        if overrides is not None:
            overrides = set(overrides)
        else:
            overrides = set()
        rules[rule_name] = SMARTSGraph(smarts_string=smarts,
                                       parser=forcefield.parser,
                                       name=rule_name,
                                       overrides=overrides)
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
    atoms = list(topology.atoms())
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for rule in rules.values():
            for match_index in rule.find_matches(topology):
                atom = atoms[match_index]
                if rule.name not in atom.whitelist:
                    atom.whitelist.add(rule.name)
                    atom.blacklist |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")


def _resolve_atomtypes(topology):
    """Determine the final atomtypes from the white- and blacklists. """
    for atom in topology.atoms():
        atomtype = [rule_name for rule_name in atom.whitelist - atom.blacklist]
        if len(atomtype) == 1:
            atom.id = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError("Found multiple types for atom {} ({}): {}.".format(
                atom.index, atom.element.name, atomtype))
        else:
            raise FoyerError("Found no types for atom {} ({}).".format(
                atom.index, atom.element.name))
