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
    typemap = {atom.index: {'whitelist': set(), 'blacklist': set(), 'atomtype': None} for atom in topology.atoms()}

    rules = _load_rules(forcefield, typemap)

    # Only consider rules for elements found in topology
    subrules = dict()
    system_elements = {a.element.symbol for a in topology.atoms()}
    for key,val in rules.items():
        atom = val.node[0]['atom']
        if len(list(atom.find_data('atom_symbol'))) == 1 and not list(atom.find_data('not_expression')):
            try:
                element = next(atom.find_data('atom_symbol')).children[0]
            except IndexError:
                try:
                    atomic_num = next(atom.find_data('atomic_num')).children[0]
                    element = pt.Element[int(atomic_num)]
                except IndexError:
                    element = None
        else:
            element = None
        if element is None or element in system_elements:
            subrules[key] = val
    rules = subrules

    _iterate_rules(rules, topology, typemap, max_iter=max_iter)
    _resolve_atomtypes(topology, typemap)

    return typemap

def _load_rules(forcefield, typemap):
    """Load atomtyping rules from a forcefield into SMARTSGraphs. """
    rules = dict()
    for rule_name, smarts in forcefield.atomTypeDefinitions.items():
        if not smarts:  # We want to skip over empty smarts definitions
            continue
        overrides = forcefield.atomTypeOverrides.get(rule_name)
        if overrides is not None:
            overrides = set(overrides)
        else:
            overrides = set()
        rules[rule_name] = SMARTSGraph(smarts_string=smarts,
                                       parser=forcefield.parser,
                                       name=rule_name,
                                       overrides=overrides,
                                       typemap=typemap)
    return rules


def _iterate_rules(rules, topology, typemap, max_iter):
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
            for match_index in rule.find_matches(topology, typemap):
                atom = typemap[match_index]
                if rule.name not in atom['whitelist']:
                    atom['whitelist'].add(rule.name)
                    atom['blacklist'] |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")

    return typemap

def _resolve_atomtypes(topology, typemap):
    """Determine the final atomtypes from the white- and blacklists. """
    atoms = list(topology.atoms())
    for atom_id, atom in typemap.items():
        atomtype = [rule_name for rule_name in atom['whitelist'] - atom['blacklist']]
        if len(atomtype) == 1:
            atom['atomtype'] = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError("Found multiple types for atom {} ({}): {}.".format(
                atom_id, atoms[atom_id].element.name, atomtype))
        else:
            raise FoyerError("Found no types for atom {} ({}).".format(
                atom_id, atoms[atom_id].element.name))
