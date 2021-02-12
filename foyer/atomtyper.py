from warnings import warn

import gmso
import parmed
import parmed.periodic_table as pt
from foyer.exceptions import FoyerError
from foyer.smarts_graph import SMARTSGraph

def find_atomtypes(topology, forcefield, max_iter=10):
    """Determine atomtypes for all atoms.

    Parameters
    ----------
    topology : parmed.Structure, or gmso.Topology
        The topology that we are trying to atomtype, must matched with backend.
    forcefield : foyer.Forcefield
        The forcefield object.
    max_iter : int, optional, default=10
        The maximum number of iterations.

    """
    # Handle Parmed Structure
    if isinstance(topology, parmed.Structure):
        typemap = {atom.idx: {'whitelist': set(), 'blacklist': set(),
            'atomtype': None} for atom in topology.atoms}

        system_elements = set()
        for a in topology.atoms:
            # First add non-element types, which are strings, then elements
            if a.name.startswith('_'):
                if a.name in forcefield.non_element_types:
                    system_elements.add(a.name)
            else:
                if 0 < a.atomic_number <= pt.KNOWN_ELEMENTS:
                    element = pt.Element[a.atomic_number]
                    system_elements.add(element)
                else:
                    raise FoyerError(
                        'Parsed atom {} as having neither an element '
                        'nor non-element type.'.format(a)
                    )

    # Handle GMSO Topology
    elif isinstance(topology, gmso.Topology):
        typemap = {topology.get_index(site): {'whitelist': set(), 'blacklist': set(),
            'atomtype': None} for site in topology.sites}

        system_elements = set()
        for site in topology.sites:
            # First add non-element types, which are strings, then elements
            if site.element:
                element = site.element.symbol
                system_elements.add(element)
            else:
                #if site.name in forcefield.non_element_types:
                if site.name:
                    system_elements.add(site.name)
                else:
                    raise FoyerError(
                        f'Parsed atom {site} as having neither an element '
                        'nor non-element type.'
                    )
    else:
        raise TypeError("Passed {}, acceptable objects are parmed structure and gmso topologies".format(
            type(topology)
            )
        )
    rules = _load_rules(forcefield, typemap)

    # Only consider rules for elements found in topology
    subrules = dict()
    for key, val in rules.items():
        atom = val.nodes[0]['atom']
        if len(list(atom.find_data('atom_symbol'))) == 1 and \
                    not list(atom.find_data('not_expression')):
            try:
                element = next(atom.find_data('atom_symbol')).children[0]
            except IndexError:
                try:
                    atomic_num = next(atom.find_data('atomic_num')).children[0]
                    element = pt.Element[atomic_num]
                except IndexError:
                    element = None
        else:
            element = None
        if element is None or element in system_elements:
            subrules[key] = val
    rules = subrules
    typemap = _iterate_rules(rules, topology, typemap, max_iter=max_iter)
    typemap = _resolve_atomtypes(topology, typemap)

    return typemap

def _load_rules(forcefield, typemap):
    """Load atomtyping rules from a forcefield into SMARTSGraphs. """
    rules = dict()
    # For every SMARTS string in the force field,
    # create a SMARTSGraph object
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
    """Iteratively run all the rules until the white- and blacklists converge.

    Parameters
    ----------
    rules : dict
        A dictionary mapping rule names (typically atomtype names) to
        SMARTSGraphs that evaluate those rules.
    topology : parmed.Structure or gmso.Topology
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
                # This conditional is not strictly necessary, but it prevents
                # redundant set addition on later iterations
                if rule.name not in atom['whitelist']:
                    atom['whitelist'].add(rule.name)
                    atom['blacklist'] |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")
    return typemap

def _resolve_atomtypes(structure, typemap):
    """Determine the final atomtypes from the white- and blacklists. """
    for atom_id, atom in typemap.items():
        atomtype = [rule_name for rule_name in
                    atom['whitelist'] - atom['blacklist']]
        if len(atomtype) == 1:
            atom['atomtype'] = atomtype[0]
        elif isinstance(structure, gmso.Topology):
            if len(atomtype) > 1:
                raise FoyerError("Found multiple types for Site {} ({}): {}.".format(
                    atom_id, structure.sites[atom_id].element, atomtype))
            else:
                raise FoyerError("Found no types for Site {} ({}).".format(
                    atom_id, structure.sites[atom_id]))
        elif isinstance(structure, parmed.Structure):
            if len(atomtype) > 1:
                raise FoyerError("Found multiple types for atom {} ({}): {}.".format(
                    atom_id, structure.atoms[atom_id].atomic_number, atomtype))
            else:
                raise FoyerError("Found no types for atom {} ({}).".format(
                    atom_id, structure.atoms[atom_id].atomic_number))

    return typemap
