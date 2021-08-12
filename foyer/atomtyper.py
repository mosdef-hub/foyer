"""Determine proper atom types for chemical systems."""
from warnings import warn

import ele
from ele.exceptions import ElementError
from gmso import Topology
from parmed import Structure

from foyer.exceptions import FoyerError
from foyer.smarts_graph import SMARTSGraph
from foyer.topology_graph import TopologyGraph


def find_atomtypes(structure, forcefield, max_iter=10):
    """Determine atomtypes for all atoms.

    Parameters
    ----------
    structure : parmed.Structure, or gmso.Topology, or TopologyGraph
        The topology that we are trying to atomtype. If a parmed.Structure or
        gmso.Topology is provided, it will be convert to a TopologyGraph before
        atomtyping.
    forcefield : foyer.Forcefield
        The forcefield object.
    max_iter : int, optional, default=10
        The maximum number of iterations.

    """
    topology_graph = structure

    if isinstance(structure, Structure):
        topology_graph = TopologyGraph.from_parmed(structure)
    elif isinstance(structure, Topology):
        topology_graph = TopologyGraph.from_gmso_topology(structure)

    typemap = {
        atom_index: {"whitelist": set(), "blacklist": set(), "atomtype": None}
        for atom_index in topology_graph.atoms(data=False)
    }

    rules = _load_rules(forcefield, typemap)

    # Only consider rules for elements found in topology
    subrules = dict()

    system_elements = set()
    for _, atom_data in topology_graph.atoms(data=True):
        # First add non-element types, which are strings, then elements
        name = atom_data.name
        if name.startswith("_"):
            if name in forcefield.non_element_types:
                system_elements.add(name)
        else:
            atomic_number = atom_data.atomic_number
            atomic_symbol = atom_data.element
            try:
                element_from_num = ele.element_from_atomic_number(
                    atomic_number
                ).symbol
                element_from_sym = ele.element_from_symbol(atomic_symbol).symbol
                assert element_from_num == element_from_sym
                system_elements.add(element_from_num)
            except ElementError:
                raise FoyerError(
                    "Parsed atom {} as having neither an element "
                    "nor non-element type.".format(name)
                )
            except AssertionError:
                raise FoyerError(
                    f"Parsed atom {name} has mismatching atom number ({atomic_number}) "
                    f"and atom symbol ({atomic_symbol})."
                )

    for key, val in rules.items():
        atom = val.nodes[0]["atom"]
        if len(list(atom.find_data("atom_symbol"))) == 1 and not list(
            atom.find_data("not_expression")
        ):
            try:
                element = next(atom.find_data("atom_symbol")).children[0]
            except IndexError:
                try:
                    atomic_num = next(atom.find_data("atomic_num")).children[0]
                    element = ele.element_from_atomic_number(atomic_num).symbol
                except IndexError:
                    element = None
        else:
            element = None
        if element is None or element in system_elements:
            subrules[key] = val
    rules = subrules

    _iterate_rules(rules, topology_graph, typemap, max_iter=max_iter)
    _resolve_atomtypes(topology_graph, typemap)

    return typemap


def _load_rules(forcefield, typemap):
    """Load atomtyping rules from a forcefield into SMARTSGraphs."""
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
        rules[rule_name] = SMARTSGraph(
            smarts_string=smarts,
            parser=forcefield.parser,
            name=rule_name,
            overrides=overrides,
            typemap=typemap,
        )
    return rules


def _iterate_rules(rules, topology_graph, typemap, max_iter):
    """Iterate through all the rules until the white- and blacklists converge.

    Parameters
    ----------
    rules : dict
        A dictionary mapping rule names (typically atomtype names) to
        SMARTSGraphs that evaluate those rules.
    topology_graph : TopologyGraph
        The topology graph that we are trying to atomtype.
    max_iter : int
        The maximum number of iterations.

    """
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for rule in rules.values():
            for match_index in rule.find_matches(topology_graph, typemap):
                atom = typemap[match_index]
                # This conditional is not strictly necessary, but it prevents
                # redundant set addition on later iterations
                if rule.name not in atom["whitelist"]:
                    atom["whitelist"].add(rule.name)
                    atom["blacklist"] |= rule.overrides
                    found_something = True
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")
    return typemap


def _resolve_atomtypes(topology_graph, typemap):
    """Determine the final atomtypes from the white- and blacklists."""
    atoms = {
        atom_idx: data for atom_idx, data in topology_graph.atoms(data=True)
    }
    for atom_id, atom in typemap.items():
        atomtype = [
            rule_name for rule_name in atom["whitelist"] - atom["blacklist"]
        ]
        if len(atomtype) == 1:
            atom["atomtype"] = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError(
                "Found multiple types for atom {} ({}): {}.".format(
                    atom_id, atoms[atom_id].atomic_number, atomtype
                )
            )
        else:

            raise FoyerError(
                "Found no types for atom {} ({}).".format(
                    atom_id, atoms[atom_id].atomic_number
                )
            )
