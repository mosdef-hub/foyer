from warnings import warn

from oset import oset as OrderedSet

from foyer.exceptions import FoyerError
from foyer.rule import Rule

RULE_NAME_TO_RULE = dict()


def find_atomtypes(atoms, forcefield, debug=False):
    """Determine atomtypes for all atoms.

    This function is fairly general in that it can function on any list of atom
    objects as long as they have a property `neighbors`, which is a list of
    other atoms that they are bonded to as well as the attributes `whitelist`
    and `blacklist` which are sets - if they are ordered sets it simplifies
    debugging a little bit.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for.
    forcefield : simtk.openmm.app.Forcefield object
        The forcefield object.
    debug : bool, default=True
        Provides debug information about the logical consistency of the
        atomtyping rules.

    See also
    --------
    _sanitize

    """

    for atom in atoms:
        atom.whitelist = OrderedSet()
        atom.blacklist = OrderedSet()
        # if atom.element:
        #     atom.element_name = pt.Element[atom.element]
        # else:
        #     # TODO: more robust element detection
        #     atom.element_name = atom.name

    _load_rules(forcefield)
    # if debug:
    #     _sanitize()
    _iterate_rules(atoms, max_iter=10)
    _resolve_atomtypes(atoms)

def _load_rules(forcefield):
    global RULE_NAME_TO_RULE
    RULE_NAME_TO_RULE = dict()
    all_names = forcefield._atomTypeDefinitions.keys()
    for rule_name, smarts in forcefield._atomTypeDefinitions.items():
        overrides = forcefield._atomTypeOverrides.get(rule_name)
        RULE_NAME_TO_RULE[rule_name] = Rule(rule_name, forcefield.parser, smarts, overrides=overrides)

# TODO: enable this for speedup
# def _build_rule_map():
#     """Build up a tree of element types-->neighbor counts-->rules. """
#     for rule_name, smarts in RULE_NAME_TO_RULE.items():
#         element_type = smarts.anchor
#         neighbor_count = smarts.n_neighbors
#
#         if element_type not in RULE_MAP:
#             RULE_MAP[element_type] = dict()
#         if neighbor_count not in RULE_MAP[element_type]:
#             RULE_MAP[element_type][neighbor_count] = []
#         RULE_MAP[element_type][neighbor_count].append(rule_name)


def _iterate_rules(atoms, max_iter=10):
    """Iteratively run all the rules until the white- and backlists converge.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for.
    max_iter : int, optional, default=10
        The maximum number of iterations.

    """
    for _ in range(max_iter):
        max_iter -= 1
        found_something = False
        for atom in atoms:
            for rule in RULE_NAME_TO_RULE.values():
                if rule.name not in atom.whitelist:
                    # # Only run rules with matching element and neighbor counts.
                    # if atom.element_name in RULE_MAP:
                    #     if len(atom.bond_partners) in RULE_MAP[atom.element_name]:
                    #         for rule in RULE_MAP[atom.element_name][len(atom.bond_partners)]:
                    if rule.matches(atom):
                        atom.whitelist.add(rule.name)
                        atom.blacklist |= rule.overrides
                        found_something = True
                    #     else:
                    #         warn("No rule for {}-neighbor '{}' atom".format(
                    #             len(atom.bond_partners), atom.element_name))
                    # else:
                    #     warn("No rule for atom name '{}'".format(atom.element_name))
        if not found_something:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")

def _resolve_atomtypes(atoms):
    """Determine the final atomtypes from the white- and blacklists."""
    for i, atom in enumerate(atoms):
        atomtype = [rule_name for rule_name in atom.whitelist - atom.blacklist]

        if len(atomtype) == 1:
            atom.id = atomtype[0]
        elif len(atomtype) > 1:
            raise FoyerError("Found multiple types for atom {0} ({1}): {2}.".format(i, atom.element.name, atomtype))
        else:
            raise FoyerError("Found no types for atom {0} ({1}).".format(i, atom.element.name))

# def neighbor_element_types(atom):
#     """Returns the number of neighbors of each element type for an `atom`.
#
#     The dict maintained is `NEIGHBOR_TYPES_MAP` and is organized as follows:
#         atom: defaultdict{element: number of neighbors of that element type}
#     E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
#         Atom: {'C': 3, 'H': 1}
#
#     If the queried `atom` is not already in `NEIGHBOR_TYPES_MAP`, its entry will
#     be added.
#
#     """
#     if atom not in NEIGHBOR_TYPES_MAP:
#         neighbors = defaultdict(int)
#         for neighbor in atom.bond_partners:
#             name = neighbor.element_name
#             neighbors[name] += 1
#         NEIGHBOR_TYPES_MAP[atom] = neighbors
#     return NEIGHBOR_TYPES_MAP[atom]
#

# def neighbor_whitelist_types(atom):
#     """Returns the number of neighbors of each type of whitelisted rule.
#
#     The dict maintained is `NEIGHBOR_WHITELIST_MAP` and is organized as follows:
#         atom: defaultdict{rule id: number of neighbors with rule id whitelisted}
#     E.g. for a carbon atom in the middle of an alkane chain:
#         Atom: {'136': 2, '140': 2}
#
#     If the queried `atom` is not already in `NEIGHBOR_WHITELIST_MAP`, an empty
#     dict is returned. This behavior differs from `neighbor_element_types()`
#     because whitelisted rules are not known a priori. New entries are only
#     added when the @Whitelist decorator is used.
#
#     See Also `increment_neighbor_whitelist()`.
#     """
#     if atom in NEIGHBOR_WHITELIST_MAP:
#         return NEIGHBOR_WHITELIST_MAP[atom]
#     else:
#         return dict()


# def increment_neighbor_whitelist(atom, whitelist_type):
#     """Increment the counts in an atom's `NEIGHBOR_WHITELIST_MAP` entry. """
#     for neighbor in atom.bond_partners:
#         if neighbor not in NEIGHBOR_WHITELIST_MAP:
#             NEIGHBOR_WHITELIST_MAP[neighbor] = defaultdict(int)
#         NEIGHBOR_WHITELIST_MAP[neighbor][whitelist_type] += 1
#

# ------------------------------------- #
# Sanitization and associated functions #
# ------------------------------------- #


# def _sanitize():
#     """Analyze all rules for possible inconsistencies.
#
#     This function serves primarily as a tool for developers who intend to add
#     new rules or modify existing ones. Ideally, it will help you identify and
#     correct logical inconsistencies as early as possible. Additionally, it
#     suggests other rules that you may want to consider blacklisting.
#
#     """
#     supported_elements = find_all_supported_elements()
#     rule_matches = find_all_rule_matches(supported_elements)
#
#     # Build directed graphs showing which rules blacklist each other.
#     for key, rules in rule_matches.items():
#         # Only consider patterns matched by multiple rules.
#         if len(rules) < 2:
#             continue
#
#         element_type, pattern = key
#         graph = nx.DiGraph()
#         for rule_number in rules:
#             graph.add_node(rule_number)
#             blacklisted_rules = set()
#             decorators = get_decorators(RULE_NAME_TO_RULE[rule_number])
#
#             for dec in decorators:
#                 if isinstance(dec, Blacklist):
#                     blacklisted_rules.update(dec.rule_numbers)
#             for blacklisted_rule in blacklisted_rules:
#                 graph.add_edge(rule_number, blacklisted_rule)
#
#         if not nx.is_connected(graph.to_undirected()):
#             draw_rule_graph('unconnected', graph, element_type, pattern)
#
#         if not nx.is_directed_acyclic_graph(graph):
#             draw_rule_graph('not_DAG', graph, element_type, pattern)
#
#         # Check for multiple sinks. This is not necessarily incorrect.
#         sinks = []
#         for node in graph.nodes():
#             if len(nx.descendants(graph, node)) == 0:
#                 sinks.append(node)
#         if len(sinks) > 1:
#             draw_rule_graph('multiple_sinks', graph, element_type, pattern, sinks=sinks)
#
#         # Check for multiple sources. This is not necessarily incorrect.
#         sources = []
#         for node in graph.nodes():
#             if len(nx.ancestors(graph, node)) == 0:
#                 sources.append(node)
#         if len(sources) > 1:
#             draw_rule_graph('multiple_sources', graph, element_type, pattern, sources=sources)
#
#
# def find_all_supported_elements():
#     """Find all elements currently supported by rules.
#
#     Returns
#     -------
#     supported_elements : list
#         A sorted list of all the supported elements.
#     """
#     supported_elements = set()
#     for rule_number, rule in RULE_NAME_TO_RULE.items():
#         decorators = get_decorators(rule)
#
#         element_type = None
#         for dec in decorators:
#             if isinstance(dec, Element):
#                 check_duplicate_element(element_type, rule_number)
#                 element_type = dec.element_type
#                 supported_elements.add(element_type)
#     supported_elements = list(supported_elements)
#     supported_elements.sort()
#     return supported_elements
#
#
# def find_all_rule_matches(supported_elements):
#     """Find all elements and combinations of neighbor types that have a rule.
#
#     Parameters
#     ----------
#     supported_elements : list
#         A sorted list of all the supported elements.
#
#     Returns
#     -------
#     rule_matches : dict
#         All patterns and rules that apply. The dictionary is structured as
#         follows:
#             key: (element, (neighbor element 1, neighbor element 2, etc..))
#             value: set(rule numbers)
#         Example entry (from time of writing this comment):
#             ('C', ('C', 'C', 'H')): set(['145', '142'])
#
#     See also
#     --------
#     find_all_supported_elements
#
#     """
#     rule_matches = dict()
#     for rule_number, rule in RULE_NAME_TO_RULE.items():
#         decorators = get_decorators(rule)
#
#         element_type = None
#         neighbor_count = None
#         for dec in decorators:
#             if isinstance(dec, Element):
#                 element_type = dec.element_type
#             if isinstance(dec, NeighborCount):
#                 check_duplicate_neighbor_count(neighbor_count, rule_number)
#                 neighbor_count = dec.count
#         # All POSSIBLE combinations of elements and neighbors.
#         all_patterns = set(combinations_with_replacement(supported_elements,
#                                                          neighbor_count))
#
#         # Remove the ones that don't actually have a rule.
#         removed_patterns = set()
#         for dec in decorators:
#             if isinstance(dec, NeighborsExactly):
#                 for pattern in all_patterns:
#                     if not pattern.count(dec.neighbor_type) == dec.count:
#                         removed_patterns.add(pattern)
#             elif isinstance(dec, NeighborsAtLeast):
#                 for pattern in all_patterns:
#                     if not pattern.count(dec.neighbor_type) >= dec.count:
#                         removed_patterns.add(pattern)
#             elif isinstance(dec, NeighborsAtMost):
#                 for pattern in all_patterns:
#                     if not pattern.count(dec.neighbor_type) <= dec.count:
#                         removed_patterns.add(pattern)
#         all_patterns.difference_update(removed_patterns)
#
#         for pattern in all_patterns:
#             if (element_type, pattern) not in rule_matches:
#                 rule_matches[(element_type, pattern)] = {rule_number}
#             else:
#                 rule_matches[(element_type, pattern)].add(rule_number)
#     return rule_matches
