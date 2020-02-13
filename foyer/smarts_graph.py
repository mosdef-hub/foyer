from collections import OrderedDict, defaultdict
import itertools

import networkx as nx
from networkx.algorithms import isomorphism
import parmed.periodic_table as pt

from foyer.smarts import SMARTS


class SMARTSGraph(nx.Graph):
    """A graph representation of a SMARTS pattern.

    Attributes
    ----------
    smarts_string : str
        The SMARTS string outlined in the force field
    parser : foyer.smarts.SMARTS
        The parser whose grammar rules convert the SMARTSstring 
        into the AST
    name : str
    overrides : set
        Rules or SMARTSGraph over which this SMARTSGraph takes precedence

    Other Parameters
    ----------
    args
    kwargs

    Attributes
    ----------
    graph_matcher : smarts_graph.SMARTSMatcher
        implementation of VF2 that handles subgraph matching
    """
    # Because the first atom in a SMARTS string is always the one we want to
    # type, the graph's nodes needs to be ordered.

    def __init__(self, smarts_string, parser=None, name=None, overrides=None,
                typemap=None,
                 *args, **kwargs):
        super(SMARTSGraph, self).__init__(*args, **kwargs)

        self.smarts_string = smarts_string
        self.name = name
        self.overrides = overrides
        self.typemap = typemap

        if parser is None:
            self.ast = SMARTS().parse(smarts_string)
        else:
            self.ast = parser.parse(smarts_string)

        self._atom_indices = OrderedDict()
        self._add_nodes()
        self._add_edges(self.ast)
        self._add_label_edges()
        self._graph_matcher = None

    def _add_nodes(self):
        """Add all atoms in the SMARTS string as nodes in the graph."""
        for n, atom in enumerate(self.ast.find_data('atom')):
            self.add_node(n, atom=atom)
            self._atom_indices[id(atom)] = n

    def _add_edges(self, ast_node, trunk=None):
        """"Add all bonds in the SMARTS string as edges in the graph."""
        atom_indices = self._atom_indices
        for ast_child in ast_node.children:
            if ast_child.data == 'atom':
                atom_idx = atom_indices[id(ast_child)]
                if trunk is not None:
                    trunk_idx = atom_indices[id(trunk)]
                    self.add_edge(atom_idx, trunk_idx)
                trunk = ast_child
            elif ast_child.data == 'branch':
                self._add_edges(ast_child, trunk)

    def _add_label_edges(self):
        """Add edges between all atoms with the same atom_label in rings."""
        # We need each individual label and atoms with multiple ring labels
        # would yield e.g. the string '12' so split those up.
        label_digits = defaultdict(list)
        for node, attr in self.nodes(data=True):
            atom = attr["atom"]
            for label in atom.find_data("atom_label"):
                digits = list(label.children[0])
                for digit in digits:
                    label_digits[digit].append(atom)

        for label, (atom1, atom2) in label_digits.items():
            atom1_idx = self._atom_indices[id(atom1)]
            atom2_idx = self._atom_indices[id(atom2)]
            self.add_edge(atom1_idx, atom2_idx)

    def _node_match(self, host, pattern):
        """ Determine if two graph nodes are equal """
        atom_expr = pattern['atom'].children[0]
        atom = host['atom']
        return self._atom_expr_matches(atom_expr, atom)

    def _atom_expr_matches(self, atom_expr, atom):
        """ Helper function for evaluating SMARTS string expressions """
        if atom_expr.data == 'not_expression':
            return not self._atom_expr_matches(atom_expr.children[0], atom)
        elif atom_expr.data in ('and_expression', 'weak_and_expression'):
            return (self._atom_expr_matches(atom_expr.children[0], atom) and
                    self._atom_expr_matches(atom_expr.children[1], atom))
        elif atom_expr.data == 'or_expression':
            return (self._atom_expr_matches(atom_expr.children[0], atom) or
                    self._atom_expr_matches(atom_expr.children[1], atom))
        elif atom_expr.data == 'atom_id':
            return self._atom_id_matches(atom_expr.children[0], atom, self.typemap)
        elif atom_expr.data == 'atom_symbol':
            return self._atom_id_matches(atom_expr, atom, self.typemap)
        else:
            raise TypeError('Expected atom_id, atom_symbol, and_expression, '
                            'or_expression, or not_expression. '
                            'Got {}'.format(atom_expr.data))

    @staticmethod
    def _atom_id_matches(atom_id, atom, typemap):
        """ Helper func for comparing atomic indices, symbols, neighbors, rings """
        atomic_num = atom.element
        if atom_id.data == 'atomic_num':
            return atomic_num == int(atom_id.children[0])
        elif atom_id.data == 'atom_symbol':
            if str(atom_id.children[0]) == '*':
                return True
            elif str(atom_id.children[0]).startswith('_'):
                # Store non-element elements in .name
                return atom.name == str(atom_id.children[0])
            else:
                return atomic_num == pt.AtomicNum[str(atom_id.children[0])]
        elif atom_id.data == 'has_label':
            label = atom_id.children[0][1:]  # Strip the % sign from the beginning.
            return label in typemap[atom.idx]['whitelist']
        elif atom_id.data == 'neighbor_count':
            return len(atom.bond_partners) == int(atom_id.children[0])
        elif atom_id.data == 'ring_size':
            cycle_len = int(atom_id.children[0])
            for cycle in typemap[atom.idx]['cycles']:
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_id.data == 'ring_count':
            n_cycles = len(typemap[atom.idx]['cycles'])
            if n_cycles == int(atom_id.children[0]):
                return True
            return False
        elif atom_id.data == 'matches_string':
            raise NotImplementedError('matches_string is not yet implemented')

    def find_matches(self, structure, typemap):
        """Return sets of atoms that match this SMARTS pattern in a topology.

        Notes:
        ------
        When this function gets used in atomtyper.py, we actively modify the
        white- and blacklists of the atoms in `topology` after finding a match.
        This means that between every successive call of
        `subgraph_isomorphisms_iter()`, the topology against which we are
        matching may have actually changed. Currently, we take advantage of this
        behavior in some edges cases (e.g. see `test_hexa_coordinated` in
        `test_smarts.py`).

        """
        # Note: Needs to be updated in sync with the grammar in `smarts.py`.
        ring_tokens = ['ring_size', 'ring_count']
        has_ring_rules = any(list(self.ast.find_data(token))
                             for token in ring_tokens)
        _prepare_atoms(structure, typemap, compute_cycles=has_ring_rules)

        top_graph = nx.Graph()
        top_graph.add_nodes_from(((a.idx, {'atom': a})
                                  for a in structure.atoms))
        top_graph.add_edges_from(((b.atom1.idx, b.atom2.idx)
                                  for b in structure.bonds))

        if self._graph_matcher is None:
            atom = nx.get_node_attributes(self, name='atom')[0]
            if len(list(atom.find_data('atom_symbol'))) == 1 and \
                        not list(atom.find_data('not_expression')):
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
            self._graph_matcher = SMARTSMatcher(top_graph, self,
                                                node_match=self._node_match,
                                                element=element,
                                                typemap=typemap)

        matched_atoms = set()
        for mapping in self._graph_matcher.subgraph_isomorphisms_iter():
            mapping = {node_id: atom_id for atom_id, node_id in mapping.items()}
            # The first node in the smarts graph always corresponds to the atom
            # that we are trying to match.
            atom_index = mapping[0]
            # Don't yield duplicate matches found via matching the pattern in a
            # different order.
            if atom_index not in matched_atoms:
                matched_atoms.add(atom_index)
                yield atom_index


class SMARTSMatcher(isomorphism.vf2userfunc.GraphMatcher):
    """ Inherits and implements VF2 for a SMARTSGraph"""
    def __init__(self, G1, G2, node_match, element, typemap):
        super(SMARTSMatcher, self).__init__(G1, G2, node_match)
        self.element = element
        # TODO: Parse out nodes containing other elements (see git history)
        self.valid_nodes = G1.nodes()

    def candidate_pairs_iter(self):
        """Iterator over candidate pairs of nodes in G1 and G2."""
        # All computations are done using the current state!
        G2_nodes = self.G2_nodes

        # First we compute the inout-terminal sets.
        T1_inout = set(self.inout_1.keys()) - set(self.core_1.keys())
        T2_inout = set(self.inout_2.keys()) - set(self.core_2.keys())

        # If T1_inout and T2_inout are both nonempty.
        # P(s) = T1_inout x {min T2_inout}
        if T1_inout and T2_inout:
            for node in T1_inout:
                yield node, min(T2_inout)
        else:
            # First we determine the candidate node for G2
            other_node = min(G2_nodes - set(self.core_2))
            host_nodes = self.valid_nodes if other_node == 0 else self.G1.nodes()
            for node in host_nodes:
                if node not in self.core_1:
                    yield node, other_node

        # For all other cases, we don't have any candidate pairs.


def _find_chordless_cycles(bond_graph, max_cycle_size):
    """Find all chordless cycles (i.e. rings) in the bond graph

    Traverses the bond graph to determine all cycles (i.e. rings) each
    atom is contained within. Algorithm has been adapted from:
    https://stackoverflow.com/questions/4022662/find-all-chordless-cycles-in-an-undirected-graph/4028855#4028855
    """
    cycles = [[] for _ in bond_graph.nodes]

    '''
    For all nodes we need to find the cycles that they are included within.
    '''
    for i, node in enumerate(bond_graph.nodes):
        neighbors = list(bond_graph.neighbors(node))
        pairs = list(itertools.combinations(neighbors, 2))
        '''
        Loop over all pairs of neighbors of the node. We will see if a ring
        exists that includes these branches.
        '''
        for pair in pairs:
            '''
            We need to store all node sequences that could be rings. We will
            update this as we traverse the graph.
            '''
            connected = False
            possible_rings = []

            last_node = pair[0]
            ring = [last_node, node, pair[1]]
            possible_rings.append(ring)

            if bond_graph.has_edge(last_node, pair[1]):
                cycles[i].append(ring)
                connected = True

            while not connected:
                '''
                Branch and create a new list of possible rings
                '''
                new_possible_rings = []
                for possible_ring in possible_rings:
                    next_neighbors = list(bond_graph.neighbors(possible_ring[-1]))
                    for next_neighbor in next_neighbors:
                        if next_neighbor != possible_ring[-2]:
                            new_possible_rings.append(possible_ring + \
                                                      [next_neighbor])
                possible_rings = new_possible_rings

                for possible_ring in possible_rings:
                    if bond_graph.has_edge(possible_ring[-1], last_node):
                        if any([bond_graph.has_edge(possible_ring[-1], 
                                internal_node)
                                for internal_node in possible_ring[1:-2]]):
                            pass
                        else:
                            cycles[i].append(possible_ring)
                            connected = True

                if not possible_rings or len(possible_rings[0]) == max_cycle_size:
                    break

    return cycles


def _prepare_atoms(structure, typemap, compute_cycles=False):
    """Compute cycles and add white-/blacklists to atoms."""
    atom1 = structure.atoms[0]#next(topology.atoms())
    has_whitelists = 'whitelist' in typemap[atom1.idx]
    has_cycles = 'cycles' in typemap[atom1.idx]
    compute_cycles = compute_cycles and not has_cycles

    if compute_cycles or not has_whitelists:
        for atom in structure.atoms:
            if compute_cycles:
                typemap[atom.idx]['cycles'] = set()
            if not has_whitelists:
                typemap[atom.idx]['whitelist'] = set()
                typemap[atom.idx]['blacklist'] = set()

    if compute_cycles:
        bond_graph = nx.Graph()
        bond_graph.add_nodes_from(structure.atoms)
        bond_graph.add_edges_from([(b.atom1, b.atom2) for b in structure.bonds])
        all_cycles = _find_chordless_cycles(bond_graph, max_cycle_size=8)
        for atom, cycles in zip(bond_graph.nodes, all_cycles):
            for cycle in cycles:
                typemap[atom.idx]['cycles'].add(tuple(cycle))
