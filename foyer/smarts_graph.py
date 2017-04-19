from collections import OrderedDict, defaultdict
import sys

import networkx as nx
from networkx.algorithms import isomorphism
from oset import oset as OrderedSet
import parmed.periodic_table as pt

from foyer.smarts import SMARTS


class SMARTSGraph(nx.Graph):
    """A graph representation of a SMARTS pattern.

    Attributes
    ----------
    smarts_string : str
    parser : foyer.smarts.SMARTS
    name : str
    overrides : set

    Other Parameters
    ----------
    args
    kwargs
    """
    # Because the first atom in a SMARTS string is always the one we want to
    # type, the graph's nodes needs to be ordered.
    node_dict_factory = OrderedDict

    def __init__(self, smarts_string, parser=None, name=None, overrides=None,
                 *args, **kwargs):
        super(SMARTSGraph, self).__init__(*args, **kwargs)

        self.smarts_string = smarts_string
        self.name = name
        self.overrides = overrides

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
        for n, atom in enumerate(self.ast.select('atom')):
            self.add_node(n, atom=atom)
            self._atom_indices[id(atom)] = n

    def _add_edges(self, ast_node, trunk=None):
        """"Add all bonds in the SMARTS string as edges in the graph."""
        atom_indices = self._atom_indices
        for atom in ast_node.tail:
            if atom.head == 'atom':
                atom_idx = atom_indices[id(atom)]
                if atom.is_first_kid and atom.parent().head == 'branch':
                    trunk_idx = atom_indices[id(trunk)]
                    self.add_edge(atom_idx, trunk_idx)
                if not atom.is_last_kid:
                    if atom.next_kid.head == 'atom':
                        next_idx = atom_indices[id(atom.next_kid)]
                        self.add_edge(atom_idx, next_idx)
                    elif atom.next_kid.head == 'branch':
                        trunk = atom
                else:  # We traveled through the whole branch.
                    return
            elif atom.head == 'branch':
                self._add_edges(atom, trunk)

    def _add_label_edges(self):
        """Add edges between all atoms with the same atom_label in rings."""
        labels = self.ast.select('atom_label')
        if not labels:
            return

        # We need each individual label and atoms with multiple ring labels
        # would yield e.g. the string '12' so split those up.
        label_digits = defaultdict(list)
        for label in labels:
            digits = list(label.tail[0])
            for digit in digits:
                label_digits[digit].append(label.parent())

        for label, (atom1, atom2) in label_digits.items():
            atom1_idx = self._atom_indices[id(atom1)]
            atom2_idx = self._atom_indices[id(atom2)]
            self.add_edge(atom1_idx, atom2_idx)

    def _node_match(self, host, pattern):
        atom_expr = pattern['atom'].tail[0]
        atom = host['atom']
        return self._atom_expr_matches(atom_expr, atom)

    def _atom_expr_matches(self, atom_expr, atom):
        if atom_expr.head == 'not_expression':
            return not self._atom_expr_matches(atom_expr.tail[0], atom)
        elif atom_expr.head in ('and_expression', 'weak_and_expression'):
            return (self._atom_expr_matches(atom_expr.tail[0], atom) and
                    self._atom_expr_matches(atom_expr.tail[1], atom))
        elif atom_expr.head == 'or_expression':
            return (self._atom_expr_matches(atom_expr.tail[0], atom) or
                    self._atom_expr_matches(atom_expr.tail[1], atom))
        elif atom_expr.head == 'atom_id':
            return self._atom_id_matches(atom_expr.tail[0], atom)
        elif atom_expr.head == 'atom_symbol':
            return self._atom_id_matches(atom_expr, atom)
        else:
            raise TypeError('Expected atom_id, atom_symbol, and_expression, '
                            'or_expression, or not_expression. '
                            'Got {}'.format(atom_expr.head))

    @staticmethod
    def _atom_id_matches(atom_id, atom):
        atomic_num = atom.element.atomic_number
        if atom_id.head == 'atomic_num':
            return atomic_num == int(atom_id.tail[0])
        elif atom_id.head == 'atom_symbol':
            if str(atom_id.tail[0]) == '*':
                return True
            elif str(atom_id.tail[0]).startswith('_'):
                return atom.element.name == str(atom_id.tail[0])
            else:
                return atomic_num == pt.AtomicNum[str(atom_id.tail[0])]
        elif atom_id.head == 'has_label':
            label = atom_id.tail[0][1:]  # Strip the % sign from the beginning.
            return label in atom.whitelist
        elif atom_id.head == 'neighbor_count':
            return len(atom.bond_partners) == int(atom_id.tail[0])
        elif atom_id.head == 'ring_size':
            cycle_len = int(atom_id.tail[0])
            for cycle in atom.cycles:
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_id.head == 'ring_count':
            n_cycles = len(atom.cycles)
            if n_cycles == int(atom_id.tail[0]):
                return True
            return False
        elif atom_id.head == 'matches_string':
            raise NotImplementedError('matches_string is not yet implemented')

    def find_matches(self, topology):
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
        has_ring_rules = any(self.ast.select(token)
                             for token in ring_tokens)
        _prepare_atoms(topology, compute_cycles=has_ring_rules)

        top_graph = nx.Graph()
        top_graph.add_nodes_from(((a.index, {'atom': a})
                                  for a in topology.atoms()))
        top_graph.add_edges_from(((b[0].index, b[1].index)
                                  for b in topology.bonds()))

        if self._graph_matcher is None:
            atom = nx.get_node_attributes(self, 'atom')[0]
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
            self._graph_matcher = SMARTSMatcher(top_graph, self,
                                                node_match=self._node_match,
                                                element=element)

        # The first node in the smarts graph always corresponds to the atom
        # that we are trying to match.
        first_atom = next(self.nodes_iter())
        matched_atoms = set()
        for mapping in self._graph_matcher.subgraph_isomorphisms_iter():
            mapping = {node_id: atom_id for atom_id, node_id in mapping.items()}
            atom_index = mapping[first_atom]
            # Don't yield duplicate matches found via matching the pattern in a
            # different order.
            if atom_index not in matched_atoms:
                matched_atoms.add(atom_index)
                yield atom_index


class SMARTSMatcher(isomorphism.vf2userfunc.GraphMatcher):
    def __init__(self, G1, G2, node_match, element):
        super(SMARTSMatcher, self).__init__(G1, G2, node_match)
        self.element = element
        if element not in [None, '*']:
            self.valid_nodes = [n for n, atom in nx.get_node_attributes(G1, 'atom').items()
                                if atom.element.symbol == element]
        else:
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


def _prepare_atoms(topology, compute_cycles=False):
    """Compute cycles and add white-/blacklists to atoms."""
    atom1 = next(topology.atoms())
    has_whitelists = hasattr(atom1, 'whitelist')
    has_cycles = hasattr(atom1, 'cycles')
    compute_cycles = compute_cycles and not has_cycles

    if compute_cycles or not has_whitelists:
        for atom in topology.atoms():
            if compute_cycles:
                atom.cycles = set()
            if not has_whitelists:
                atom.whitelist = OrderedSet()
                atom.blacklist = OrderedSet()

    if compute_cycles:
        bond_graph = nx.Graph()
        bond_graph.add_nodes_from(topology.atoms())
        bond_graph.add_edges_from(topology.bonds())
        cycles = nx.cycle_basis(bond_graph)
        for cycle in cycles:
            for atom in cycle:
                atom.cycles.add(tuple(cycle))
