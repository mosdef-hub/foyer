import itertools as it
from collections import OrderedDict, defaultdict

import networkx as nx
from networkx.algorithms import isomorphism
import parmed.periodic_table as pt

from foyer.exceptions import FoyerError
from foyer.smarts import SMARTS


class SMARTSGraph(nx.Graph):
    """

    Attributes
    ----------
    smarts_string : str
    parser : foyer.smarts.SMARTS.PARSER
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

        self._add_nodes()
        self._add_edges(self.ast)
        self._add_label_edges()

    def _add_nodes(self):
        """Add all atoms in the SMARTS string as nodes in the graph. """
        for atom in self.ast.select('atom'):
            self.add_node(id(atom), atom=atom)

    def _add_edges(self, ast_node, trunk=None):
        """"Add all bonds in the SMARTS string as edges in the graph. """
        for atom in ast_node.tail:
            if atom.head == 'atom':
                if atom.is_first_kid and atom.parent().head == 'branch':
                    if trunk is None:
                        raise FoyerError("Can't add branch without a trunk")
                    self.add_edge(id(atom), id(trunk))
                if not atom.is_last_kid:
                    if atom.next_kid.head == 'atom':
                        self.add_edge(id(atom), id(atom.next_kid))
                    elif atom.next_kid.head == 'branch':
                        trunk = atom
                else:  # We traveled through the whole branch.
                    return
            elif atom.head == 'branch':
                self._add_edges(atom, trunk)

    def _add_label_edges(self):
        """Add edges between all atoms with the same atom_label in rings. """
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

        for label1, label2 in it.permutations(label_digits, 2):
            if label1 == label2:
                atom1, atom2 = label_digits[label1], label_digits[label2]
                self.add_edge(atom1, atom2)

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
            raise TypeError('Expected and_expression, or_expression,'
                            ' or atom_id, got {}'.format(atom_expr.head))

    def _atom_id_matches(self, atom_id, atom):
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
            return label in (rule_name for rule_name in atom.whitelist)
        elif atom_id.head == 'neighbor_count':
            return len(atom.bond_partners) == int(atom_id.tail[0])
        elif atom_id.head == 'ring_size':
            cycle_len = int(atom_id.tail[0])
            for cycle in atom.cycles:
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_id.head == 'matches_string':
            raise NotImplementedError('matches_string feature is not yet implemented')

    def find_matches(self, topology):
        """Return sets of atoms that match this SMARTS pattern in a topology.

        Notes:
        ------
        When this function gets used in atomtyper.py, we actively modify the
        white- and blacklists of the atoms in `topology` after finding a match.
        What this means is that between every successive call of
        `subgraph_isomorphisms_iter()`, the topology against which we are
        matching may have actually changed. Currently, we take advantage of this
        behavior in some edges cases (e.g. see `test_hexa_coordinated` in
        `test_smarts.py`).
        """
        if topology is None:
            return False

        g = nx.Graph()
        g.add_nodes_from(((a.index, {'atom': a})
                          for a in topology.atoms()))
        g.add_edges_from(((b.atom1.index, b.atom2.index)
                          for b in topology.bonds()))

        gm = isomorphism.GraphMatcher(g, self, node_match=self._node_match)
        for mapping in gm.subgraph_isomorphisms_iter():
            # The first node in the smarts graph always corresponds to the atom
            # that we are trying to match.
            mapping = {node_id: atom_id for atom_id, node_id in mapping.items()}
            atom_index = mapping[next(self.nodes_iter())]
            yield atom_index


if __name__ == '__main__':
    from foyer.tests.utils import get_fn
    import parmed as pmd
    from foyer.forcefield import generate_topology

    mol2 = pmd.load_file('opls_validation/methanol/methanol.top', xyz='opls_validation/methanol/methanol.gro')
    top, _ = generate_topology(mol2)

    from oset import oset as OrderedSet
    atoms = list(top.atoms())
    for atom in atoms:
        atom.whitelist = OrderedSet()
        # if atom.element.symbol == 'N':
        #     atom.whitelist.add('opls_760')
        # if atom.element.symbol == 'H':
        #     atom.whitelist.add('opls_140')
        atom.blacklist = OrderedSet()

    # opls_763 = SMARTSGraph('H[C;X4][N;%opls_760]')
    graph = SMARTSGraph('HC(H)(H)OH')

    for _ in range(20):
        graph = SMARTSGraph('HC(H)(H)OH')
        print(graph.node.__class__)
        print([graph.node[x]['atom'].tail[0].tail[0] for x in graph.nodes()])
        # for i, mapping in enumerate(graph.find_matches(top)):
        #     print(i, mapping, atoms[mapping])

