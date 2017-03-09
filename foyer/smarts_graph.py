import itertools as it
from collections import OrderedDict

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
    def __init__(self, smarts_string, parser=None, name=None, overrides=None,
                 *args, **kwargs):
        super(SMARTSGraph, self).__init__(*args, **kwargs)
        self.node_dict_factory = OrderedDict

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
        atoms = self.ast.select('atom')
        for atom in atoms:
            self.add_node(id(atom), atom=atom)

    def _add_edges(self, ast_node, trunk=None):
        """"Add all bonds in the SMARTS string as edges in the graph. """
        for atom in ast_node.tail:
            if atom.head == 'atom':
                if atom.is_first_kid and atom.parent().head == 'branch':
                    if trunk is None:
                        raise FoyerError("Can't add branch without a trunk")
                    self.add_edge(id(atom), id(trunk))
                elif not atom.is_last_kid:
                    if atom.next_kid.head == 'atom':
                        self.add_edge(id(atom), id(atom.next_kid),)
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
        label_digits = []
        for label in labels:
            label_digits.extend(*label)

        for label1, label2 in it.permutations(label_digits, 2):
            if label1.tail[0] == label2.tail[0]:
                atom1, atom2 = label1.parent(), label2.parent()
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
        if topology is None:
            return False

        g = nx.Graph()
        g.add_nodes_from(((a.index, {'atom': a})
                          for a in topology.atoms()))
        g.add_edges_from(((b.atom1.index, b.atom2.index)
                          for b in topology.bonds()))

        gm = isomorphism.GraphMatcher(g, self, node_match=self._node_match)

        for mapping in gm.subgraph_isomorphisms_iter():
            ordered_mapping = []
            mapping = {id(v): k for k, v in mapping.items()}
            for node in self.nodes():
                smarts_atom = self.node[node]['atom']
                atom_index = mapping[id(node)]
                ordered_mapping.append((smarts_atom, atom_index))
            yield ordered_mapping


if __name__ == '__main__':
    from foyer.tests.utils import get_fn
    import parmed as pmd
    from foyer.forcefield import generate_topology

    mol2 = pmd.load_file(get_fn('ring.mol2'), structure=True)
    top, _ = generate_topology(mol2)

    pattern = SMARTSGraph('[C;X2;r6][#6;X2][#6;X2]')
    pattern = SMARTSGraph('[N;X3]([O;X1])([O;X1])[C;X3;r6]')
    for i, mapping in enumerate(pattern.find_matches(top)):
        print(i, mapping)

