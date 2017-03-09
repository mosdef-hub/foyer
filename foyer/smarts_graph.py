import itertools as it
from collections import OrderedDict

import networkx as nx
from networkx.algorithms import isomorphism
import parmed.periodic_table as pt

from foyer.smarts import SMARTS


class SMARTSGraph(nx.Graph):
    def __init__(self, smarts_string, parser=None, name=None, overrides=None, *args, **kwargs):
        super(SMARTSGraph, self).__init__(*args, **kwargs)
        self.smarts_string = smarts_string
        if parser is None:
            self.parser = SMARTS().PARSER
        else:
            self.parser = parser
        self.ast = self.parser.parse(smarts_string)

        self._add_nodes()
        self._add_edges(self.ast)
        self._add_label_edges()
        self.name = name
        self.overrides = overrides
        self.node_dict_factory = OrderedDict

    def _add_nodes(self):
        atoms = self.ast.select('atom')
        for atom in atoms:
            self.add_node(id(atom), atom=atom)

    def _add_edges(self, ast_node, trunk=None):
        for atom in ast_node.tail:
            if atom.head == 'atom':
                if atom.is_first_kid and atom.parent().head == 'branch':
                    assert trunk is not None, "Can't add branch without a trunk"
                    self.add_edge(id(atom), id(trunk))
                elif not atom.is_last_kid:
                    if atom.next_kid.head == 'atom':
                        self.add_edge(id(atom), id(atom.next_kid),)
                    elif atom.next_kid.head == 'branch':
                        trunk = atom
                else:  # We traveled through the whole branch.
                    return
            elif atom.head == 'branch':
                self._add_edges(id(atom), id(trunk))

    def _add_label_edges(self):
        """Connect all atoms with the same atom_label in rings. """
        labels = self.ast.select('atom_label')

        if not labels:
            return

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
        atomic_num = atom.element._atomic_number
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

        from simtk.openmm.app import Topology
        assert isinstance(topology, Topology)

        g = nx.Graph()
        g.add_nodes_from(((a.index, { 'atom': a }) for a in topology.atoms()), )
        g.add_edges_from(((b.atom1.index, b.atom2.index) for b in topology.bonds()))

        gm = isomorphism.GraphMatcher(g, self, node_match=self._node_match)

        for mapping in gm.subgraph_isomorphisms_iter():
            ordered_mapping = []
            mapping = {id(v): k for k, v in mapping.items()}
            # import pdb
            # pdb.set_trace()
            for node in self.nodes():
                ordered_mapping.append((self.node[node]['atom'], mapping[id(node)]))

            yield ordered_mapping

if __name__ == '__main__':
    #S1 = SMARTSGraph('O([H&X1])(H)')

    # S1 = SMARTSGraph('[#6]1[#6][#6][#6][#6][#6]1')
    from foyer.tests.utils import get_fn
    import parmed as pmd
    mol2 = pmd.load_file(get_fn('ring.mol2'), structure=True)
    from foyer.forcefield import generate_topology
    top,_ = generate_topology(mol2)


    pattern = SMARTSGraph('[C;X2][#6;X2][#6;X2]')


    # print('===============')
    # print(len(S1.nodes()))
    # for n in S1.nodes():
    #     print(n)
    # print(len(S1.edges()))
    # for e in S1.edges():
    #     print(e)

    for i,mapping in enumerate(pattern.find_matches(top)):
        print(i, mapping)

