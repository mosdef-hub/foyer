import itertools as it

import networkx as nx
from networkx.algorithms import isomorphism


from foyer.smarts import SMARTS


class SMARTSGraph(nx.Graph):
    def __init__(self, smarts_string, parser=None, *args, **kwargs):
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

    def _add_nodes(self):
        atoms = self.ast.select('atom')
        for atom in atoms:
            # Plyplus STree objects by default  evaluate the below but we care
            # about uniqueness of ast objects, not just their content.

            # From plyplus.strees.py:
            #
            # def __hash__(self):
            #     return hash((self.head, tuple(self.tail)))
            #
            # def __eq__(self, other):
            #     try:
            #         return self.head == other.head and self.tail == other.tail
            #     except AttributeError:
            #         return False

            def _node_equality(self, other):
                return id(self) == id(other)

            def _node_hash(self):
                return id(self)

            def _node_lt(self, other):
                return id(self) < id(other)

            atom.__class__.__eq__ = _node_equality
            atom.__class__.__hash__ = _node_hash
            atom.__class__.__lt__ = _node_lt
            self.add_node(atom, attr=atom)

    def _add_edges(self, ast_node, trunk=None):
        for atom in ast_node.tail:
            if atom.head == 'atom':
                if atom.is_first_kid and atom.parent().head == 'branch':
                    assert trunk is not None, "Can't add branch without a trunk"
                    self.add_edge(atom, trunk)
                elif not atom.is_last_kid:
                    if atom.next_kid.head == 'atom':
                        self.add_edge(atom, atom.next_kid,)
                    elif atom.next_kid.head == 'branch':
                        trunk = atom
                else:  # We traveled through the whole branch.
                    return
            elif atom.head == 'branch':
                self._add_edges(atom, trunk)

    def _add_label_edges(self):
        """Connect all atoms with the same atom_label in rings. """
        labels = self.ast.select('atom_label')
        if not labels:
            return
        for label1, label2 in it.permutations(labels, 2):
            if label1.tail[0] == label2.tail[0]:
                atom1, atom2 = label1.parent(), label2.parent()
                self.add_edge(atom1, atom2)

    def _node_match(self, node1, node2):
        #return node1['attr'] == node2['attr']
        return True

    def __eq__(self, other):
        gm = isomorphism.GraphMatcher(self, other, node_match=self._node_match)
        return isomorphism.GraphMatcher.is_isomorphic(gm)

    def __contains__(self, item):
        if item is None:
            return False
        gm = isomorphism.GraphMatcher(self, item, node_match=self._node_match)
        return isomorphism.GraphMatcher.subgraph_is_isomorphic(gm)


if __name__ == '__main__':
    #S1 = SMARTSGraph('O([H&X1])(H)')

    S1 = SMARTSGraph('[#6]1[#6][#6][#6][#6][#6]1')

    print('===============')
    print(len(S1.nodes()))
    for n in S1.nodes():
        print(n)
    print(len(S1.edges()))
    for e in S1.edges():
        print(e)
