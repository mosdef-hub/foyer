from networkx.algorithms import isomorphism
from foyer.smarts import SMARTS
from plyplus.strees import STree

import networkx as nx


class Ast2Nx:
    """
    A class first to parse the SMARTS string to AST tree and finally convert to NetworkX graph.
    * Subgraph Isomorphism
        Support subgraph isomorphism by VF2 algorithm from NetworkX.
    \ref{VF2 algorithm}
    http://networkx.readthedocs.io/en/stable/reference/algorithms.isomorphism.vf2.html?highlight=VF2#subgraph-isomorphism
    Actually, we found that subgraph isomorphism can be achieved in linear time for planer graph, however, it is really
    hard to implement the algorithm described in the following paper which require rich mathematics knowledge.
    \cite{Subgraph Isomorphism in Planar Graphs and Related Problems}
    https://arxiv.org/abs/cs/9911003
    * Visualization
        Support graph visualization by Matplotlib. Another way to visualize is to generate the .graphml file and adopt
    some professional tools such as Cytoscape. Interestingly, Cytoscape has already had a plug-in called NetMatchStar
    to achieve subgraph isomorphism which can also highlight all the occurrences of the subgraph.
    \ref{Cytoscape}
    http://www.cytoscape.org/
    \ref{NetMatchStar}
    http://apps.cytoscape.org/apps/netmatchstar
    """

    def __init__(self, smarts_str):
        # variables
        self.PARSER = SMARTS()  # initial SMARTS parser
        self.AST = self.PARSER.parse(smarts_str)  # initial AST tree structure which is parsed from SMARTS string
        self.AST.select('start')  # in order to generate the order inside of the AST such as next_kid
        self.atom_with_id = {}  # all the atoms(nodes) added to the graph should have a unique id
        self.atom_name = {}  # subgraph isomorphism match based on the name attribute of the nodes. if two nodes have the same name, then they will be considered equal. Same mechanism for the name of edge
        self.atom_with_label = {}  # to address the ring issue, we introduce atom_label to the AST such that if two atoms own the same label number, they will be connected together
        self.NetworkX = nx.Graph()  # wait to convert

        # invoke initial functions to initialize
        self._assign_id()  # assign a unique id and name to each of the atoms, which will initialize self.atom_id & self.atom_name
        self._set_atoms_with_label()  # get all the atoms that have the atom_label attribute

        # invoke converting functions to generate graph
        self._add_nodes()
        self._add_edges(self.AST)
        self._add_label_edges()

    # ===== initialize function ===== #
    def _assign_id(self):
        # assign a unique id to each atom
        atom_id = 0
        atoms = self.AST.select('atom')
        for atom in atoms:
            self.atom_with_id[id(atom)] = (atom, atom_id)
            self.atom_name[id(atom)] = self._set_atom_name(atom)
            atom_id += 1

    def _set_atoms_with_label(self):
        # set all atoms with atom_label attribute
        atoms_with_label = self.AST.select('atom_label')
        for atom_ in atoms_with_label:
            assert atom_.parent().head == 'atom', "the parent of atom_label has to be atom."
            self.atom_with_label[id(atom_.parent())] = list(atom_.tail[0])

    def _set_atom_name(self, atom):
        # ToDo: deal with the any name situation, i.e., (*)
        atom_name_list = []
        for atom_name_ in atom.tail:
            if atom_name_.head == 'atom_label':
                # atom_name_list.append('labelled')  # not sure whether we need the atom_label information or not
                pass
            else:
                atom_name_list.append(str(atom_name_))
        return 'atom(' + ', '.join(atom_name_list) + ')'

    def get_atom_with_id(self, atom):
        # given atom return atom_with_id
        if id(atom) in self.atom_with_id:
            return self.atom_with_id[id(atom)]
        return self.atom_with_id[atom]

    def get_atom_name(self, atom):
        # given atom return atom_name
        if id(atom) in self.atom_name:
            return self.atom_name[id(atom)]
        return self.atom_name[atom]

    def get_edge_name(self, atom1, atom2):
        # given two atoms return the edge name between them
        return self.get_atom_name(atom1) + '-' + self.get_atom_name(atom2)

    # ===== graph generation function ===== #
    def _add_nodes(self):
        # add all nodes to the graph
        atoms = self.AST.select('atom')
        for atom in atoms:
            atom_name = self.get_atom_name(atom)
            self.NetworkX.add_node(self.get_atom_with_id(atom), name=atom_name)

    def _add_edges(self, ASTtree, trunk=None):
        # add all edges to the graph
        for atom in ASTtree.tail:
            if atom.head == 'atom':
                # atom is the type that want to add to the graph
                if atom.is_first_kid and (atom.parent().head == 'branch'):
                    # if this atom is the first one in its branch then it should connect to the trunk
                    assert trunk is not None, "can't add branch to a None root!"
                    self.NetworkX.add_edge(self.get_atom_with_id(atom), self.get_atom_with_id(trunk),
                                           name=self.get_edge_name(id(atom), id(trunk)))
                if not atom.is_last_kid:
                    # if this atom is not the last one, it should connect to the its next atom
                    if atom.next_kid.head == 'atom':
                        self.NetworkX.add_edge(self.get_atom_with_id(atom), self.get_atom_with_id(atom.next_kid),
                                               name=self.get_edge_name(atom, atom.next_kid))
                    elif atom.next_kid.head == 'branch':
                        # if the next atom is a new branch then this atom should be the trunk for the new branch
                        trunk = atom
                else:
                    return  # we already travel through the whole branch
            elif atom.head == 'branch':
                # a new branch appeared, so we recursively travel to the new branch
                self._add_edges(atom, trunk)

    def _add_label_edges(self):
        # connect all the atoms with same atom_label together
        for atom_id, labels in self.atom_with_label.items():
            for atom_id_inner, labels_inner in self.atom_with_label.items():
                if atom_id_inner is not atom_id:
                    for label_inner in labels_inner:
                        if label_inner in labels:
                            self.NetworkX.add_edge(self.get_atom_with_id(atom_id_inner), self.get_atom_with_id(atom_id),
                                                   name=self.get_edge_name(atom_id_inner, atom_id))

    def to_file(self, name_graphml='AST2NX.graphml'):
        # write to a graphml file which can be read by a lot of professional visualization tools such as Cytoscape
        if name_graphml.endswith('.graphml'):
            nx.write_graphml(self.NetworkX, name_graphml)
        else:
            nx.write_graphml(self.NetworkX, name_graphml + '.graphml')

    # ===== isomorphism functions ===== #
    def foyer_node_match(self, G1_node, G2_node):
        # ToDo: need a completed atom matching rules here.
        # ToDo: For example, atomic_num(6) ?= atom(C), how to deal with and_expression, and how to deal with wildcard (*)?
        if G1_node == G2_node:
            return True
        return False

    def foyer_edge_match(self, G1_edge, G2_edge):
        # ToDo: need a completed edge matching rules here.
        if G1_edge == G2_edge:
            return True
        return False

    def __eq__(self, other):
        # entire graph isomorphism
        GM = isomorphism.GraphMatcher(self.NetworkX, other.NetworkX,
                                      node_match=self.foyer_node_match, edge_match=self.foyer_edge_match)
        return isomorphism.GraphMatcher.is_isomorphic(GM)

    def __contains__(self, item):
        # subgraph isomorphic
        # ToDO: iterator for subgraph isomorphism \ref{subgraph_isomorphisms_iter} http://networkx.readthedocs.io/en/stable/reference/generated/networkx.algorithms.isomorphism.GraphMatcher.subgraph_isomorphisms_iter.html#subgraph-isomorphisms-iter
        GM = isomorphism.GraphMatcher(self.NetworkX, item.NetworkX,
                                      node_match=self.foyer_node_match, edge_match=self.foyer_edge_match)
        return isomorphism.GraphMatcher.subgraph_is_isomorphic(GM)


# override the "less or equal to" operator for STree class for the graph isomorphism
def __lt__(self, other):
    return str(self) < str(other)
STree.__lt__ = __lt__
