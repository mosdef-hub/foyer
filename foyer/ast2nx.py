from networkx.algorithms import isomorphism
from foyer.smarts import SMARTS
from plyplus.strees import STree

import networkx as nx


class Ast2Nx:
    """ A class first to parse the SMARTS string to AST tree and finally convert to NetworkX graph.
    Support graph/subgraph isomorphism using VF2 algorithm.

    Parameters
    ----------
    smarts_str : SMARTS string
        The SMARTS string used to generate Networkx graph.

    Attributes
    ----------
    PARSER : SMARTS
        The SMARTS parser.
    AST : STree
        The abstract syntax tree generated from SMARTS.
    atom_with_id : dict, {id(atom): (atom, unique_id)}
        Assign each atom with a unique id.
    atom_name : dict, {id(atom): name}
        Assign each atom with a name. The name can be repeated in order to make grpah/sub-graph isomorphism.
    atom_with_label : dict, {id(atom}: atom_label}
        In SMARTS, atom_label is used to mark the jointed point in rings
    referrers : set
        Other compounds that reference this part with labels.
    """

    def __init__(self, smarts_str):

        self.PARSER = SMARTS()
        self.AST = self.PARSER.parse(smarts_str)
        self.AST.select('start')
        self.atom_with_id = {}
        self.atom_name = {}
        self.atom_with_label = {}
        self.NetworkX = nx.Graph()

        # invoke initial functions to initialize.
        self._assign_id()
        self._set_atoms_with_label()

        # invoke converting functions to generate graph.
        self._add_nodes()
        self._add_edges(self.AST)
        self._add_label_edges()

    def _assign_id(self):
        """ assign a unique id to each atom.
        """
        atom_id = 0
        atoms = self.AST.select('atom')
        for atom in atoms:
            self.atom_with_id[id(atom)] = (atom, atom_id)
            self.atom_name[id(atom)] = self._set_atom_name(atom)
            atom_id += 1

    def _set_atoms_with_label(self):
        """ assign an atom_label to each labelled atom.
        """
        atoms_with_label = self.AST.select('atom_label')
        for atom_ in atoms_with_label:
            assert atom_.parent().head == 'atom', "the parent of atom_label has to be atom."
            self.atom_with_label[id(atom_.parent())] = list(atom_.tail[0])

    def _set_atom_name(self, atom):
        """ the name assigned to each atom.
        The name here is particularly designed for isomorphism.
        The graph isomorphism algorithm will first check the syntactics, i.e., structure,
        and then check the identity of each node. The name here can be used for identity checking.
        """
        # ToDo: Need more details to specifically handle the atom name.
        atom_name_list = []
        for atom_name_ in atom.tail:
            if atom_name_.head == 'atom_label':
                # ignore atom_label information.
                pass
            else:
                atom_name_list.append(str(atom_name_))
        return 'atom(' + ', '.join(atom_name_list) + ')'

    def _get_atom_with_id(self, atom):
        """ given atom return atom_with_id.
        """
        if id(atom) in self.atom_with_id:
            return self.atom_with_id[id(atom)]
        return self.atom_with_id[atom]

    def _get_atom_name(self, atom):
        """ given atom return atom_name.
        """
        if id(atom) in self.atom_name:
            return self.atom_name[id(atom)]
        return self.atom_name[atom]

    def _get_edge_name(self, atom1, atom2):
        """ given two atoms return the edge name between them.
        """
        return self._get_atom_name(atom1) + '-' + self._get_atom_name(atom2)

    def _add_nodes(self):
        """ add all nodes to the graph.
        """
        atoms = self.AST.select('atom')
        for atom in atoms:
            atom_name = self._get_atom_name(atom)
            self.NetworkX.add_node(self._get_atom_with_id(atom), name=atom_name)

    def _add_edges(self, ASTtree, trunk=None):
        """ add all edges to the graph.
        """
        for atom in ASTtree.tail:
            if atom.head == 'atom':
                # atom is the type that want to add to the graph
                if atom.is_first_kid and (atom.parent().head == 'branch'):
                    # if this atom is the first one in its branch then it should connect to the trunk
                    assert trunk is not None, "can't add branch to a None root!"
                    self.NetworkX.add_edge(self._get_atom_with_id(atom), self._get_atom_with_id(trunk),
                                           name=self._get_edge_name(id(atom), id(trunk)))
                if not atom.is_last_kid:
                    # if this atom is not the last one, it should connect to the its next atom
                    if atom.next_kid.head == 'atom':
                        self.NetworkX.add_edge(self._get_atom_with_id(atom), self._get_atom_with_id(atom.next_kid),
                                               name=self._get_edge_name(atom, atom.next_kid))
                    elif atom.next_kid.head == 'branch':
                        # if the next atom is a new branch then this atom should be the trunk for the new branch
                        trunk = atom
                else:
                    return  # we already travel through the whole branch
            elif atom.head == 'branch':
                # a new branch appeared, so we recursively travel to the new branch
                self._add_edges(atom, trunk)

    def _add_label_edges(self):
        """ Connect all atoms with the same atom_label.
        """
        for atom_id, labels in self.atom_with_label.items():
            for atom_id_inner, labels_inner in self.atom_with_label.items():
                if atom_id_inner is not atom_id:
                    for label_inner in labels_inner:
                        if label_inner in labels:
                            self.NetworkX.add_edge(self._get_atom_with_id(atom_id_inner),
                                                   self._get_atom_with_id(atom_id),
                                                   name=self._get_edge_name(atom_id_inner, atom_id))

    def to_file(self, name_graphml='AST2NX.graphml'):
        """ write to a graphml file which can be read by a lot of professional visualization tools such as Cytoscape.
        """
        if name_graphml.endswith('.graphml'):
            nx.write_graphml(self.NetworkX, name_graphml)
        else:
            nx.write_graphml(self.NetworkX, name_graphml + '.graphml')

    # ===== isomorphism functions ===== #
    def foyer_node_match(self, G1_node, G2_node):
        """ the matching rule for node/atom.
        For example, atomic_num(6) ?= atom(C), how to deal with and_expression, and how to deal with wildcard (*)?
        """
        # ToDo: need a completed atom matching rules here.
        if G1_node == G2_node:
            return True
        return False

    def foyer_edge_match(self, G1_edge, G2_edge):
        """ the matching rule for edge/bond.
        """
        # ToDo: need a completed atom matching rules here.
        if G1_edge == G2_edge:
            return True
        return False

    def __eq__(self, other):
        """ whole graph ismorphism.
        """
        GM = isomorphism.GraphMatcher(self.NetworkX, other.NetworkX,
                                      node_match=self.foyer_node_match, edge_match=self.foyer_edge_match)
        return isomorphism.GraphMatcher.is_isomorphic(GM)

    def __contains__(self, item):
        """ subgraph ismorphism.
        """
        # ToDo: iterator for sub-graph isomorphism.
        GM = isomorphism.GraphMatcher(self.NetworkX, item.NetworkX,
                                      node_match=self.foyer_node_match, edge_match=self.foyer_edge_match)
        return isomorphism.GraphMatcher.subgraph_is_isomorphic(GM)


def __lt__(self, other):
    """ override the "less or equal to" operator for STree class for the graph isomorphism
    """
    return str(self) < str(other)
STree.__lt__ = __lt__
