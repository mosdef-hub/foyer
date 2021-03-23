import networkx as nx
from parmed import Structure


class AtomData:
    """Stores atom data necessary for atom typing

    Parameters
    ----------
    index: int
        The index of the atom in the topology graph
    name: str
        The name of the atom
    atomic_number: int
        The atomic number of the atom
    element: str
        The element associated with the atom
    bond_partners: list of int
        The list of length two containing bond partner atom indices for this atom
    **kwargs
        Any extra information for this atom

    ToDo: After 3.6 support removed, convert this to a dataclass
    """

    def __init__(self, index, name, atomic_number, element, bond_partners, **kwargs):
        self.index = index
        self.name = name
        self.atomic_number = atomic_number
        self.element = element
        self.bond_partners = bond_partners
        for key, value in kwargs.items():
            setattr(self, key, value)


class TopologyGraph(nx.Graph):
    """A general TopologyGraph

    This class subclasses nx.Graph to provide a general
    topology Graph for atom typing in foyer. Each node in
    this graph is identified by atom indices and edge in this
    graph represents a bond between the atoms.
    """

    def __init__(self, *args, **kwargs):
        super(TopologyGraph, self).__init__(*args, **kwargs)

    def add_atom(self, index, name, atomic_number, element, bond_partners, **kwargs):
        """Add an atom to the topology graph
        Parameters
        ----------
        index: int
            The index of the atom in the topology graph
        name: str
            The name of the atom
        atomic_number: int
            The atomic number of the atom
        element: str
            The element associated with the atom
        bond_partners: list of int
            The list of length two containing bond partner atom indices for this atom
        **kwargs
            Any extra information for this atom

        See Also
        --------
        foyer.topology_graph.AtomData
            The class used to store atom data.
        """
        atom_data = AtomData(
            index, name, atomic_number, element, bond_partners, **kwargs
        )
        self.add_node(index, atom_data=atom_data)

    def add_bond(self, atom_1_index, atom_2_index):
        """Add a bond(edge) between two atoms in this TopologyGraph

        Parameters
        ----------
        atom_1_index: int
            The index of the first atom that forms this bond
        atom_2_index: int
            The index of the second atom that forms this bond
        """
        self.add_edge(atom_1_index, atom_2_index)

    def atoms(self, data=False):
        if data:
            for idx, data in self.nodes(data=data):
                yield idx, data["atom_data"]
        else:
            for idx in self.nodes(data=data):
                yield idx

    @classmethod
    def from_parmed(cls, structure: Structure) -> nx.Graph:
        """Return a TopologyGraph with relevant attributes from a parmed Structure

        Parameters
        ----------
        structure: Structure
            The parmed structure

        Returns
        -------
        TopologyGraph
            The equivelent TopologyGraph of the parmed Structure `structure`
        """
        topology_graph = cls()
        for atom in structure.atoms:
            topology_graph.add_atom(
                name=atom.name,
                index=atom.idx,
                atomic_number=atom.atomic_number,
                element=atom.element,
                bond_partners=[bonded_atom.idx for bonded_atom in atom.bond_partners],
            )

        for bond in structure.bonds:
            topology_graph.add_bond(bond.atom1.idx, bond.atom2.idx)

        return topology_graph
