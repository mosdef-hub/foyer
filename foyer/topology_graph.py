import networkx as nx
from parmed import Structure
from gmso import Topology

from foyer.exceptions import FoyerError


class AtomData:
    """Stores atom data necessary for atom typing

    Parameters
    ----------
    index: int
        The index of the atom in the topology graph
    name: str
        The name of the atom
    atomic_number: optional, int, default=None
        The atomic number if the atom represents an element
    element: optional, str, default=None
        The element symbol associated with the atom
        if it represents an element
    **kwargs
        Any extra information for this atom

    ToDo: After 3.6 support removed, convert this to a dataclass
    """

    def __init__(self, index, name, atomic_number=None, element=None, **kwargs):
        self.index = index
        self.name = name
        self.atomic_number = atomic_number
        self.element = element
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

    def add_atom(self, index, name, atomic_number=None, element=None, **kwargs):
        """Add an atom to the topology graph
        Parameters
        ----------
        index: int
            The index of the atom in the topology graph
        name: str
            The name of the atom. For non-element type,
            the name must start with an underscore (_)
        atomic_number: optional, int, default=None
            The atomic number if the atom represents an element
        element: optional, str, default=None
            The element symbol associated with the atom
            if it represents an element
        **kwargs
            Any extra information for this atom

        See Also
        --------
        foyer.topology_graph.AtomData
            The class used to store atom data
        """
        if not name.startswith("_") and not (atomic_number and element):
            raise FoyerError(
                "For atoms representing an element, please include "
                "either the atomic_number or element symbol for the atom"
            )

        atom_data = AtomData(index, name, atomic_number, element, **kwargs)
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

    def add_bond_partners(self):
        for atom_idx, data in self.nodes(data=True):
            data["bond_partners"] = list(self.neighbors(atom_idx))

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
            if atom.name.startswith('_'):
                atomic_number = None
                element = None
            else:
                atomic_number = atom.atomic_number
                element = atom.element_name

            topology_graph.add_atom(
                name=atom.name,
                index=atom.idx,
                atomic_number=atomic_number,
                element=element,
            )

        for bond in structure.bonds:
            topology_graph.add_bond(bond.atom1.idx, bond.atom2.idx)

        return topology_graph

    @classmethod
    def from_gmso(cls, topology: Topology) -> nx.Graph:
        """Return a TopologyGraph with relevant attributes from a gmso Topology

        Parameters
        ----------
        topology : Topology
            The GMSO Topology

        Returns
        -------
        TopologyGraph
            The equivalent TopologyGraph of the GMSO Topology `topology`
        """
        topology_graph = cls()
        for atom in topology.sites:
            if atom.name.startswith('_'):
                atomic_number = None
                element = None
            else:
                atomic_number = atom.element.atomic_number
                element = atom.element.symbol

            topology_graph.add_atom(
                name=atom.name,
                index=topology.get_index(atom),
                atomic_number=atomic_number,
                element=element
            )
        for bond in topology.bonds:
            topology_graph.add_bond(topology.get_index(bond.connection_members[0]),
                                    topology.get_index(bond.connection_members[1]))

        return topology_graph

