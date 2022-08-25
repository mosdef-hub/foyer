"""Module to represent chemical systems as graph structures."""
from typing import TYPE_CHECKING, Optional

import networkx as nx
from parmed import Structure

from foyer.exceptions import FoyerError

if TYPE_CHECKING:
    from openff.toolkit.topology import Topology as OpenFFTopology

    from foyer.utils.io import import_

    gmso = import_("gmso")


class AtomData:
    """Store atom data necessary for atom typing.

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
    """A general TopologyGraph.

    This class subclasses nx.Graph to provide a general
    topology Graph for atom typing in foyer. Each node in
    this graph is identified by atom indices and edge in this
    graph represents a bond between the atoms.
    """

    def __init__(self, *args, **kwargs):
        super(TopologyGraph, self).__init__(*args, **kwargs)

    def add_atom(
        self,
        index: int,
        name: str,
        atomic_number: Optional[int] = None,
        # consider renaming to `symbol`, `element` implies a custom object
        element: Optional[str] = None,
        **kwargs,
    ):
        """Add an atom to the topology graph.

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
                "both the atomic_number or element symbol for the atom"
            )

        atom_data = AtomData(index, name, atomic_number, element, **kwargs)
        self.add_node(index, atom_data=atom_data)

    def add_bond(self, atom_1_index, atom_2_index):
        """Add a bond(edge) between two atoms in this TopologyGraph.

        Parameters
        ----------
        atom_1_index: int
            The index of the first atom that forms this bond
        atom_2_index: int
            The index of the second atom that forms this bond
        """
        self.add_edge(atom_1_index, atom_2_index)

    def atoms(self, data=False):
        """Iterate through atoms in the TopologyGraph."""
        if data:
            for idx, data in self.nodes(data=data):
                yield idx, data["atom_data"]
        else:
            for idx in self.nodes(data=data):
                yield idx

    def add_bond_partners(self):
        """Add atom indices for atoms involved in a bond."""
        for atom_idx, data in self.nodes(data=True):
            data["bond_partners"] = list(self.neighbors(atom_idx))

    @classmethod
    def from_parmed(cls, structure: Structure):
        """Return a TopologyGraph with relevant attributes from a parmed Structure.

        Parameters
        ----------
        structure: Structure
            The parmed structure

        Returns
        -------
        TopologyGraph
            The equivalent TopologyGraph of the parmed Structure `structure`
        """
        if not isinstance(structure, Structure):
            raise TypeError(
                f"Expected `structure` to be of type {Structure}. "
                f"Got {type(structure).__name__} instead"
            )

        topology_graph = cls()
        for atom in structure.atoms:
            if atom.name.startswith("_"):
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
    def from_openff_topology(cls, openff_topology: "OpenFFTopology"):
        """Return a TopologyGraph with relevant attributes from an openForceField topology.

        Parameters
        ----------
        openff_topology: openff.toolkit.Topology
            The openFF Topology

        Returns
        -------
        TopologyGraph

            The equivalent TopologyGraph of the openFF Topology `openff_topology`
        """
        try:
            from openff.toolkit.topology import Topology
        except ImportError as e:
            raise ImportError(
                "`TopologyGraph.from_openff_topology` requires that "
                "the OpenFF Toolkit is installed."
            ) from e

        if not isinstance(openff_topology, Topology):
            raise TypeError(
                f"Expected `openff_topology` to be of type {Topology}. "
                f"Got {type(openff_topology).__name__} instead"
            )

        uses_old_api = hasattr(Topology(), "_topology_molecules")

        top_graph = cls()

        if uses_old_api:
            from parmed import periodic_table as pt

            for top_atom in openff_topology.topology_atoms:
                atom = top_atom.atom
                element_symbol = pt.Element[atom.atomic_number]
                top_graph.add_atom(
                    name=atom.name,
                    index=top_atom.topology_atom_index,
                    atomic_number=atom.atomic_number,
                    element=element_symbol,
                )

            for top_bond in openff_topology.topology_bonds:
                atoms_indices = [
                    atom.topology_atom_index for atom in top_bond.atoms
                ]
                top_graph.add_bond(atoms_indices[0], atoms_indices[1])

            return top_graph

        else:
            from openff.units.elements import SYMBOLS

            for atom in openff_topology.atoms:
                atom_index = openff_topology.atom_index(atom)
                element_symbol = SYMBOLS[atomic_number]
                top_graph.add_atom(
                    name=atom.name,
                    index=atom_index,
                    atomic_number=atom.atomic_number,
                    element=element_symbol,
                )

            for bond in openff_topology.bonds:
                atoms_indices = [
                    openff_topology.atom_index(atom) for atom in bond.atoms
                ]
                top_graph.add_bond(*atoms_indices)

            return top_graph

    @classmethod
    def from_gmso_topology(cls, gmso_topology: "gmso.Topology"):
        """Return a TopologyGraph with relevant attributes from an GMSO topology.

        Parameters
        ----------
        gmso_topology: gmso.Topology
            The GMSO Topology

        Returns
        -------
        TopologyGraph
            The equivalent TopologyGraph of the openFF Topology `openff_topology`
        """
        try:
            import gmso
        except ImportError as e:
            raise ImportError(
                "`TopologyGraph.from_gmso_topology` requires that GMSO is installed."
            ) from e

        if not isinstance(gmso_topology, gmso.Topology):
            raise TypeError(
                f"Expected `gmso_topology` to be of type {gmso.Topology}. "
                f"Got {type(gmso_topology).__name__} instead"
            )

        top_graph = cls()
        for atom in gmso_topology.sites:
            if isinstance(atom, gmso.Atom):
                if atom.name.startswith("_"):
                    top_graph.add_atom(
                        name=atom.name,
                        index=gmso_topology.get_index(atom),
                        atomic_number=None,
                        element=atom.name,
                    )

                else:
                    top_graph.add_atom(
                        name=atom.name,
                        index=gmso_topology.get_index(atom),
                        atomic_number=atom.element.atomic_number,
                        element=atom.element.symbol,
                    )

        for top_bond in gmso_topology.bonds:
            atoms_indices = [
                gmso_topology.get_index(atom)
                for atom in top_bond.connection_members
            ]
            top_graph.add_bond(atoms_indices[0], atoms_indices[1])

        return top_graph
