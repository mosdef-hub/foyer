"""Module to represent chemical systems as graph structures."""
import enum
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


class BondOrder(enum.Enum):
    """Enum to represent various bond orders from multiple sources."""

    SINGLE = "1"
    DOUBLE = "2"
    TRIPLE = "3"
    AMIDE = "am"
    AROMATIC = "ar"
    DUMMY = "du"
    NOTCONNECTED = "nc"
    UNKNOWN = "un"

    @classmethod
    def _missing_(cls, value):
        """If value cannot be found, default to unknown."""
        return cls.UNKNOWN


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
        symbol: Optional[str] = None,
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
        symbol: optional, str, default=None
            The element symbol associated with the atom
            if it represents an element
        **kwargs
            Any extra information for this atom

        See Also
        --------
        foyer.topology_graph.AtomData
            The class used to store atom data
        """
        if not name.startswith("_") and not (atomic_number and symbol):
            raise FoyerError(
                "For atoms representing an element, please include "
                "both the atomic_number or symbol for the atom"
            )

        atom_data = AtomData(index, name, atomic_number, symbol, **kwargs)
        self.add_node(index, atom_data=atom_data)

    def add_bond(self, atom_1_index, atom_2_index, bond_type=BondOrder.SINGLE):
        """Add a bond(edge) between two atoms in this TopologyGraph.

        Parameters
        ----------
        atom_1_index: int
            The index of the first atom that forms this bond
        atom_2_index: int
            The index of the second atom that forms this bond
        bond_type: str, default = BondOrder.SINGLE
            The type of bond being added, can be:
            "1": single,
            "2": double,
            "3": triple,
            "am": amide,
            "ar": aromatic,
            "un": unknown,
            "du": dummy
            "nc": not connected
        """
        self.add_edge(atom_1_index, atom_2_index, {"bond_type": bond_type})

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
                symbol = None
            else:
                atomic_number = atom.atomic_number
                symbol = atom.element_name

            topology_graph.add_atom(
                name=atom.name,
                index=atom.idx,
                atomic_number=atomic_number,
                symbol=symbol,
            )
        bond_type_dict = {
            round(1, 2): BondOrder.SINGLE,
            round(2, 2): BondOrder.DOUBLE,
            round(3, 2): BondOrder.TRIPLE,
            round(1.5, 2): BondOrder.AROMATIC,
            round(1.25, 2): BondOrder.AMIDE,
        }
        for bond in structure.bonds:
            topology_graph.add_bond(
                bond.atom1.idx,
                bond.atom2.idx,
                bond_type=bond_type_dict.get(
                    round(bond.order, 2), BondOrder.UNKNOWN
                ),
            )

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

        top_graph = cls()

        from openff.units.elements import SYMBOLS

        for atom in openff_topology.atoms:
            atom_index = openff_topology.atom_index(atom)
            element_symbol = SYMBOLS[atom.atomic_number]
            top_graph.add_atom(
                name=atom.name,
                index=atom_index,
                atomic_number=atom.atomic_number,
                symbol=element_symbol,
            )

        bond_type_dict = {
            "Single": BondOrder.SINGLE,
            "Double": BondOrder.DOUBLE,
            "Triple": BondOrder.TRIPLE,
            "Aromatic": BondOrder.AROMATIC,
            "Amide": BondOrder.AMIDE,
        }

        for bond in openff_topology.bonds:
            atoms_indices = [
                openff_topology.atom_index(atom) for atom in bond.atoms
            ]
            top_graph.add_bond(
                atoms_indices[0],
                atoms_indices[1],
                bond_type=bond_type_dict.get(top_bond.type, BondOrder.UNKNOWN),
            )

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
                        symbol=atom.name,
                    )

                else:
                    top_graph.add_atom(
                        name=atom.name,
                        index=gmso_topology.get_index(atom),
                        atomic_number=atom.element.atomic_number,
                        symbol=atom.element.symbol,
                    )
        # until gmso has bond orders, assume single bonds
        for top_bond in gmso_topology.bonds:
            atoms_indices = [
                gmso_topology.get_index(atom)
                for atom in top_bond.connection_members
            ]
            top_graph.add_bond(
                atoms_indices[0], atoms_indices[1], bond_type=BondOrder.SINGLE
            )

        return top_graph
