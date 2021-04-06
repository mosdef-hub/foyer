import networkx as nx
import pytest

from foyer.atomtyper import find_atomtypes
from foyer.forcefield import Forcefield
from foyer.topology_graph import TopologyGraph
from foyer.utils.io import has_gmso, has_openff_toolkit


@pytest.mark.skipif(
    condition=not (has_gmso or has_openff_toolkit),
    reason="openff-toolkit and gmso not installed",
)
class TestTopologyGraph:
    @pytest.fixture(scope="session")
    def openff_topology_graph(self):
        from openff.toolkit.topology import Molecule, Topology

        openff_ethane = Molecule.from_smiles("CC")
        openff_ethane = Topology.from_molecules(openff_ethane)
        return TopologyGraph.from_openff_topology(openff_ethane)

    @pytest.fixture(scope="session")
    def gmso_topology_graph(self):
        import mbuild as mb
        from gmso.external import from_mbuild

        ethane = mb.load("CC", smiles=True)
        return TopologyGraph.from_gmso_topology(from_mbuild(ethane))

    @pytest.fixture(scope="session")
    def parmed_topology_graph(self):
        # Got errors for parmed.rdkit.from_smiles
        # [23:55:10] Molecule does not have explicit Hs. Consider calling AddHs
        import mbuild as mb

        ethane = mb.conversion.to_parmed(mb.load("CC", smiles=True))
        return TopologyGraph.from_parmed(ethane)

    def test_graph_equivalence(
        self, openff_topology_graph, gmso_topology_graph, parmed_topology_graph
    ):
        assert nx.is_isomorphic(openff_topology_graph, gmso_topology_graph)
        assert nx.is_isomorphic(gmso_topology_graph, parmed_topology_graph)
        assert nx.is_isomorphic(openff_topology_graph, parmed_topology_graph)

    def test_graph_atomdata_equivalence(
        self, openff_topology_graph, gmso_topology_graph, parmed_topology_graph
    ):
        atom_data_gmso_dict = {}
        atom_data_openff_dict = {}
        atom_data_parmed_dict = {}
        for (
            (openff_idx, atom_data_openff),
            (gmso_idx, atom_data_gmso),
            (parmed_idx, atom_data_parmed),
        ) in zip(
            openff_topology_graph.atoms(data=True),
            gmso_topology_graph.atoms(data=True),
            parmed_topology_graph.atoms(data=True),
        ):
            atom_data_openff_dict[openff_idx] = {
                "index": atom_data_openff.index,
                "element": atom_data_openff.element,
                "atomic_number": atom_data_openff.atomic_number,
            }

            atom_data_gmso_dict[gmso_idx] = {
                "index": atom_data_gmso.index,
                "element": atom_data_gmso.element,
                "atomic_number": atom_data_gmso.atomic_number,
            }

            atom_data_parmed_dict[parmed_idx] = {
                "index": atom_data_parmed.index,
                "element": atom_data_parmed.element,
                "atomic_number": atom_data_parmed.atomic_number,
            }

        idx = 0
        while True:
            try:
                assert (
                    atom_data_openff_dict[idx]
                    == atom_data_openff_dict[idx]
                    == atom_data_parmed_dict[idx]
                )
            except KeyError:
                break
            idx += 1

    def test_atom_typing(
        self, openff_topology_graph, gmso_topology_graph, parmed_topology_graph
    ):
        # ToDo: More robust testing for atomtyping
        opls = Forcefield(name="oplsaa")
        openff_typemap = find_atomtypes(openff_topology_graph, forcefield=opls)
        gmso_typemap = find_atomtypes(gmso_topology_graph, forcefield=opls)
        parmed_typemap = find_atomtypes(parmed_topology_graph, forcefield=opls)
        assert openff_typemap
        assert gmso_typemap
        assert parmed_typemap

    def test_from_type_error(self):
        with pytest.raises(TypeError):
            TopologyGraph.from_openff_topology('NonOpenFFTopology')

        with pytest.raises(TypeError):
            TopologyGraph.from_gmso_topology('NonGMSOTopology')

        with pytest.raises(TypeError):
            TopologyGraph.from_parmed('NonParmedStructure')
