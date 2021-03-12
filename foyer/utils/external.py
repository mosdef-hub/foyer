import requests
import networkx as nx
from parmed import Structure


def get_ref(ref_url, headers):
    bibtex_ref = requests.get(ref_url, headers=headers)
    if bibtex_ref.ok:
        return bibtex_ref
    else:
        return None


def networkx_from_parmed(structure: Structure) -> nx.Graph:
    """Return a networkx graph with relevant attributes from a parmed Structure"""
    topology_graph = nx.Graph()
    for atom in structure.atoms:
        topology_graph.add_node(
            atom.idx,
            name=atom.name,
            index=atom.idx,
            atomic_number=atom.atomic_number,
            atom=atom
        )
    for bond in structure.bonds:
        topology_graph.add_edge(
            bond.atom1.idx,
            bond.atom2.idx
        )
    return topology_graph

