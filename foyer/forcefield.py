import itertools
import os

import networkx as nx
import parmed.gromacs as gmx
from parmed.gromacs import GromacsTopologyFile
import parmed as pmd
from six import string_types

from foyer.atomtyper import find_atomtypes


def apply_forcefield(structure, forcefield, debug=True):
    """Apply a forcefield to a Topology. """
    if isinstance(forcefield, string_types):
        if forcefield.lower() in ['opls-aa', 'oplsaa', 'opls']:
            ff_path = os.path.join(gmx.GROMACS_TOPDIR, 'oplsaa.ff/forcefield.itp')
        elif forcefield.lower() in ['trappeua']:
            ff_path = os.path.join(gmx.GROMACS_TOPDIR, 'trappeua.ff/forcefield.itp')
        else:
            raise ValueError("Unsupported forcefield: '{0}'".format(forcefield))
        ff = GromacsTopologyFile(ff_path, parametrize=False)
    elif isinstance(forcefield, pmd.Structure):
        ff = forcefield

    find_atomtypes(structure.atoms, forcefield, debug=debug)

    ff.atoms = structure.atoms
    ff.bonds = structure.bonds
    create_bonded_forces(ff)
    ff.parametrize()
    return ff


def create_bonded_forces(structure, angles=True, dihedrals=True,
                         impropers=False):
    """Convert Bonds to ForcefieldBonds and find angles and dihedrals. """
    bondgraph = nx.Graph()
    bondgraph.add_edges_from(((b.atom1, b.atom2) for b in structure.bonds))

    if any([angles, dihedrals, impropers]):
        for node_1 in bondgraph.nodes_iter():
            neighbors_1 = bondgraph.neighbors(node_1)
            if len(neighbors_1) > 1:
                if angles:
                    create_angles(structure, node_1, neighbors_1)
                if dihedrals:
                    for node_2 in neighbors_1:
                        if node_2.idx > node_1.idx:
                            neighbors_2 = bondgraph.neighbors(node_2)
                            if len(neighbors_2) > 1:
                                create_dihedrals(structure, node_1, neighbors_1, node_2, neighbors_2)
                if impropers and len(neighbors_1) >= 3:
                    create_impropers(structure, node_1, neighbors_1)


def create_angles(structure, node, neighbors):
    """Add all possible angles around a node to a structure. """
    for pair in itertools.combinations(neighbors, 2):
        angle = pmd.Angle(pair[0], node, pair[1])
        structure.angles.append(angle)


def create_dihedrals(structure, node_1, neighbors_1, node_2, neighbors_2):
    """Add all possible dihedrals around a pair of nodes to a structure. """
    # We need to make sure we don't remove the node from the neighbor lists
    # that we will be re-using in the following iterations.
    neighbors_1 = set(neighbors_1) - {node_2}
    neighbors_2.remove(node_1)

    for pair in itertools.product(neighbors_1, neighbors_2):
        if pair[0] != pair[1]:
            dihedral = pmd.Dihedral(pair[0], node_1, node_2, pair[1])
            structure.dihedrals.append(dihedral)
            structure.rb_torsions.append(dihedral)


def create_impropers(structure, node, neighbors):
    """Add all possible impropers around a node to a structure. """
    for triplet in itertools.combinations(neighbors, 3):
        improper = pmd.Improper(node, triplet[0], triplet[1], triplet[2])
        structure.impropers.append(improper)
