import itertools
import os
from warnings import warn
import glob
try:
    from functools import lru_cache
except:
    from functools32 import lru_cache

from parmed.gromacs.gromacstop import GromacsTopologyFile
from parmed.parameters import ParameterSet
from parmed.topologyobjects import AtomType, BondType
from pkg_resources import resource_filename

import parmed as pmd
from simtk.openmm.app.forcefield import _convertParameterToNumber, HarmonicBondGenerator, HarmonicAngleGenerator, \
    NonbondedGenerator
import simtk.openmm.app.element as elem
from foyer.atomtyper import find_atomtypes
from simtk.openmm import app
import networkx as nx
import parmed.periodic_table as pt


class Forcefield(app.ForceField):
    def __init__(self, *files):
        super(Forcefield, self).__init__(*files)

    @classmethod
    @lru_cache(maxsize=32)
    def _available(cls):
        resource_extension = '.xml'
        resource_dir = resource_filename('foyer', 'ff')
        resource_files = set(glob.glob(os.path.join(resource_dir, '*'+resource_extension)))

        resource_dict = {}
        for resource_path in resource_files:
            base_path, resource_fn = os.path.split(resource_path)
            basename = resource_fn[:-len(resource_extension)]
            resource_dict[basename] = resource_path

        return resource_dict

    @classmethod
    def available(cls):
        return Forcefield._available().keys()

    @classmethod
    def by_name(cls, forcefield_name):
        if forcefield_name in cls.available():
            return Forcefield(cls._available()[forcefield_name])
        else:
            raise IOError('Forcefield {} cannot be found'.format(forcefield_name))

    def registerAtomType(self, parameters):
        """Register a new atom type."""
        super(Forcefield, self).registerAtomType(parameters)

        # foyer requires that ForceField should have a _atomTypeDefinitions property
        if not hasattr(self, '_atomTypeDefinitions'):
            self._atomTypeDefinitions = dict()
        if not hasattr(self, '_atomTypeOverrides'):
            self._atomTypeOverrides = dict()
        if not hasattr(self, '_atomTypeDesc'):
            self._atomTypeDesc = dict()

        name = parameters['name']
        if 'def' in parameters:
            self._atomTypeDefinitions[name] = parameters['def']

        if 'overrides' in parameters:
            self._atomTypeOverrides[name] = parameters['overrides']

        if 'des' in parameters:
            self._atomTypeDesc[name] = parameters['desc']

    def apply(self, structure, in_place=True, debug=False):
        """Apply a forcefield to a Topology. """

        if not structure.bonds:
            warn('Structure contains no bonds: \n{}\n'.format(structure))

        if not in_place:
            structure = structure.copy(structure.__class__)

        # cast structure to a gromacs topology
        structure = GromacsTopologyFile.from_structure(structure)

        # set atomic numbers if needed
        for atom in structure.atoms:
            if atom.atomic_number == 0:
                if atom.name in pt.AtomicNum:
                    atom.atomic_number = pt.AtomicNum[str(atom.name)]
                else:
                    try:
                        atom.atomic_number = int(atom.name)
                    except ValueError:
                        warn('Cannot set atomic number. Element name: {}\n'.format(atom.name))

        self.find_atomtypes(structure.atoms, debug=debug)
        self.create_bonded_interactions(structure)
        self.parameterize(structure)

        return structure

    def find_atomtypes(self, atoms, debug=False):
        # call the atomtyper
        return find_atomtypes(atoms, self, debug=debug)

    def create_bonded_interactions(self, structure, angles=True, dihedrals=True,
                                   impropers=False, pairs=True):
        """Generate all possible angles, dihedrals and 1-4 pairs. """
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
                                    create_dihedrals(structure, node_1, neighbors_1,
                                                     node_2, neighbors_2, pairs)
                    if impropers and len(neighbors_1) >= 3:
                        create_impropers(structure, node_1, neighbors_1)

    def parameterize(self, structure):
        pass
        # structure.parametrize()


def create_angles(structure, node, neighbors):
    """Add all possible angles around a node to a structure. """
    for pair in itertools.combinations(neighbors, 2):
        angle = pmd.Angle(pair[0], node, pair[1])
        structure.angles.append(angle)


def create_dihedrals(structure, node_1, neighbors_1, node_2, neighbors_2,
                     pairs=True):
    """Add all possible dihedrals around a pair of nodes to a structure. """
    neighbors_1 = set(neighbors_1) - {node_2}
    neighbors_2 = set(neighbors_2) - {node_1}

    for pair in itertools.product(neighbors_1, neighbors_2):
        if pair[0] != pair[1]:
            dihedral = pmd.Dihedral(pair[0], node_1, node_2, pair[1])
            # if hasattr(structure, 'parameterset'):
            #     if structure.parameterset.dihedral_types:
            #         structure.dihedrals.append(dihedral)
            #     if structure.parameterset.rb_torsion_types:
            #         structure.rb_torsions.append(dihedral)
            # else:
            structure.dihedrals.append(dihedral)
            # if pairs:
            #     pair = pmd.NonbondedException(pair[0], pair[1])
            #     structure.adjusts.append(pair)


def create_impropers(structure, node, neighbors):
    """Add all possible impropers around a node to a structure. """
    for triplet in itertools.combinations(neighbors, 3):
        improper = pmd.Improper(node, triplet[0], triplet[1], triplet[2])
        structure.impropers.append(improper)

