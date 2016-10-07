import itertools
import os
from warnings import warn
import glob
from functools import lru_cache
from pkg_resources import resource_filename

# import networkx as nx
import parmed.gromacs as gmx
from parmed.gromacs import GromacsTopologyFile
import parmed as pmd
from simtk.openmm.app.forcefield import _convertParameterToNumber
import simtk.openmm.app.element as elem

from six import string_types

from foyer.atomtyper import find_atomtypes
from simtk.openmm import app
import networkx as nx


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
        name = parameters['name']
        if name in self._atomTypes:
            raise ValueError('Found multiple definitions for atom type: '+name)
        atomClass = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element = parameters['element']
            if not isinstance(element, elem.Element):
                element = elem.get_by_symbol(element)
        self._atomTypes[name] = (atomClass, mass, element)
        if atomClass in self._atomClasses:
            typeSet = self._atomClasses[atomClass]
        else:
            typeSet = set()
            self._atomClasses[atomClass] = typeSet
        typeSet.add(name)
        self._atomClasses[''].add(name)

        # foyer requires that ForceField should have a _atomTypeDefinitions property
        if not hasattr(self, '_atomTypeDefinitions'):
            self._atomTypeDefinitions = dict()
        if not hasattr(self, '_atomTypeOverrides'):
            self._atomTypeOverrides = dict()
        if not hasattr(self, '_atomTypeDesc'):
            self._atomTypeDesc = dict()


        if 'def' in parameters:
            self._atomTypeDefinitions[name] = parameters['def']

        if 'overrides' in parameters:
            self._atomTypeOverrides[name] = parameters['overrides']

        if 'des' in parameters:
            self._atomTypeDesc[name] = parameters['desc']

    def apply(self, structure, debug=False):
        """Apply a forcefield to a Topology. """

        if not structure.bonds:
            warn('Structure contains no bonds: \n{}\n'.format(structure))

        self.find_atomtypes(structure.atoms, debug=debug)
        self.create_forces(structure)
        # self.parametrize(structure)

    def find_atomtypes(self, atoms, debug=False):
        return find_atomtypes(atoms, self, debug=debug)

    def create_forces(self, structure, angles=True, dihedrals=True,
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


def create_angles(structure, node, neighbors):
    """Add all possible angles around a node to a structure. """
    for pair in itertools.combinations(neighbors, 2):
        angle = pmd.Angle(pair[0], node, pair[1])
        structure.angles.append(angle)


def create_dihedrals(structure, node_1, neighbors_1, node_2, neighbors_2,
                     pairs=True):
    """Add all possible dihedrals around a pair of nodes to a structure. """
    # We need to make sure we don't remove the node from the neighbor lists
    # that we will be re-using in the following iterations.
    neighbors_1 = set(neighbors_1) - {node_2}
    neighbors_2.remove(node_1)

    for pair in itertools.product(neighbors_1, neighbors_2):
        if pair[0] != pair[1]:
            dihedral = pmd.Dihedral(pair[0], node_1, node_2, pair[1])
            if structure.parameterset.dihedral_types:
                structure.dihedrals.append(dihedral)
            if structure.parameterset.rb_torsion_types:
                structure.rb_torsions.append(dihedral)
            if pairs:
                pair = pmd.NonbondedException(pair[0], pair[1])
                structure.adjusts.append(pair)


def create_impropers(structure, node, neighbors):
    """Add all possible impropers around a node to a structure. """
    for triplet in itertools.combinations(neighbors, 3):
        improper = pmd.Improper(node, triplet[0], triplet[1], triplet[2])
        structure.impropers.append(improper)




def _loadFile(*args, **kwargs):
    slf = args[0]
    args = args[1:]
    return slf.orig_loadFile(*args, **kwargs)


# if __name__ == '__main__':
#     ff = load(os.path.join(os.path.split(__file__)[0], 'ff', 'ff.xml'))
#     print(ff)
