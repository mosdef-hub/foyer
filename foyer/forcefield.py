import glob
import itertools
import os
from pkg_resources import resource_filename
import warnings

import mbuild as mb
import networkx as nx
import numpy as np
import parmed as pmd
import simtk.openmm.app.element as elem
import simtk.unit as u
from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app.forcefield import (NoCutoff, CutoffNonPeriodic, HBonds,
                                         AllBonds, HAngles, NonbondedGenerator,
                                         _convertParameterToNumber)

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError
from foyer import smarts


def generate_topology(non_omm_topology, non_element_types=None):
    """Create an OpenMM Topology from another supported topology structure. """
    if non_element_types is None:
        non_element_types = set()

    if isinstance(non_omm_topology, pmd.Structure):
        return _topology_from_parmed(non_omm_topology, non_element_types)
    elif isinstance(non_omm_topology, mb.Compound):
        return _topology_from_mbuild(non_omm_topology, non_element_types)
    else:
        raise FoyerError('Unknown topology format: {}\n'
                         'Supported formats are: '
                         '"parmed.Structure", '
                         '"mbuild.Compound", '
                         '"openmm.app.Topology"'.format(non_omm_topology))


def _topology_from_parmed(structure, non_element_types):
    """Convert a ParmEd Structure to an OpenMM Topology. """
    topology = app.Topology()
    chain = topology.addChain()

    residue = topology.addResidue(structure.title, chain)
    atoms = dict()  # pmd.Atom: omm.Atom

    for pmd_atom in structure.atoms:
        name = pmd_atom.name
        if pmd_atom.element in non_element_types:
            element = non_element_types[pmd_atom.element]
        else:
            if pmd_atom.atomic_number != 0:
                element = elem.Element.getByAtomicNumber(pmd_atom.atomic_number)
            else:
                element = elem.Element.getBySymbol(pmd_atom.name)

        omm_atom = topology.addAtom(name, element, residue)
        atoms[pmd_atom] = omm_atom
        omm_atom.bond_partners = []

    for bond in structure.bonds:
        atom1 = atoms[bond.atom1]
        atom2 = atoms[bond.atom2]
        topology.addBond(atom1, atom2)
        atom1.bond_partners.append(atom2)
        atom2.bond_partners.append(atom1)
    if structure.box_vectors and np.any([x._value for x in structure.box_vectors]):
        topology.setPeriodicBoxVectors(structure.box_vectors)

    positions = structure.positions
    return topology, positions


def _topology_from_mbuild(compound, non_element_types):
    """Convert an mBuild Compound to an OpenMM Topology. """
    topology = app.Topology()
    chain = topology.addChain()

    residue = topology.addResidue(compound.name, chain)
    atoms = dict()  # mb.Particle: omm.Atom

    for mb_particle in compound.particles():
        name = mb_particle.name
        if mb_particle.name in non_element_types:
            element = non_element_types[mb_particle.name]
        else:
            element = elem.Element.getBySymbol(mb_particle.name)

        omm_atom = topology.addAtom(name, element, residue)
        atoms[mb_particle] = omm_atom
        omm_atom.bond_partners = []

    for bond in compound.bonds():
        atom1 = atoms[bond[0]]
        atom2 = atoms[bond[1]]
        topology.addBond(atom1, atom2)
        atom1.bond_partners.append(atom2)
        atom2.bond_partners.append(atom1)
    bounding_box = compound.boundingbox
    box_vectors = [[0, 0, 0] * u.nanometer,
                   [0, 0, 0] * u.nanometer,
                   [0, 0, 0] * u.nanometer]
    box_vectors[0][0] = (compound.periodicity[1] or
                         bounding_box.lengths[0]) * u.nanometer
    box_vectors[1][1] = (compound.periodicity[1] or
                         bounding_box.lengths[1]) * u.nanometer
    box_vectors[2][2] = (compound.periodicity[2] or
                         bounding_box.lengths[2]) * u.nanometer
    topology.setPeriodicBoxVectors(box_vectors)

    positions = compound.xyz
    return topology, positions


class Forcefield(app.ForceField):
    """

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, default=None
        Name of a forcefield to load that is packaged within foyer.


    """
    def __init__(self, forcefield_files=None, name=None):
        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self._included_forcefields = dict()
        self.non_element_types = dict()

        all_files_to_load = []
        if forcefield_files is not None:
            if isinstance(forcefield_files, (list, tuple, set)):
                for file in forcefield_files:
                    all_files_to_load.append(file)
            else:
                all_files_to_load.append(forcefield_files)

        if name is not None:
            try:
                file = self.included_forcefields[name]
            except KeyError:
                raise IOError('Forcefield {} cannot be found'.format(name))
            else:
                all_files_to_load.append(file)

        super(Forcefield, self).__init__(*all_files_to_load)

        self.parser = smarts.SMARTS(self.non_element_types)

    @property
    def included_forcefields(self):
        if any(self._included_forcefields):
            return self._included_forcefields

        resource_extension = '.xml'
        ff_dir = resource_filename('foyer', 'forcefields')
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, '*'+resource_extension)))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename = ff_file[:-len(resource_extension)]
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

    def _create_element(self, element, mass):
        if not isinstance(element, elem.Element):
            try:
                element = elem.get_by_symbol(element)
            except KeyError:
                # Enables support for non-atomistic "element types"
                if element not in self.non_element_types:
                    warnings.warn('Non-atomistic element type detected. '
                                  'Creating custom element for {}'.format(element))
                element = elem.Element(number=0,
                                       mass=mass,
                                       name=element,
                                       symbol=element)
        return element

    def registerAtomType(self, parameters):
        """Register a new atom type. """
        name = parameters['name']
        if name in self._atomTypes:
            raise ValueError('Found multiple definitions for atom type: '+name)
        atomClass = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element = self._create_element(parameters['element'], mass)
            self.non_element_types[element.name] = element

        self._atomTypes[name] = self.__class__._AtomType(name, atomClass, mass, element)
        if atomClass in self._atomClasses:
            typeSet = self._atomClasses[atomClass]
        else:
            typeSet = set()
            self._atomClasses[atomClass] = typeSet
        typeSet.add(name)
        self._atomClasses[''].add(name)

        name = parameters['name']
        if 'def' in parameters:
            self.atomTypeDefinitions[name] = parameters['def']
        if 'overrides' in parameters:
            overrides = set(parameters['overrides'].split(","))
            if overrides:
                self.atomTypeOverrides[name] = overrides
        if 'des' in parameters:
            self.atomTypeDesc[name] = parameters['desc']

    def apply(self, topology, *args, **kwargs):
        if not isinstance(topology, app.Topology):
            topology, positions = generate_topology(topology, self.non_element_types)
        else:
            pass
        box_vectors = topology.getPeriodicBoxVectors()
        system = self.createSystem(topology, *args, **kwargs)
        structure = pmd.openmm.load_topology(topology=topology, system=system)
        structure.positions = positions
        if box_vectors is not None:
            structure.box_vectors = box_vectors
        return structure

    def createSystem(self, topology, nonbondedMethod=NoCutoff,
                     nonbondedCutoff=1.0*u.nanometer, constraints=None,
                     rigidWater=True, removeCMMotion=True, hydrogenMass=None,
                     **args):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """

        # Atomtype the system.
        G = nx.Graph()
        G.add_nodes_from(topology.atoms())
        G.add_edges_from(topology.bonds())
        cycles = nx.cycle_basis(G)

        for atom in topology.atoms():
            atom.cycles = set()

        for cycle in cycles:
            for atom in cycle:
                atom.cycles.add(tuple(cycle))

        find_atomtypes(atoms=list(topology.atoms()), forcefield=self)

        data = app.ForceField._SystemData()
        data.atoms = list(topology.atoms())
        for atom in data.atoms:
            data.excludeAtomWith.append([])

        # Make a list of all bonds
        for bond in topology.bonds():
            data.bonds.append(app.ForceField._BondData(bond[0].index, bond[1].index))

        # Record which atoms are bonded to each other atom
        bondedToAtom = []
        for i in range(len(data.atoms)):
            bondedToAtom.append(set())
            data.atomBonds.append([])
        for i in range(len(data.bonds)):
            bond = data.bonds[i]
            bondedToAtom[bond.atom1].add(bond.atom2)
            bondedToAtom[bond.atom2].add(bond.atom1)
            data.atomBonds[bond.atom1].append(i)
            data.atomBonds[bond.atom2].append(i)

        # TODO: Better way to lookup nonbonded parameters...?
        nonbonded_params = None
        for generator in self.getGenerators():
            if isinstance(generator, NonbondedGenerator):
                nonbonded_params = generator.params.paramsForType
                break

        for chain in topology.chains():
            for res in chain.residues():
                for atom in res.atoms():
                    data.atomType[atom] = atom.id
                    if nonbonded_params:
                        params = nonbonded_params[atom.id]
                        data.atomParameters[atom] = params

        # Create the System and add atoms
        sys = mm.System()
        for atom in topology.atoms():
            # Look up the atom type name, returning a helpful error message if it cannot be found.
            if atom not in data.atomType:
                raise Exception("Could not identify atom type for atom '%s'." % str(atom))
            typename = data.atomType[atom]

            # Look up the type name in the list of registered atom types, returning a helpful error message if it cannot be found.
            if typename not in self._atomTypes:
                msg  = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
                msg += "Known atom types are: %s" % str(self._atomTypes.keys())
                raise Exception(msg)

            # Add the particle to the OpenMM system.
            mass = self._atomTypes[typename].mass
            sys.addParticle(mass)

        # Adjust hydrogen masses if requested.
        if hydrogenMass is not None:
            if not u.is_quantity(hydrogenMass):
                hydrogenMass *= u.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.
        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Make a list of all unique angles
        uniqueAngles = set()
        for bond in data.bonds:
            for atom in bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions
        uniquePropers = set()
        for angle in data.angles:
            for atom in bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions
        for atom in range(len(bondedToAtom)):
            bondedTo = bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    data.impropers.append((atom, subset[0], subset[1], subset[2]))

        # Identify bonds that should be implemented with constraints
        if constraints == AllBonds or constraints == HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.name.startswith('H') or atom2.name.startswith('H')
        if rigidWater:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH':
                    bond.isConstrained = True

        # Identify angles that should be implemented with constraints
        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.name.startswith('H'):
                    numH += 1
                if atom3.name.startswith('H'):
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.name.startswith('O')))
        else:
            data.isAngleConstrained = len(data.angles)*[False]
        if rigidWater:
            for i in range(len(data.angles)):
                angle = data.angles[i]
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                    data.isAngleConstrained[i] = True

        # Add virtual sites
        for atom in data.virtualSites:
            (site, atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(atoms[0], atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                sys.setVirtualSite(index, mm.LocalCoordinatesSite(atoms[0], atoms[1], atoms[2],
                                                                  mm.Vec3(site.originWeights[0], site.originWeights[1], site.originWeights[2]),
                                                                  mm.Vec3(site.xWeights[0], site.xWeights[1], site.xWeights[2]),
                                                                  mm.Vec3(site.yWeights[0], site.yWeights[1], site.yWeights[2]),
                                                                  mm.Vec3(site.localPos[0], site.localPos[1], site.localPos[2])))

        # Add forces to the System
        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.
        for script in self._scripts:
            exec(script, locals())
        return sys
