import collections
import glob
import itertools
import os
from tempfile import mktemp, mkstemp
import xml.etree.ElementTree as ET

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
from pkg_resources import resource_filename
import requests
import warnings
import re

import mbuild as mb
import numpy as np
import parmed as pmd
import simtk.openmm.app.element as elem
import foyer.element as custom_elem
import simtk.unit as u
from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app.forcefield import (NoCutoff, CutoffNonPeriodic, HBonds,
                                         AllBonds, HAngles, NonbondedGenerator,
                                         _convertParameterToNumber)

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError
from foyer import smarts
from foyer.validator import Validator


def preprocess_forcefield_files(forcefield_files=None):
    if forcefield_files is None:
        return None

    preprocessed_files = []

    for xml_file in forcefield_files:
        if not hasattr(xml_file,'read'):
            f = open(xml_file)
            _,suffix = os.path.split(xml_file)
        else:
            f = xml_file
            suffix=""

        # read and preprocess
        xml_contents = f.read()
        f.close()
        xml_contents = re.sub(r"(def\w*=\w*[\"\'])(.*)([\"\'])", lambda m: m.group(1) + re.sub(r"&(?!amp;)", r"&amp;", m.group(2)) + m.group(3),
                              xml_contents)

        try:
            '''
            Sort topology objects by precedence, defined by the number of
            `type` attributes specified, where a `type` attribute indicates
            increased specificity as opposed to use of `class`
            '''
            root = ET.fromstring(xml_contents)
            for element in root:
                if 'Force' in element.tag:
                    element[:] = sorted(element, key=lambda child: (
                        -1 * len([attr_name for attr_name in child.keys()
                                    if 'type' in attr_name])))
            xml_contents = ET.tostring(root, method='xml').decode()
        except ET.ParseError:
            '''
            Provide the user with a warning if sorting could not be performed.
            This indicates a bad XML file, which will be passed on to the
            Validator to yield a more descriptive error message.
            '''
            warnings.warn('Invalid XML detected. Could not auto-sort topology '
                          'objects by precedence.')

        # write to temp file
        _, temp_file_name = mkstemp(suffix=suffix)
        with open(temp_file_name, 'w') as temp_f:
            temp_f.write(xml_contents)

        # append temp file name to list
        preprocessed_files.append(temp_file_name)

    return preprocessed_files


def generate_topology(non_omm_topology, non_element_types=None,
        residues=None):
    """Create an OpenMM Topology from another supported topology structure."""
    if non_element_types is None:
        non_element_types = set()

    if isinstance(non_omm_topology, pmd.Structure):
        return _topology_from_parmed(non_omm_topology, non_element_types)
    elif isinstance(non_omm_topology, mb.Compound):
        pmdCompoundStructure = non_omm_topology.to_parmed(residues=residues)
        return _topology_from_parmed(pmdCompoundStructure, non_element_types)
    else:
        raise FoyerError('Unknown topology format: {}\n'
                         'Supported formats are: '
                         '"parmed.Structure", '
                         '"mbuild.Compound", '
                         '"openmm.app.Topology"'.format(non_omm_topology))


def _topology_from_parmed(structure, non_element_types):
    """Convert a ParmEd Structure to an OpenMM Topology."""
    topology = app.Topology()
    residues = dict()
    for pmd_residue in structure.residues:
        chain = topology.addChain()
        omm_residue = topology.addResidue(pmd_residue.name, chain)
        residues[pmd_residue] = omm_residue
    atoms = dict()  # pmd.Atom: omm.Atom

    for pmd_atom in structure.atoms:
        name = pmd_atom.name
        if pmd_atom.name in non_element_types:
            element = non_element_types[pmd_atom.name]
        else:
            if (isinstance(pmd_atom.atomic_number, int) and
                    pmd_atom.atomic_number != 0):
                element = elem.Element.getByAtomicNumber(pmd_atom.atomic_number)
            else:
                element = elem.Element.getBySymbol(pmd_atom.name)

        omm_atom = topology.addAtom(name, element, residues[pmd_atom.residue])
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


def _topology_from_residue(res):
    """Converts a openmm.app.Topology.Residue to openmm.app.Topology.

    Parameters
    ----------
    res : openmm.app.Topology.Residue
        An individual residue in an openmm.app.Topology

    Returns
    -------
    topology : openmm.app.Topology
        The generated topology

    """
    topology = app.Topology()
    chain = topology.addChain()
    new_res = topology.addResidue(res.name, chain)

    atoms = dict()  # { omm.Atom in res : omm.Atom in *new* topology }

    for res_atom in res.atoms():
        topology_atom = topology.addAtom(name=res_atom.name,
                         element=res_atom.element,
                         residue=new_res)
        atoms[res_atom] = topology_atom
        topology_atom.bond_partners = []

    for bond in res.bonds():
        atom1 = atoms[bond.atom1]
        atom2 = atoms[bond.atom2]
        topology.addBond(atom1, atom2)
        atom1.bond_partners.append(atom2)
        atom2.bond_partners.append(atom1)

    return topology


def _check_independent_residues(topology):
    """Check to see if residues will constitute independent graphs."""
    for res in topology.residues():
        atoms_in_residue = set([atom for atom in res.atoms()])
        bond_partners_in_residue = [item for sublist in [atom.bond_partners for atom in res.atoms()] for item in sublist]
        # Handle the case of a 'residue' with no neighbors
        if not bond_partners_in_residue:
            continue
        if set(atoms_in_residue) != set(bond_partners_in_residue):
            return False
    return True


def _update_atomtypes(unatomtyped_topology, res_name, prototype):
    """Update atomtypes in residues in a topology using a prototype topology.

    Atomtypes are updated when residues in each topology have matching names.

    Parameters
    ----------
    unatomtyped_topology : openmm.app.Topology
        Topology lacking atomtypes defined by `find_atomtypes`.
    prototype : openmm.app.Topology
        Prototype topology with atomtypes defined by `find_atomtypes`.

    """
    for res in unatomtyped_topology.residues():
        if res.name == res_name:
            for old_atom, new_atom_id in zip([atom for atom in res.atoms()], [atom.id for atom in prototype.atoms()]):
                old_atom.id = new_atom_id

def _error_or_warn(error, msg):
    """Raise an error or warning if topology objects are not fully parameterized.

    Parameters
    ----------
    error : bool
        If True, raise an error, else raise a warning
    msg : str
        The message to be provided with the error or warning
    """
    if error:
        raise Exception(msg)
    else:
        warnings.warn(msg)


class Forcefield(app.ForceField):
    """Specialization of OpenMM's Forcefield allowing SMARTS based atomtyping.

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, default=None
        Name of a forcefield to load that is packaged within foyer.


    """
    def __init__(self, forcefield_files=None, name=None, validation=True, debug=False):
        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self.atomTypeRefs = dict()
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

        preprocessed_files = preprocess_forcefield_files(all_files_to_load)
        if validation:
            for ff_file_name in preprocessed_files:
                Validator(ff_file_name, debug)

        super(Forcefield, self).__init__(*preprocessed_files)
        self.parser = smarts.SMARTS(self.non_element_types)
        self._SystemData = self._SystemData()

    @property
    def included_forcefields(self):
        if any(self._included_forcefields):
            return self._included_forcefields

        ff_dir = resource_filename('foyer', 'forcefields')
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, '*.xml')))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
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
                element = custom_elem.Element(number=0,
                                       mass=mass,
                                       name=element,
                                       symbol=element)
            else:
                return element, False

        return element, True

    def registerAtomType(self, parameters):
        """Register a new atom type. """
        name = parameters['name']
        if name in self._atomTypes:
            raise ValueError('Found multiple definitions for atom type: ' + name)
        atom_class = parameters['class']
        mass = _convertParameterToNumber(parameters['mass'])
        element = None
        if 'element' in parameters:
            element, custom = self._create_element(parameters['element'], mass)
            if custom:
                self.non_element_types[element.symbol] = element

        self._atomTypes[name] = self.__class__._AtomType(name, atom_class, mass, element)
        if atom_class in self._atomClasses:
            type_set = self._atomClasses[atom_class]
        else:
            type_set = set()
            self._atomClasses[atom_class] = type_set
        type_set.add(name)
        self._atomClasses[''].add(name)

        name = parameters['name']
        if 'def' in parameters:
            if parameters['def']:
                self.atomTypeDefinitions[name] = parameters['def']
        if 'overrides' in parameters:
            overrides = set(atype.strip() for atype
                            in parameters['overrides'].split(","))
            if overrides:
                self.atomTypeOverrides[name] = overrides
        if 'des' in parameters:
            self.atomTypeDesc[name] = parameters['desc']
        if 'doi' in parameters:
            dois = set(doi.strip() for doi in parameters['doi'].split(','))
            self.atomTypeRefs[name] = dois

    def apply(self, topology, references_file=None, use_residue_map=True,
              assert_angle_params=True, assert_dihedral_params=True,
              assert_improper_params=False, *args, **kwargs):
        """Apply the force field to a molecular structure

        Parameters
        ----------
        topology : openmm.app.Topology or parmed.Structure or mbuild.Compound
            Molecular structure to apply the force field to
        references_file : str, optional, default=None
            Name of file where force field references will be written (in Bibtex
            format)
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            proper dihedrals.
        assert_improper_params : bool, optional, default=False
            If True, Foyer will exit if parameters are not found for all system
            improper dihedrals.
        """
        if not isinstance(topology, app.Topology):
            residues = kwargs.get('residues')
            topology, positions = generate_topology(topology,
                    self.non_element_types, residues=residues)
        else:
            positions = np.empty(shape=(topology.getNumAtoms(), 3))
            positions[:] = np.nan
        box_vectors = topology.getPeriodicBoxVectors()
        topology = self.run_atomtyping(topology, use_residue_map=use_residue_map)
        system = self.createSystem(topology, *args, **kwargs)

        structure = pmd.openmm.load_topology(topology=topology, system=system)

        '''
        Check that all topology objects (angles, dihedrals, and impropers)
        have parameters assigned. OpenMM will generate an error if bond parameters
        are not assigned.
        '''
        data = self._SystemData

        if data.angles and (len(data.angles) != len(structure.angles)):
            msg = ("Parameters have not been assigned to all angles. Total "
                   "system angles: {}, Parameterized angles: {}"
                   "".format(len(data.angles), len(structure.angles)))
            _error_or_warn(assert_angle_params, msg)

        proper_dihedrals = [dihedral for dihedral in structure.dihedrals
                            if not dihedral.improper]
        if data.propers and len(data.propers) != \
                len(proper_dihedrals) + len(structure.rb_torsions):
            msg = ("Parameters have not been assigned to all proper dihedrals. "
                   "Total system dihedrals: {}, Parameterized dihedrals: {}. "
                   "Note that if your system contains torsions of Ryckaert-"
                   "Bellemans functional form, all of these torsions are "
                   "processed as propers.".format(len(data.propers),
                   len(proper_dihedrals) + len(structure.rb_torsions)))
            _error_or_warn(assert_dihedral_params, msg)

        improper_dihedrals = [dihedral for dihedral in structure.dihedrals
                              if dihedral.improper]
        if data.impropers and len(data.impropers) != \
                len(improper_dihedrals) + len(structure.impropers):
            msg = ("Parameters have not been assigned to all impropers. Total "
                   "system impropers: {}, Parameterized impropers: {}. "
                   "Note that if your system contains torsions of Ryckaert-"
                   "Bellemans functional form, all of these torsions are "
                   "processed as propers".format(len(data.impropers),
                   len(improper_dihedrals) + len(structure.impropers)))
            _error_or_warn(assert_improper_params, msg)

        structure.bonds.sort(key=lambda x: x.atom1.idx)
        structure.positions = positions
        if box_vectors is not None:
            structure.box_vectors = box_vectors
        if references_file:
            atom_types = set(atom.type for atom in structure.atoms)
            self._write_references_to_file(atom_types, references_file)

        return structure

    def run_atomtyping(self, topology, use_residue_map=True):
        """Atomtype the topology

        Parameters
        ----------
        topology : openmm.app.Topology
            Molecular structure to find atom types of
        use_residue_map : boolean, optional, default=True
            Store atomtyped topologies of residues to a dictionary that maps
            them to residue names.  Each topology, including atomtypes, will be
            copied to other residues with the same name. This avoids repeatedly
            calling the subgraph isomorphism on idential residues and should
            result in better performance for systems with many identical
            residues, i.e. a box of water. Note that for this to be applied to
            independent molecules, they must each be saved as different
            residues in the topology.
        """
        if use_residue_map:
            independent_residues = _check_independent_residues(topology)

            if independent_residues:
                residue_map = dict()

                for res in topology.residues():
                    if res.name not in residue_map.keys():
                        residue = _topology_from_residue(res)
                        find_atomtypes(residue, forcefield=self)
                        residue_map[res.name] = residue

                for key, val in residue_map.items():
                    _update_atomtypes(topology, key, val)

            else:
                find_atomtypes(topology, forcefield=self)

        else:
            find_atomtypes(topology, forcefield=self)

        if not all([a.id for a in topology.atoms()][0]):
            raise ValueError('Not all atoms in topology have atom types')

        return topology

    def createSystem(self, topology, nonbondedMethod=NoCutoff,
                     nonbondedCutoff=1.0 * u.nanometer, constraints=None,
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

        # Overwrite previous _SystemData object
        self._SystemData = app.ForceField._SystemData()

        data = self._SystemData
        data.atoms = list(topology.atoms())
        for atom in data.atoms:
            data.excludeAtomWith.append([])

        # Make a list of all bonds
        for bond in topology.bonds():
            data.bonds.append(app.ForceField._BondData(bond[0].index, bond[1].index))

        # Record which atoms are bonded to each other atom
        bonded_to_atom = []
        for i in range(len(data.atoms)):
            bonded_to_atom.append(set())
            data.atomBonds.append([])
        for i in range(len(data.bonds)):
            bond = data.bonds[i]
            bonded_to_atom[bond.atom1].add(bond.atom2)
            bonded_to_atom[bond.atom2].add(bond.atom1)
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
                    transfer_mass = hydrogenMass - sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    mass = sys.getParticleMass(atom1.index) - transfer_mass
                    sys.setParticleMass(atom1.index, mass)

        # Set periodic boundary conditions.
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sys.setDefaultPeriodicBoxVectors(box_vectors[0],
                                             box_vectors[1],
                                             box_vectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a '
                             'Topology that does not specify periodic box '
                             'dimensions')

        # Make a list of all unique angles
        unique_angles = set()
        for bond in data.bonds:
            for atom in bonded_to_atom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        unique_angles.add((atom, bond.atom1, bond.atom2))
                    else:
                        unique_angles.add((bond.atom2, bond.atom1, atom))
            for atom in bonded_to_atom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        unique_angles.add((bond.atom1, bond.atom2, atom))
                    else:
                        unique_angles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(unique_angles))

        # Make a list of all unique proper torsions
        unique_propers = set()
        for angle in data.angles:
            for atom in bonded_to_atom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        unique_propers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        unique_propers.add((angle[2], angle[1], angle[0], atom))
            for atom in bonded_to_atom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        unique_propers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        unique_propers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(unique_propers))

        # Make a list of all unique improper torsions
        for atom in range(len(bonded_to_atom)):
            bonded_to = bonded_to_atom[atom]
            if len(bonded_to) > 2:
                for subset in itertools.combinations(bonded_to, 3):
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
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(
                    atoms[0], atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(
                    atoms[0], atoms[1], atoms[2],
                    site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(
                    atoms[0], atoms[1], atoms[2],
                    site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                local_coord_site = mm.LocalCoordinatesSite(
                    atoms[0], atoms[1], atoms[2],
                    mm.Vec3(site.originWeights[0], site.originWeights[1], site.originWeights[2]),
                    mm.Vec3(site.xWeights[0], site.xWeights[1], site.xWeights[2]),
                    mm.Vec3(site.yWeights[0], site.yWeights[1], site.yWeights[2]),
                    mm.Vec3(site.localPos[0], site.localPos[1], site.localPos[2]))
                sys.setVirtualSite(index, local_coord_site)

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

    def _write_references_to_file(self, atom_types, references_file):
        atomtype_references = {}
        for atype in atom_types:
            try:
                atomtype_references[atype] = self.atomTypeRefs[atype]
            except KeyError:
                warnings.warn("Reference not found for atom type '{}'."
                              "".format(atype))
        unique_references = collections.defaultdict(list)
        for atomtype, dois in atomtype_references.items():
            for doi in dois:
                unique_references[doi].append(atomtype)
        unique_references = collections.OrderedDict(sorted(unique_references.items()))
        with open(references_file, 'w') as f:
            for doi, atomtypes in unique_references.items():
                url = "http://dx.doi.org/" + doi
                headers = {"accept": "application/x-bibtex"}
                bibtex_ref = requests.get(url, headers=headers).text
                note = (',\n\tnote = {Parameters for atom types: ' +
                        ', '.join(sorted(atomtypes)) + '}')
                bibtex_ref = bibtex_ref[:-2] + note + bibtex_ref[-2:]
                f.write('{}\n'.format(bibtex_ref))
