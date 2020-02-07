import collections
import glob
import itertools
import os
from tempfile import NamedTemporaryFile
import xml.etree.ElementTree as ET

from pkg_resources import resource_filename
import warnings
import re

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
from foyer.xml_writer import write_foyer
from foyer.utils.io import import_, has_mbuild
from foyer.utils.external import get_ref


def preprocess_forcefield_files(forcefield_files=None):
    """Pre-process foyer Forcefield XML files"""
    if forcefield_files is None:
        return None

    preprocessed_files = []

    for xml_file in forcefield_files:
        if not hasattr(xml_file, 'read'):
            f = open(xml_file)
            _, suffix = os.path.split(xml_file)
        else:
            f = xml_file
            suffix = ""

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
        temp_file = NamedTemporaryFile(suffix=suffix, delete=False)
        with open(temp_file.name, 'w') as temp_f:
            temp_f.write(xml_contents)

        # append temp file name to list
        preprocessed_files.append(temp_file.name)

    return preprocessed_files


def generate_topology(non_omm_topology, non_element_types=None, residues=None):
    """Create an OpenMM Topology from another supported topology structure."""
    if non_element_types is None:
        non_element_types = set()

    if isinstance(non_omm_topology, pmd.Structure):
        return _topology_from_parmed(non_omm_topology, non_element_types)
    elif has_mbuild:
        mb = import_('mbuild')
        if (non_omm_topology, mb.Compound):
            pmd_comp_struct = non_omm_topology.to_parmed(residues=residues)
            return _topology_from_parmed(pmd_comp_struct, non_element_types)
    else:
        raise FoyerError('Unknown topology format: {}\n'
                         'Supported formats are: '
                         '"parmed.Structure", '
                         '"mbuild.Compound", '
                         '"openmm.app.Topology"'.format(non_omm_topology))


def _structure_from_residue(residue, parent_structure):
    """Convert a ParmEd Residue to an equivalent Structure."""
    structure = pmd.Structure()
    for atom in residue.atoms:
        structure.add_atom(atom, resname=residue.name, resnum=residue.number)

    for bond in parent_structure.bonds:
        if bond.atom1 in residue.atoms and bond.atom2 in residue.atoms:
            structure.bonds.append(bond)

    idx_offset = min([a.idx for a in structure])
    for atom in structure.atoms:
        atom._idx -= idx_offset

    return structure


def _topology_from_parmed(structure, non_element_types):
    """Convert a ParmEd Structure to an OpenMM Topology."""
    topology = app.Topology()
    residues = dict()
    for pmd_residue in structure.residues:
        chain = topology.addChain()
        omm_residue = topology.addResidue(pmd_residue.name, chain)
        # Index ParmEd residues on name & number, no other info i.e. chain
        residues[(pmd_residue.name, pmd_residue.idx)] = omm_residue
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

        omm_atom = topology.addAtom(name, element, residues[(pmd_atom.residue.name, pmd_atom.residue.idx)])
        omm_atom.id = pmd_atom.id
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


def _check_independent_residues(structure):
    """Check to see if residues will constitute independent graphs."""
    for res in structure.residues:
        atoms_in_residue = set([*res.atoms])
        bond_partners_in_residue = [item for sublist in [atom.bond_partners for atom in res.atoms] for item in sublist]
        # Handle the case of a 'residue' with no neighbors
        if not bond_partners_in_residue:
            continue
        if set(atoms_in_residue) != set(bond_partners_in_residue):
            return False
    return True


def _unwrap_typemap(structure, residue_map):
    master_typemap = {atom.idx: {'whitelist': set(), 'blacklist': set(), 'atomtype': None} for atom in structure.atoms}
    for res in structure.residues:
        for res_ref, val in residue_map.items():
            if id(res.name) == id(res_ref):
                for i, atom in enumerate(res.atoms):
                    master_typemap[int(atom.idx)]['atomtype'] = val[i]['atomtype']
    return master_typemap


def _separate_urey_bradleys(system, topology):
    """ Separate urey bradley bonds from harmonic bonds in OpenMM System

    Parameters
    ---------
    topology : openmm.app.Topology
        Molecular structure to find atom types of
    system : openmm System

    """
    atoms = [a for a in topology.atoms()]
    bonds = [b for b in topology.bonds()]
    ub_force = mm.HarmonicBondForce()
    harmonic_bond_force = mm.HarmonicBondForce()
    for force_idx, force in enumerate(system.getForces()):
        if isinstance(force, mm.HarmonicBondForce):
            for bond_idx in range(force.getNumBonds()):
                if ((atoms[force.getBondParameters(bond_idx)[0]],
                    atoms[force.getBondParameters(bond_idx)[1]]) not in bonds and
                    (atoms[force.getBondParameters(bond_idx)[1]],
                     atoms[force.getBondParameters(bond_idx)[0]]) not in bonds):
                        ub_force.addBond(*force.getBondParameters(bond_idx))
                else:
                    harmonic_bond_force.addBond(
                        *force.getBondParameters(bond_idx))
            system.removeForce(force_idx)

    system.addForce(harmonic_bond_force)
    system.addForce(ub_force)


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


def _check_bonds(data, structure, assert_bond_params):
    """Check if any bonds lack paramters."""
    if data.bonds:
        missing = [b for b in structure.bonds
                   if b.type is None]
        if missing:
            nmissing = len(structure.bonds) - len(missing)
            msg = ("Parameters have not been assigned to all bonds. "
                   "Total system bonds: {}, Parametrized bonds: {}"
                   "".format(len(structure.bonds), nmissing))
            _error_or_warn(assert_bond_params, msg)


def _check_angles(data, structure, verbose, assert_angle_params):
    """Check if all angles were found and parametrized."""
    if verbose:
        for omm_ids in data.angles:
            missing_angle = True
            for pmd_angle in structure.angles:
                pmd_ids = (pmd_angle.atom1.idx, pmd_angle.atom2.idx, pmd_angle.atom3.idx)
                if pmd_ids == omm_ids:
                    missing_angle = False
            if missing_angle:
                print("Missing angle with ids {} and types {}.".format(
                    omm_ids, [structure.atoms[idx].type for idx in omm_ids]))

    if data.angles and (len(data.angles) != len(structure.angles)):
        msg = ("Parameters have not been assigned to all angles. Total "
               "system angles: {}, Parameterized angles: {}"
               "".format(len(data.angles), len(structure.angles)))
        _error_or_warn(assert_angle_params, msg)



def _check_dihedrals(data, structure, verbose,
                     assert_dihedral_params, assert_improper_params):
    """Check if all dihedrals, including impropers, were found and parametrized."""
    proper_dihedrals = [dihedral for dihedral in structure.dihedrals
                        if not dihedral.improper]

    if verbose:
        for omm_ids in data.propers:
            missing_dihedral = True
            for pmd_proper in structure.rb_torsions:
                pmd_ids = (pmd_proper.atom1.idx, pmd_proper.atom2.idx, pmd_proper.atom3.idx, pmd_proper.atom4.idx)
                if pmd_ids == omm_ids:
                    missing_dihedral = False
            if missing_dihedral:
                print('missing improper with ids {}'.format(pmd_ids))

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
        self.atomTypeClasses = dict()
        self.atomTypeElements = dict()
        self._included_forcefields = dict()
        self.non_element_types = dict()
        self._version = None
        self._name = None

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
        try:
            super(Forcefield, self).__init__(*preprocessed_files)
        finally:
            for ff_file_name in preprocessed_files:
                os.remove(ff_file_name)

        if isinstance(forcefield_files, str):
            self._version = self._parse_version_number(forcefield_files)
            self._name = self._parse_name(forcefield_files)
        elif isinstance(forcefield_files, list):
            self._version = [self._parse_version_number(f) for f in forcefield_files]
            self._name = [self._parse_name(f) for f in forcefield_files]

        self.parser = smarts.SMARTS(self.non_element_types)
        self._SystemData = self._SystemData()

    @property
    def version(self):
        return self._version

    @property
    def name(self):
        return self._name

    @property
    def included_forcefields(self):
        if any(self._included_forcefields):
            return self._included_forcefields

        ff_dir = resource_filename('foyer', 'forcefields')
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, 'xml/*.xml')))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

    def _parse_version_number(self, forcefield_file):
        with open(forcefield_file, 'r') as f:
            tree = ET.parse(f)
            root = tree.getroot()
            try:
                return root.attrib['version']
            except KeyError:
                warnings.warn(
                    'No force field version number found in force field XML file.'
                )
                return None

    def _parse_name(self, forcefield_file):
        with open(forcefield_file, 'r') as f:
            tree = ET.parse(f)
            root = tree.getroot()
            try:
                return root.attrib['name']
            except KeyError:
                warnings.warn(
                    'No force field name found in force field XML file.'
                )
                return None

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
        if 'element' in parameters:
            self.atomTypeElements[name] = parameters['element']
        if 'class' in parameters:
            self.atomTypeClasses[name] = parameters['class']

    def apply(self, structure, references_file=None, use_residue_map=True,
              assert_bond_params=True, assert_angle_params=True,
              assert_dihedral_params=True, assert_improper_params=False,
              combining_rule='geometric', verbose=False, *args, **kwargs):
        """Apply the force field to a molecular structure

        Parameters
        ----------
        structure : parmed.Structure or mbuild.Compound
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
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system
            proper dihedrals.
        assert_improper_params : bool, optional, default=False
            If True, Foyer will exit if parameters are not found for all system
            improper dihedrals.
        combining_rule : str, optional, default='geometric'
            The combining rule of the system, stored as an attribute of the
            ParmEd structure. Accepted arguments are `geometric` and `lorentz`.
        verbose : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potentially problematic details it encounters.
        """
        if self.atomTypeDefinitions == {}:
            raise FoyerError('Attempting to atom-type using a force field '
                    'with no atom type defitions.')

        if not isinstance(structure, pmd.Structure):
            mb = import_('mbuild')
            if isinstance(structure, mb.Compound):
                structure = structure.to_parmed(**kwargs)

        typemap = self.run_atomtyping(structure, use_residue_map=use_residue_map, **kwargs)

        self._apply_typemap(structure, typemap)

        return self.parametrize_system(structure=structure,
            references_file=references_file, assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params, assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params, combining_rule=combining_rule,
            verbose=verbose, *args, **kwargs)

    def run_atomtyping(self, structure, use_residue_map=True, **kwargs):
        """Atomtype the topology

        Parameters
        ----------
        structure : parmed.structure.Structure
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
            independent_residues = _check_independent_residues(structure)

            if independent_residues:
                residue_map = dict()

                # Need to call this only once and store results for later id() comparisons
                for res_id, res in enumerate(structure.residues):
                    if structure.residues[res_id].name not in residue_map.keys():
                        tmp_res = _structure_from_residue(res, structure)
                        typemap = find_atomtypes(tmp_res, forcefield=self)
                        residue_map[res.name] = typemap

                typemap = _unwrap_typemap(structure, residue_map)

            else:
                typemap = find_atomtypes(structure, forcefield=self)

        else:
            typemap = find_atomtypes(structure, forcefield=self)

        return typemap

    def parametrize_system(self, structure=None,
                           references_file=None, assert_bond_params=True,
                           assert_angle_params=True,
                           assert_dihedral_params=True,
                           assert_improper_params=False,
                           combining_rule='geometric', verbose=False,
                           *args, **kwargs):

        topology, positions = _topology_from_parmed(structure, self.non_element_types)

        system = self.createSystem(topology, *args, **kwargs)

        _separate_urey_bradleys(system, topology)

        data = self._SystemData

        structure = pmd.openmm.load_topology(topology=topology, system=system)
        structure.bonds.sort(key=lambda x: x.atom1.idx)
        structure.positions = positions
        box_vectors = topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            structure.box_vectors = box_vectors

        _check_bonds(data, structure, assert_bond_params)
        _check_angles(data, structure, verbose, assert_angle_params)
        _check_dihedrals(data, structure, verbose,
                              assert_dihedral_params, assert_improper_params)

        if references_file:
            atom_types = set(atom.type for atom in structure.atoms)
            self._write_references_to_file(atom_types, references_file)

        # TODO: Check against the name of the force field and/or store
        # combining rule directly in XML, i.e.
        # if self.name == 'oplsaa':
        structure.combining_rule = combining_rule

        return structure

    def createSystem(self, topology, nonbondedMethod=NoCutoff,
                     nonbondedCutoff=1.0 * u.nanometer, constraints=None,
                     rigidWater=True, removeCMMotion=True, hydrogenMass=None,
                     switchDistance=None,
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
        switchDistance : float=None
            The distance at which the potential energy switching function is turned on for
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        args['switchDistance'] = switchDistance
        # Overwrite previous _SystemData object
        self._SystemData = app.ForceField._SystemData()

        data = self._SystemData
        data.atoms = list(topology.atoms())
        for _ in data.atoms:
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

    def _apply_typemap(self, structure, typemap):
        """Add atomtypes to the topology according to the typemap"""
        for atom in structure.atoms:
            atom.id = typemap[atom.idx]['atomtype']

        if not all([a.id for a in structure.atoms][0]):
            raise ValueError('Not all atoms in topology have atom types')

    def _prepare_structure(self, topology, **kwargs):
        """Separate positions and other topological information"""
        if not isinstance(topology, app.Topology):
            residues = kwargs.get('residues')
            topology, positions = generate_topology(topology,
                                                    self.non_element_types, residues=residues)
        else:
            positions = np.empty(shape=(topology.getNumAtoms(), 3))
            positions[:] = np.nan

        return topology, positions

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
                url = "http://api.crossref.org/works/{}/transform/application/x-bibtex".format(doi)
                headers = {"accept": "application/x-bibtex"}
                bibtex_ref = get_ref(url, headers=headers)
                if bibtex_ref is None:
                    warnings.warn('Could not get ref for doi'.format(doi))
                    continue
                else:
                    bibtex_text = bibtex_ref.text
                note = (',\n\tnote = {Parameters for atom types: ' +
                        ', '.join(sorted(atomtypes)) + '}')
                bibtex_text = bibtex_text[:-2] + note + bibtex_text[-2:]
                f.write('{}\n'.format(bibtex_text))


pmd.Structure.write_foyer = write_foyer
