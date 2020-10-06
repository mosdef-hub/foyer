import collections
import glob
import itertools
import os
from tempfile import NamedTemporaryFile
import xml.etree.ElementTree as ET
from lxml import etree
from pkg_resources import resource_filename
import warnings
import re

import numpy as np
import foyer.element as custom_elem

import gmso
#from gmso.external.convert_foyer import from_foyer
from gmso.external import from_mbuild

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


# Copy from original forcefield.py
def preprocess_forcefield_files(forcefield_files=None, backend='gmso'):
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

    if backend == 'openmm':
        return preprocessed_files
    elif backend == 'gmso':
        # Run through the forcefield XML conversion
        return forcefield_files
    else:
        raise FoyerError('Backend not supported')

class Forcefield(object):
    """General Forcefield object that can be created by either GMSO Forcefield or OpenMM Forcefield

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, None
        Name of a forcefield to load that is packaged within foyer.
    backend : str, optional, default='openmm'
        Name of the backend used to store all the Types' information.
        Can choose between 'openmm' and 'gmso'.

    """
    def __init__(self, forcefield_files=None, name=None,
                       validation=True, backend='gmso',
                       debug=False):
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

        if forcefield_files is not None:
            if isinstance(forcefield_files, (list, tuple, set)):
                all_files_to_load = list(forcefield_files)
            else:
                all_files_to_load = [forcefield_files]

        if name is not None:
            try:
                file = self.included_forcefields[name]
            except KeyError:
                raise IOError('Forcefild {} cannot be found.'.format(name))
            else:
                all_files_to_load = [file]

        # Preprocessed the input files
        preprocessed_files = preprocess_forcefield_files(all_files_to_load, backend=backend)
        if validation:
            for ff_file_name in preprocessed_files:
                Validator(ff_file_name, debug)

        # Load in an internal forcefield object depends on given backend
        if backend == 'gmso':
            self._parse_gmso(*preprocessed_files)
        elif backend == 'openmm':
            self._parse_mm(*preprocessed_files)
        elif backend == 'openff':
            self._parse_ff(*preprocessed_files)
        else:
            raise FoyerError("Unsupported backend")

        #Remove the temporary files afterward
            for ff_file_name in preprocessed_files:
                os.remove(ff_file_name)

        self.parser = smarts.SMARTS(self.non_element_types)

    @property
    def version(self):
        return self._version

    @property
    def name(self):
        return self._name

    # Parse forcefield meta information
    def _parse_gmso(self, forcefield_files):
        """Parse meta fata information when using GMSO as backend
        """
        self.ff = gmso.ForceField(forcefield_files)
        self._version = self.ff.version
        self._name = self.ff.name
        for name, atype in self.ff.atom_types.items():
            self.atomTypeDefinitions[name] = atype.definition
            self.atomTypeOverrides[name] = atype.overrides
            self.atomTypeDesc[name] = atype.description
            self.atomTypeRefs[name] = atype.doi
            self.atomTypeClasses[name] = atype.atomclass
            #self.atomTypeElements[name] = atype.element

    def _parse_mm(self, forcefield_files):
        """Parse meta data information when using OpenMM as backend
        """
        self.ff = app.ForceField(forcefield_files)
        tree = ET.parse(forcefield_files)
        root = tree.getroot()
        self._version = root.attrib.get('version')
        self._name = root.attrib.get('name')
        
        for atypes_group in root.findall('AtomTypes'):
            for atype in atypes_group:
                name = atype.attrib['name']
                if 'def' in atype.attrib:
                    self.atomTypeDefinitions[name] = atype.attrib['def']
                if 'overrides' in atype.attrib:
                    overrides = set(atype_name.strip() for atype_name in
                                    atype.attrib['overrides'].split(','))
                    if overrides:
                        self.atomTypeOverrides[name] = overrides
                if 'desc' in atype.attrib:
                    self.atomTypeDesc[name] = atype.attrib['desc']
                if 'doi' in atype.attrib:
                    dois = set(doi.strip() for doi in
                               atype.attrib['doi'].split(','))
                    self.atomTypeRefs[name] = dois
                if 'element' in atype.attrib:
                    # Could potentially use ele here instead of just a string
                    self.atomTypeElements[name] = atype.attrib['element']
                if 'class' in atype.attrib:
                    self.atomTypeClasses[name] = atype.attrib['class']
        return None

    def _parse_ff(self, forcefield_files):
        """Parse meta data information when using OpenFF as backend
        """
        self.ff = app.ForceField(forcefield.files)
        return None


    def apply(self, top, references_file=None, use_residue_map=True,
    assert_bond_params=True, assert_angle_params=True,
    assert_dihedral_params=True, assert_improper_params=True,
    verbose=False, *args, **kwargs):
        """Apply the force field to a molecular topology

        Parameters
        ----------
        top : gmso.Topology, or mb.Compound
            Molecular Topology to apply the force field to
        references_file : str, optional, defaut=None
            Name of file where force field refrences will be written (in Bibtex format).
        use_residue_map : bool, optional, default=True
            Options to speed up if there are a lot of repeated
            subtopology within the topology (assuming they all
            have the same name).
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all system dihedrals.
        assert_improper_params : bool, optional, default=True
            If True, Foyer will exit if parametes are not found for all system impropers.
        verbose : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or potential problematic details it encounters.
        """
        if self.atomTypeDefinitions == {}:
            raise FoyerError('Attempting to atom-type using a forcefield '
                'with no atom type definitions.')

        typemap = self._run_atomtyping(top,
                    use_residue_map=use_residue_map,
                    **kwargs)

        return self._parametrize(top=top, typemap=typemap,
            references_file=references_file,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
            verbose=verbose,
            *args, **kwargs)

    def _run_atomtyping(self, top, use_residue_map=True, **kwargs):
        """Atomtype the topology

        Parameters
        ----------
        top : gmso.Topology
            Molecular Topology to find atom types for.
        use_residue_map : bool, optional, default=True
            Spped up options for duplicates subtopology.
        """

        if use_residue_map:
            # Detect duplicates subtopology/residues
            # (do matching by name, assert same number
            # of atoms)
            typemap = find_atomtypes(top, forcefield=self)
        else:
            typemap = find_atomtypes(top, forcefield=self)

        return typemap

    def _parametrize(self, top=None,
                        typemap=dict(),
                        references_file=None,
                        assert_bond_params=True,
                        assert_angle_params=True,
                        assert_dihedral_params=True,
                        assert_improper_params=True,
                        verbose=False,
                        **kwargs):
        """Parametrize the Topology from the typemap provided

        Assign AtomTypes and BondTypes to Atoms and Bonds, respectively.
        Creat Angles, Dihedrals, Impropers and assing corresponding
        AngleType, DihedralTypes, and ImproperTypes.
        Parameters
        ----------

        Returns
        -------
        """
        # Generate missing angles a
        top.identify_connections()

        # Assign AtomTypes
        for atom in top.sites:
            atom.atom_type = self.ff.atom_types.get(typemap[top.get_index(atom)]['atomtype'])
        if not all([a.atom_type for a in top.sites]):
            raise ValueError('Not all atoms in topology have atom types')

        # Assign BondTypes
        for bond in top.bonds:
            # Still need to deal with permutation
            member_types = sorted([atom.atom_type.name for atom in bond.connection_members])
            btype_name = '~'.join(member_types)
            bond.bond_type = self.ff.bond_types.get(btype_name)

        # Assign AngleTypes
        for angle in top.angles:
            # Still need to deal with permutation
            member_types = [atom.atom_type.name for atom in angle.connection_members]
            agtype_name = '~'.join(member_types)
            angle.angle_type = self.ff.angle_types.get(agtype_name)

        # Assign DihedralTypes
        for dihedral in top.dihedrals:
            # Still need to deal with permutation
            member_types = [atom.atom_type.name for atom in dihedral.connection_members]
            dtype_name = '~'.join(member_types)
            dihedral.dihedral_type = self.ff.dihedral_types.get(dtype_name)

        # Assing ImproperTypes
        for improper in top.impropers:
            # Still need to deal with permutation
            member_types = [atom.atom_type.name for atom in improper.connection_members]
            itype_name = '~'.join(member_types)
            improper.improper_type = self.ff.improper_types.get(itype_name)

        #check_paramters(top, assert_bond_params,
        #                     assert_angle_params,
        #                     assert_dihedral_params,
        #                     assert_improper_params,
        #                     debug)

    def _check_parameters(self, top,
                            assert_bond_params=True,
                            assert_angle_params=True,
                            assert_dihedral_params=True,
                            assert_improper_params=True,
                            debug=False):
        """Check if the parameters are fulling filled

        Parameteres
        -----------
        assert_bond_params : bool, optional, default=True
            Check if all bonds have params
        assert_angle_params : bool, optional, default=True
            Check if all angles have params
        assert_dihedral_params : bool, optional, default=True
            Check if all dihedrals have params
        assert_improper_params : bool, optional, default=True
            Check if all impropers have params
        """
        missing_bond_params = dict()
        missing_angle_params = dict()
        missing_dihedral_params = dict()
        missing_improper_params = dict()

        for bond in top.bonds:
            if not bond.bond_type:
                missing_bond_params[bond.name] = [a.atom_type.name
                for a in bond.connection_members]
        for angle in top.angles:
            if not angle.angle_type:
                missing_angle_params[angle.name] = [a.atom_type.name
                for a in angle.connection_members]
        for dihedral in top.dihedrals:
            if not dihedral.dihedral_type:
                missing_dihedral_params[dihedral.name] = [a.atom_type.name
                for a in dihedral.connection_members]
        for improper in top.impropers:
            if not improper.improper_type:
                missing_improper_params[improper.name] = [a.atom_type.name
                for a in improper.connection_members]

        if debug:
            from pprint import pprint
            if missing_bond_params:
                print('Bonds with missing parameters: ')
                pprint(missing_bond_params)
            if missing_angle_params:
                print('Angles with missing parameters: ')
                pprint(missing_angle_params)
            if missing_dihedral_params:
                print('Dihedral with missing parameters: ')
                pprint(missing_dihedral_params)
            if missing_improper_params:
                print('Improper with missing parameters: ')
                pprint(missing_improper_params)

        if assert_bond_params and missing_bond_params:
            raise FoyerError('Some bonds are missing parameters. '
                             'Change debug=True for more information')
        if assert_angle_params and missing_angle_params:
            raise FoyerError('Some angles are missing parameters. '
                             'Change debug=True for more information')
        if assert_dihedral_params and missing_dihedral_params:
            raise FoyerError('Some dihedrals are missing parameters. '
                             'Change debug=True for more information')
        if assert_improper_params and missing_improper_params:
            raise FoyerError('Some impropers are missing parameters. '
                             'Change debug=True for more information')

