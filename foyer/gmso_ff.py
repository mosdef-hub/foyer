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

from simtk import openmm as mm

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError
from foyer import smarts
from foyer.validator import Validator
from foyer.xml_writer import write_foyer
from foyer.utils.io import import_, has_mbuild
from foyer.utils.external import get_ref

def preprocessed_forcefield_files(forcefield_files=None):
    """Preprocess foyer Forcefield XML files

    Make sure the provided forcefield XML files have the correct
    format and convert them to GMSO forcefield XML format.
    """
    # Validation can be copied from the old code
    # Converting will rely on Ray's XML conversion PR
    return forcefield_files


class Forcefield(gmso.ForceField):
    """Specialization of GMSO ForceField allowing SMARTS based atomtyping.

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, default=None
        Name of a forcefield to load that is packaged within foyer.
    validation : bool, optional, default=True
        Validate the given forcefield files
    backend : bool, optional, default='gmso'
        Backend used to store ForceField types information.
        Can be picked from either 'gmso' or 'openmm'.
    debug : bool, optional, default=False
        Prompt to print our extra messages for debugging purpose.
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
        preprocessed_files = preprocessed_forcefield_files(all_files_to_load)
        if validation:
            for ff_file_name in preprocessed_files:
                Validator(ff_file_name, debug)

        # Load in an internal forcefield object depends on given backend
        super(Forcefield, self).__init__(*preprocessed_files)
        '''
        if backend=='gmso':
            self.ff = gmso.ForceField(*preprocessed_files)
        elif backend=='openmm':
            self.ff = mm.ForceField(*preprocessed_files)
        else:
            raise Exception("Unsupported backend")
        '''
        # Remove the file
        """
        Remove the converted files afterward
        finally:
            for ff_file_name in preprocessed_files:
                os.remove(ff_file_name)
        """
        for name, atype in self.atom_types.items():
            self.atomTypeDefinitions[name] = atype.definition
            self.atomTypeOverrides[name] = atype.overrides
            self.atomTypeDesc[name] = atype.description
            self.atomTypeRefs[name] = atype.doi
            self.atomTypeClasses[name] = atype.atomclass
            #self.atomTypeElements[name] = atype.element

        self.parser = smarts.SMARTS(self.non_element_types)


    def apply(self, top, reference_file, use_residue_map,
              assert_atom_params=True,
              assert_bond_params=True,
              assert_angle_params=True,
              assert_dihedral_params=True,
              assert_improper_params=True):
        """Apply the forcefield to a molecular structure

        Parameters
        ----------
        top : mbuild.Compound or gmso.Topology
            Molecular structure to apply the forcefield to.
        reference_file : str, optional, default=None
            Name of the file where forcefield refrences (doi)
            will be written to (in Bibtex format).
        use_residue_map :
            Speed up measure for systems with a lot of repeated units.
        assert_atom_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all atoms.
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all proper dihedrals.
        asser_improper_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all improper dihedrals.
        debug : bool, optional, default=False
            If True, Foyer will print out missing atoms, bonds, anlges, and dihedrals
            that are missing types.
        Return
        ------
        topology : gmso.Topology
            Return typed gmso.Topology object
        """
        if self.atomTypeDefinitions == {}:
            raise FoyerError('Attempting to atom-type usinga force field with no atom type definition')

        if not isinstance(top, gmso.Topology):
            mb = import_('mbuild')
            top = from_mbuild(mb.load(top))

        # Generate typemap
        typemap = self.run_atomtyping(top)

        # Apply typemap
        self._apply_typemap(top, typemap)

        # Parameterize
        self.parameterize_system(top)
        return top

    def run_atomtyping(self, top):
        """ Atomtype the topology

        Parameteres
        -----------

        Returns
        -------
        """
        # TO DO: speed up options (similar to use_residue_map)

        typemap = find_atomtypes(top, forcefield=self)
        return typemap

    def parameterize_system(self, top, typemap,
                            assert_atom_params=True,
                            assert_bond_params=True,
                            assert_angle_params=True,
                            assert_dihedral_params=True,
                            assert_improper_params=True,
                            debug=True):

        """ Parameterize systems

        Assign AtomTypes, BondTypes to atoms, bonds. Create
        angles and dihedral and assign corresponding AngleType
        and DihedralTypes

        Parameters
        ----------

        Returns
        -------
        """
        # Generate missing angles, dihedrals, and improper
        top.identify_connections()

        # Assigns AtomTypes
        for idx in typemap:
            top.sites[idx].atom_type = self.atom_types.typemap[idx]['atomtype']

        # Assigns BondTypes
        for bond in top.bonds:
            # Note: still need to deal with permutation
            btype_name = '{}~{}'.format(
                            bond.connection_members[0],
                            bond.connection_members[1])
            bond.bond_type = self.bond_types.get(btype_name)

        # Assign angle type
        for angle in top.angles:
            # Note: still need to deal with permutation
            agtype_name = '{}~{}~{}'.format(
                            angle.connection_members[0],
                            angle.connection_members[1],
                            angle.connection_members[2])
            angle.angle_type = self.angle_types.get(agtype_name)

        # Assign proper dihedral type
        for dihedral in top.dihedrals:
            dtype_name = '{}~{}~{}~{}'.format(
                            dihedral.connection_members[0],
                            dihedral.connection_members[1],
                            dihedral.connection_members[2],
                            dihedral.connection_members[3])
            dihedral.dihedral_type = self.dihedral_types.get(dtype_name)

        # Assign improper type
        for improper in top.impropers:
            itype_name = '{}~{}~{}~{}'.format(
                            improper.connection_members[0],
                            improper.connection_members[1],
                            improper.connection_members[2],
                            improper.connection_members[3])
            improper.dihedral_type = self.improper_types.get(itype_name)

        _check_parameters(top,
            assert_atom_params=assert_atom_params,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
            debug=debug)

        return top

    def _check_parameters(self, top,
                            assert_atom_params=True,
                            assert_bond_params=True,
                            assert_angle_params=True,
                            assert_dihedral_params=True,
                            debug=True):
        """Check if the params are fulling filled

        Parameters
        ----------

        """
        missing_types = dict()
        missing_types['atom'] = [atom for atom in top.sites
                                      if not atom.atom_type]
        missing_types['bond'] = [bond for bond in top.bonds
                                      if not bond.bond_type]
        missing_types['angle'] = [angle for angle in top.angles
                            if not angle.angle_type]
        missing_types['dihedral'] = [dihedral for dihedral in top.dihedral
                            if not dihedral.dihedral_type]
        missing_types['improper'] = [improper for improper in top.improper
                            if not improper.improper_type]
        if debug:
            for component in missing_types:
                if len(missing_types[component]) > 0:
                    print('Detect missing {0} with {0} type: '.format(component))
                    print(missing_types[component])

        # Need to figure out a way to raise a collective Error
        # messages
        if assert_atom_params and missing_types['atom']:
            raise FoyerError
        if assert_bond_params and missing_types['bond']:
            raise FoyerError
        if assert_angle_params and missing_types['angle']:
            raise FoyerError
        if assert_dihedral_params and missing_types['dihedral']:
            raise FoyerError
        if assert_improper_params and missing_types['improper']:
            raise FoyerError

        return None


