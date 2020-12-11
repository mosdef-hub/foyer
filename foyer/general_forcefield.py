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
from pathlib import Path

import numpy as np
import foyer.element as custom_elem

import gmso
#from gmso.external.convert_foyer import from_foyer
from gmso.external import from_mbuild

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
        preprocessed_files.append(Path(temp_file.name))

    if backend == 'gmso':
        # Run through the forcefield XML conversion
        from gmso.external.convert_foyer import convert_foyer
        for idx, file in enumerate(preprocessed_files):
            tmp_gmso = f'tmp_gmso_{idx}.xml'
            convert_foyer(foyer_xml=file,
                          gmso_xml=tmp_gmso)
            os.remove(file)
            tmp_processed_files.append(tmp_gmso)

    else:
        raise FoyerError('Backend not supported')

    return tmp_processed_files


class Forcefield(object):
    """General Forcefield object that can be created by either a GMSO XML
    forcefield file of Foyer XML forcefield file

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, None
        Name of a forcefield to load that is packaged within foyer.
    backend : str, optional, default='gmso'
        Name of the backend used to store all the Types' information.
        At this point, 'gmso' is the only valid backend, but this set up
        allow future backend to be implemented easier.
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
                raise IOError('Forcefield {} cannot be found.'.format(name))
            else:
                all_files_to_load = [file]

        # Preprocessed the input files

        preprocessed_files = preprocess_forcefield_files(
                                all_files_to_load,
                                backend=backend)

        # Load in an internal forcefield object depends on given backend
        if backend == 'gmso':
            self._parse_gmso(preprocessed_files)
        else:
            raise FoyerError('Backend not supported')

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
        return None

    def apply(self, top, references_file=None, use_residue_map=True,
              assert_bond_params=True, assert_angle_params=True,
              assert_dihedral_params=True, assert_improper_params=True,
              verbose=False, backend='gmso', *args, **kwargs):
        """Apply the force field to a molecular topology

        Parameters
        ----------
        top : gmso.Topology, or mb.Compound
            Molecular Topology to apply the force field to
        references_file : str, optional, defaut=None
            Name of file where force field refrences will be written
            (in Bibtex format).
        use_residue_map : bool, optional, default=True
            Options to speed up if there are a lot of repeated
            subtopology within the topology (assuming they all
            have the same name).
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system dihedrals.
        assert_improper_params : bool, optional, default=True
            If True, Foyer will exit if parametes are not found for all
            system impropers.
        verbose : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potential problematic details it encounters.
        backend : str, optional, default='gmso'
            Name of the backend used to store all the Types' information.
            At this point, 'gmso' is the only valid backend, but this set up
            allow future backend to be implemented easier.
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
            backend=backend,
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

        # TO DO in another PR
        if use_residue_map:
            # Detect duplicates subtopology/residues
            # (do matching by name, assert same number
            # of atoms)
            # Not implemented yet
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
                        backend='gmso',
                        debug=False,
                        **kwargs):
        """Parametrize the Topology from the typemap provided

        Assign AtomTypes and BondTypes to Atoms and Bonds, respectively.
        Creat Angles, Dihedrals, Impropers and assing corresponding
        AngleType, DihedralTypes, and ImproperTypes.

        Parameters
        ----------
        top : gmso.Topology
            gmso.Topology object that needed to be parametrized
        typemap : dict
            typemap generated by the atomtyper
        references_file : str, optional, defaut=None
            Name of file where force field refrences will be written
            (in Bibtex format).
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system dihedrals.
        assert_improper_params : bool, optional, default=True
            If True, Foyer will exit if parametes are not found for all
            system impropers.
        verbose : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potential problematic details it encounters.
        Returns
        -------
        """
        # Generate missing connections (angles, dihedrals, and impropers)
        top.identify_connections()

        if backend=='gmso':
            self._parametrize_gmsoFF(top=top,
                                typemap=typemap)
        else:
            raise FoyerError('Backend not supported')

        self._check_parameters(top, assert_bond_params,
                             assert_angle_params,
                             assert_dihedral_params,
                             assert_improper_params,
                             debug)
        return top

    def _parametrize_gmsoFF(self, top, typemap):
        """Parametrize a Topology with gmso.ForceField"""
        # Assign AtomTypes
        for atom in top.sites:
            atom.atom_type = self.ff.atom_types.get(
                                typemap[top.get_index(atom)]['atomtype'])
        if not all([a.atom_type for a in top.sites]):
            raise ValueError('Not all atoms in topology have atom types')

        # Assign BondTypes
        for bond in top.bonds:
            equiv_bmembers = bond.equivalent_members()
            btype_names = ['{}~{}'.format(
                                    bmem[0].atom_type.name,
                                    bmem[1].atom_type.name)
                                    for bmem in equiv_bmembers]
            for btype_name in btype_names:
                bond.bond_type = self.ff.bond_types.get(btype_name)
                if bond.bond_type:
                    # Grab the first match and then break
                    break

        # Assign AngleTypes
        for angle in top.angles:
            equiv_agmembers = angle.equivalent_members()
            agtype_names = ['{}~{}~{}'.format(
                                agmem[0].atom_type.name,
                                agmem[1].atom_type.name,
                                agmem[2].atom_type.name)
                                for agmem in equiv_agmembers]
            for agtype_name in agtype_names:
                angle.angle_type = self.ff.angle_types.get(agtype_name)
                if angle.angle_type:
                    # Grab the first match and then break
                    break

        # Assign DihedralTypes
        for dihedral in top.dihedrals:
            equiv_dmembers = dihedral.equivalent_members()
            dtype_names = ['{}~{}~{}~{}'.format(
                                dmem[0].atom_type.name,
                                dmem[1].atom_type.name,
                                dmem[2].atom_type.name,
                                dmem[3].atom_type.name)
                                for dmem in equiv_dmembers]
            for dtype_name in dtype_names:
                dihedral.dihedral_type = self.ff.dihedral_types.get(dtype_name)
                if dihedral.dihedral_type:
                    # Grab the first match and then break
                    break

        # Assing ImproperTypes
        for improper in top.impropers:
            equiv_imembers = improper.equivalent_members()
            itype_names = ['{}~{}~{}~{}'.format(
                                imem[0].atom_type.name,
                                imem[1].atom_type.name,
                                imem[2].atom_type.name,
                                imem[3].atom_type.name)
                           for imem in equiv_imembers]
            for itype_name in itype_names:
                improper.improper_type = self.ff.improper_types.get(itype_name)
                if improper.improper_type:
                    # Grab the first match and then break
                    break

        return top

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
