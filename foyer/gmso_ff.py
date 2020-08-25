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
import foyer.element as custom_elem

import gmso
from gmso import ForceField

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError
from foyer import smarts
from foyer.validator import Validator
from foyer.xml_writer import write_foyer
from foyer.utils.io import import_, has_mbuild
from foyer.utils.external import get_ref

# Copy from original forcefield.py
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

class Forcefield(gmso.ForceField):
    """General Forcefield object that can be created by either GMSO Forcefield or OpenMM Forcefield

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load
    name : str, optional, None
        Name of a forcefield to load that is packaged within foyer
    backend : str, optional, default='openmm'
        Name of the backend used to store all the Types' information.
        Can choose between 'openmm' and 'gmso'

    """
    def __init__(self, forcefield_files=None,
                 name=None, validation=True, debug=False):
        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self.atomTypeRefs = dict()
        self.atomTypeClasses = dict()
        self.atomTypeElements = dict()
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
        try:
            super(Forcefield, self).__init__(*preprocessed_files)
        finally:
            for ff_file_name in preprocessed_files:
                os.remove(ff_file_name)

        self.parser = smarts.SMARTS(self.non_element_types)
        # Need to find a way to update the foyer specific dicts

    @classmethod
    def from_xml(cls, xml_locs):
        """Overwrite the original gmso.Forcefield.from_xml

        This class method creates a Forcefield object from a foyer XML files,
        overwritten the classmethod in GMSO Forcefield which createa the object
        from the GMSO XML file format.

        Paramaters
        ----------
        xml_locs : str or iterable of str
            string or iterable of strings containing the forcefield XML locations

        Returns
        -------
        forcefield : gmso.Forcefield
            A gmso.Forcefield objection existed inside the foyer.Forcefield object
        """
        return None
        

