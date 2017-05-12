from collections import Counter
from os.path import join, split, abspath
from warnings import warn

from lxml import etree
from lxml.etree import DocumentInvalid
import networkx as nx
from plyplus.common import ParseError

from foyer.exceptions import (ValidationError, ValidationWarning,
                              raise_collected)
from foyer.smarts_graph import SMARTSGraph


class Validator(object):
    def __init__(self, ff_file_name):
        from foyer.forcefield import preprocess_forcefield_files
        preprocessed_ff_file_name = preprocess_forcefield_files([ff_file_name])

        ff_tree = etree.parse(preprocessed_ff_file_name[0])
        self.validate_xsd(ff_tree)

        self.atom_type_names = ff_tree.xpath('/ForceField/AtomTypes/Type/@name')
        self.atom_types = ff_tree.xpath('/ForceField/AtomTypes/Type')

        self.validate_class_type_exclusivity(ff_tree)

        # Loading forcefield should succeed, because XML can be parsed and
        # basics have been validated.
        from foyer.forcefield import Forcefield
        self.smarts_parser = Forcefield(preprocessed_ff_file_name, validation=False).parser

        self.validate_smarts()
        self.validate_overrides()

    @staticmethod
    def validate_xsd(ff_tree, xsd_file=None):
        if xsd_file is None:
            xsd_file = join(split(abspath(__file__))[0], 'forcefields', 'ff.xsd')

        xmlschema_doc = etree.parse(xsd_file)
        xmlschema = etree.XMLSchema(xmlschema_doc)

        error_texts = {'missing_atom_type_in_nonbonded':
                           ("Atom type {} is found in NonbondedForce at line {}"
                            " but undefined in AtomTypes"),
                       'nonunique_atomtype_name':
                           "Atom type {} is defined a second time at line {}",
                       'atomtype_name_key':
                           "Atom type {} is defined a second time at line {}"
        }

        def create_error(keyword, message, line):
            atomtype = message[message.find("[") + 1:message.find("]")]
            error_text = error_texts[keyword].format(atomtype, line)
            return ValidationError(error_text, ex, line)

        try:
            xmlschema.assertValid(ff_tree)
        except DocumentInvalid as ex:
            message = ex.error_log.last_error.message
            line = ex.error_log.last_error.line
            # rewrite error message for constraint violation
            if ex.error_log.last_error.type_name == "SCHEMAV_CVC_IDC":
                for keyword in error_texts:
                    if keyword in message:
                        raise create_error(keyword, message, line)
                else:
                    raise ValidationError('Unhandled XML validation error. '
                                          'Please consider submitting a bug report.', ex, line)
            raise

    def validate_class_type_exclusivity(self, ff_tree):
        sections = {'HarmonicBondForce/Bond': 2,
                    'HarmonicAngleForce/Angle': 3,
                    'RBTorsionForce/Proper': 4}

        errors = []
        for element, num_atoms in sections.items():
            valid_attribs = set()
            for n in range(1, num_atoms + 1):
                valid_attribs.add('class{}'.format(n))
                valid_attribs.add('type{}'.format(n))

            for entry in ff_tree.xpath('/ForceField/{}'.format(element)):
                attribs = [valid for valid in valid_attribs
                           if entry.attrib.get(valid) is not None]
                if num_atoms != len(attribs):
                    error = ValidationError(
                        'Invalid number of "class" and/or "type" attributes for'
                        ' {} at line {}'.format(element, entry.sourceline),
                        None, entry.sourceline)
                    errors.append(error)
                number_endings = Counter([a[-1] for a in attribs])
                if not all(1 == x for x in number_endings.values()):
                    error = ValidationError(
                        'Only one "class" or "type" attribute may be defined'
                        ' for each atom in a bonded force. See line {}'.format(entry.sourceline),
                        None, entry.sourceline)
                    errors.append(error)

                referenced_types = []
                for valid in valid_attribs:
                    if valid.startswith('type'):
                        atomtype = entry.attrib.get(valid)
                        if atomtype:
                            referenced_types.append(atomtype)

                for atomtype in referenced_types:
                    if atomtype not in self.atom_type_names:
                        error = ValidationError(
                            'Atom type {} is found in {} at line {} but'
                            ' undefined in AtomTypes'.format(atomtype, element.split('/')[0], entry.sourceline),
                            None, entry.sourceline)
                        errors.append(error)
        raise_collected(errors)

    def validate_smarts(self):
        missing_smarts = []
        errors = []
        for entry in self.atom_types:
            smarts_string = entry.attrib.get('def')
            name = entry.attrib['name']
            if smarts_string is None:
                missing_smarts.append(name)
                continue
            # make sure smarts string can be parsed
            try:
                self.smarts_parser.parse(smarts_string)
            except ParseError as ex:
                if " col " in ex.args[0]:
                    column = ex.args[0][ex.args[0].find(" col ") + 5:].strip()
                    column = " at character {} of {}".format(column, smarts_string)
                else:
                    column = ""

                malformed = ValidationError(
                    "Malformed SMARTS string{} on line {}".format(column, entry.sourceline),
                    ex, entry.sourceline)
                errors.append(malformed)
                continue

            # make sure referenced labels exist
            smarts_graph = SMARTSGraph(smarts_string, parser=self.smarts_parser,
                                       name=name, overrides=entry.attrib.get('overrides'))
            for atom_expr in nx.get_node_attributes(smarts_graph, 'atom').values():
                labels = atom_expr.select('has_label')
                for label in labels:
                    atom_type = label.tail[0][1:]
                    if atom_type not in self.atom_type_names:
                        undefined = ValidationError(
                            "Reference to undefined atomtype '{}' in SMARTS "
                            "string '{}' at line {}".format(atom_type, entry.attrib['def'], entry.sourceline),
                            None, entry.sourceline)
                        errors.append(undefined)
        raise_collected(errors)
        if missing_smarts:
            warn("The following atom types do not have smarts definitions: {}".format(
                ', '.join(missing_smarts)), ValidationWarning)

    def validate_overrides(self):
        errors = []
        for entry in self.atom_types:
            overrides = entry.attrib.get('overrides')
            if not overrides:
                continue
            overridden_types = [at.strip() for at in overrides.split(',') if at]
            for atom_type in overridden_types:
                if atom_type not in self.atom_type_names:
                    undefined = ValidationError(
                        "Reference to undefined atomtype '{}' in 'overrides' "
                        "'{}' at line {}".format(atom_type, entry.attrib['overrides'], entry.sourceline),
                        None, entry.sourceline)
                    errors.append(undefined)
        raise_collected(errors)


if __name__ == '__main__':
    from foyer.tests.utils import get_fn
    v = Validator(get_fn('validate_types.xml'))
    v = Validator(get_fn('validationerror_validate_types.xml'))
