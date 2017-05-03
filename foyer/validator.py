from os.path import join, split, abspath
from warnings import warn

from lxml import etree
from lxml.etree import DocumentInvalid
import networkx as nx
from plyplus.common import ParseError

from foyer.exceptions import ValidationError, ValidationWarning, MultipleValidationError
from foyer.smarts_graph import SMARTSGraph


class Validator(object):
    def __init__(self, ff_file_name):
        from foyer.forcefield import preprocess_forcefield_files
        preprocessed_ff_file_name = preprocess_forcefield_files([ff_file_name])

        ff_tree = etree.parse(preprocessed_ff_file_name[0])
        self.validate_xsd(ff_tree)

        # Loading forcefield should succeed, because XML can be parsed and
        # basics have been validated.
        from foyer.forcefield import Forcefield
        self.smarts_parser = Forcefield(preprocessed_ff_file_name, validation=False).parser

        self.atom_type_names = ff_tree.xpath('/ForceField/AtomTypes/Type/@name')
        self.atom_types = ff_tree.xpath('/ForceField/AtomTypes/Type')
        self.validate_smarts()
        self.validate_overrides()

    def validate_xsd(self, ff_tree, xsd_file=None):
        if xsd_file is None:
            xsd_file = join(split(abspath(__file__))[0], 'forcefields', 'ff.xsd')

        xmlschema_doc = etree.parse(xsd_file)
        xmlschema = etree.XMLSchema(xmlschema_doc)

        try:
            xmlschema.assertValid(ff_tree)
        except DocumentInvalid as ex:
            message = ex.error_log.last_error.message
            line = ex.error_log.last_error.line
            # rewrite error message for constraint violation
            if ex.error_log.last_error.type_name == "SCHEMAV_CVC_IDC":
                if "missing_atom_type_in_nonbonded" in message:
                    atomtype = message[message.find("[") + 1:message.find("]")]
                    raise ValidationError(
                        "Atom type {} is found in NonbondedForce at line {} but"
                        " undefined in AtomTypes".format(atomtype, line), ex, line)
                elif "nonunique_atomtype_name" in message:
                    atomtype = message[message.find("[") + 1:message.find("]")]
                    raise ValidationError(
                        "Atom type {} is defined a second time at line {}".format(
                            atomtype, line), ex, line)
                elif "atomtype_name_key" in message:
                    atomtype = message[message.find("[") + 1:message.find("]")]
                    raise ValidationError(
                        "Atom type {} is defined a second time at line {}".format(
                            atomtype, line), ex, line)
                else:
                    raise ValidationError('Unhandled XML validation error. '
                                          'Please consider submitting a bug report.', ex, line)
            raise

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

                malformed = ValidationError("Malformed SMARTS string{} on line {}".format(
                    column, entry.sourceline), ex, entry.sourceline)
                errors.append(malformed)
                print(malformed)
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
                            "Reference to undefined atomtype '{}' in SMARTS string '{}' at line {}".format(
                                atom_type, entry.attrib['def'], entry.sourceline), None, entry.sourceline)
                        errors.append(undefined)
        if len(errors) > 1:
            raise MultipleValidationError(errors)
        elif len(errors) == 1:
            raise errors[0]
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
                        "Reference to undefined atomtype '{}' in 'overrides' '{}' at line {}".format(
                            atom_type, entry.attrib['overrides'], entry.sourceline), None, entry.sourceline)
                    errors.append(undefined)
        if len(errors) > 1:
            raise MultipleValidationError(errors)
        elif len(errors) == 1:
            raise errors[0]
