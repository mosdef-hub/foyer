from os.path import join, split, abspath
from warnings import warn

from lxml import etree, objectify
from lxml.etree import XMLSyntaxError, Element, DocumentInvalid
import networkx as nx
from plyplus.common import ParseError

from foyer.exceptions import ValidationError, ValidationWarning
from foyer.smarts_graph import SMARTSGraph


class Validator(object):
    def __init__(self, ff_file_name):
        ff_tree = etree.parse(ff_file_name)
        self.validate_xsd(ff_tree)

        # Loading forcefield should succeed, because XML can be parsed and
        # basics have been validated.
        from foyer.forcefield import Forcefield
        self.smarts_parser = Forcefield(ff_file_name, validation=False).parser

        self.validate_smarts(ff_tree)
        # TODO: figure out what exceptions are raised, and provide good error messages

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
            raise

    def validate_smarts(self, ff_tree):
        results = ff_tree.xpath('/ForceField/AtomTypes/Type')
        atom_types = ff_tree.xpath('/ForceField/AtomTypes/Type/@name')

        missing_smarts = []
        for r in results:
            smarts_string = r.attrib.get('def')
            name = r.attrib['name']
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

                raise ValidationError("Malformed SMARTS string{} on line {}".format(
                    column, r.sourceline), ex, r.sourceline)

            # make sure referenced labels exist
            smarts_graph = SMARTSGraph(smarts_string, parser=self.smarts_parser,
                                       name=name, overrides=r.attrib.get('overrides'))
            for atom_expr in nx.get_node_attributes(smarts_graph, 'atom').values():
                labels = atom_expr.select('has_label')
                for label in labels:
                    atom_type = label.tail[0][1:]
                    if atom_type not in atom_types:
                        raise ValidationError(
                            "Reference to undefined atomtype {} in SMARTS string" " '{}' at line {}".format(
                                atom_type, r.attrib['def'], r.sourceline), None, r.sourceline)

        warn("The following atom types do not have smarts definitions: {}".format(', '.join(missing_smarts)),
             ValidationWarning)

    def validate_overrides(self, ff_tree):
        results = ff_tree.xpath('/ForceField/AtomTypes/Type[@name]')

        for r in results:
            smarts_string = r.attrib['def']
            print(smarts_string)
            print(r.sourceline)
            # make sure smarts string can be parsed
            self.smarts.parse(smarts_string)


if __name__ == '__main__':
    # ff_file_name = join(split(abspath(__file__))[0], 'forcefields', 'oplsaa.xml')

    from foyer.tests.utils import get_fn

    v = Validator(get_fn('bad_ff_missingsmarts.xml'))
