from io import StringIO
from lxml import etree, objectify
from lxml.etree import XMLSyntaxError, Element
from os.path import join, split, abspath

from foyer.smarts import SMARTS


class Validator(object):
    def __init__(self, optional_names=''):
        self.smarts = SMARTS(optional_names)

    def validate(self, ff_file_name):
        # parse XML
        ff_tree = etree.parse(ff_file_name)
        # TODO: figure out what exceptions are raised, and provide good error messages

        # validate tree against schema
        self.validate_xsd(ff_tree)
        # TODO: figure out what exceptions are raised, and provide good error messages

        # validate SMARTS strings
        self.validate_smarts(ff_tree)
        # TODO: figure out what exceptions are raised, and provide good error messages

    def validate_xsd(self, ff_tree, xsd_file=None):
        if xsd_file is None:
            xsd_file = join(split(abspath(__file__))[0], 'forcefields', 'ff.xsd')

        xmlschema_doc = etree.parse(xsd_file)
        xmlschema = etree.XMLSchema(xmlschema_doc)

        xmlschema.assertValid(ff_tree)

    def validate_smarts(self, ff_tree):
        results = ff_tree.xpath('/ForceField/AtomTypes/Type[@def]')
        for r in results:
            smarts_string = r.attrib['def']
            print(smarts_string)
            print(r.sourceline)
            # make sure smarts string can be parsed
            self.smarts.parse(smarts_string)

            # make sure referenced labels exist

    def validate_overrides(self, ff_tree):
        results = ff_tree.xpath('/ForceField/AtomTypes/Type[@name]')

        for r in results:
            smarts_string = r.attrib['def']
            print(smarts_string)
            print(r.sourceline)
            # make sure smarts string can be parsed
            self.smarts.parse(smarts_string)




if __name__ == '__main__':
    ff_file_name = join(split(abspath(__file__))[0], 'forcefields', 'oplsaa.xml')

    v = Validator()
    v.validate(ff_file_name)