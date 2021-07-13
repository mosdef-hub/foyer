import glob
import os

import pytest
from lxml.etree import DocumentInvalid, XMLSyntaxError
from pkg_resources import resource_filename

from foyer.exceptions import (
    MultipleValidationError,
    ValidationError,
    ValidationWarning,
)
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import glob_fn
from foyer.validator import Validator

XMLS = glob_fn("*.xml")
ERRORS = {
    "validationerror": (ValidationError, MultipleValidationError),
    "xmlsyntaxerror": XMLSyntaxError,
    "documentinvalid": DocumentInvalid,
}

FF_DIR = resource_filename("foyer", "forcefields")
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, "xml/*.xml"))


class TestValidator(BaseTest):
    @pytest.mark.parametrize("ff_file", XMLS)
    def test_xmls(self, ff_file):
        file_name = os.path.split(ff_file)[1]
        if "error" in file_name:
            error_type = ERRORS[file_name.split("_")[0]]
            with pytest.raises(error_type):
                Validator(ff_file)
        elif file_name.startswith("warning"):
            with pytest.warns(ValidationWarning):
                Validator(ff_file)
        else:
            Validator(ff_file)

    @pytest.mark.parametrize("ff_file", FORCEFIELDS)
    def test_forcefields(self, ff_file):
        Validator(ff_file)
