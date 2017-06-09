import glob
import os
from pkg_resources import resource_filename

from lxml.etree import XMLSyntaxError, DocumentInvalid
import pytest

from foyer.tests.utils import glob_fn
from foyer.exceptions import (ValidationError, ValidationWarning,
                             MultipleValidationError)
from foyer.validator import Validator

XMLS = glob_fn('*.xml')
ERRORS = {'validationerror': (ValidationError, MultipleValidationError),
          'xmlsyntaxerror': XMLSyntaxError,
          'documentinvalid': DocumentInvalid,
}

FF_DIR = resource_filename('foyer', 'forcefields')
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, '*.xml'))


@pytest.mark.parametrize('ff_file', XMLS)
def test_xmls(ff_file):
    file_name = os.path.split(ff_file)[1]
    if 'error' in file_name:
        error_type = ERRORS[file_name.split('_')[0]]
        with pytest.raises(error_type):
            Validator(ff_file)
    elif file_name.startswith('warning'):
        with pytest.warns(ValidationWarning):
            Validator(ff_file)
    else:
        Validator(ff_file)


@pytest.mark.parametrize('ff_file', FORCEFIELDS)
def test_forcefields(ff_file):
    Validator(ff_file)



