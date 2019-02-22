from __future__ import division, print_function

import glob
import sys
from os.path import join, split, abspath
import os
import importlib
from unittest import SkipTest
import inspect
import textwrap

import numpy as np


class DelayImportError(ImportError, SkipTest):
    pass


MESSAGES = dict()
MESSAGES['mbuild'] = '''

The code at {filename}:{line_number} requires the "mbuild" package

mbuild can be installed using:

# conda install -c mosdef mbuild

or

# pip install mbuild

'''


def import_(module):
    """Import a module, and issue a nice message to stderr if the module isn't installed.
    Parameters
    ----------
    module : str
        The module you'd like to import, as a string
    Returns
    -------
    module : {module, object}
        The module object
    Examples
    --------
    >>> # the following two lines are equivalent. the difference is that the
    >>> # second will check for an ImportError and print you a very nice
    >>> # user-facing message about what's wrong (where you can install the
    >>> # module from, etc) if the import fails
    >>> import tables
    >>> tables = import_('tables')
    """
    try:
        return importlib.import_module(module)
    except ImportError as e:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = 'The code at {filename}:{line_number} requires the ' + module + ' package'
            e = ImportError('No module named %s' % module)

        frame, filename, line_number, function_name, lines, index = \
            inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(filename=os.path.basename(filename), line_number=line_number)
        m = textwrap.dedent(m)

        bar = '\033[91m' + '#' * max(len(line) for line in m.split(os.linesep)) + '\033[0m'

        print('', file=sys.stderr)
        print(bar, file=sys.stderr)
        print(m, file=sys.stderr)
        print(bar, file=sys.stderr)
        raise DelayImportError(m)


try:
    import mbuild
    has_mbuild = True
    del mbuild
except ImportError:
    has_mbuild = False
