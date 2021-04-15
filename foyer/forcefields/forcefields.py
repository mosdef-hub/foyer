import os
import glob
from pkg_resources import resource_filename
import warnings

from foyer import Forcefield
from foyer.validator import ValidationWarning


def get_ff_path():
    return [resource_filename('foyer', 'forcefields')]


def get_forcefield_paths(forcefield_name=None):
    for dir_path in get_ff_path():
        file_pattern = os.path.join(dir_path, 'xml/*.xml')
        file_paths = [file_path for file_path in glob.glob(file_pattern)]
    return file_paths


def get_forcefield(name=None):
    if name is None:
        raise ValueError('Need a force field name')
    file_paths = get_forcefield_paths()
    try:
        ff_path = next(val for val in file_paths if name in val)
    except StopIteration:
        raise ValueError(
            f"Could not find forcefield named {name} in path {get_ff_path()}"
            )
    return Forcefield(ff_path)


def load_OPLSAA():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ValidationWarning)
        return get_forcefield(name='oplsaa')


def load_TRAPPE_UA():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        return get_forcefield(name='trappe-ua')
