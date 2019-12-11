import os
import glob
from pkg_resources import resource_filename

from foyer import Forcefield


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
        raise ValueError('Could not find force field with name {}'
                ' in path {}'.format(name, get_ff_path()))
    return Forcefield(ff_path)


def load_OPLSAA():
    return get_forcefield(name='oplsaa')


def load_TRAPPE_UA():
    return get_forcefield(name='trappe-ua')
