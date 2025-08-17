"""Support user-created forcefield XML files."""

import glob
import importlib.resources as resources
import os

from foyer import Forcefield


def get_ff_path():
    """Return path to forcefield locations."""
    ff_dir = resources.files("foyer").joinpath("forcefields")
    return [ff_dir]


def get_forcefield_paths(forcefield_name=None):
    """Return file paths to forcefield XMLs."""
    for dir_path in get_ff_path():
        file_pattern = os.path.join(dir_path, "xml/*.xml")
        file_paths = [file_path for file_path in glob.glob(file_pattern)]
    return file_paths


def get_forcefield(name=None):
    """Find forcefield based on name."""
    if name is None:
        raise ValueError("Need a force field name")
    file_paths = get_forcefield_paths()
    try:
        ff_path = next(val for val in file_paths if name in val)
    except StopIteration:
        raise ValueError(
            f"Could not find forcefield named {name} in path {get_ff_path()}"
        )
    return Forcefield(ff_path)


def load_OPLSAA():
    """Load internal forcefield XML for the packaged OPLS-AA forcefield."""
    return get_forcefield(name="oplsaa")


def load_TRAPPE_UA():
    """Load internal forcefield XML for the packaged TRAPPE-UA forcefield."""
    return get_forcefield(name="trappe-ua")
