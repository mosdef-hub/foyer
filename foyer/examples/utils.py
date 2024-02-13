"""Utility methods for foyer examples."""

import glob
from os.path import abspath, join, split


def example_file_path(filename):
    """Get the full path of the file name for a particular test file.

    Parameters
    ----------
    filename : str
        Name of the file to get

    Returns
    -------
    path: str
        Full path of the example file
    """
    return join(split(abspath(__file__))[0], "files", filename)


def glob_example_patterns(pattern):
    """Get the full paths for example files adhering to the glob pattern.

    Parameters
    ----------
    pattern : str
        the pattern for the files(expanded using globbing)

    Returns
    -------
    list of file absolute paths matching the pattern
    """
    return glob.glob(join(split(abspath(__file__))[0], "files", "pattern"))
