"""Foyer: Atomtyping and forcefield applying. """

from __future__ import print_function

from setuptools import find_packages, setup

#####################################
VERSION = "0.9.2"
ISRELEASED = True
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + ".dev0"
#####################################


setup(
    name="foyer",
    version=__version__,
    description=__doc__.split("\n")[0],
    long_description=__doc__,
    author="Janos Sallai, Christoph Klein",
    author_email="janos.sallai@vanderbilt.edu, christoph.klein@vanderbilt.edu",
    url="https://github.com/mosdef-hub/foyer",
    download_url="https://github.com/mosdef-hub/foyer/tarball/{}".format(
        __version__
    ),
    packages=find_packages(),
    package_data={
        "foyer": [
            "foyer/tests/*.txt",
            "foyer/tests/files/*.mol2",
            "foyer/tests/files/*.pdb",
            "foyer/tests/files/*.xml",
            "foyer/forcefields/*.xml",
            "opls_validation/*/*.top",
            "opls_validation/*/*.gro",
            "opls_validation/*/*.mol2",
            "opls_validation/oplsaa.ff/*",
            "examples/files/*",
        ]
    },
    entry_points={
        "foyer.forcefields": [
            "load_OPLSAA = foyer.forcefields.forcefields:load_OPLSAA",
            "load_TRAPPE_UA = foyer.forcefields.forcefields:load_TRAPPE_UA",
        ]
    },
    package_dir={"foyer": "foyer"},
    include_package_data=True,
    license="MIT",
    zip_safe=False,
    keywords="foyer",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS",
    ],
)
