import glob
import os
from pkg_resources import resource_filename

import mbuild as mb
import parmed as pmd
import parmed.unit as u
import pytest

from foyer import Forcefield
from foyer.tests.utils import get_fn


FF_DIR = resource_filename('foyer', 'forcefields')
FORCEFIELDS = glob.glob(os.path.join(FF_DIR, '*.xml'))


def test_load_files():
    for ff_file in FORCEFIELDS:
        ff1 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) > 0

        ff2 = Forcefield(forcefield_files=ff_file)
        assert len(ff1._atomTypes) == len(ff2._atomTypes)

        ff3 = Forcefield(name=ff_file.split('/')[-1].strip('.xml'))
        assert len(ff1._atomTypes) == len(ff3._atomTypes)


def test_duplicate_type_definitions():
    with pytest.raises(ValueError):
        ff4 = Forcefield(name='oplsaa', forcefield_files=FORCEFIELDS)


def test_from_parmed():
    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    oplsaa = Forcefield(FORCEFIELDS[0])
    ethane = oplsaa.apply(mol2)

    assert sum((1 for at in ethane.atoms if at.type == 'opls_135')) == 2
    assert sum((1 for at in ethane.atoms if at.type == 'opls_140')) == 6
    assert len(ethane.bonds) == 7
    assert all(x.type for x in ethane.bonds)
    assert len(ethane.angles) == 12
    assert all(x.type for x in ethane.angles)
    assert len(ethane.rb_torsions) == 9
    assert all(x.type for x in ethane.dihedrals)

    mol2 = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    mol2.box_vectors = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    oplsaa = Forcefield(FORCEFIELDS[0])
    ethane = oplsaa.apply(mol2)

    assert ethane.box_vectors == mol2.box_vectors


def test_from_mbuild():
    mol2 = mb.load(get_fn('ethane.mol2'))
    oplsaa = Forcefield(FORCEFIELDS[0])
    ethane = oplsaa.apply(mol2)

    assert sum((1 for at in ethane.atoms if at.type == 'opls_135')) == 2
    assert sum((1 for at in ethane.atoms if at.type == 'opls_140')) == 6
    assert len(ethane.bonds) == 7
    assert all(x.type for x in ethane.bonds)
    assert len(ethane.angles) == 12
    assert all(x.type for x in ethane.angles)
    assert len(ethane.rb_torsions) == 9
    assert all(x.type for x in ethane.dihedrals)

    boundingbox = mol2.boundingbox
    assert ethane.box_vectors[0][0].value_in_unit(u.nanometers) == boundingbox.lengths[0]
    assert ethane.box_vectors[1][1].value_in_unit(u.nanometers) == boundingbox.lengths[1]
    assert ethane.box_vectors[2][2].value_in_unit(u.nanometers) == boundingbox.lengths[2]
