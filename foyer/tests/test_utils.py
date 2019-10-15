import numpy as np
import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.tests.utils import get_fn
from foyer.utils.nbfixes import apply_nbfix


def test_apply_nbfix():
    opls = Forcefield(name='oplsaa')
    ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    ethane = opls.apply(ethane)
    ethane_tweaked = apply_nbfix(
        struct=ethane,
        atom_type1='opls_135',
        atom_type2='opls_140',
        sigma=0.4,
        epsilon=50.0,
    )

    assert not ethane.has_NBFIX()
    assert ethane_tweaked.has_NBFIX()

    # 0.44898.... is rmin, which parmed uses internally in place of sigma
    for atom in ethane_tweaked:
        if atom.atom_type.name == 'opls_135':
            assert np.allclose(
                atom.atom_type.nbfix['opls_140'][:2],
                [0.44898481932374923, 50.0]
            )
        elif atom.atom_type.name == 'opls_140':
            assert np.allclose(
                atom.atom_type.nbfix['opls_135'][:2],
                [0.44898481932374923, 50.0]
            )

def test_apply_nbfix_bad_atom_type():
    opls = Forcefield(name='oplsaa')
    ethane = pmd.load_file(get_fn('ethane.mol2'), structure=True)
    ethane = opls.apply(ethane)
    with pytest.raises(ValueError):
        apply_nbfix(
            struct=ethane,
            atom_type1='opls__typo_135',
            atom_type2='opls_140',
            sigma=0.4,
            epsilon=50.0,
        )
    with pytest.raises(ValueError):
        apply_nbfix(
            struct=ethane,
            atom_type1='opls_135',
            atom_type2='opls_141',
            sigma=0.4,
            epsilon=50.0,
        )
