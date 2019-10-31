import pytest
import foyer


def test_basic_import():
    assert 'forcefields' in dir(foyer)


def test_loading_forcefields():
    #funcs = [func for func in dir(foyer.forcefields) if 'load' in func and '__' not in func]
    for func in dir(foyer.forcefields):
        if 'load_' in func and '__' not in func:
            eval('foyer.forcefields.' + func)()


def test_load_forcefield():
    OPLSAA = foyer.forcefields.get_forcefield(name='oplsaa')
    TRAPPE_UA = foyer.forcefields.get_forcefield(name='trappe-ua')
    with pytest.raises(ValueError):
        foyer.forcefields.get_forcefield('bogus_name')
