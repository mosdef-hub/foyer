import pytest
import foyer


def test_basic_import():
    assert 'forcefields' in dir(foyer)


def test_loading_forcefields():
    """Test that the forcefield loader functions run without error"""
    available_loaders = foyer.forcefield.get_available_forcefield_loaders()
    for loader in available_loaders:
        loader()


def test_load_forcefield():
    OPLSAA = foyer.forcefields.get_forcefield(name='oplsaa')
    TRAPPE_UA = foyer.forcefields.get_forcefield(name='trappe-ua')
    with pytest.raises(ValueError):
        foyer.forcefields.get_forcefield('bogus_name')
