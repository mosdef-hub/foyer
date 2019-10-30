import pytest
import foyer


def test_basic_import():
    assert 'forcefields' in dir(foyer)


@pytest.mark.parametrize('ff_loader', ['load_OPLSAA', 'load_TRAPPE_UA'])
def test_forcefields_exist(ff_loader):
    assert ff_loader in dir(foyer.forcefields)


def test_load_forcefield():
    OPLSAA = foyer.forcefields.get_forcefield(name='oplsaa')
    TRAPPE_UA = foyer.forcefields.get_forcefield(name='trappe-ua')
    with pytest.raises(ValueError):
        foyer.forcefields.get_forcefield('bogus_name')
