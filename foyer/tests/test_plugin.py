import pytest
import foyer


def test_basic_import():
    assert 'forcefields' in dir(foyer)


@pytest.mark.parametrize('ff_name', ['OPLSAA', 'TRAPPE_UA'])
def test_forcefields_exist(ff_name):
    ff_name in dir(foyer.forcefields)


def test_load_forcefield():
    OPLSAA = foyer.forcefields.get_forcefield(name='oplsaa')
    TRAPPE_UA = foyer.forcefields.get_forcefield(name='trappe-ua')
    with pytest.raises(ValueError):
        foyer.forcefields.get_forcefield('bogus_name')
