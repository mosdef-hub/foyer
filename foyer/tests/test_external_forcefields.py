import pytest


def test_basic_import():
    import foyer
    assert 'external_forcefields' in dir(foyer)
    import foyer.external_forcefields
