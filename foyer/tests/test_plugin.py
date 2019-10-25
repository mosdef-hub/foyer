import pytest


def test_basic_import():
    import foyer
    assert 'forcefields' in dir(foyer)
    import foyer.forcefields.forcefields
