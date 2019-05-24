import pytest

@pytest.fixture(scope="session")
def initdir(tmpdir):
    tmpdir.chdir()
