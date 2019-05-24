import pytest

@pytest.fixture(autouse=True)
def initdir(tmpdir):
    tmpdir.chdir()
