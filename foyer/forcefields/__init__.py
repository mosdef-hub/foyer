"""Plugin support for forcefields."""

from importlib.metadata import entry_points


class ForceFields:
    """Container class for entrypoint-based forcefield files."""

    pass


forcefields = ForceFields()

for entry_point in entry_points(group="foyer.forcefields"):
    setattr(forcefields, entry_point.name, entry_point.load())
