"""Plugin support for forcefields."""
from pkg_resources import iter_entry_points


class ForceFields:
    """Container class for entrypoint-based forcefield files."""

    pass


forcefields = ForceFields()

for entry_point in iter_entry_points(group="foyer.forcefields", name=None):
    setattr(forcefields, entry_point.name, entry_point.load())
