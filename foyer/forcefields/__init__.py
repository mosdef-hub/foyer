from pkg_resources import iter_entry_points


class ForceFields:
    pass

forcefields = ForceFields()

for entry_point in iter_entry_points(group='foyer.forcefields', name=None):
    setattr(forcefields, entry_point.name, entry_point.load())
