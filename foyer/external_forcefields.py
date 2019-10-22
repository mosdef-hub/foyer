from pkg_resources import iter_entry_points


class ExternalForceFields:
    pass

external_forcefields = ExternalForceFields()

for entry_point in iter_entry_points(group='foyer.external_forcefields', name=None):
    setattr(external_forcefields, entry_point.name, entry_point.load())
