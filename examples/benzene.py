import os

import parmed as pmd
from pkg_resources import resource_filename

from foyer.forcefield import apply_forcefield


resource_dir = resource_filename('foyer', '../opls_validation')
top_filename = os.path.join(resource_dir, 'benzene.top')
gro_filename = os.path.join(resource_dir, 'benzene.gro')
ff_filename = os.path.join(resource_dir, 'oplsaa.ff/forcefield.itp')

structure = pmd.load_file(top_filename, xyz=gro_filename)
parametrized = apply_forcefield(structure, forcefield=ff_filename, debug=False)

parametrized.save('benzene.gro', overwrite=True)
parametrized.save('benzene.top', overwrite=True)
