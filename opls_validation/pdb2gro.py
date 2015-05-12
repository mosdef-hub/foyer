from glob import glob
import os

pdb_files = glob('*.pdb')
for pdb_file in pdb_files:
    name, ext = os.path.splitext(pdb_file)
    gro_file = '{0}.gro'.format(name)
    os.system('gmx editconf -f {0} -o {1}'.format(pdb_file, name))

