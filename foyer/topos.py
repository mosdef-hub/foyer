import os
import glob
import parmed as pmd
from pkg_resources import resource_filename


class Topos(object):

    def __init__(self):
        self.topo_dict = self.all_topos()

    def all_topos(self):

        resource_dir = resource_filename('foyer', '../opls_validation')
        top_files = set(glob.glob(os.path.join(resource_dir, '*.top')))

        topo_dict = {}

        for top_path in top_files:
            base_path, top_filename = os.path.split(top_path)
            basename = top_filename[:-4]
            gro_file = '{}-gas.gro'.format(basename)
            gro_path = os.path.join(base_path, gro_file)
            topo_dict[basename] = (top_path, gro_path)

        return topo_dict

    def load_topo(self, basename):

        top_path, gro_path = self.topo_dict[basename]
        structure = pmd.gromacs.GromacsTopologyFile(top_path, xyz=gro_path,
                                                    parametrize=False)
        structure.title = structure.title.replace(' GAS', '')
        return structure
