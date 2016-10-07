import os
import glob
from functools import lru_cache
from pkg_resources import resource_filename
import parmed as pmd

class Topology(pmd.Structure):

    def __init__(self):
        super(Topology, self).__init__()

    @classmethod
    @lru_cache(maxsize=32)
    def available(cls):

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

    @classmethod
    def by_name(cls, basename):

        top_path, gro_path = cls.available()[basename]
        structure = pmd.gromacs.GromacsTopologyFile(top_path, xyz=gro_path,
                                                    parametrize=False)
        structure.title = structure.title.replace(' GAS', '')
        return structure
