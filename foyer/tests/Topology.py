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
            mol2_file = '{}.mol2'.format(basename)
            mol2_path = os.path.join(base_path, mol2_file)
            topo_dict[basename] = (mol2_path, top_path, gro_path)

        return topo_dict

    @classmethod
    def by_name(cls, basename, parameterized=False):

        mol2_path, top_path, gro_path = cls.available()[basename]
        if parameterized:

            structure = pmd.gromacs.GromacsTopologyFile(top_path, xyz=gro_path,
                                                    parametrize=False)
            structure.title = structure.title.replace(' GAS', '')
        else:
            # structure = pmd.gromacs.GromacsTopologyFile(pdb_path, xyz=pdb_path,
            #                                         parametrize=False)
            structure = pmd.load_file(mol2_path, structure=True)
            structure.title = structure.title.replace(' GAS', '')

        return structure

    def strip_parameterization(self):
        for atom in self.atoms:
            atom.type = None
            atom.charge = None
            atom.atom_type = None
        # keep bonds
        self.bond_types.clear()
        self.angles.clear()
        self.angle_types.clear()
        self.dihedrals.clear()
        self.dihedral_types.clear()
        self.impropers.clear()
        self.improper_types.clear()
