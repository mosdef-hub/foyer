import itertools
import os
import warnings

from intermol.forces import *
from intermol.gromacs.gromacs_parser import GromacsParser, default_gromacs_include_dir
from oset import oset as OrderedSet

from foyer.atomtyper import find_atomtypes


def apply_forcefield(intermol_system, forcefield, debug=True):
    """Apply a forcefield to a Topology. """
    if forcefield.lower() in ['opls-aa', 'oplsaa', 'opls']:
        ff = Forcefield('oplsaa')
    elif forcefield.lower() in ['trappeua']:
        ff = Forcefield('trappeua')
    else:
        raise ValueError("Unsupported forcefield: '{0}'".format(forcefield))

    bondgraph = prepare_atoms(intermol_system)
    # Nodes are tuples of (atom, moleculetype).
    atoms = [atom for atom, _ in bondgraph.nodes()]
    find_atomtypes(atoms, forcefield, debug=debug)
    ff.resolve_bondingtypes(bondgraph)
    propogate_atomtyping(intermol_system)
    enumerate_forcefield_terms(intermol_system, bondgraph, ff)

    # Copy over defaults.
    intermol_system.nonbonded_function = ff.system.nonbonded_function
    intermol_system.combination_rule = ff.system.combination_rule
    intermol_system.genpairs = ff.system.genpairs
    intermol_system.lj_correction = ff.system.lj_correction
    intermol_system.coulomb_correction = ff.system.coulomb_correction

    #intermol_system.gen_pairs(n_excl=4)
    #ff = ff.prune()


def propogate_atomtyping(intermol_system):
    """Copy atomtype and bondingtype info to each atom instance in the system.

    intermol_system:
    """
    for moltype in intermol_system.molecule_types.values():
        first_mol = None
        for molecule in moltype.molecules:
            if first_mol is None:
                first_mol = molecule
                continue
            for typed_atom, untyped_atom in zip(first_mol.atoms, molecule.atoms):
                untyped_atom._atomtype = typed_atom._atomtype
                untyped_atom.bondingtype = typed_atom.bondingtype


def prepare_atoms(intermol_system):
    """Add neighbors and white- and blacklists to each atom.

    Note
    ----
    The use of ordered sets is not strictly necessary but it helps when
    debugging because it shows the order in which rules are added.

    Parameters
    ----------
    atoms : list of Atom objects
        The atoms whose atomtypes you are looking for. Atoms must have a
        property `neighbors` which is a list of other atoms that they are
        bonded to.

    """
    bondgraph = intermol_system.bondgraph_per_moleculetype
    for atom, mol_type in bondgraph.nodes():
        atom.neighbors = [atom for atom, _ in bondgraph.neighbors((atom, mol_type))]
        atom.whitelist = OrderedSet()
        atom.blacklist = OrderedSet()
    return bondgraph


def enumerate_forcefield_terms(intermol_system, bondgraph, forcefield, angles=True,
                               dihedrals=True, impropers=False):
    """Convert Bonds to ForcefieldBonds and find angles and dihedrals. """
    create_bonds(intermol_system, forcefield)

    if any([angles, dihedrals, impropers]):
        for node_1 in bondgraph.nodes_iter():
            neighbors_1 = bondgraph.neighbors(node_1)
            if len(neighbors_1) > 1:
                if angles:
                    create_angles(intermol_system, forcefield, node_1, neighbors_1)
                if dihedrals:
                    for node_2 in neighbors_1:
                        if node_2[0].index > node_1[0].index:
                            neighbors_2 = bondgraph.neighbors(node_2)
                            if len(neighbors_2) > 1:
                                create_dihedrals(intermol_system, forcefield, node_1, neighbors_1, node_2, neighbors_2)
                if impropers and len(neighbors_1) >= 3:
                    create_impropers(intermol_system, node_1, neighbors_1)


def create_bonds(intermol_system, forcefield):
    """Convert from tuples of (Atom1, Atom2) to ForcefieldBonds. """

    for mol_type in intermol_system.molecule_types.values():
        for molecule in mol_type.molecules:
            for bond in mol_type.bonds:
                atom1 = molecule.atoms[bond.atom1 - 1]
                atom2 = molecule.atoms[bond.atom2 - 1]
                bondingtypes = tuple([atom1.bondingtype, atom2.bondingtype])
                # TODO: Hide lookup logic.
                try:
                    bondtype = forcefield.bondtypes[bondingtypes]
                except KeyError:
                    try:
                        bondtype = forcefield.bondtypes[bondingtypes[::-1]]
                    except KeyError:
                        raise ValueError('No bondtype exists for bondingtypes {0}'.format(bondingtypes))
                bond.forcetype = bondtype
            break  # Only loop through one of the molecules.


def create_angles(intermol_system, forcefield, node, neighbors):
    """Find all angles around a node. """
    atom2 = node[0]
    mol_type = node[1]
    neighbor_atoms = [atom for atom, _ in neighbors]

    for pair in itertools.combinations(neighbor_atoms, 2):
        atom1 = pair[0]
        atom3 = pair[1]
        bondingtypes = tuple([atom1.bondingtype, atom2.bondingtype, atom3.bondingtype])
        # TODO: Hide lookup logic.
        try:
            angletype = forcefield.angletypes[bondingtypes]
        except KeyError:
            try:
                angletype = forcefield.angletypes[bondingtypes[::-1]]
            except KeyError:
                raise ValueError('No angletype exists for bondingtypes {0}'.format(bondingtypes))
        angle = Angle(atom1.index, atom2.index, atom3.index)
        angle.forcetype = angletype

        intermol_system.angletypes[bondingtypes] = angletype
        mol_type.angles.add(angle)


def create_dihedrals(intermol_system, forcefield, node_1, neighbors_1, node_2, neighbors_2):
    """Find all dihedrals around a pair of nodes. """
    # We need to make sure we don't remove the node from the neighbor lists
    # that we will be re-using in the following iterations.
    neighbors_1 = set(neighbors_1) - {node_2}
    neighbors_2.remove(node_1)

    atom2 = node_1[0]
    atom3 = node_2[0]
    mol_type = node_1[1]
    neighbor1_atoms = [atom for atom, _ in neighbors_1]
    neighbor2_atoms = [atom for atom, _ in neighbors_2]
    for pair in itertools.product(neighbor1_atoms, neighbor2_atoms):
        if pair[0] != pair[1]:
            atom1 = pair[0]
            atom4 = pair[1]
            # TODO: Hide lookup logic.
            try:
                bondingtypes = tuple([atom1.bondingtype, atom2.bondingtype, atom3.bondingtype, atom4.bondingtype, False])
                dihedraltype = forcefield.dihedraltypes[bondingtypes]
            except KeyError:
                try:
                    bondingtypes = tuple([atom4.bondingtype, atom3.bondingtype, atom2.bondingtype, atom1.bondingtype, False])
                    dihedraltype = forcefield.dihedraltypes[bondingtypes]
                except KeyError:
                    warnings.warn('No dihedraltype exists for bondingtypes {0}'.format(bondingtypes))
                    continue
            dihedral = Dihedral(atom1.index, atom2.index, atom3.index, atom4.index)
            if len(dihedraltype) == 1:
                dihedral.forcetype = next(iter(dihedraltype))
            else:
                raise NotImplementedError

            intermol_system.dihedraltypes[bondingtypes] = dihedraltype
            mol_type.dihedrals.add(dihedral)


def create_impropers(topology, node, neighbors):
    """Find all impropers around a node. """
    for triplet in itertools.combinations(neighbors, 3):
        atoms = sort_atoms_alphabetically([node, triplet[2]])
        if atoms[0] == node:
            topology.add_ff_improper(ForcefieldImproper(
                node, triplet[0], triplet[1], triplet[2]))
        elif atoms[1] == node:
            topology.add_ff_improper(ForcefieldImproper(
                triplet[2], triplet[1], triplet[0], node))


class Forcefield(GromacsParser):
    """A container class for the OPLS forcefield."""

    def __init__(self, forcefield_string):
        """Populate the database using files bundled with GROMACS."""
        ff_file = os.path.join(default_gromacs_include_dir(),
                               '{0}.ff/forcefield.itp'.format(forcefield_string))
        super(Forcefield, self).__init__(ff_file, None)
        self.read()

    def read(self):
        """ """
        self.current_directive = None
        self.if_stack = list()
        self.else_stack = list()

        self.atomtypes = self.system.atomtypes
        self.bondtypes = self.system.bondtypes
        self.angletypes = self.system.angletypes
        self.dihedraltypes = self.system.dihedraltypes
        self.implicittypes = dict()
        self.pairtypes = dict()
        self.cmaptypes = dict()
        self.nonbondedtypes = dict()

        # Parse the itp file into a set of plain text, intermediate
        # TopMoleculeType objects.
        self.process_file(self.top_filename)

    def resolve_bondingtypes(self, bondgraph):
        """ """
        for n, (atom, _) in enumerate(bondgraph.nodes()):
            if atom.atomtype[0] not in self.atomtypes:
                print("Could not find atomtype: '{0}' in forcefield.".format(atom.atomtype[0]))
            else:
                atom.bondingtype = self.atomtypes[atom.atomtype[0]].bondtype
                atom.cgnr = n
                if not atom.charge:
                    atom.charge = (0, self.atomtypes[atom.atomtype[0]].charge)
                if not atom.mass:
                    atom.mass = (0, self.atomtypes[atom.atomtype[0]].mass)

    # def find_atom_types(self, bondtype):
    #     """If the id is the atom type, return the AtomType object. """
    #     matching_atom_types = []
    #     bondtype = str(bondtype)
    #     for kind, atomtype in self.atomtypes.items():
    #         if bondtype.endswith('*'):
    #             # id is a wildcard ending in *
    #             prefix = bondtype[:-1]
    #
    #             if atomtype.bondtype.startswith(prefix):
    #                 matching_atom_types.append(kind)
    #             elif kind.startswith(prefix):
    #                 matching_atom_types.append(kind)
    #         else:
    #             # id is not a wildcard
    #             if bondtype == atomtype.bondtype:
    #                 matching_atom_types.append(kind)
    #     return matching_atom_types
    #
    # def prune(self):
    #     """Create force field with only information relevant to topology. """
    #
    #     bonds_to_remove = set()
    #     for bond in self.top.ff_bonds:
    #         if bond.kind not in self.bondtypes:
    #             print("Could not find bondtype: '{0}' in forcefield.".format(bond.kind))
    #             bonds_to_remove.add(bond)
    #     self.top._ff_bonds = self.top._ff_bonds - bonds_to_remove
    #
    #     angles_to_remove = set()
    #     for angle in self.top.ff_angles:
    #         if angle.kind not in self.angletypes:
    #             print("Could not find angletype: '{0}' in forcefield.".format(angle.kind))
    #             angles_to_remove.add(angle)
    #     self.top._ff_angles = self.top._ff_angles - angles_to_remove
    #
    #     dihedrals_to_remove = set()
    #     for dihedral in self.top.ff_dihedrals:
    #         if dihedral.kind not in self.dihedraltypes:
    #             print("Could not find dihedraltype: '{0}' in forcefield.".format(dihedral.kind))
    #             dihedrals_to_remove.add(dihedral)
    #     self.top._ff_dihedrals = self.top._ff_dihedrals - dihedrals_to_remove
    #
    #     impropers_to_remove = set()
    #     for improper in self.top.ff_impropers:
    #         if improper.kind == ('O_2', 'C_2', 'OS', 'CT'):
    #             print("Keeping explicitcly defined improper: '{0}'".format(improper.kind))
    #         elif improper.kind not in self.impropertypes:
    #             print("Could not find impropertype: '{0}' in forcefield.".format(improper.kind))
    #             impropers_to_remove.add(improper)
    #     self.top._ff_impropers = self.top._ff_impropers - impropers_to_remove
    #
    #
    #     # ff = Forcefield(self.topology)
    # ff.atomtypes.update(self.atomtypes)
    #
    # retained_types = [atom.atomtype for atom in self.topology.atoms]
    #
    # # Prune the atom types
    # for atomtype in list(ff.atomtypes.keys()):
    #     if atomtype not in retained_types:
    #         del ff.atomtypes[atomtype]
    #
    #
    # # Prune the bond types, resolving wildcards.
    # for (bondtype1, bondtype2), bond_type in self.bondtypes.items():
    #     atom_types1 = ff.find_atom_types(bondtype1)
    #     atom_types2 = ff.find_atom_types(bondtype2)
    #
    #     # For every combination of matching atom kinds, create a bond type.
    #     for (atom_type1, atom_type2) in itertools.product(atom_types1, atom_types2):
    #         pair = (atom_type1, atom_type2)
    #         ff.bondtypes[pair] = bond_type
    #
    # # Prune the angle types, resolving wildcards.
    # for (bondtype1, bondtype2, bondtype3), angle_type in self.angletypes.items():
    #     atom_types1 = ff.find_atom_types(bondtype1)
    #     atom_types2 = ff.find_atom_types(bondtype2)
    #     atom_types3 = ff.find_atom_types(bondtype3)
    #
    #     # For every combination of matching atom kinds, create an angle type.
    #     for (atom_type1, atom_type2, atom_type3) in itertools.product(atom_types1, atom_types2, atom_types3):
    #         triplet = (atom_type1, atom_type2, atom_type3)
    #         ff.angletypes[triplet] = angle_type
    # return ff
