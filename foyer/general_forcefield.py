"""Foyer general Forcefield class and utility method."""
import collections
import glob
import os
import warnings
from copy import deepcopy
from tempfile import NamedTemporaryFile

import gmso
import mbuild as mb
from gmso.core import element
from gmso.external import from_mbuild
from pkg_resources import resource_filename

from foyer import smarts
from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError, ValidationError
from foyer.utils.external import get_ref
from foyer.utils.io import import_
from foyer.validator import Validator


# Copy from original forcefield.py
def preprocess_forcefield_files(forcefield_files=None, backend="gmso"):
    """Pre-process foyer Forcefield XML files."""
    if forcefield_files is None:
        return None

    tmp_processed_files = list()
    if backend == "gmso":
        # Run through the forcefield XML conversion
        from gmso.external.convert_foyer_xml import from_foyer_xml

        for idx, file in enumerate(forcefield_files):
            _, suffix = os.path.split(file)
            tempfile = NamedTemporaryFile(suffix=suffix, delete=False)
            try:
                Validator(ff_file_name=file, debug=False)
                from_foyer_xml(
                    foyer_xml=str(file),
                    gmso_xml=str(tempfile.name),
                    overwrite=True,
                )
            except:
                warnings.warn(
                    f"Could not convert {str(file)}, attempt to read in as is."
                )
                os.popen(f"cp {file} {tempfile.name}")
            tmp_processed_files.append(tempfile.name)
    else:
        raise FoyerError("Backend not supported." 'Supports backend: "gmso".')

    return tmp_processed_files


class Forcefield(object):
    """General Forcefield object that can be created by either a GMSO XML forcefield file of Foyer XML forcefield file.

    Parameters
    ----------
    forcefield_files : list of str, optional, default=None
        List of forcefield files to load.
    name : str, optional, None
        Name of a forcefield to load that is packaged within foyer.
    backend : str, optional, default='gmso'
        Name of the backend used to store potential types information.
        At this point, 'gmso' is the only valid backend, but this set up
        allow future backend to be implemented more easily.
    """

    def __init__(
        self, forcefield_files=None, name=None, backend="gmso", **kwargs
    ):
        self.atomTypeDefinitions = dict()
        self.atomTypeOverrides = dict()
        self.atomTypeDesc = dict()
        self.atomTypeRefs = dict()
        self.atomTypeClasses = dict()
        self.atomTypeElements = dict()
        self._included_forcefields = dict()
        self.backend = backend
        self.non_element_types = dict()
        self._version = None
        self._name = None

        all_files_to_load = []
        if forcefield_files is not None:
            if isinstance(forcefield_files, (list, tuple, set)):
                all_files_to_load = list(forcefield_files)
            else:
                all_files_to_load = [forcefield_files]

        if name is not None:
            try:
                file = self.included_forcefields[name]
            except KeyError:
                raise IOError("Forcefield {} cannot be found.".format(name))
            else:
                all_files_to_load = [file]

        # Preprocessed the input files
        preprocessed_files = preprocess_forcefield_files(
            all_files_to_load, backend=backend
        )

        # Load in an internal forcefield object depends on given backend
        if backend == "gmso":
            self._parse_gmso(preprocessed_files, **kwargs)
        else:
            raise FoyerError(
                "Backend not supported." 'Supoprts backend: "gmso".'
            )

        # Remove the temporary files afterward
        for ff_file_name in preprocessed_files:
            os.remove(ff_file_name)

        self.parser = smarts.SMARTS(self.non_element_types)

    @property
    def version(self):
        """Return version of the loaded Forcefield."""
        return self._version

    @property
    def name(self):
        """Return name of the loaded Forcefield."""
        return self._name

    @property
    def included_forcefields(self):
        """Access forcefield included with foyer."""
        if any(self._included_forcefields):
            return self._included_forcefields

        ff_dir = resource_filename("foyer", "forcefields")
        ff_filepaths = set(glob.glob(os.path.join(ff_dir, "xml/*.xml")))

        for ff_filepath in ff_filepaths:
            _, ff_file = os.path.split(ff_filepath)
            basename, _ = os.path.splitext(ff_file)
            self._included_forcefields[basename] = ff_filepath
        return self._included_forcefields

    # Parse forcefield meta information
    def _parse_gmso(self, forcefield_files, **kwargs):
        """Parse metadata information when using GMSO as backend."""
        if forcefield_files:
            self.ff = gmso.ForceField(forcefield_files, **kwargs)
            self._version = self.ff.version
            self._name = self.ff.name
            for name, atype in self.ff.atom_types.items():
                self.atomTypeDefinitions[name] = atype.definition
                self.atomTypeOverrides[name] = atype.overrides
                self.atomTypeDesc[name] = atype.description
                self.atomTypeRefs[name] = set(atype.doi.split(","))
                self.atomTypeClasses[name] = atype.atomclass
                if atype.tags.get("element"):
                    ele = atype.tags["element"]
                    if element.element_by_symbol(ele):
                        self.atomTypeElements[atype.name] = ele
                    else:
                        self.non_element_types[ele] = None
                else:
                    # Register atomtype with missing atomtype as atomistic (empty string)
                    self.atomTypeElements[name] = ""
        return None

    def apply(
        self,
        top,
        references_file=None,
        use_residue_map=True,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=True,
        combining_rule="geometric",
        debug=False,
        *args,
        **kwargs,
    ):
        """Apply the force field to a molecular topology.

        Parameters
        ----------
        top : gmso.Topology, and structures/files that can be loaded by mbuild
            Molecular Topology to apply the force field to.
        references_file : str, optional, defaut=None
            Name of file where force field references will be written
            (in Bibtex format).
        use_residue_map : bool, optional, default=True
            Options to speed up if there are a lot of repeated
            subtopology within the topology (assuming they all
            have the same name).
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system dihedrals.
        assert_improper_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system impropers.
        combining_rule : str, optional, default='geometric'
            The combining rule of the system, stored as an attribute of the
            GMSO topology. Available options: 'geometric', 'lorentz'.
        debug : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potential problematic details it encounters.
        """
        if self.atomTypeDefinitions == {}:
            raise FoyerError(
                "Attempting to atom-type using a forcefield "
                "with no atom type definitions."
            )

        if self.backend == "gmso":
            if not isinstance(top, gmso.Topology):
                mb = import_("mbuild")
                tmp_top = mb.load(top)
                top = gmso.external.from_mbuild(tmp_top)

            assert isinstance(top, gmso.Topology)
            typemap = self._run_atomtyping(
                top, use_residue_map=use_residue_map, **kwargs
            )

        return self._parametrize(
            top=top,
            typemap=typemap,
            references_file=references_file,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
            debug=debug,
            combining_rule=combining_rule,
            *args,
            **kwargs,
        )

    def _run_atomtyping(self, top, use_residue_map=True, **kwargs):
        """Atomtype the topology.

        Parameters
        ----------
        top : gmso.Topology
            Molecular Topology to be atomtyped.
        use_residue_map : bool, optional, default=True
            A speed-up option that utilizes previously atom typed
            molecules as a template for future atom typing of
            identical molecules, instead of reevaluating them each
            individually. To be implemented
        """
        if isinstance(top, mb.Compound):
            top = from_mbuild(top)
        # TO DO in another PR
        if use_residue_map:
            # Detect duplicates subtopology/residues
            # (do matching by name, assert same number
            # of atoms)
            # Not implemented yet
            typemap = find_atomtypes(top, forcefield=self)
        else:
            typemap = find_atomtypes(top, forcefield=self)

        return typemap

    def _parametrize(
        self,
        top=None,
        typemap=dict(),
        references_file=None,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=True,
        combining_rule="geometric",
        debug=False,
        **kwargs,
    ):
        """Parametrize the Topology from the typemap provided.

        Assign AtomTypes and BondTypes to Atoms and Bonds, respectively.
        Create Angles, Dihedrals, Impropers and assign corresponding
        AngleTypes, DihedralTypes, and ImproperTypes.

        Parameters
        ----------
        top : gmso.Topology
            gmso.Topology object that needed to be parametrized
        typemap : dict
            typemap generated by the atomtyper
        references_file : str, optional, default=None
            Name of the file where force field references will be written
            (in Bibtex format).
        assert_bond_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system bonds.
        assert_angle_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system angles.
        assert_dihedral_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system dihedrals.
        assert_improper_params : bool, optional, default=True
            If True, Foyer will exit if parameters are not found for all
            system impropers.
        combining_rule : str, optional, default='geometric'
            The combining rule of the system, stored as an attribute of the
            GMSO topology. Available options: 'geometric', 'lorentz'.
        debug : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potential problematic details it encounters.

        Returns
        -------
        top : gmso.Topology
            A parametrized Topology.
        """
        # Generate missing connection objects (angles, dihedrals, and impropers)
        if isinstance(top, mb.Compound):
            top = gmso.external.from_mbuild(top)
        top.identify_connections()

        if self.backend == "gmso":
            self._parametrize_gmsoFF(
                top=top, typemap=typemap, combining_rule=combining_rule
            )
        else:
            raise FoyerError("Backend not supported")

        top.update_topology()
        self._check_parameters(
            top,
            assert_bond_params,
            assert_angle_params,
            assert_dihedral_params,
            assert_improper_params,
            debug,
        )

        if references_file:
            atom_types = set(site.atom_type for site in top.sites)
            self._write_references_to_file(atom_types, references_file)

        top.typed = True
        return top

    def _parametrize_gmsoFF(self, top, typemap, combining_rule):
        """Parametrize a Topology with gmso.ForceField."""
        # Assign AtomTypes
        for atom in top.sites:
            atom.atom_type = deepcopy(
                self.ff.get_potential(
                    "atom_type", typemap[top.get_index(atom)]["atomtype"]
                )
            )
        if not all(a.atom_type for a in top.sites):
            raise ValueError("Not all atoms in topology have atom types")

        # Assign BondTypes
        for bond in top.bonds:
            self._connection_type_lookup(bond)
        # Assign AngleTypes
        for angle in top.angles:
            self._connection_type_lookup(angle)
        # Assign DihedralTypes
        for dihedral in top.dihedrals:
            self._connection_type_lookup(dihedral)
        # Assign ImproperTypes
        for improper in top.impropers:
            self._connection_type_lookup(improper)
        # Assign combining rules
        top.combining_rule = combining_rule
        return top

    def _connection_type_lookup(self, connection):
        if isinstance(connection, gmso.Bond):
            bmem = connection.connection_members
            btype_name = [bmem[i].atom_type.name for i in range(2)]
            btype_class = [bmem[i].atom_type.atomclass for i in range(2)]
            for name in [btype_name, btype_class]:
                connection.bond_type = deepcopy(
                    self.ff.get_potential("bond_type", name, warn=True)
                )
                if connection.bond_type:
                    break

        elif isinstance(connection, gmso.Angle):
            agmem = connection.connection_members
            agtype_name = [agmem[i].atom_type.name for i in range(3)]
            agtype_class = [agmem[i].atom_type.atomclass for i in range(3)]

            for name in [agtype_name, agtype_class]:
                connection.angle_type = deepcopy(
                    self.ff.get_potential("angle_type", name, warn=True)
                )
                if connection.angle_type:
                    break

        elif isinstance(connection, gmso.Dihedral):
            dmem = connection.connection_members
            dtype_name = [dmem[i].atom_type.name for i in range(4)]
            dtype_class = [dmem[i].atom_type.atomclass for i in range(4)]

            for name in [dtype_name, dtype_class]:
                connection.dihedral_type = deepcopy(
                    self.ff.get_potential("dihedral_type", name, warn=True)
                )
                if connection.dihedral_type:
                    break

        elif isinstance(connection, gmso.Improper):
            imem = connection.connection_members
            itype_name = [imem[i].atom_type.name for i in range(4)]
            itype_class = [imem[i].atom_type.atomclass for i in range(4)]

            for name in [itype_name, itype_class]:
                connection.improper_type = deepcopy(
                    self.ff.get_potential("improper_type", name, warn=True)
                )
                if connection.improper_type:
                    break

    def _check_parameters(
        self,
        top,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=True,
        debug=False,
    ):
        """Check if the parameters are fulfilled.

        Parameters
        ----------
        top : gmso.Topology
            Molecular Topology whose parameters are being checked.
        assert_bond_params : bool, optional, default=True
            Check if all bonds have params
        assert_angle_params : bool, optional, default=True
            Check if all angles have params
        assert_dihedral_params : bool, optional, default=True
            Check if all dihedrals have params
        assert_improper_params : bool, optional, default=True
            Check if all impropers have params
        debug : bool, optional, default=False
            If True, Foyer will print debug-level information about notable or
            potential problematic details it encounters.
        """
        missing_bond_params = dict()
        missing_angle_params = dict()
        missing_dihedral_params = dict()
        missing_improper_params = dict()

        for bond in top.bonds:
            if not bond.bond_type:
                missing_bond_params[bond.name] = [
                    a.atom_type.name for a in bond.connection_members
                ]
        for angle in top.angles:
            if not angle.angle_type:
                missing_angle_params[angle.name] = [
                    a.atom_type.name for a in angle.connection_members
                ]
        for dihedral in top.dihedrals:
            if not dihedral.dihedral_type:
                missing_dihedral_params[dihedral.name] = [
                    a.atom_type.name for a in dihedral.connection_members
                ]
        for improper in top.impropers:
            if not improper.improper_type:
                missing_improper_params[improper.name] = [
                    a.atom_type.name for a in improper.connection_members
                ]

        if debug:
            from pprint import pprint

            if missing_bond_params:
                print("Bonds with missing parameters: ")
                pprint(missing_bond_params)
            if missing_angle_params:
                print("Angles with missing parameters: ")
                pprint(missing_angle_params)
            if missing_dihedral_params:
                print("Dihedral with missing parameters: ")
                pprint(missing_dihedral_params)
            if missing_improper_params:
                print("Improper with missing parameters: ")
                pprint(missing_improper_params)

        if assert_bond_params and missing_bond_params:
            raise FoyerError(
                "Some bonds are missing parameters. "
                "Change debug=True for more information"
            )
        if assert_angle_params and missing_angle_params:
            raise FoyerError(
                "Some angles are missing parameters. "
                "Change debug=True for more information"
            )
        if assert_dihedral_params and missing_dihedral_params:
            raise FoyerError(
                "Some dihedrals are missing parameters. "
                "Change debug=True for more information"
            )
        if assert_improper_params and missing_improper_params:
            raise FoyerError(
                "Some impropers are missing parameters. "
                "Change debug=True for more information"
            )
        return (
            missing_bond_params,
            missing_angle_params,
            missing_dihedral_params,
            missing_improper_params,
        )

    def _write_references_to_file(self, atom_types, references_file):
        atomtype_references = {}
        for atype in atom_types:
            try:
                atomtype_references[atype.name] = self.atomTypeRefs[atype.name]
            except KeyError:
                warnings.warn(
                    "Reference not found for atom type '{}'." "".format(atype)
                )
        unique_references = collections.defaultdict(list)
        for atomtype, dois in atomtype_references.items():
            for doi in dois:
                unique_references[doi].append(atomtype)
        unique_references = collections.OrderedDict(
            sorted(unique_references.items())
        )
        with open(references_file, "w") as f:
            for doi, atomtypes in unique_references.items():
                url = "http://api.crossref.org/works/{}/transform/application/x-bibtex".format(
                    doi
                )
                headers = {"accept": "application/x-bibtex"}
                bibtex_ref = get_ref(url, headers=headers)
                if bibtex_ref is None:
                    warnings.warn("Could not get ref for doi {}".format(doi))
                    continue
                else:
                    bibtex_text = bibtex_ref.text
                note = (
                    ",\n\tnote = {Parameters for atom types: "
                    + ", ".join(sorted(atomtypes))
                    + "}"
                )
                bibtex_text = bibtex_text[:-2] + note + bibtex_text[-2:]
                f.write("{}\n".format(bibtex_text))


# TO DO
# gmso.Topology.write_foyer = write_foyer
