try:
    from functools import lru_cache
except:
    from functools32 import lru_cache
import glob
import itertools
import os
from pkg_resources import resource_filename

import mbuild as mb
import parmed as pmd
import simtk.openmm.app.element as elem
import simtk.unit as unit
from simtk import openmm as mm
from simtk.openmm import app
from simtk.openmm.app.forcefield import (NoCutoff, CutoffNonPeriodic,
                                         CutoffPeriodic, Ewald, HBonds,
                                         AllBonds, HAngles, NonbondedGenerator)

from foyer.atomtyper import find_atomtypes
from foyer.exceptions import FoyerError


def generate_topology(non_omm_topology):
    topology = app.Topology()
    chain = topology.addChain()
    if isinstance(non_omm_topology, pmd.Structure):
        residue = topology.addResidue(non_omm_topology.title, chain)
        atoms = dict()  # pmd.Atom: omm.Atom

        # Create atoms in the residue.
        for pmd_atom in non_omm_topology.atoms:
            name = pmd_atom.name
            if pmd_atom.atomic_number != 0:
                element = elem.Element.getByAtomicNumber(pmd_atom.atomic_number)
            else:  # TODO: more robust element detection or enforcement of symbols
                element = elem.Element.getByAtomicNumber(pmd_atom.atomic_number)

            omm_atom = topology.addAtom(name, element, residue)
            atoms[pmd_atom] = omm_atom
            omm_atom.bond_partners = []

        # Create bonds.
        for bond in non_omm_topology.bonds:
            atom1 = atoms[bond.atom1]
            atom2 = atoms[bond.atom2]
            topology.addBond(atom1, atom2)
            atom1.bond_partners.append(atom2)
            atom2.bond_partners.append(atom1)

    elif isinstance(non_omm_topology, mb.Compound):
        residue = topology.addResidue(non_omm_topology.name, chain)
        atoms = dict()  # mb.Particle: omm.Atom

        # Create atoms in the residue.
        for mb_particle in non_omm_topology.particles():
            name = mb_particle.name
            element = elem.Element.getBySymbol(mb_particle.name)

            omm_atom = topology.addAtom(name, element, residue)
            atoms[mb_particle] = omm_atom
            omm_atom.bond_partners = []

        # Create bonds.
        for atom1, atom2 in non_omm_topology.bonds():
            topology.addBond(atom1, atom2)
            atom1.bond_partners.append(atom2)
            atom2.bond_partners.append(atom1)
    else:
        raise FoyerError('Unknown topology format: {}\n'
                         'Supported formats are: '
                         '"parmed.Structure", '
                         '"openmm.app.Topology"'.format(topology))
    return topology


class Forcefield(app.ForceField):
    def __init__(self, *files):
        self._atomTypeDefinitions = dict()
        self._atomTypeOverrides = dict()
        self._atomTypeDesc = dict()
        super(Forcefield, self).__init__(*files)

    @classmethod
    @lru_cache(maxsize=32)
    def _available(cls):
        resource_extension = '.xml'
        resource_dir = resource_filename('foyer', 'forcefields')
        resource_files = set(glob.glob(os.path.join(resource_dir, '*'+resource_extension)))

        resource_dict = {}
        for resource_path in resource_files:
            base_path, resource_fn = os.path.split(resource_path)
            basename = resource_fn[:-len(resource_extension)]
            resource_dict[basename] = resource_path

        return resource_dict

    @classmethod
    def available(cls):
        return Forcefield._available().keys()

    @classmethod
    def by_name(cls, forcefield_name):
        if forcefield_name in cls.available():
            return Forcefield(cls._available()[forcefield_name])
        else:
            raise IOError('Forcefield {} cannot be found'.format(forcefield_name))

    def registerAtomType(self, parameters):
        """Register a new atom type. """
        super(Forcefield, self).registerAtomType(parameters)

        name = parameters['name']
        if 'def' in parameters:
            self._atomTypeDefinitions[name] = parameters['def']
        if 'overrides' in parameters:
            self._atomTypeOverrides[name] = parameters['overrides']
        if 'des' in parameters:
            self._atomTypeDesc[name] = parameters['desc']

    def apply(self, topology, *args, **kwargs):
        if not isinstance(topology, app.Topology):
            topology = generate_topology(topology)
        else:
            pass

        system = self.createSystem(topology, *args, **kwargs)
        return pmd.openmm.load_topology(topology=topology, system=system)

    def createSystem(self, topology, nonbondedMethod=NoCutoff,
                     nonbondedCutoff=1.0*unit.nanometer, constraints=None,
                     rigidWater=True, removeCMMotion=True, hydrogenMass=None,
                     residueTemplates=dict(), verbose=False, **args):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        residueTemplates : dict=dict()
           Key: Topology Residue object
           Value: string, name of _TemplateData residue template object to use for
                  (Key) residue
           This allows user to specify which template to apply to particular Residues
           in the event that multiple matching templates are available (e.g Fe2+ and Fe3+
           templates in the ForceField for a monoatomic iron ion in the topology).
        args
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        if isinstance(topology, pmd.Structure):
            topology = generate_topology(topology)
        elif isinstance(topology, app.Topology):
            pass
        else:
            raise FoyerError('Unknown topology format: {}\n'
                             'Supported formats are: '
                             '"parmed.Structure", '
                             '"openmm.app.Topology"'.format(topology))

        # Atomtype the system.
        find_atomtypes(atoms=list(topology.atoms()),
                       forcefield=self,
                       debug=verbose)

        data = app.ForceField._SystemData()
        data.atoms = list(topology.atoms())
        for atom in data.atoms:
            data.excludeAtomWith.append([])

        # Make a list of all bonds
        for bond in topology.bonds():
            data.bonds.append(app.ForceField._BondData(bond[0].index, bond[1].index))

        # Record which atoms are bonded to each other atom
        bondedToAtom = []
        for i in range(len(data.atoms)):
            bondedToAtom.append(set())
            data.atomBonds.append([])
        for i in range(len(data.bonds)):
            bond = data.bonds[i]
            bondedToAtom[bond.atom1].add(bond.atom2)
            bondedToAtom[bond.atom2].add(bond.atom1)
            data.atomBonds[bond.atom1].append(i)
            data.atomBonds[bond.atom2].append(i)

        # TODO: Better way to lookup nonbonded parameters...?
        for generator in self.getGenerators():
            if isinstance(generator, NonbondedGenerator):
                nonbonded_params = generator.params.paramsForType
                break

        for chain in topology.chains():
            for res in chain.residues():
                for atom in res.atoms():
                    data.atomType[atom] = atom.type
                    params = nonbonded_params[atom.type]
                    params['charge'] = 0
                    data.atomParameters[atom] = params

        # Create the System and add atoms
        sys = mm.System()
        for atom in topology.atoms():
            # Look up the atom type name, returning a helpful error message if it cannot be found.
            if atom not in data.atomType:
                raise Exception("Could not identify atom type for atom '%s'." % str(atom))
            typename = data.atomType[atom]

            # Look up the type name in the list of registered atom types, returning a helpful error message if it cannot be found.
            if typename not in self._atomTypes:
                msg  = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
                msg += "Known atom types are: %s" % str(self._atomTypes.keys())
                raise Exception(msg)

            # Add the particle to the OpenMM system.
            mass = self._atomTypes[typename].mass
            sys.addParticle(mass)

        # Adjust hydrogen masses if requested.
        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.
        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Make a list of all unique angles
        uniqueAngles = set()
        for bond in data.bonds:
            for atom in bondedToAtom[bond.atom1]:
                if atom != bond.atom2:
                    if atom < bond.atom2:
                        uniqueAngles.add((atom, bond.atom1, bond.atom2))
                    else:
                        uniqueAngles.add((bond.atom2, bond.atom1, atom))
            for atom in bondedToAtom[bond.atom2]:
                if atom != bond.atom1:
                    if atom > bond.atom1:
                        uniqueAngles.add((bond.atom1, bond.atom2, atom))
                    else:
                        uniqueAngles.add((atom, bond.atom2, bond.atom1))
        data.angles = sorted(list(uniqueAngles))

        # Make a list of all unique proper torsions
        uniquePropers = set()
        for angle in data.angles:
            for atom in bondedToAtom[angle[0]]:
                if atom not in angle:
                    if atom < angle[2]:
                        uniquePropers.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        uniquePropers.add((angle[2], angle[1], angle[0], atom))
            for atom in bondedToAtom[angle[2]]:
                if atom not in angle:
                    if atom > angle[0]:
                        uniquePropers.add((angle[0], angle[1], angle[2], atom))
                    else:
                        uniquePropers.add((atom, angle[2], angle[1], angle[0]))
        data.propers = sorted(list(uniquePropers))

        # Make a list of all unique improper torsions
        for atom in range(len(bondedToAtom)):
            bondedTo = bondedToAtom[atom]
            if len(bondedTo) > 2:
                for subset in itertools.combinations(bondedTo, 3):
                    data.impropers.append((atom, subset[0], subset[1], subset[2]))

        # Identify bonds that should be implemented with constraints
        if constraints == AllBonds or constraints == HAngles:
            for bond in data.bonds:
                bond.isConstrained = True
        elif constraints == HBonds:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                bond.isConstrained = atom1.name.startswith('H') or atom2.name.startswith('H')
        if rigidWater:
            for bond in data.bonds:
                atom1 = data.atoms[bond.atom1]
                atom2 = data.atoms[bond.atom2]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH':
                    bond.isConstrained = True

        # Identify angles that should be implemented with constraints
        if constraints == HAngles:
            for angle in data.angles:
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                numH = 0
                if atom1.name.startswith('H'):
                    numH += 1
                if atom3.name.startswith('H'):
                    numH += 1
                data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.name.startswith('O')))
        else:
            data.isAngleConstrained = len(data.angles)*[False]
        if rigidWater:
            for i in range(len(data.angles)):
                angle = data.angles[i]
                atom1 = data.atoms[angle[0]]
                atom2 = data.atoms[angle[1]]
                atom3 = data.atoms[angle[2]]
                if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                    data.isAngleConstrained[i] = True

        # Add virtual sites
        for atom in data.virtualSites:
            (site, atoms, excludeWith) = data.virtualSites[atom]
            index = atom.index
            data.excludeAtomWith[excludeWith].append(index)
            if site.type == 'average2':
                sys.setVirtualSite(index, mm.TwoParticleAverageSite(atoms[0], atoms[1], site.weights[0], site.weights[1]))
            elif site.type == 'average3':
                sys.setVirtualSite(index, mm.ThreeParticleAverageSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'outOfPlane':
                sys.setVirtualSite(index, mm.OutOfPlaneSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1], site.weights[2]))
            elif site.type == 'localCoords':
                sys.setVirtualSite(index, mm.LocalCoordinatesSite(atoms[0], atoms[1], atoms[2],
                                                                  mm.Vec3(site.originWeights[0], site.originWeights[1], site.originWeights[2]),
                                                                  mm.Vec3(site.xWeights[0], site.xWeights[1], site.xWeights[2]),
                                                                  mm.Vec3(site.yWeights[0], site.yWeights[1], site.yWeights[2]),
                                                                  mm.Vec3(site.localPos[0], site.localPos[1], site.localPos[2])))

        # Add forces to the System
        for force in self._forces:
            force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(sys, data, args)

        # Execute scripts found in the XML files.
        for script in self._scripts:
            exec(script, locals())
        return sys
    #
    # def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0 * unit.nanometer,
    #                  constraints=None, rigidWater=True, removeCMMotion=True, hydrogenMass=None, **args):
    #     """Construct an OpenMM System representing a Topology with this force field.
    #
    #     Parameters
    #     ----------
    #     topology : Topology
    #         The Topology for which to create a System
    #     nonbondedMethod : object=NoCutoff
    #         The method to use for nonbonded interactions.  Allowed values are
    #         NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
    #     nonbondedCutoff : distance=1*nanometer
    #         The cutoff distance to use for nonbonded interactions
    #     constraints : object=None
    #         Specifies which bonds and angles should be implemented with constraints.
    #         Allowed values are None, HBonds, AllBonds, or HAngles.
    #     rigidWater : boolean=True
    #         If true, water molecules will be fully rigid regardless of the value
    #         passed for the constraints argument
    #     removeCMMotion : boolean=True
    #         If true, a CMMotionRemover will be added to the System
    #     hydrogenMass : mass=None
    #         The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
    #         added to a hydrogen is subtracted from the heavy atom to keep
    #         their total mass the same.
    #     args
    #          Arbitrary additional keyword arguments may also be specified.
    #          This allows extra parameters to be specified that are specific to
    #          particular force fields.
    #
    #     Returns
    #     -------
    #     system
    #         the newly created System
    #     """
    #     data = app.ForceField._SystemData()
    #     data.atoms = list(topology.atoms())
    #     for atom in data.atoms:
    #         data.excludeAtomWith.append([])
    #
    #     # Make a list of all bonds
    #
    #     for bond in topology.bonds():
    #         data.bonds.append(app.ForceField._BondData(bond[0].index, bond[1].index))
    #
    #     # Record which atoms are bonded to each other atom
    #
    #     bondedToAtom = []
    #     for i in range(len(data.atoms)):
    #         bondedToAtom.append(set())
    #         data.atomBonds.append([])
    #     for i in range(len(data.bonds)):
    #         bond = data.bonds[i]
    #         bondedToAtom[bond.atom1].add(bond.atom2)
    #         bondedToAtom[bond.atom2].add(bond.atom1)
    #         data.atomBonds[bond.atom1].append(i)
    #         data.atomBonds[bond.atom2].append(i)
    #
    #     # Find the template matching each residue and assign atom types.
    #     # If no templates are found, attempt to use residue template generators to create new templates (and potentially atom types/parameters).
    #
    #     for chain in topology.chains():
    #         for res in chain.residues():
    #             # # Attempt to match one of the existing templates.
    #             # [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
    #             # if matches is None:
    #             #     # No existing templates match.  Try any registered residue template generators.
    #             #     for generator in self._templateGenerators:
    #             #         if generator(self, res):
    #             #             # This generator has registered a new residue template that should match.
    #             #             [template, matches] = self._getResidueTemplateMatches(res, bondedToAtom)
    #             #             if matches is None:
    #             #                 # Something went wrong because the generated template does not match the residue signature.
    #             #                 raise Exception(
    #             #                     'The residue handler %s indicated it had correctly parameterized residue %s, but the generated template did not match the residue signature.' % (
    #             #                     generator.__class__.__name__, str(res)))
    #             #             else:
    #             #                 # We successfully generated a residue template.  Break out of the for loop.
    #             #                 break
    #             #
    #             # # Raise an exception if we have found no templates that match.
    #             # if matches is None:
    #             #     raise ValueError('No template found for residue %d (%s).  %s' % (
    #             #     res.index + 1, res.name, _findMatchErrors(self, res)))
    #
    #             # Store parameters for the matched residue template.
    #             # matchAtoms = dict(zip(matches, res.atoms()))
    #             for atom in res.atoms():
    #                 data.atomType[atom] = atom.type
    #                 data.atomParameters[atom]
    #                 for site in template.virtualSites:
    #                     if match == site.index:
    #                         data.virtualSites[atom] = (
    #                         site, [matchAtoms[i].index for i in site.atoms], matchAtoms[site.excludeWith].index)
    #
    #     # Create the System and add atoms
    #
    #     sys = mm.System()
    #     for atom in topology.atoms():
    #         # Look up the atom type name, returning a helpful error message if it cannot be found.
    #         if atom not in data.atomType:
    #             raise Exception("Could not identify atom type for atom '%s'." % str(atom))
    #         typename = data.atomType[atom]
    #
    #         # Look up the type name in the list of registered atom types, returning a helpful error message if it cannot be found.
    #         if typename not in self._atomTypes:
    #             msg = "Could not find typename '%s' for atom '%s' in list of known atom types.\n" % (typename, str(atom))
    #             msg += "Known atom types are: %s" % str(self._atomTypes.keys())
    #             raise Exception(msg)
    #
    #         # Add the particle to the OpenMM system.
    #         mass = self._atomTypes[typename].mass
    #         sys.addParticle(mass)
    #
    #     # Adjust hydrogen masses if requested.
    #
    #     if hydrogenMass is not None:
    #         if not unit.is_quantity(hydrogenMass):
    #             hydrogenMass *= unit.dalton
    #         for atom1, atom2 in topology.bonds():
    #             if atom1.element == elem.hydrogen:
    #                 (atom1, atom2) = (atom2, atom1)
    #             if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
    #                 transferMass = hydrogenMass - sys.getParticleMass(atom2.index)
    #                 sys.setParticleMass(atom2.index, hydrogenMass)
    #                 sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index) - transferMass)
    #
    #     # Set periodic boundary conditions.
    #
    #     boxVectors = topology.getPeriodicBoxVectors()
    #     if boxVectors is not None:
    #         sys.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
    #     elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
    #         raise ValueError(
    #             'Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')
    #
    #     # Make a list of all unique angles
    #
    #     uniqueAngles = set()
    #     for bond in data.bonds:
    #         for atom in bondedToAtom[bond.atom1]:
    #             if atom != bond.atom2:
    #                 if atom < bond.atom2:
    #                     uniqueAngles.add((atom, bond.atom1, bond.atom2))
    #                 else:
    #                     uniqueAngles.add((bond.atom2, bond.atom1, atom))
    #         for atom in bondedToAtom[bond.atom2]:
    #             if atom != bond.atom1:
    #                 if atom > bond.atom1:
    #                     uniqueAngles.add((bond.atom1, bond.atom2, atom))
    #                 else:
    #                     uniqueAngles.add((atom, bond.atom2, bond.atom1))
    #     data.angles = sorted(list(uniqueAngles))
    #
    #     # Make a list of all unique proper torsions
    #
    #     uniquePropers = set()
    #     for angle in data.angles:
    #         for atom in bondedToAtom[angle[0]]:
    #             if atom not in angle:
    #                 if atom < angle[2]:
    #                     uniquePropers.add((atom, angle[0], angle[1], angle[2]))
    #                 else:
    #                     uniquePropers.add((angle[2], angle[1], angle[0], atom))
    #         for atom in bondedToAtom[angle[2]]:
    #             if atom not in angle:
    #                 if atom > angle[0]:
    #                     uniquePropers.add((angle[0], angle[1], angle[2], atom))
    #                 else:
    #                     uniquePropers.add((atom, angle[2], angle[1], angle[0]))
    #     data.propers = sorted(list(uniquePropers))
    #
    #     # Make a list of all unique improper torsions
    #
    #     for atom in range(len(bondedToAtom)):
    #         bondedTo = bondedToAtom[atom]
    #         if len(bondedTo) > 2:
    #             for subset in itertools.combinations(bondedTo, 3):
    #                 data.impropers.append((atom, subset[0], subset[1], subset[2]))
    #
    #     # Identify bonds that should be implemented with constraints
    #
    #     if constraints == AllBonds or constraints == HAngles:
    #         for bond in data.bonds:
    #             bond.isConstrained = True
    #     elif constraints == HBonds:
    #         for bond in data.bonds:
    #             atom1 = data.atoms[bond.atom1]
    #             atom2 = data.atoms[bond.atom2]
    #             bond.isConstrained = atom1.name.startswith('H') or atom2.name.startswith('H')
    #     if rigidWater:
    #         for bond in data.bonds:
    #             atom1 = data.atoms[bond.atom1]
    #             atom2 = data.atoms[bond.atom2]
    #             if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH':
    #                 bond.isConstrained = True
    #
    #     # Identify angles that should be implemented with constraints
    #
    #     if constraints == HAngles:
    #         for angle in data.angles:
    #             atom1 = data.atoms[angle[0]]
    #             atom2 = data.atoms[angle[1]]
    #             atom3 = data.atoms[angle[2]]
    #             numH = 0
    #             if atom1.name.startswith('H'):
    #                 numH += 1
    #             if atom3.name.startswith('H'):
    #                 numH += 1
    #             data.isAngleConstrained.append(numH == 2 or (numH == 1 and atom2.name.startswith('O')))
    #     else:
    #         data.isAngleConstrained = len(data.angles) * [False]
    #     if rigidWater:
    #         for i in range(len(data.angles)):
    #             angle = data.angles[i]
    #             atom1 = data.atoms[angle[0]]
    #             atom2 = data.atoms[angle[1]]
    #             atom3 = data.atoms[angle[2]]
    #             if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
    #                 data.isAngleConstrained[i] = True
    #
    #     # Add virtual sites
    #
    #     for atom in data.virtualSites:
    #         (site, atoms, excludeWith) = data.virtualSites[atom]
    #         index = atom.index
    #         data.excludeAtomWith[excludeWith].append(index)
    #         if site.type == 'average2':
    #             sys.setVirtualSite(index, mm.TwoParticleAverageSite(atoms[0], atoms[1], site.weights[0], site.weights[1]))
    #         elif site.type == 'average3':
    #             sys.setVirtualSite(index, mm.ThreeParticleAverageSite(atoms[0], atoms[1], atoms[2], site.weights[0],
    #                                                                   site.weights[1], site.weights[2]))
    #         elif site.type == 'outOfPlane':
    #             sys.setVirtualSite(index, mm.OutOfPlaneSite(atoms[0], atoms[1], atoms[2], site.weights[0], site.weights[1],
    #                                                         site.weights[2]))
    #         elif site.type == 'localCoords':
    #             sys.setVirtualSite(index, mm.LocalCoordinatesSite(atoms[0], atoms[1], atoms[2],
    #                                                               mm.Vec3(site.originWeights[0], site.originWeights[1],
    #                                                                       site.originWeights[2]),
    #                                                               mm.Vec3(site.xWeights[0], site.xWeights[1],
    #                                                                       site.xWeights[2]),
    #                                                               mm.Vec3(site.yWeights[0], site.yWeights[1],
    #                                                                       site.yWeights[2]),
    #                                                               mm.Vec3(site.localPos[0], site.localPos[1],
    #                                                                       site.localPos[2])))
    #
    #     # Add forces to the System
    #
    #     for force in self._forces:
    #         force.createForce(sys, data, nonbondedMethod, nonbondedCutoff, args)
    #     if removeCMMotion:
    #         sys.addForce(mm.CMMotionRemover())
    #
    #     # Let force generators do postprocessing
    #
    #     for force in self._forces:
    #         if 'postprocessSystem' in dir(force):
    #             force.postprocessSystem(sys, data, args)
    #
    #     # Execute scripts found in the XML files.
    #
    #     for script in self._scripts:
    #         exec(script, locals())
    #     return sys


#
# def create_angles(structure, node, neighbors):
#     """Add all possible angles around a node to a structure. """
#     for pair in itertools.combinations(neighbors, 2):
#         angle = pmd.Angle(pair[0], node, pair[1])
#         structure.angles.append(angle)
#
#
# def create_dihedrals(structure, node_1, neighbors_1, node_2, neighbors_2,
#                      pairs=True):
#     """Add all possible dihedrals around a pair of nodes to a structure. """
#     neighbors_1 = set(neighbors_1) - {node_2}
#     neighbors_2 = set(neighbors_2) - {node_1}
#
#     for pair in itertools.product(neighbors_1, neighbors_2):
#         if pair[0] != pair[1]:
#             dihedral = pmd.Dihedral(pair[0], node_1, node_2, pair[1])
#             # if hasattr(structure, 'parameterset'):
#             #     if structure.parameterset.dihedral_types:
#             #         structure.dihedrals.append(dihedral)
#             #     if structure.parameterset.rb_torsion_types:
#             #         structure.rb_torsions.append(dihedral)
#             # else:
#             structure.dihedrals.append(dihedral)
#             # if pairs:
#             #     pair = pmd.NonbondedException(pair[0], pair[1])
#             #     structure.adjusts.append(pair)
#
#
# def create_impropers(structure, node, neighbors):
#     """Add all possible impropers around a node to a structure. """
#     for triplet in itertools.combinations(neighbors, 3):
#         improper = pmd.Improper(node, triplet[0], triplet[1], triplet[2])
#         structure.impropers.append(improper)
