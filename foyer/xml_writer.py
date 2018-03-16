import collections
from lxml import etree as ET

import numpy as np


def write_foyer(self, filename, forcefield=None, unique=True):
    """
    To-dos:
    - Possibly standardize, and definitely document, the rounding for each
      parameter.
    - Break writers for each section into separate helper functions.
    - Add `coulomb14scale` and `lj14scale` attributes to `NonbondedForce`
    - Add a function that searches for class equivalences in force definitions
      with equivalent parameters. This would allow the generated XML to write
      parameter definitions with classes instead of types, reducing the XML size.
    - Sort by atom ids when `unique=False` (currently no sorting is performed)
    """
    if not (len(self.bonds) > 0 and self.bonds[0].type is not None):
        raise Exception("Cannot write Foyer XML from an unparametrized "
                        "Structure.")
    root = ET.Element('ForceField')
    _write_atoms(root, self.atoms, forcefield, unique)
    if len(self.bonds) > 0 and self.bonds[0].type is not None:
        _write_bonds(root, self.bonds, unique)
    if len(self.angles) > 0 and self.angles[0].type is not None:
        _write_angles(root, self.angles, unique)
    if len(self.dihedrals) > 0 and self.dihedrals[0].type is not None:
        _write_periodic_torsions(root, self.dihedrals, unique)
    if len(self.rb_torsions) > 0 and self.rb_torsions[0].type is not None:
        _write_rb_torsions(root, self.rb_torsions, unique)

    _remove_duplicate_elements(root, unique)

    ET.ElementTree(root).write(filename, pretty_print=True)

def _write_atoms(root, atoms, forcefield, unique):
    atomtypes = ET.SubElement(root, 'AtomTypes')
    nonbonded = ET.SubElement(root, 'NonbondedForce')
    atomtype_attrs = collections.OrderedDict([
        ('name', 'name'),
        ('class', 'forcefield.atomTypeClasses[name]'),
        ('element', 'forcefield.atomTypeElements[name]'),
        ('mass', 'atom.mass'),
        ('def', 'forcefield.atomTypeDefinitions[name]'),
        ('desc', 'forcefield.atomTypeDesc[name]'),
        ('doi', 'forcefield.atomTypeRefs[name]')
        ])
    nonbonded_attrs = collections.OrderedDict([
        ('type', 'name'),
        ('charge', 'round(atom.charge, 4)'),
        ('sigma', 'round(atom.atom_type.sigma, 4)'),
        ('epsilon', 'round(atom.atom_type.epsilon * 4.184, 6)')
        ])
    for atom in atoms:
        atomtype = ET.SubElement(atomtypes, 'Type')
        nb_force = ET.SubElement(nonbonded, 'Atom')

        name = atom.atom_type.name
        for key, val in atomtype_attrs.items():
            try:
                label = str(eval(val))
            except (AttributeError, KeyError):
                label = ''
            atomtype.set(key, label)

        if not unique:
            nb_force.set('id', str(atom.idx))
        for key, val in nonbonded_attrs.items():
            label = str(eval(val))
            nb_force.set(key, label)

def _write_bonds(root, bonds, unique):
    bond_forces = ET.SubElement(root, 'HarmonicBondForce')
    bond_attrs = collections.OrderedDict([
        ('type1', 'atypes[0]'),
        ('type2', 'atypes[1]'),
        ('length', 'round(bond.type.req / 10, 4)'),
        ('k', 'round(bond.type.k * 4.184 * 200, 1)')
        ])
    for bond in bonds:
        bond_force = ET.SubElement(bond_forces, 'Bond')
        atypes = [atype for atype in [bond.atom1.type, bond.atom2.type]]
        if unique:
            atypes = sorted(atypes)
        else:
            bond_force.set('id1', str(bond.atom1.idx))
            bond_force.set('id2', str(bond.atom2.idx))
        for key, val in bond_attrs.items():
            label = str(eval(val))
            bond_force.set(key, label)

def _write_angles(root, angles, unique):
    angle_forces = ET.SubElement(root, 'HarmonicAngleForce')
    angle_attrs = collections.OrderedDict([
        ('type1', 'atypes[0]'),
        ('type2', 'atypes[1]'),
        ('type3', 'atypes[2]'),
        ('angle', 'round(angle.type.theteq * (np.pi / 180.), 10)'),
        ('k', 'round(angle.type.k * 4.184 * 2, 3)')
        ])
    for angle in angles:
        angle_force = ET.SubElement(angle_forces, 'Angle')
        atypes = [atype for atype in [angle.atom1.type,
                                      angle.atom2.type,
                                      angle.atom3.type]]
        if unique:
            # Sort the first and last atom types
            atypes[::2] = sorted(atypes[::2])
        else:
            angle_force.set('id1', str(angle.atom1.idx))
            angle_force.set('id2', str(angle.atom2.idx))
            angle_force.set('id3', str(angle.atom3.idx))
        for key, val in angle_attrs.items():
            label = str(eval(val))
            angle_force.set(key, label)

def _write_periodic_torsions(root, dihedrals, unique):
    periodic_torsion_forces = ET.SubElement(root, 'PeriodicTorsionForce')
    dihedral_attrs = collections.OrderedDict([
        ('type1', 'atypes[0]'),
        ('type2', 'atypes[1]'),
        ('type3', 'atypes[2]'),
        ('type4', 'atypes[3]'),
        ('periodicity1', 'dihedral.type.per'),
        ('phase1', 'round(dihedral.type.phase * (np.pi / 180.), 8)')
        ('k1', 'round(dihedral.type.phi_k * 4.184, 3)')
        ])
    for dihedral in dihedrals:
        dihedral_force = ET.SubElement(periodic_torsion_forces, dihedral_type)
        atypes = [atype for atype in [dihedral.atom1.type,
                                      dihedral.atom2.type,
                                      dihedral.atom3.type,
                                      dihedral.atom4.type]]
        if dihedral.improper:
            dihedral_type = 'Improper'
            # We want the central atom listed first and then sort the
            # remaining atom types.
            atypes[0], atypes[2] = atypes[2], atypes[0]
            if unique:
                atypes[1:] = sorted(atypes[1:])
            else:
                dihedral_force.set('id1', str(dihedral.atom3.idx))
                dihedral_force.set('id2', str(dihedral.atom2.idx))
                dihedral_force.set('id3', str(dihedral.atom1.idx))
                dihedral_force.set('id4', str(dihedral.atom4.idx))
        else:
            dihedral_type = 'Proper'
            if unique:
                if atypes[0] > atypes[-1]:
                    atypes = atypes[::-1]
            else:
                dihedral_force.set('id1', str(dihedral.atom1.idx))
                dihedral_force.set('id2', str(dihedral.atom2.idx))
                dihedral_force.set('id3', str(dihedral.atom3.idx))
                dihedral_force.set('id4', str(dihedral.atom4.idx))
        for key, val in dihedral_attrs.items():
            label = str(eval(val))
            dihedral_force.set(key, label)

def _write_rb_torsions(root, rb_torsions, unique):
    rb_torsion_forces = ET.SubElement(root, 'RBTorsionForce')
    for rb_torsion in rb_torsions:
        rb_torsion_force = ET.SubElement(rb_torsion_forces, 'Proper')
        atypes = [atype for atype in [rb_torsion.atom1.type,
                                      rb_torsion.atom2.type,
                                      rb_torsion.atom3.type,
                                      rb_torsion.atom4.type]]
        if unique:
            if atypes[0] > atypes[-1]:
                atypes = atypes[::-1]
        else:
            rb_torsion_force.set('id1', str(rb_torsion.atom1.idx))
            rb_torsion_force.set('id2', str(rb_torsion.atom2.idx))
            rb_torsion_force.set('id3', str(rb_torsion.atom3.idx))
            rb_torsion_force.set('id4', str(rb_torsion.atom4.idx))
        rb_torsion_force.set('type1', atypes[0])
        rb_torsion_force.set('type2', atypes[1])
        rb_torsion_force.set('type3', atypes[2])
        rb_torsion_force.set('type4', atypes[3])
        for c_id in range(6):
            rb_torsion_force.set('c{}'.format(c_id),
                str(round(getattr(rb_torsion.type, 'c{}'.format(c_id)) * 4.184,
                          4)))

def _remove_duplicate_elements(root, unique):
    sortby = {'AtomTypes': ['name'],
              'HarmonicBondForce': ['type1', 'type2'],
              'HarmonicAngleForce': ['type1', 'type2', 'type3'],
              'PeriodicTorsionForce': ['type1', 'type2', 'type3', 'type4'],
              'RBTorsionForce': ['type1', 'type2', 'type3', 'type4'],
              'NonbondedForce': ['type']}
    prev = None
    for child in root:
        if not unique and child.tag != 'AtomTypes':
            continue
        child[:] = sorted(child, key=lambda elem:
            tuple(elem.get(id) for id in sortby[child.tag]))
        elems_to_remove = []
        for elem in child:
            if _elements_equal(elem, prev):
                elems_to_remove.append(elem)
                continue
            prev = elem
        for elem_to_remove in elems_to_remove:
            child.remove(elem_to_remove)

def _elements_equal(e1, e2):
    if type(e1) != type(e2): return False
    if e1.tag != e1.tag: return False
    if e1.text != e2.text: return False
    if e1.tail != e2.tail: return False
    if e1.attrib != e2.attrib: return False
    if len(e1) != len(e2): return False
    return all([_elements_equal(c1, c2) for c1, c2 in zip(e1, e2)])
