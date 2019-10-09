from __future__ import division

import collections
from lxml import etree as ET
from foyer.smarts_graph import SMARTSGraph
import networkx as nx
import warnings

import numpy as np


def write_foyer(self, filename, forcefield=None, unique=True):
    """Outputs a Foyer XML from a ParmEd Structure

    Information from a ParmEd Structure is used to create a Foyer force
    field XML. If the Forcefield used to parameterize the Structure is
    passed to this function through the `forcefield` argument then the
    resulting XML should be able to be used to exacty reproduce the
    parameterization. Otherwise, all topological information will be
    written, but the `AtomTypes` section will be missing information
    (such as SMARTS definitions).

    This function can also be used to output the complete topology for a
    system if the `unique` argument is set to `False`. This can be useful
    for small molecules to be used to test the implementation of a Foyer
    force field XML.

    Additional features to be added include:
    - Documentation (and better standardization) of parameter rounding.
    - Grouping of duplicate parameter definitions by atom class (the current
      implementation considers parameters to be unique if atom `type`
      definitions differ)
    - Sort parameters by atom ids when `unique=False` (currently no sorting
      is performed, this is minor but would make files slightly more readable)

    Parameters
    ----------
    filename : str
        Name of the Foyer XML file to be written
    forcefield : foyer.Forcefield, optional, default=None
        Foyer Forcefield used to parameterize the ParmEd Structure. This
        is used to obtain additional information that is not available
        from the Structure itself, e.g. atomtype SMARTS definitions.
    unique : boolean, optional, default=True
        Write only unique elements. If False, elements are written for each
        atom, bond, etc. in the system. `unique=False` is primarily used
        for storing the topology of "test" molecules for a Foyer forcefield.

    """
    # Assume if a Structure has a bond and bond type that the Structure is
    # parameterized. ParmEd uses the same logic to denote parameterization.
    if not (len(self.bonds) > 0 and self.bonds[0].type is not None):
        raise Exception("Cannot write Foyer XML from an unparametrized "
                        "Structure.")

    root = ET.Element('ForceField')
    _write_atoms(self, root, self.atoms, forcefield, unique)
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


def _write_atoms(self, root, atoms, forcefield, unique):
    atomtypes = ET.SubElement(root, 'AtomTypes')
    nonbonded = ET.SubElement(root, 'NonbondedForce')
    nonbonded.set('coulomb14scale', str(_infer_coulomb14scale(self)))
    nonbonded.set('lj14scale', str(_infer_lj14scale(self)))
    atomtype_attrs = collections.OrderedDict([
        ('name', 'name'),
        ('class', 'forcefield.atomTypeClasses[name]'),
        ('element', 'forcefield.atomTypeElements[name]'),
        ('mass', 'atom.mass'),
        ('def', 'forcefield.atomTypeDefinitions[name]'),
        ('desc', 'forcefield.atomTypeDesc[name]'),
        ('doi', 'forcefield.atomTypeRefs[name]'),
        ('overrides', 'forcefield.atomTypeOverrides[name]')
        ])
    atom_type_set = set([atom.atom_type.name for atom in atoms])
    for atom in atoms:
        atomtype = ET.SubElement(atomtypes, 'Type')
        nb_force = ET.SubElement(nonbonded, 'Atom')

        name = atom.atom_type.name
        for key, val in atomtype_attrs.items():
            '''
            I'm not crazy about using `eval` here, but this is where we want
            to handle `AttributeError` and `KeyError` exceptions to pass a
            blank string for these attributes.
            '''
            try:
                if key == 'doi':
                    label = eval(val)#[a for a in label]
                    label = ','.join([a for a in label])
                elif key == 'overrides':
                # Only write overrides atomtypes if they are in atom_type_set
                    label = []
                    original_label = []
                    for item in eval(val):
                        original_label.append(item)
                        if item in atom_type_set:
                            label.append(item)
                    if len(label) == 0:
                        label = ''
                    else:
                        label = ','.join([a for a in label])
                # Write out the original overrides atomtypes as a comment
                    atomtype.append(ET.Comment('Note: original overrides=\"{}\"'.format(','.join([a for a in original_label]))))
                else:
                    label = str(eval(val))
            except (AttributeError, KeyError):
                label = ''
            atomtype.set(key, label)

        if not unique:
            nb_force.set('id', str(atom.idx))
        nb_force.set('type', name)
        nb_force.set('charge', str(round(atom.charge, 4)))
        nb_force.set('sigma', str(round(atom.atom_type.sigma/10, 4)))
        nb_force.set('epsilon', str(round(atom.atom_type.epsilon * 4.184, 6)))

    _update_defs(atomtypes, nonbonded, forcefield)

def _update_defs(atomtypes, nonbonded, forcefield):
    def_list = [i.get('def') for i in atomtypes.iterchildren()]
    name_list = [i.get('name') for i in atomtypes.iterchildren()]
    smarts_list = list()
    smarts_parser = forcefield.parser
    for smarts_string, name in zip(def_list, name_list):
        smarts_graph = SMARTSGraph(smarts_string, parser=smarts_parser,
                                   name=name)
        for atom_expr in nx.get_node_attributes(smarts_graph, name='atom').values():
            labels = atom_expr.find_data('has_label')
            for label in labels:
                atom_type = label.children[0][1:]
                smarts_list.append(atom_type)
    smarts_list = list(set(smarts_list))
    extra_types = [i for i in smarts_list if i not in name_list]

    for extra in extra_types:
        for i, definition in enumerate(def_list):
            if extra in definition:
                warnings.warn('Removing undefined atom type `{}`'
                        ' from SMARTS string `{}`'.format(
                            extra, definition))
                extra_edit = '%' + extra
                extra_index = definition.find(extra_edit)
                if definition[extra_index-1] == ';':
                    new_def = definition.replace(extra_edit + ',' , '')
                else:
                    new_def = definition.replace(',' + extra_edit, '')
                atomtypes[i].set('def', new_def)

def _write_bonds(root, bonds, unique):
    bond_forces = ET.SubElement(root, 'HarmonicBondForce')
    for bond in bonds:
        bond_force = ET.SubElement(bond_forces, 'Bond')
        atypes = [atype for atype in [bond.atom1.type, bond.atom2.type]]
        if unique:
            atypes = sorted(atypes)
        else:
            bond_force.set('id1', str(bond.atom1.idx))
            bond_force.set('id2', str(bond.atom2.idx))
        for id in range(2):
            bond_force.set('type{}'.format(id+1), atypes[id])
        bond_force.set('length', str(round(bond.type.req / 10, 4)))
        bond_force.set('k', str(round(bond.type.k * 4.184 * 200, 1)))


def _write_angles(root, angles, unique):
    angle_forces = ET.SubElement(root, 'HarmonicAngleForce')
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
        for id in range(3):
            angle_force.set('type{}'.format(id+1), atypes[id])
        angle_force.set('angle', str(round(angle.type.theteq * (np.pi / 180), 10)))
        angle_force.set('k', str(round(angle.type.k * 4.184 * 2, 3)))


def _write_periodic_torsions(root, dihedrals, unique):
    periodic_torsion_forces = ET.SubElement(root, 'PeriodicTorsionForce')
    last_dihedral_force = None
    for dihedral in dihedrals:
        if dihedral.improper:
            dihedral_type = 'Improper'
        else:
            dihedral_type = 'Proper'
        dihedral_force = ET.SubElement(periodic_torsion_forces, dihedral_type)
        atypes = [atype for atype in [dihedral.atom1.type,
                                      dihedral.atom2.type,
                                      dihedral.atom3.type,
                                      dihedral.atom4.type]]
        if dihedral.improper:
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
            if unique:
                if atypes[0] > atypes[-1]:
                    atypes = atypes[::-1]
            else:
                dihedral_force.set('id1', str(dihedral.atom1.idx))
                dihedral_force.set('id2', str(dihedral.atom2.idx))
                dihedral_force.set('id3', str(dihedral.atom3.idx))
                dihedral_force.set('id4', str(dihedral.atom4.idx))
        for id in range(4):
            dihedral_force.set('type{}'.format(id+1), atypes[id])
        dihedral_force.set('periodicity1', str(dihedral.type.per))
        dihedral_force.set('phase1',
                           str(round(dihedral.type.phase * (np.pi / 180), 8)))
        dihedral_force.set('k1', str(round(dihedral.type.phi_k * 4.184, 3)))
        if last_dihedral_force is not None:
            # Check to see if this current dihedral force needs to be
            # "merged" into the last dihedral force
            last_dihedral_tuple = (last_dihedral_force.attrib['type1'],
                                    last_dihedral_force.attrib['type2'],
                                    last_dihedral_force.attrib['type3'],
                                    last_dihedral_force.attrib['type4'])
            current_dihedral_tuple = (dihedral_force.attrib['type1'],
                                        dihedral_force.attrib['type2'],
                                        dihedral_force.attrib['type3'],
                                        dihedral_force.attrib['type4'])
            if (last_dihedral_tuple == current_dihedral_tuple and 
                _unique_periodictorsion_parameters(last_dihedral_force,
                    dihedral_force)):
                # Merge the last and current dihedral forces
                # Find the nth periodicity we can set
                n = 1
                while 'periodicity{}'.format(n) in last_dihedral_force.attrib:
                    n +=1
                last_dihedral_force.attrib['periodicity{}'.format(n)] = \
                        dihedral_force.attrib['periodicity1']
                last_dihedral_force.attrib['phase{}'.format(n)] = \
                        dihedral_force.attrib['phase1']
                last_dihedral_force.attrib['k{}'.format(n)] = \
                        dihedral_force.attrib['k1']
                periodic_torsion_forces.remove(dihedral_force)
            else:
                last_dihedral_force = dihedral_force

        else:
            last_dihedral_force = dihedral_force


def _unique_periodictorsion_parameters(dihedral1, dihedral2):
    """ Return true if dihedral1 contains the parameters of dihedral2

    Parameters
    ---------
    dihedral1: ET.subelement
        This is the "larger" dihedral ETelement that is collecting multiple
        periodicities
    dihedral2: ET.subelement
        This should only contain periodicity1, phase1, k1 attributes
    """
    n = 1
    param_tuples = set()
    while 'periodicity{}'.format(n) in dihedral1.attrib:
        param_tuples.add((dihedral1.attrib['periodicity{}'.format(n)],
                            dihedral1.attrib['phase{}'.format(n)],
                            dihedral1.attrib['k{}'.format(n)]))
        n+=1
    if (dihedral2.attrib['periodicity1'], dihedral2.attrib['phase1'], dihedral2.attrib['k1']) in param_tuples:
        return False
    else:
        return True


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
        for id in range(4):
            rb_torsion_force.set('type{}'.format(id+1), atypes[id])
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
    """
    Note: This was grabbed, basically verbatim, from:
    https://stackoverflow.com/questions/7905380/testing-equivalence-of-xml-etree-elementtree
    """
    if type(e1) != type(e2): return False
    if e1.tag != e2.tag: return False
    if e1.text != e2.text: return False
    if e1.tail != e2.tail: return False
    if e1.attrib != e2.attrib: return False
    if len(e1) != len(e2): return False
    return all([_elements_equal(c1, c2) for c1, c2 in zip(e1, e2)])


def _infer_coulomb14scale(struct):
    """Attempt to infer the coulombic 1-4 scaling factor by parsing the
    list of adjusts in the structure"""

    coul14 = [t.type.chgscale for t in struct.adjusts]

    if len(set(coul14)) == 1:
        return coul14[0]
    else:
        raise ValueError(
            'Structure has inconsistent 1-4 coulomb scaling factors. This is '
            'currently not supported'
        )


def _infer_lj14scale(struct):
    """Attempt to infer the Lennard-Jones 1-4 scaling factor by parsing the
    list of adjusts in the structure"""

    lj14scale = list()

    for adj in struct.adjusts:
        type1 = adj.atom1.atom_type
        type2 = adj.atom2.atom_type
        expected_sigma = (type1.sigma + type2.sigma) * 0.5
        expected_epsilon = (type1.epsilon * type2.epsilon) ** 0.5

        # We expect sigmas to be the same but epsilons to be scaled by a factor
        if not np.isclose(adj.type.sigma, expected_sigma):
            raise ValueError(
                'Unexpected 1-4 sigma value found in adj {}. Expected {}'
                'and found {}'.format(adj, adj.type.sigma, expected_sigma)
            )

        lj14scale.append(adj.type.epsilon/expected_epsilon)

    unique_lj14_scales = np.unique(np.array(lj14scale).round(8))
    if len(unique_lj14_scales) == 1:
        return lj14scale[0]
    else:
        raise ValueError(
            'Structure has inconsistent 1-4 LJ scaling factors. This is '
            'currently not supported. Found these factors: '
            '{}'.format(unique_lj14_scales)
        )
