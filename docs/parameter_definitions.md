### Parameter definitions

Parameter definitions within force field XMLs follow the same conventions
as defined in the [OpenMM documentation](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields).
Currently, only certain functional forms for molecular forces are supported,
while future developments are expected to allow Foyer to support any desired
functional form, including reactive and tabulated potentials.
The currently supported functional forms for molecular forces are:

 - **Nonbonded**
   - [Lennard-Jones (12-6)](http://docs.openmm.org/7.0.0/userguide/application.html#nonbondedforce)
 - **Bonds**
   - [Harmonic](http://docs.openmm.org/7.0.0/userguide/application.html#harmonicbondforce)
 - **Angles**
   - [Harmonic](http://docs.openmm.org/7.0.0/userguide/application.html#harmonicangleforce)
 - **Torsions (proper)**
   - [Periodic](http://docs.openmm.org/7.0.0/userguide/application.html#periodictorsionforce)
   - [Ryckaert-Bellemans](http://docs.openmm.org/7.0.0/userguide/application.html#rbtorsionforce)
 - **Torsions (improper)**
   - [Periodic](http://docs.openmm.org/7.0.0/userguide/application.html#periodictorsionforce)

Definitions for each molecular force follow the OpenMM standard.

#### Classes vs. Types
OpenMM allows users to specify either a
[`class` or a `type`](http://docs.openmm.org/7.0.0/userguide/application.html#atom-types-and-atom-classes),
to define each particle within the force definition.
Here, `type` refers to a specific atom type (as defined in the
`<AtomTypes>` section), while `class` refers to a more general
description that can apply to multiple atomtypes (i.e. multiple atomtypes
may share the same class). This aids in limiting the number of force
definitions required in a force field XML, as many similar atom types may
share force parameters.

#### Assigning parameters by specificity
Foyer deviates from OpenMM's convention when matching force definitions in
a force field XML to instances of these forces in a molecular system.
In OpenMM, forces are assigned according to the first matching definition
in a force field XML, even if multiple matching definitions exist.
In contrast, Foyer assigns force parameters based on definition
specificity, where definitions containing more `type` attributes are
considered to be more specific.

**Example:**
```
<RBTorsionForce>
  <Proper class1="CT" class2="CT" class3="CT" class4="CT" c0="2.9288" c1="-1.4644" c2="0.2092" c3="-1.6736" c4="0.0" c5="0.0"/>
  <Proper type1="opls_961" type2="opls_136" type3="opls_136" type4="opls_136" c0="-0.987424" c1="0.08363" c2="-0.08368" c3="-0.401664" c4="1.389088" c5="0.0"/>
</RBTorsionForce>
```
Above, two proper torsions are defined, both describing a torsional force between
four tetrahedral carbons. However, the first definition features four `class`
attributes and zero `type` attributes, as this describes a general dihedral for
all tetrahedral carbons. The second definition features zero `class` attributes
and four `type` attributes, and describes a more specific dihedral for the case
where one end carbon is of type `'opls_961'` (a fluorinated carbon) and the
remaining three carbons are of type `'opls_136'` (alkane carbons). Now consider
we want to use a force field containing the above torsion definitions to assign
parameters to some molecular system that features partially fluorinated alkanes.
When assigning torsion parameters to a quartet of atoms where one end carbon is
fluorinated (`'opls_961'`) and the remaining three are hydrogenated (`'opls_136'`),
if using the OpenMM procedure for parameter assignment the more general
`'CT-CT-CT-CT'` torsion parameters (the first definition above) would be assigned
because this would be the first matching definition in the force field XML.
However, with Foyer, the second definition will be auto-detected as more specific,
due to the greater number of `type` attributes (4 vs. 0) and those parameters will
be assigned instead.

It should be noted that if two definitions feature the same specificity level
(i.e. the same number of `type` definitions) then automatic detection of precedence
is not possible and parameter assignment will follow the OpenMM procedure
whereby parameters from the first matching force definition in the XML will
be assigned.
