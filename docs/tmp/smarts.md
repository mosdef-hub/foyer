### SMARTS-based atomtyping

Foyer allows users to describe atomtypes using a modified version of 
[SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
You may already be familiar with
[SMILES](https://www.wikiwand.com/en/Simplified_molecular-input_line-entry_system)
representations for describing chemical structures. SMARTS is a straightforward
extension of this notation.

#### Basic usage
Consider the following example defining the OPLS-AA atomtypes for a methyl group
carbon and its hydrogen atoms:
```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
  <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
 </AtomTypes>
</ForceField>
```

This `.xml` format is an extension of the [OpenMM force field format](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields)
The above example utilizes two additional `.xml` attributes supported by foyer:
`def` and `desc`. The atomtype that we are attempting to match is always the
__first__ token in the SMARTS string, in the above example, `[C;X4]` and `H`.
The `opls_135` (methyl group carbon) is defined by a SMARTS
string indicated a carbon with 4 bonds, a carbon neighbor and 3
hydrogen neighbors. The `opls_140` (alkane hydrogen) is defined simply as a
hydrogen atom bonded to a carbon with 4 bonds.


#### Overriding atomtypes
When multiple atomtype definitions can apply to a given atom, we must establish
precedence between those definitions. Many other atomtypers determine rule
precedence by placing more specific rules first and evaluate those in sequence,
breaking out of the loop as soon as a match is found.

While this approach works, it becomes more challenging to maintain the correct
ordering of rules as the number of atomtypes grows. Foyer iteratively runs all
rules on all atoms and each atom maintains a whitelist (rules that apply) and a
blacklist (rules that have been superceded by another rule). The set difference
between the white- and blacklists yields the correct atomtype if the force field
is implemented correctly.

We can add a rule to a blacklist using the `overrides` attribute in the `.xml`
file. To illustrate an example where overriding can be used consider the
following types describing alkenes and benzene:

```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_141" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)C" desc="alkene C (R2-C=)"/>
  <Type name="opls_142" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)H" desc="alkene C (RH-C=)"/>
  <Type name="opls_144" class="HC" element="H" mass="1.00800"  def="[H][C;X3]" desc="alkene H"/>
  <Type name="opls_145" class="CA" element="C" mass="12.01100" def="[C;X3;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1" overrides="opls_141,opls_142"/>
  <Type name="opls_146" class="HA" element="H" mass="1.00800"  def="[H][C;%opls_145]" overrides="opls_144" desc="benzene H"/>
 </AtomTypes>
</ForceField>
```

If we're atomtyping a benzene molecule, the carbon atoms will match the SMARTS
patterns for both `opls_142` and `opls_145`. Without the `overrides` attribute,
foyer will notify you that multiple atomtypes were found for each carbon.
Providing the `overrides` indicates that if the `opls_145` pattern matches, it
should supercede the specified rules.

#### Current Grammar Supported
We currently do not (yet) support all of [SMARTS' features](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html). [Here](https://github.com/mosdef-hub/foyer/issues/63) we're keeping track of which portions are supported.
