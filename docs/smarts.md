### SMARTS-based atomtyping

Foyer allows users to describe atomtypes using a modified version of 
[SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html).

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


