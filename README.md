### Foyer: A package for atom-typing as well as applying and disseminating forcefields

[![Linux Build Status](https://travis-ci.org/mosdef-hub/foyer.svg?branch=master)](https://travis-ci.org/mosdef-hub/foyer)
[![Windows Build status](https://ci.appveyor.com/api/projects/status/r6b2ny2hjo1t1ulb/branch/master?svg=true)](https://ci.appveyor.com/project/ctk3b/foyer/branch/master)
[![PyPI Version](https://badge.fury.io/py/foyer.svg)](https://pypi.python.org/pypi/foyer)
[![Anaconda Badge](https://anaconda.org/mosdef/foyer/badges/version.svg)](https://anaconda.org/mosdef/foyer)
[![Coverage Status](https://coveralls.io/repos/github/mosdef-hub/foyer/badge.svg?branch=master)](https://coveralls.io/github/mosdef-hub/foyer?branch=master)

Annotate an [OpenMM .xml force field](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields)
file with SMARTS-based atomtypes:

```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
  <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
 </AtomTypes>
</ForceField>
```

Apply the forcefield to arbitrary chemical topologies. We currently support:

* [OpenMM.Topology](http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.app.topology.Topology.html#)
* [ParmEd.Structure](http://parmed.github.io/ParmEd/html/structure.html)
* [mBuild.Compound](http://mosdef-hub.github.io/mbuild/data_structures.html)

```python
from foyer import Forcefield
import parmed as pmd

untyped_ethane = pmd.load_file('ethane.mol2', structure=True)
oplsaa = Forcefield(forcefield_files='oplsaa.xml')
ethane = oplsaa.apply(untyped_ethane)

# Save to any format supported by ParmEd
ethane.save('ethane.top')
ethane.save('ethane.gro')
```
### Getting started?
Check out our example template for disseminating force fields:
https://github.com/mosdef-hub/forcefield_template

### [Installation instructions](docs/installation.md)

### [SMARTS-based atomtyping](docs/smarts.md)

### [Force field validation](docs/validation.md)

### [Defining force field parameters](docs/parameter_definitions.md)

### [Usage examples](docs/usage_examples.md)


#### [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

This material is based upon work supported by the National Science Foundation under grants NSF ACI-1047828 and NSF ACI-1535150. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
