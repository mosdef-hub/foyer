### Usage Examples

Foyer supports atomtyping of both all-atom and coarse-grained molecular systems, and
also allows for separate force fields to be used to atom-type separate portions of a
molecular system.

#### Creating a box of ethane
Here we use mBuild to construct a box filled with ethane molecules and use Foyer to
atom-type the system, applying the OPLS force field, and save to run-able GROMACS 
files.
```python
import mbuild as mb
from mbuild.examples import Ethane

ethane_box = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
ethane_box.save('ethane-box.gro')
ethane_box.save('ethane-box.top', forcefield_name='oplsaa')
```
----

#### Creating a box of coarse-grained ethane
Again we will use mBuild to construct a box filled with ethane molecules.  However,
now we will model ethane using a united-atom description and apply the TraPPE force
field during atom-typing.  Note how particle names are prefixed with an underscore so
that Foyer knows these particles are non-atomistic.
```python
import mbuild as mb

ethane_UA = mb.Compound()
ch3_1 = mb.Particle(name='_CH3', pos=[0, 0, 0])
ch3_2 = mb.Particle(name='_CH3', pos=[0.15, 0, 0])
ethane_UA.add([ch3_1, ch3_2])
ethane_UA.add_bond((ch3_1, ch3_2))

ethane_UA_box = mb.fill_box(ethane_UA, 100, box=[2, 2, 2])
ethane_UA_box.save('ethane-UA-box.gro')
ethane_UA_box.save('ethane-UA-box.top', forcefield_name='trappe-ua')
```
----
#### Combining force fields
In some instances, the use of multiple force fields may be desired to describe a
molecular system.  For example, the user may want to use one force field for a
surface and another for a fluid in the same system.  Foyer supports this
functionality by allowing the user to separately atom-type parts of a system.  In
this example, we take a system featuring bulk united atom ethane above a silica 
surface and apply the OPLS force field to the surface and the TraPPE force field to
the ethane.  The two atomtyped Parmed structures are then combined using a simple
'\+' operator and can be saved to Gromacs files.
```python
from foyer import Forcefield
from foyer.tests.utils import get_fn
import mbuild as mb
from mbuild.examples import Ethane
from mbuild.lib.atoms import H
from mbuild.lib.bulk_materials import AmorphousSilica

interface = mb.SilicaInterface(bulk_silica=AmorphousSilica())
interface = mb.Monolayer(surface=interface, chains=H(), guest_port_name='up')

box = mb.Box(mins=[0, 0, max(interface.xyz[:,2])],
             maxs=interface.periodicity + [0, 0, 4]) 

ethane_box = mb.fill_box(compound=Ethane(), n_compounds=200, box=box)

opls = Forcefield(name='oplsaa')
opls_silica = Forcefield(forcefield_files=get_fn('opls-silica.xml'))
ethane_box = opls.apply(ethane_box)
interface = opls_silica.apply(interface)

system = interface + ethane_box

system.save('ethane-silica.gro')
system.save('ethane-silica.top')
```
