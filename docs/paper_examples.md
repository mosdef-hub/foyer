### Paper Examples

Contained below are the toy examples from the *Usage Examples* section of the [`foyer` paper](https://arxiv.org/pdf/1812.06779.pdf).
The source code selections are listed below on this page, there are [Jupyter Notebooks](https://github.com/mosdef-hub/foyer/tree/master/docs/examples) where you can try these examples yourself.
Note that these examples are meant to showcase the abilities of `foyer` through rather simplistic examples.
If the user would like to examine more in-depth and complex examples, please refer to the [tutorial repository](https://github.com/mosdef-hub/mosdef_tutorials).


Below is *Listing 6* from the paper, a python script to fill a $$2x2x2 nm$$ box with 100 ethane molecules.
The system is then atomtyped using the OPLS-AA forcefield.
There are two approaches to the same problem detailed below in this listing, the first approach uses the `forcefield_files` function argument from [`mBuild`](https://github.com/mosdef-hub/mbuild) to atomptype the system (using foyer under the hood).
While the second approach creates a `foyer` `Forcefield` object, which then calls its `apply` function, operating on the `mBuild` `Compound` to return the properly atomtyped structure.

```python
import mbuild as mb
from mbuild.examples import Ethane
from foyer.test.utils import get_fn
from foyer import Forcefield

### Approach 1 ###
# Create the chemical topology
ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
# Apply and save the topology
ethane_fluid.save(’ethane-box.top’, forcefield_files=get_fn(’oplsaa_alkane.xml’))
ethane_fluid.save(’ethane-box.gro’)

### Approach 2 ###
# Create the chemical topology
ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
# Load the forcefield
opls_alkane = Forcefield(forcefield_files=get_fn(’oplsaa_alkane.xml’))
# Apply the forcefield to atom-type
ethane_fluid = opls_alkane.apply(ethane_fluid)
# Save the atom-typedsystem
ethane_fluid.save(’ethane-box.top’, overwrite=True)
ethane_fluid.save(’ethane-box.gro’, overwrite=True)
```

The other example listing from the text showcases the ability to create two separate chemical topologies and applying different forcefield files to each.
The two parameterized systems that are generated are then combined into a single `ParmEd` `Structure` and saved to disk.

```python
from foyer import Forcefield
from foyer.test.utils import get_fn
import mbuild as mb
from mbuild.examples import Ethane
from mbuild.lib.atoms import H
from mbuild.lib.bulk_materials import AmorphousSilica

# Create a silica substrate, capping surface oxygens with hydrogen
silica=mb.SilicaInterface(bulk_silica=AmorphousSilica())
silica_substrate=mb.Monolayer(surface=silica,chains=H(),guest_port_name=’up’)
# Determine the box dimensions dictated by the silica substrate
box=mb.Box(mins=[0, 0,max(silica.xyz[:,2])],maxs=silica.periodicity+ [0, 0, 4])
# Fill the box with ethane
ethane_fluid=mb.fill_box(compound=Ethane(),n_compounds=200,box=box)
# Load the forcefields
opls_silica=Forcefield(forcefield_files=get_fn(’opls-silica.xml’))
opls_alkane=Forcefield(forcefield_files=get_fn(’oplsaa_alkane.xml’))
# Apply the forcefields
silica_substrate=opls_silica.apply(silica_substrate)
ethane_fluid=opls_alkane.apply(ethane_fluid)
# Merge the two topologies
system=silica_substrate+ethane_fluid
# Save the atom-typed system
system.save(’ethane-silica.top’)
system.save(’ethane-silica.gro’)
```


