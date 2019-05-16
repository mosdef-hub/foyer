Paper Examples
~~~~~~~~~~~~~~

Contained below are the toy examples from the *Usage Examples* section of the `foyer paper <https://arxiv.org/pdf/1812.06779.pdf>`__. The source code selections are listed below on this page, there are `Jupyter
Notebooks <https://github.com/mosdef-hub/foyer/tree/master/docs/examples>`__
where you can try these examples yourself. Note that these examples are
meant to showcase the abilities of ``foyer`` through simple examples. If
the user would like to examine more in-depth examples using ``foyer``
with ``mBuild``, please refer to the `tutorial
repository <https://github.com/mosdef-hub/mosdef_tutorials>`__.

Below is *Listing 6* from the paper, a python script to fill a :math:`2x2x2 nm` 
box with 100 ethane molecules. The system is then atomtyped using the
OPLS-AA forcefield. There are two approaches to the same problem
detailed below in this listing, the first approach uses the
``forcefield_files`` function argument from
`mBuild <https://github.com/mosdef-hub/mbuild>`__ to atomptype the
system (using foyer under the hood). While the second approach creates a
``foyer`` ``Forcefield`` object, which then calls its ``apply``
function, operating on the ``mBuild`` ``Compound`` to return the
properly atomtyped structure. Note that in all instances when using
``foyer``, the chemical system of interest is converted into a
``ParmEd`` ``Structure``. Even the ``mBuild`` ``Compounds``, when
calling the ``save`` routine, are converted into a ``ParmEd``
``Structure`` before ``foyer`` can atomtype them. The object returned by
``foyer`` after the atomtypes have been applied are ``ParmEd``
``Structures``. This is subject to change in later iterations of
``foyer``.

Example 1
^^^^^^^^^

.. code:: python
    
    import mbuild as mb
    from mbuild.examples import Ethane
    from foyer.tests.utils import get_fn
    from foyer import Forcefield

    """ Applying a force field while saving from mBuild """
    # Create the chemical topology
    ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
    # Apply and save the topology
    ethane_fluid.save("ethane-box.top", forcefield_files=get_fn("oplsaa_alkane.xml"))
    ethane_fluid.save("ethane-box.gro")

    """ Applying a force field directly with foyer """
    # Create the chemical topology
    ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
    # Load the forcefield
    opls_alkane = Forcefield(forcefield_files=get_fn("oplsaa_alkane.xml"))
    # Apply the forcefield to atom-type
    ethane_fluid = opls_alkane.apply(ethane_fluid)
    # Save the atom-typed system
    ethane_fluid.save("ethane-box.top", overwrite=True)
    ethane_fluid.save("ethane-box.gro", overwrite=True)


---------------------------------------

Example 2
^^^^^^^^^

The other example listing from the text showcases the ability to create
two separate chemical topologies and applying different forcefield files
to each. The two parameterized systems that are generated are then
combined into a single ``ParmEd`` ``Structure`` and saved to disk.

.. code:: python

    from foyer import Forcefield
    from foyer.tests.utils import get_fn
    import mbuild as mb
    from mbuild.examples import Ethane
    from mbuild.lib.atoms import H
    from mbuild.lib.bulk_materials import AmorphousSilica

    # Create a silica substrate, capping surface oxygens with hydrogen
    silica=mb.recipes.SilicaInterface(bulk_silica=AmorphousSilica())
    silica_substrate=mb.recipes.Monolayer(surface=silica,chains=H(),guest_port_name="up")
    # Determine the box dimensions dictated by the silica substrate
    box=mb.Box(mins=[0, 0,max(silica.xyz[:,2])],maxs=silica.periodicity+ [0, 0, 4])
    # Fill the box with ethane
    ethane_fluid=mb.fill_box(compound=Ethane(),n_compounds=200,box=box)
    # Load the forcefields
    opls_silica=Forcefield(forcefield_files=get_fn("oplsaa_with_silica.xml"))
    opls_alkane=Forcefield(forcefield_files=get_fn("oplsaa_alkane.xml"))
    # Apply the forcefields
    silica_substrate=opls_silica.apply(silica_substrate)
    ethane_fluid=opls_alkane.apply(ethane_fluid)
    # Merge the two topologies
    system=silica_substrate+ethane_fluid
    # Save the atom-typed system
    system.save("ethane-silica.top")
    system.save("ethane-silica.gro")
