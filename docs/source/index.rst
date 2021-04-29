Foyer
======
*A package for atom-typing as well as applying and disseminating force fields*

|License|
|Citing|
|Anaconda|
|CodeCov|
|Azure|

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
.. |Citing| image:: https://img.shields.io/badge/cite-foyer-green
   :target: reference/citing.html
.. |Anaconda| image:: https://anaconda.org/conda-forge/foyer/badges/version.svg
   :target: https://anaconda.org/conda-forge/foyer
.. |Codecov| image:: https://codecov.io/gh/mosdef-hub/foyer/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/mosdef-hub/foyer
.. |Azure| image:: https://dev.azure.com/mosdef/mosdef/_apis/build/status/mosdef-hub.foyer?branchName=master
   :target: https://dev.azure.com/mosdef/mosdef/_build/latest?definitionId=2&branchName=master


Overview
~~~~~~~~

**Foyer** is an open-source Python tool that provides a framework for the application
and dissemination of classical molecular modeling force fields. Importantly,
it enables users to define and apply atom-typing rules in a format that is
simultaneously human- and machine-readable. A primary goal of **foyer**
is to eliminate ambiguity in the atom-typing and force field application steps of
molecular simulations in order to improve reproducibility. Foyer force fields are
defined in an XML format derived from the OpenMM XML.
`SMARTS strings <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_
are used to define the chemical context of each atom type and "overrides" are used
to define clear precedence of different atom types. **Foyer** is designed to
be compatible with the other tools in the
`Molecular Simulation Design Framework (MoSDeF) ecosystem <https://mosdef.org>`_.


Resources
~~~~~~~~~

* :doc:`Installation guide <getting_started/install>`: Instructions for installing foyer.
* :doc:`Quickstart <getting_started/quickstart>`: A brief introduction to foyer.
* `MoSDeF  <https://mosdef.org>`_: Learn more about the **Mo**\ lecular **S**\ imulation **De**\ sign **F**\ ramework.
* `Foyer paper <https://www.sciencedirect.com/science/article/pii/S0927025619303040>`_: The journal article describing foyer.
* `GitHub repository <https://github.com/mosdef-hub/foyer>`_: Download the source code or contribute to the development of foyer.
* `Issue Tracker <https://github.com/mosdef-hub/foyer/issues>`_: Report issues and request features.


Citation
~~~~~~~~

If you use foyer in your research, please cite the
`foyer paper <https://www.sciencedirect.com/science/article/pii/S0927025619303040>`_.
See :doc:`here <reference/citing>` for details.


Installation
~~~~~~~~~~~~

Complete installation instructions are :doc:`here <getting_started/install>`.
A conda installation is available:

.. code-block:: bash

    conda create --name foyer foyer -c conda-forge

Example
~~~~~~~

Annotate an `OpenMM .xml force
field <http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields>`__
file with SMARTS-based atomtypes:

.. code:: xml

    <ForceField>
     <AtomTypes>
      <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
      <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
     </AtomTypes>
    </ForceField>

Apply the forcefield to arbitrary chemical topologies. We currently
support:

-  `OpenMM.Topology <http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.app.topology.Topology.html#>`__
-  `ParmEd.Structure <http://parmed.github.io/ParmEd/html/structure.html>`__
-  `mBuild.Compound <http://mosdef-hub.github.io/mbuild/data_structures.html>`__

.. code:: python

    from foyer import Forcefield
    import parmed as pmd

    untyped_ethane = pmd.load_file('ethane.mol2', structure=True)
    oplsaa = Forcefield(forcefield_files='oplsaa.xml')
    ethane = oplsaa.apply(untyped_ethane)

    # Save to any format supported by ParmEd
    ethane.save('ethane.top')
    ethane.save('ethane.gro')

Getting started?
~~~~~~~~~~~~~~~~

Check out our example template for disseminating force fields:
https://github.com/mosdef-hub/forcefield_template


Credits
~~~~~~~

This material is based upon work supported by the National Science
Foundation under grants NSF ACI-1047828 and NSF ACI-1535150. Any
opinions, findings, and conclusions or recommendations expressed in this
material are those of the author(s) and do not necessarily reflect the
views of the National Science Foundation.

Table of Contents
~~~~~~~~~~~~~~~~~

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    getting_started/install
    getting_started/docker
    getting_started/quickstart

.. toctree::
    :maxdepth: 2
    :caption: Topic Guides

    topic_guides/smarts
    topic_guides/parameter_definitions
    topic_guides/ffapply
    topic_guides/usage_examples
    topic_guides/paper_examples

.. toctree::
    :maxdepth: 2
    :caption: Reference

    reference/units
    reference/validation
    reference/ffclass
    reference/citing
    reference/license
