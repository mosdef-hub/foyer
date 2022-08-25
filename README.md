### Foyer: A package for atom-typing as well as applying and disseminating forcefields

[![Gitter chat](https://badges.gitter.im/mosdef-hub/gitter.svg)](https://gitter.im/mosdef-hub/Lobby)
[![CI](https://github.com/mosdef-hub/foyer/actions/workflows/CI.yaml/badge.svg)](https://github.com/mosdef-hub/foyer/actions/workflows/CI.yaml)
[![Anaconda Badge](https://anaconda.org/conda-forge/foyer/badges/version.svg)](https://anaconda.org/conda-forge/foyer)
[![codecov](https://codecov.io/gh/mosdef-hub/foyer/branch/master/graph/badge.svg)](https://codecov.io/gh/mosdef-hub/foyer)
[![DOI](https://zenodo.org/badge/34077879.svg)](https://zenodo.org/badge/latestdoi/34077879)


## Overview
Foyer is an open-source Python tool for defining and applying force field atom-typing
rules in a format that is both human- and machine-readable.  It parametrizes chemical topologies,
generating, syntactically correct input files for various simulation engines. Foyer provides a framework for force field
dissemination, helping to eliminate ambiguity in atom-typing and improving reproducibility
(for more information, see [our paper](https://www.sciencedirect.com/science/article/pii/S0927025619303040) or its corresponding [pre-print](https://arxiv.org/pdf/1812.06779.pdf)).

#### Foyer within the MoSDeF Ecosystem
<p align="center">
  <img src="docs/images/mosdef_graphic_foyer.svg?raw=true" alt="Foyer within the MoSDeF Ecosystem" width="500" height="500"/>
</p>

Foyer defines force fields in an XML format, where SMARTS strings are used to define the chemical context
of a particular atom type and “overrides” are used to set rule precedence, rather than a rigid hierarchical scheme.
Foyer builds upon the [OpenMM .xml force field](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields)
file, annotated with SMARTS-based atomtypes, e.g.:

```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
  <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
 </AtomTypes>
</ForceField>
```

Foyer can apply the forcefield to arbitrary chemical topologies. We currently support:

* [OpenMM.Topology](https://github.com/openmm/openmm/blob/7.6.0/wrappers/python/openmm/app/topology.py#L70)
* [ParmEd.Structure](http://parmed.github.io/ParmEd/html/structure.html)
* [mBuild.Compound](http://mosdef-hub.github.io/mbuild/data_structures.html)
* [gmso.Topology](https://gmso.mosdef.org/en/stable/data_structures.html#gmso.Topology)
* [openff.tookit.topology.Topology](https://open-forcefield-toolkit.readthedocs.io/en/0.9.2/api/generated/openff.toolkit.topology.Topology.html#openff-toolkit-topology-topology)

Application of a force field can be as simple as:
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

The `Foyer` package is part of the [Molecular Simulation Design Framework (MoSDeF) project](http://mosdef.org/).
Libraries in the MoSDeF ecosystem are designed to provide utilities neccessary to streamline
a researcher's simulation workflow. When setting up simulation studies,
we also recommend users to follow the [TRUE](https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938)
(Transparent, Reproducible, Usable-by-others, and Extensible) standard, which is a set of common
practices meant to improve the reproducibility of computational simulation research.

## Getting started

#### Quick Start with Docker
To use `foyer` in a jupyter-notebook that runs from a docker container with all the dependencies installed use the following command:

```sh
$ docker pull mosdef/foyer:latest
$ docker run -it --name foyer -p 8888:8888 mosdef/foyer:latest\
  /opt/conda/envs/foyer-docker/bin/jupyter notebook --ip="*"
```

Alternatively, you can also start a Bourne shell directly:
```sh
$ docker run -it --name foyer mosdef/foyer:latest
```

To learn more about using `foyer` with docker, please refer to the documentation [here](https://foyer.mosdef.org/en/latest/docker.html) .


#### Getting started with SMARTS-based atom-typing
* [SMARTS-based atomtyping](docs/source/topic_guides/smarts.rst)
* [Supported SMARTS Grammar](https://github.com/mosdef-hub/foyer/issues/63)

#### Defining force fields:
* [Defining force field parameters](docs/source/topic_guides/parameter_definitions.rst)
* [Force field file validation](docs/source/reference/validation.rst)


#### Example foyer force field files:
Foyer currently includes a subset of the OPLS AA and TraPPE forcefields, currently part of the source distribution:
* https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields

Additional example force field XML files:
* https://github.com/chrisiacovella/OPLSaa_perfluoroalkanes
* https://github.com/mosdef-hub/forcefield_perfluoroethers
* https://github.com/summeraz/OPLSaa_alkylsilanes

Example template for disseminating force fields:
* https://github.com/mosdef-hub/forcefield_template


#### Using Foyer to perform atom typing:
* [Basic usage examples](docs/source/topic_guides/usage_examples.rst)
* [Detailed Jupyter notebook tutorials, including integration with mBuild](https://github.com/mosdef-hub/foyer_tutorials)
* [Jupyter notebook tutorials](https://github.com/mosdef-hub/foyer_tutorials/tree/master), from [our paper](https://arxiv.org/abs/1812.06779)

### Documentation:
* Documentation website: http://foyer.mosdef.org

### Installation instructions
* [Installation instructions](docs/source/getting_started/install.rst)

### Citing Foyer:
* If you use this package, please cite [our paper](https://www.sciencedirect.com/science/article/pii/S0927025619303040) published in [Computational Materials Science](https://www.journals.elsevier.com/computational-materials-science).
* This manuscript is also available in its pre-print form on [arxiv](https://arxiv.org/pdf/1812.06779.pdf)
* The paper and examples in this work were developed for tag [paper_COMMAT_2019](https://github.com/mosdef-hub/foyer/tree/paper_COMMAT_2019)


* Please also cite the github repository, https://github.com/mosdef-hub/foyer

#### [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

This material is based upon work supported by the National Science Foundation under grants NSF ACI-1047828 and NSF ACI-1535150. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
