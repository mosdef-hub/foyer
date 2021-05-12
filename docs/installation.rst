Installation
==============

Install from conda:

.. code:: bash

    conda install -c conda-forge -c omnia -c mosdef foyer

Install from pip:

.. code:: bash

    pip install foyer

Install an editable version from source:

.. code:: bash

    git clone https://github.com/mosdef-hub/foyer.git
    cd foyer
    pip install -e .


Install pre-commit
------------------

We use [pre-commit](https://pre-commit.com/) to automatically handle our code formatting and this package is included in the dev environment.
With the ``foyer-dev`` conda environment active, pre-commit can be installed locally as a git hook by running::

    $ pre-commit install

And (optional) all files can be checked by running::

    $ pre-commit run --all-files


Supported Python Versions
-------------------------

Python 3.6 and 3.7 are officially supported, including testing during
development and packaging. Other Python versions, such as 3.8 and 3.5 and
older, may successfully build and function but no guarantee is made.

Testing your installation
-------------------------

foyer uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

Building the documentation
--------------------------

foyer uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::

    $ conda env create -f docs-env.yml
    $ conda activate foyer-docs
    $ make html
