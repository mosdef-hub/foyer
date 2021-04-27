Installation
==============

For most users we recommend a conda installation:

.. code:: bash

    conda install -c conda-forge -c omnia foyer


If you wish to install from source, you can use the following commands:

.. code:: bash

    git clone https://github.com/mosdef-hub/foyer.git
    cd foyer
    conda env create -f environment.yml
    conda activate foyer
    pip install .

If you are using windows, you should use ``environment-win.yml`` rather than
``environment.yml``.


If you plan on contributing to the development of foyer, we recommend
you create an editable installation with all the required dependencies:

.. code:: bash

    git clone https://github.com/mosdef-hub/foyer.git
    cd foyer
    conda env create -f environment-dev.yml
    conda activate foyer-dev
    pip install -e .


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
