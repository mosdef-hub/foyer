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

