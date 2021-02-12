Installation
==============

Install with `conda <https://repo.anaconda.com/miniconda/>`_
-----------------------------------------------------
::

    $ conda install -c conda-forge foyer

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply install without specifying the channels::

    $ conda config --add channels conda-forge
    $ conda install foyer
Install with `pip <https://pypi.org/project/pip/>`_
---------------------------------------------------
::

    $ pip install foyer

Install an editable version from source
---------------------------------------
::

    $ git clone https://github.com/mosdef-hub/foyer
    $ cd foyer
    $ pip install -e .

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

