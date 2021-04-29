Forcefield Class
----------------

The primary data structure in foyer is the ``Forcefield`` class, which inherits
from the OpenMM class of the same name. The primary operation on this class is
the ``.apply()`` function, which takes a chemical topology and returns a
parametrized ``ParmEd`` ``Structure``. The user may pass some otions to this
function based on a particular use case.

.. autoclass:: foyer.forcefield.Forcefield
    :members:

