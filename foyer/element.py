import simtk.openmm.app.element as elem

class Element(elem.Element):
    """An Element represents a chemical element.

    The simtk.openmm.app.element module contains objects for all the standard chemical elements,
    such as element.hydrogen or element.carbon.  You can also call the static method Element.getBySymbol() to
    look up the Element with a particular chemical symbol.

    Element objects should be considered immutable
    """
    def __init__(self, number, name, symbol, mass):
        """Create a new element

        Parameters
        ----------
        number : int
            The atomic number of the element
        name : string
            The name of the element
        symbol : string
            The chemical symbol of the element
        mass : float
            The atomic mass of the element
        """
        ## The atomic number of the element
        self._atomic_number = number
        ## The name of the element
        self._name = name
        ## The chemical symbol of the element
        self._symbol = symbol
        ## The atomic mass of the element
        self._mass = mass
        # Index this element in a global table
        s = symbol.strip().upper()
        ## If we add a new element, we need to re-hash elements by mass
        Element._elements_by_mass = None

        if s in Element._elements_by_symbol:
            raise ValueError('Duplicate element symbol %s' % s)
