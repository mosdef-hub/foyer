Validation of force field files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Foyer performs several validation steps to help prevent malformed force
field files and SMARTS strings from making it into production code. Our
goal is to provide human readable error messages that users who may not
be intimately familiar with XML or our SMARTS parsing grammar can easily
act upon.

However, if you receive any unclear error messages or warnings we
strongly encourage you to `submit an issue <https://github.com/mosdef-hub/foyer/issues/new>`_
detailing the error message you received and, if possible, attach a minimal
example of the force field file that created the problem.

XML schema
^^^^^^^^^^

As a first line of defense, any force field files loaded by foyer is
validated by this `XML schema
definition <../foyer/forcefields/ff.xsd>`__. Here we enforce which
elements (e.g. ``HarmonicBondForce``) are valid and how their attributes
should be formatted. Additionally, the schema ensures that atomtypes are
not 1) defined more than once and that 2) atomtypes referenced in other
sections are actually defined in the ``<AtomTypes>`` element.

SMARTS validation
^^^^^^^^^^^^^^^^^

All SMARTS strings used to define atomtypes are parsed. Parsing errors
are captured and re-raised with error messages that allow you to pin
point the location of the problem in the XML file and within the SMARTS
string. Wherever possible, we attempt to provide helpful hints and we
welcome any contributions that help improve the clarity of our error
messages.

Additionally, we ensure that any atomtypes referenced using the
``%type`` or ``overrides`` `syntax <smarts.html>`__ are actually defined
in the ``<AtomTypes>`` element.
