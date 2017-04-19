class FoyerError(Exception):
    """Base class for all non-trivial errors raised by Foyer """


class FoyerWarning(Warning):
    """Base class for all non-trivial warnings raised by Foyer """


class ValidationError(FoyerError):
    """Raised when validating .xml forcefield files """
    def __init__(self, message, source, line):
        super(ValidationError, self).__init__(message)
        self.source = source
        self.line = line


class ValidationWarning(FoyerWarning):
    """Raised when validating .xml forcefield files """
    pass

