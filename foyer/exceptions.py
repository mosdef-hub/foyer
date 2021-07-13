"""Handle custom foyer exceptions."""


class FoyerError(Exception):
    """Base class for all non-trivial errors raised by Foyer."""


class FoyerWarning(Warning):
    """Base class for all non-trivial warnings raised by Foyer."""


class ValidationError(FoyerError):
    """Raised when validating .xml forcefield files."""

    def __init__(self, message, source, line):
        super(ValidationError, self).__init__(message)
        self.source = source
        self.line = line


class MultipleValidationError(FoyerError):
    """Used for grouping and raising multiple ValidationErrors of one type."""

    def __init__(self, validation_errors):
        self.validation_errors = validation_errors

    def __str__(self):
        """Represent multiple foyer ValidationErrors."""
        message = ["\n"]
        for err in self.validation_errors:
            message.append("\t" + str(err))
        return "\n".join(message)


class ValidationWarning(FoyerWarning):
    """Raised when validating .xml forcefield files."""

    pass


class MissingParametersError(FoyerError):
    """Error to be raised when parameters are missing from the forcefield."""


class MissingForceError(FoyerError):
    """Error to be raised when a particular force doesn't exist in the forcefield."""


class UnimplementedCombinationRuleError(FoyerError):
    """Error to be raised when a combination rule is requested but not implemented."""


def raise_collected(errors):
    """Return errors for all collected exceptions."""
    if len(errors) > 1:
        raise MultipleValidationError(errors)
    elif len(errors) == 1:
        raise errors[0]
