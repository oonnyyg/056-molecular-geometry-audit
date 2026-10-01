"""Input validation errors for the molecular geometry audit."""


class GeometryInputError(ValueError):
    """Raised when an input file fails validation.

    The CLI treats this as a user-facing input error: it prints the message
    and exits with a non-zero status *without* writing or overwriting any
    report file.
    """
