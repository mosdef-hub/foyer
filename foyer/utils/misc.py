
def validate_type(iter, type_):
    """Validate all the elements of the iterable are of a particular type"""
    for item in iter:
        if not isinstance(item, type_):
            raise TypeError(
                f"Expected {item} to be of type {type_.__name__} but got"
                f" {type(item).__name__} instead."
            )
