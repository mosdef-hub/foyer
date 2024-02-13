"""Methods that require external access to function."""

import requests


def get_ref(ref_url, headers):
    """Return bibtex reference for used atom-types."""
    bibtex_ref = requests.get(ref_url, headers=headers)
    if bibtex_ref.ok:
        return bibtex_ref
    else:
        return None
