import requests
import networkx as nx
from parmed import Structure


def get_ref(ref_url, headers):
    bibtex_ref = requests.get(ref_url, headers=headers)
    if bibtex_ref.ok:
        return bibtex_ref
    else:
        return None
