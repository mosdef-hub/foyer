import parmed as pmd
import pytest
from pkg_resources import resource_filename

from foyer import Forcefield
from foyer.exceptions import FoyerError
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn


class TestRunAtomTyping(BaseTest):
    @pytest.fixture(scope="session")
    def missing_overrides_ff(self):
        return Forcefield(get_fn("missing_overrides.xml"))

    def test_missing_overrides(
        self, opls_validation_benzene, missing_overrides_ff
    ):
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(opls_validation_benzene)

    def test_missing_definition(self, missing_overrides_ff):
        structure = pmd.load_file(
            get_fn("silly_chemistry.mol2"), structure=True
        )
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(structure)

    def test_element_not_found(self, oplsaa):
        from foyer.topology_graph import TopologyGraph
        #Create a TopologyGraph object for a methane molecule.
        #The central C has bad element info, which will trigger an error
        #during the atomtyping step
        top_graph = TopologyGraph()
        top_graph.add_atom(name="C",
                           index=0,
                           atomic_number=200,
                           element='Boo')
        for i in range(1, 5):
            top_graph.add_atom(name='H',
                               index=i,
                               atomic_number=1,
                               element='H')
            top_graph.add_bond(0, i)

        with pytest.raises(FoyerError):
            oplsaa.run_atomtyping(top_graph, use_residue_map=False)
