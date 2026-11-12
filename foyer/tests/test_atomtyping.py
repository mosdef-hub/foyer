import parmed as pmd
import pytest

from foyer import Forcefield
from foyer.exceptions import FoyerError
from foyer.tests.base_test import BaseTest
from foyer.tests.utils import get_fn


class TestRunAtomTyping(BaseTest):
    @pytest.fixture(scope="session")
    def missing_overrides_ff(self):
        return Forcefield(get_fn("missing_overrides.xml"))

    def test_missing_overrides(self, opls_validation_benzene, missing_overrides_ff):
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(opls_validation_benzene)

    def test_missing_definition(self, missing_overrides_ff):
        structure = pmd.load_file(get_fn("silly_chemistry.mol2"), structure=True)
        with pytest.raises(FoyerError):
            missing_overrides_ff.apply(structure)

    @pytest.mark.parametrize(
        "atomic_num, symbol", [(6, "Boo"), (200, "C"), (200, "Boo")]
    )
    def test_element_not_found(self, oplsaa, atomic_num, symbol):
        from foyer.topology_graph import TopologyGraph

        # Create a TopologyGraph object for a methane molecule.
        # The central C has bad element info, which will trigger an error
        # during the atomtyping step
        top_graph = TopologyGraph()
        top_graph.add_atom(name="C", index=0, atomic_number=atomic_num, symbol=symbol)
        for i in range(1, 5):
            top_graph.add_atom(name="H", index=i, atomic_number=1, symbol="H")
            top_graph.add_bond(0, i)

        with pytest.raises(FoyerError):
            oplsaa.run_atomtyping(top_graph, use_residue_map=False)


class TestIterateRules(BaseTest):
    @staticmethod
    def _chain(n):
        """Build a bare carbon chain (no hydrogens) as a TopologyGraph."""
        from foyer.topology_graph import TopologyGraph

        top_graph = TopologyGraph()
        for i in range(n):
            top_graph.add_atom(name="C", index=i, atomic_number=6, symbol="C")
        for i in range(n - 1):
            top_graph.add_bond(i, i + 1)
        return top_graph

    @staticmethod
    def _rules(typemap):
        from foyer.smarts import SMARTS
        from foyer.smarts_graph import SMARTSGraph

        parser = SMARTS()
        # Ordered so each %label rule runs before the rule it depends on,
        # forcing typing to take several passes to converge.
        smarts = {
            "inner": "[C;X2][C;%next]",
            "next": "[C;X2][C;%end]",
            "end": "[C;X1]",
        }
        return {
            name: SMARTSGraph(
                s, parser=parser, name=name, overrides=set(), typemap=typemap
            )
            for name, s in smarts.items()
        }

    def test_multi_pass_convergence(self, monkeypatch):
        from foyer.atomtyper import _iterate_rules
        from foyer.smarts_graph import SMARTSGraph

        top_graph = self._chain(4)
        typemap = {
            i: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for i in top_graph.atoms(data=False)
        }
        rules = self._rules(typemap)

        calls = {name: 0 for name in rules}
        original = SMARTSGraph.find_matches

        def counting_find_matches(self, topology_graph, typemap):
            calls[self.name] += 1
            return original(self, topology_graph, typemap)

        monkeypatch.setattr(SMARTSGraph, "find_matches", counting_find_matches)

        _iterate_rules(rules, top_graph, typemap, max_iter=10)

        assert typemap[0]["whitelist"] == {"end"}
        assert typemap[3]["whitelist"] == {"end"}
        assert typemap[1]["whitelist"] == {"next", "inner"}
        assert typemap[2]["whitelist"] == {"next", "inner"}

        # Rules without %label run only on the first pass; %label rules run
        # on every pass: pass 1 finds "end", pass 2 "next", pass 3 "inner",
        # pass 4 finds nothing and stops.
        assert calls["end"] == 1
        assert calls["next"] == 4
        assert calls["inner"] == 4

    def test_max_iter_warning(self, caplog):
        from foyer.atomtyper import _iterate_rules

        top_graph = self._chain(4)
        typemap = {
            i: {"whitelist": set(), "blacklist": set(), "atomtype": None}
            for i in top_graph.atoms(data=False)
        }
        rules = self._rules(typemap)
        _iterate_rules(rules, top_graph, typemap, max_iter=2)

        assert "Reached maximum iterations" in caplog.text
        assert typemap[1]["whitelist"] == {"next"}
