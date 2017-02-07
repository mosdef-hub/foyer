import itertools as it

import parmed.periodic_table as pt
from oset import oset as OrderedSet


class Rule(object):
    def __init__(self, name, parser, smarts_string, overrides=None):
        self.name = name
        self.smarts_string = smarts_string
        self.ast = parser.parse(smarts_string)
        self.ast.calc_parents()
        if overrides:
            self.overrides = set(overrides)
        else:
            self.overrides = set()

    def __repr__(self):
        return 'Rule({},{},{})'.format(self.name, self.smarts_string, self.overrides)

    def start_atom_pattern(self):
        return self.ast.tail[0]

    def _neighbor_atom_patterns(self, atom_expr):
        # the parent is an 'atom', which may have siblings
        atom = atom_expr.parent()
        if atom.is_last_kid:
            # we have reached the end of a branch: no neighbors here
            return []
        else:
            if atom.next_kid.head == 'atom':
                # we only have one neighbor, which is an atom
                # let's return the atom_expr it contains
                return [atom.next_kid]
            else:
                # we may have multiple neighbors: one or more branches, possibly followed by an atom
                assert atom.next_kid.head == 'branch'
                current_neighbor = atom.next_kid
                neighbor_atom_exprs = []
                while current_neighbor:
                    if current_neighbor.head == 'branch':
                        # add the expression of the first atom of the branch to the list
                        neighbor_atom_exprs.append(current_neighbor.tail[0])
                    if current_neighbor.head == 'atom':
                        # add the expression of the atom to the list
                        neighbor_atom_exprs.append(current_neighbor)
                        # this is an atom after the last branch (or there was no branch), so there are no more neighbors
                        break
                    if current_neighbor.is_last_kid:
                        # no more neighbors
                        break

                    # there are more neighbors
                    current_neighbor = current_neighbor.next_kid

                return neighbor_atom_exprs

    def _split_pattern(self, pattern):
        expr = pattern.tail[0]
        if len(pattern.tail) == 2:
            label = pattern.tail[1].tail[0]
            label = [int(digit) for digit in label]
        else:
            label = []
        return expr, label

    def matches(self, atom):
        return self._matches(atom,
                             atom_pattern=self.start_atom_pattern(),
                             visited_atoms=OrderedSet(),
                             labeled_atoms=dict())

    def _matches(self, atom, atom_pattern, visited_atoms, labeled_atoms):
        # Extract label from SMARTS string if present.
        atom_expr, label = self._split_pattern(atom_pattern)

        for digit in label:
            if digit in labeled_atoms:
                if labeled_atoms[digit] not in atom.bond_partners:
                     return False

        if label:
            for digit in label:
                if digit not in labeled_atoms:
                    labeled_atoms[digit] = atom

        # check if atom matches atom_expr
        if self._atom_expr_matches(atom_expr, atom):
            visited_atoms.add(atom)

            # get all neighbors in rule
            neighbor_atom_patterns = self._neighbor_atom_patterns(atom_expr)
            # print(neighbor_atom_patterns)
            if not neighbor_atom_patterns:
                # no expressions given for neighbor atoms: it's a match
                return True

            # get all neighbors in graph
            neighbor_atoms = atom.bond_partners

            # compute all combinations of rule-atom to graph-atom pairings
            n_expr_neighbors = len(neighbor_atom_patterns)
            if not len(neighbor_atoms) >= n_expr_neighbors:
                self._pop_atom(atom, visited_atoms, labeled_atoms)
                return False

            unvisited_neighbors = set(neighbor_atoms) - visited_atoms
            possible_match_sets = [zip(x, neighbor_atom_patterns)
                                   for x in it.permutations(unvisited_neighbors,
                                                            n_expr_neighbors)]

            # for all possible matchings of neighbor atoms to neighbor expressions
            for possible_match_set in possible_match_sets:
                # for all pair in a match set we check if all can be satisfied
                for neighbor_atom, neighbor_atom_pattern in possible_match_set:
                    # Check recursively if the expressions match.
                    if not self._matches(neighbor_atom,
                                         atom_pattern=neighbor_atom_pattern,
                                         visited_atoms=visited_atoms,
                                         labeled_atoms=labeled_atoms):
                        break
                else:
                    # we get here if we did not break in the loop
                    return True

            # none of the matchings work
            self._pop_atom(atom, visited_atoms, labeled_atoms)
            return False

        else:
            # the current atom did not match
            return False

    def _pop_atom(self, atom, visited_atoms, labeled_atoms):
        visited_atoms.pop(atom)
        for dig, labeled_atom in list(labeled_atoms.items()):
            if labeled_atom == atom:
                del labeled_atoms[dig]

    def _atom_expr_matches(self, atom_expr, atom):
        if atom_expr.head == 'and_expression':
            return self._atom_expr_matches(atom_expr.tail[0], atom) and \
                   self._atom_expr_matches(atom_expr.tail[1], atom)
        elif atom_expr.head == 'or_expression':
            return self._atom_expr_matches(atom_expr.tail[0], atom) or \
                   self._atom_expr_matches(atom_expr.tail[1], atom)
        elif atom_expr.head == 'atom_id':
            return self._atom_id_matches(atom_expr.tail[0], atom)
        elif atom_expr.head == 'atom_symbol':
            return self._atom_id_matches(atom_expr, atom)
        else:
            raise TypeError('Expected and_expression, or_expression,'
                            ' or atom_id, got {}'.format(atom_expr.head))

    def _atom_id_matches(self, atom_id, atom):
        if atom_id.head == 'atomic_num':
            return atom.element._atomic_number == int(atom_id.tail[0])
        elif atom_id.head == 'atom_symbol':
            if str(atom_id.tail[0]) == '*':
                return True
            elif str(atom_id.tail[0]).startswith('_'):
                return atom.element.name == str(atom_id.tail[0])
            else:
                return atom.element._atomic_number == pt.AtomicNum[str(atom_id.tail[0])]
        elif atom_id.head == 'has_label':
            label = atom_id.tail[0][1:] # cut the % sign from the beginning
            return label in (rule_name for rule_name in atom.whitelist)
        elif atom_id.head == 'neighbor_count':
            return len(atom.bond_partners) == int(atom_id.tail[0])
        elif atom_id.head == 'ring_size':
            cycle_len = int(atom_id.tail[0])
            for cycle in atom.cycles:
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_id.head == 'matches_string':
            raise NotImplementedError('matches_string feature is not yet implemented')
