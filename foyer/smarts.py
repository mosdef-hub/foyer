from collections import Counter
from parmed.periodic_table import Element
from oset import oset as OrderedSet

SYMBOLS = ['(', ')', '*', '|']


class Smarts(object):
    """

    Attributes
    ----------
    string : str

    symbols : list of str

    identifiers : list of str

    tokens : list of str

    """
    def __init__(self, smarts_string, name="", identifiers=None, tokenize=True, overrides=None):
        self.string = smarts_string
        self.name = name

        if not identifiers:
            self.identifiers = []
        else:
            self.identifiers = identifiers

        if not overrides:
            self.overrides = OrderedSet()
        else:
            self.overrides = OrderedSet(overrides)

        self.symbols = Element + SYMBOLS + identifiers
        self.symbols.sort()
        self.symbols.reverse()

        self._n_neighbors = None

        if tokenize:
            self.tokens = self.tokenize()

    @property
    def anchor(self):
        # TODO: enforce that anchor is an element
        return self.tokens[0]

    @property
    def n_neighbors(self):
        if not self._n_neighbors:
            self._n_neighbors = len(self.neighbors)
        return self._n_neighbors

    @property
    def neighbors(self):
        """Return the number of neighbors of the anchor atom."""
        position = 1
        neighbors = []
        while True:
            branch, branch_token_cnt = self.next_branch(self.tokens[position:])
            if branch:
                neighbors.append(self._branch_head(branch))
                position += branch_token_cnt
            else:
                assert position == len(self.tokens)
                break
        return neighbors

    def next_branch(self, tokens):
        """Find the next branch in a list of tokens.

        Parameters
        ----------
        tokens : list of str
            The tokens to find the next branch in.

        Returns
        -------
        branch : list of str
            The next branch with parentheses stripped if necessary.
        len_branch : int
            The number of tokens in the branch.
        """
        if not tokens:
            return [], 0

        if tokens[0] != '(':
            return tokens, len(tokens)
        else:
            in_branch = 1
            branch = []
            for t in tokens[1:]:
                if t == '(':
                    in_branch += 1
                elif t == ')':
                    in_branch -= 1
                    if in_branch == 0:
                        return branch, len(branch) + 2
                else:
                    branch.append(t)
        raise SyntaxError("Unbalanced parentheses")

    def _branch_head(self, branch):
        # return the alternatives for the first element of the branch
        head = []
        for i,token in enumerate(branch):
            if token == '|':
                continue
            head.append(token)
            if len(branch) > i+1 and branch[i+1] != '|':
                break
        return head


    def tokenize(self):
        """Tokenize a SMARTS string.

        Returns
        -------
        tokens : list of str
        """
        tokens = []
        position = 0
        len_smarts = len(self.string)
        while position < len_smarts:
            # Try to find an identifier in curly braces.
            token, char_cnt = self.next_curlybraced(position)
            if token:
                tokens.append(token)
                position += char_cnt
                continue

            # Try to find an identifier in brackets.
            token, char_cnt = self.next_bracketed(position)
            if token:
                tokens.append(token)
                position += char_cnt
                continue

            # Try to get a symbol.
            token, char_cnt = self.next_symbol(position)
            if token:
                tokens.append(token)
                position += char_cnt
                continue

            raise SyntaxError("Error tokenizing {} at index {} ([...]{})".format(
                self.string, position, self.string[position:]))
        return tokens

    def next_symbol(self, position):
        """Find the next symbol (identifier, element name, parenthesis).

        Identifiers may be followed by an arbitrary number of digits.

        Returns
        -------
        symbol : str
            The next symbol or '' if nothing was found.
        len_characters : int
            The number of characters processed or 0 if nothing was found.
        """
        symbol = self.string[position:]
        for id in self.symbols:
            if symbol.startswith(id):
                digit_cnt = 0
                if id not in SYMBOLS:
                    # Look for a number following the ID.
                    for ch in symbol[len(id):]:
                        if ch.isdigit():
                            digit_cnt += 1
                        else:
                            break
                return symbol[0:len(id)+digit_cnt], len(id) + digit_cnt
        return '', 0

    def next_bracketed(self, position):
        """Find the next bracketed token.

        Returns
        -------
        token : str
            The token with brackets stripped or '' if nothing was found.
        len_characters : int
            The number of characters processed or 0 if nothing was found.
        """
        sub_string = self.string[position:]
        if sub_string.startswith('['):
            start_idx = 1
            end_idx = sub_string.find(']')
            if end_idx == -1:
                raise SyntaxError("Cannot find closing bracket in {} )".format(sub_string))
            token = sub_string[start_idx:end_idx]
            return token, end_idx + 1
        return '', 0

    def next_curlybraced(self, position):
        """Find the next token in curly braces.

        Returns
        -------
        token : str
            The token with brackets stripped or '' if nothing was found.
        len_characters : int
            The number of characters processed or 0 if nothing was found.
        """
        sub_string = self.string[position:]
        if sub_string.startswith('{'):
            start_idx = 1
            end_idx = sub_string.find('}')
            if end_idx == -1:
                raise SyntaxError("Cannot find closing curly brace in {} )".format(sub_string))
            token = sub_string[start_idx:end_idx]

            return token, end_idx + 1
        return '', 0


    def _match_neighbors(self, neighbor_label_alternatives, neighbor_atoms):

        if not neighbor_label_alternatives and not neighbor_atoms:
            return True

        label_alternatives = neighbor_label_alternatives[0]
        neighbor_label_alternatives = neighbor_label_alternatives[1:]
        for atom in neighbor_atoms:
            if '*' in label_alternatives or atom.element_name in label_alternatives or set(label_alternatives) & set(atom.whitelist):
                nonmatched_neighbor_atoms = neighbor_atoms
                nonmatched_neighbor_atoms.remove(atom)
                success = self._match_neighbors(neighbor_label_alternatives, nonmatched_neighbor_atoms)
                if success:
                    return True
        return False


    def match(self, atom):
        # atom is a parmed atom object

        # check atom type
        if self.anchor != atom.element_name:
            return False

        # check number of neighbors
        if self.n_neighbors != len(atom.bond_partners):
            return False

        if not self._match_neighbors(self.neighbors, atom.bond_partners):
            return False

        atom.whitelist.add(self.name)
        atom.blacklist |= self.overrides

        return True


if __name__ == '__main__':
    # s = Smarts('C(H)(H)(C)C', identifiers=['Cat'])
    # s = Smarts('CatCa12345C(H)', identifiers=['Cat'])
    #s = Smarts('C(OH)(H)(H)Opls123', identifiers=[])
    # s = Smarts('C(OH)(H)(H)123', identifiers=[])
    # s = Smarts('C(OH)(H)(H)[123]C456', identifiers=[])
    # s = Smarts('C(OH)(H)(H){ 123 }{_44}', identifiers=[])
    # s = Smarts('H(141)', identifiers=['141'])
    s = Smarts('H({141}|{122}CH)*', identifiers=['141'])
    print("Tokens: {}".format(s.tokens))
    print("Anchor: {}".format(s.anchor))
    print("Neighbors of the anchor: {}".format(s.neighbors))
