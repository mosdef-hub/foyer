from parmed.periodic_table import Element

SYMBOLS = ['(', ')']


class Smarts(object):
    """

    Attributes
    ----------
    string : str

    symbols : list of str

    identifiers : list of str

    tokens : list of str

    """
    def __init__(self, smarts_string, identifiers=None, tokenize=True):
        self.string = smarts_string

        if not identifiers:
            self.identifiers = []
        else:
            self.identifiers = identifiers

        self.symbols = Element + SYMBOLS + identifiers
        self.symbols.sort()
        self.symbols.reverse()

        if tokenize:
            self.tokens = self.tokenize()

    @property
    def anchor(self):
        return self.tokens[0]

    @property
    def neighbors(self):
        """Return the number of neighbors of the anchor atom."""
        position = 1
        neighbors = []
        while True:
            branch, branch_token_cnt = self.next_branch(self.tokens[position:])
            if branch:
                neighbors.append(branch[0])
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
                    branch += t
        raise SyntaxError("Unbalanced parentheses")

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



if __name__ == '__main__':
    # s = Smarts('C(H)(H)(C)C', identifiers=['Cat'])
    s = Smarts('CatCa12345C(H)', identifiers=['Cat'])
    #s = Smarts('C(OH)(H)(H)Opls123', identifiers=[])
    # s = Smarts('C(OH)(H)(H)123', identifiers=[])
    #s = Smarts('C(OH)(H)(H)[123]C456', identifiers=[])
    s = Smarts('C(OH)(H)(H){ 123 }{_44}', identifiers=[])
    print("Tokens: {}".format(s.tokens))
    print("Anchor: {}".format(s.anchor))
    print("Neighbors of the anchor: {}".format(s.neighbors))
