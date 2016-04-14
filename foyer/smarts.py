ElementNames = [ 'EP',
            'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne','Na','Mg',
            'Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca','Sc','Ti','V' ,'Cr',
            'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
            'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
            'In','Sn','Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
            'Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',
            'At','Rn','Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm',
            'Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs',
            'Mt','Ds','Rg','Cn','Uut','Uuq','Uup','Uuh','Uus','Uuo' ]

Symbols = ['(', ')']

class Smarts(object):
    def __init__(self, smarts_string, identifiers=None):
        self.s = smarts_string

        if not identifiers:
            self.identifiers = []
        else:
            self.identifiers = identifiers

        self.symbols = ElementNames + Symbols + identifiers
        self.symbols.sort()
        self.symbols.reverse()

        self.tokenize()

    @property
    def anchor(self):
        return self.tokens[0]

    def next_branch(self, tokens):
        """Given a list of tokens, find the next branch.

        Returns a tuple of the next branch (parentheses stripped if necessary), and the number of tokens in that branch
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

    @property
    def neighbors(self):
        """Return the number of neighbors of the anchor atom."""
        i = 1
        neighbors = []
        while True:
            branch, branch_token_cnt = self.next_branch(self.tokens[i:])
            if branch:
                neighbors += [ branch[0] ]
                i += branch_token_cnt
            else:
                assert i == len(self.tokens)
                break
        return neighbors

    def next_symbol(self, i):
        """Find the next symbol (identifier, element name, parenthesis).

        Identifiers may be followed by an arbitrary number of digits.

        Returns a tuple of the string found, and the number of characters processed, or (None, None) if nothing was
        found"""
        s = self.s[i:]
        for id in self.symbols:
            if s.startswith(id):
                digit_cnt = 0
                if id not in Symbols:
                    # look for a number following the ID
                    for ch in s[len(id):]:
                        if ch.isdigit():
                            digit_cnt += 1
                        else:
                            break
                return s[0:len(id)+digit_cnt], len(id) + digit_cnt
        return None, None

    def next_bracketed(self, i):
        """Find the bracketed string.

        Returns a tuple of the string found (brackets stripped), and the number of characters processed, or
        (None, None) if nothing was found"""
        s = self.s[i:]
        if s.startswith('['):
            start_idx = 1
            end_idx = s.rfind(']')
            if end_idx == -1:
                raise SyntaxError("Cannot find closing bracket in {} )".format(s))
            token = s[start_idx:end_idx]
            return token, end_idx + 1
        return None, None

    def tokenize(self):
        """Tokenize a smarts string.

        Returns a list of tokens."""
        self.tokens = []
        i=0
        s=self.s
        while i < len(s):
            # try to find an identifier in brackets
            token, char_cnt = self.next_bracketed(i)
            if token:
                self.tokens += [token]
                i += char_cnt
                continue

            # try to get a symbol
            token, char_cnt = self.next_symbol(i)
            if token:
                self.tokens += [token]
                i += char_cnt
                continue

            # couldn't find anything syntactically correct
            raise SyntaxError("Error tokenizing {} at index {} ([...]{})".format(self.s, i, self.s[i:]))

if __name__ == '__main__':
    # s = Smiles('C(H)(H)(C)C', identifiers=['Cat'])
    # s = Smiles('CatCa12345C(H)', identifiers=['Cat'])
    #s = Smiles('C(OH)(H)(H)Opls123', identifiers=[])
    # s = Smarts('C(OH)(H)(H)123', identifiers=[])
    s = Smarts('C(OH)(H)(H)[123]C456', identifiers=[])
    print("Tokens: {}".format(s.tokens))
    print("Anchor: {}".format(s.anchor))
    print("Neighbors of the anchor: {}".format(s.neighbors))
