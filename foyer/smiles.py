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

class Smiles(object):
    def __init__(self, smiles_string, identifiers=None):
        self.s = smiles_string

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

    def next_symbol(self, s):
        for id in self.symbols:
            if s.startswith(id):
                # look for a number following the ID
                digit_cnt = 0
                for ch in s[len(id):]:
                    if ch.isdigit():
                        digit_cnt += 1
                    else:
                        break
                return s[0:len(id)+digit_cnt]
        return None

    def tokenize(self):
        self.tokens = []
        i=0
        s=self.s
        while i < len(s):
            token = self.next_symbol(s[i:])
            if token:
                self.tokens += [token]
                i += len(token)
                continue
            else:
                raise SyntaxError("Error tokenizing {} at index {}".format(self.s, i))

if __name__ == '__main__':
    # s = Smiles('C(H)(H)(C)C', identifiers=['Cat'])
    # s = Smiles('CatCa12345C(H)', identifiers=['Cat'])
    s = Smiles('C(OH)(H)(H)Opls123', identifiers=[])
    print("Tokens: {}".format(s.tokens))
    print("Anchor: {}".format(s.anchor))
    print("Neighbors of the anchor: {}".format(s.neighbors))
