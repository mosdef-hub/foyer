
import plyplus
import parmed.periodic_table as pt


smarts_grammar = plyplus.Grammar("""
    start: string;

    // Rules
    @string: chain nonlastbranch* lastbranch?;
    @nonlastbranch: LPAR branch RPAR;
    @lastbranch: branch;
    branch: string;
    @chain: atom chain | atom;
    atom: LBRACKET or_expression RBRACKET atom_label?;
    atom_label: NUM;
    ?or_expression: (or_expression or_symbol)? and_expression;
    ?and_expression: (and_expression and_symbol)? atom_id;
    @and_symbol: SEMI | AMP;
    @or_symbol: COMMA;
    //atom_id: atom_symbol | HASH atomic_num | any_atom | DOLLAR LPAR matches_string RPAR | has_label | 'D' neighbor_count | 'R' ring_size;
    atom_id: atom_symbol | HASH atomic_num | any_atom | has_label | 'D' neighbor_count | 'R' ring_size;
    atom_symbol: SYMBOL;
    atomic_num: NUM;
    any_atom: STAR;
    matches_string: string ;
    has_label: LABEL ;
    neighbor_count: NUM;
    ring_size: NUM;

    // Tokens
    HASH: '\#';
    LBRACKET: '\[';
    RBRACKET: '\]';
    LPAR: '\(';
    RPAR: '\)';
    COMMA: '\,';
    SEMI: '\;';
    AMP: '\&';
    STAR: '\*';
    DOLLAR: '\$';
    NUM: '[\d]+';
    LABEL: '\%[a-z_]+([0-9][a-z_]?)*' ;
    // Tokens for chemical elements
    SYMBOL: 'C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|H[eofgas]?|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|U|V|W|Xe|Yb?|Z[nr]';

""")


def parse(expr):
    tree = smarts_grammar.parse(expr)
    assert tree.tail
    return tree

def show_result(s):
    t = parse(s)
    print("parsing: {}\nresult :{}\n".format(s,t))


def find_atomtypes(atoms, rules):
    for atom in atoms:
        atom.whitelist = set()
        atom.blacklist = set()

    found_something = True
    while(found_something):
        for atom in atoms:
            for rule in rules:
                if rule not in atom.whitelist or rule not in atom.blacklist:
                    if rule.matches(atom):
                        atom.whitelist.add(rule)
                        atom.blacklist |= rule.overrides
                        found_something = True
        found_something = False


class Rule(object):
    def __init__(self, name, smarts_string, overrides=None):
        self.name = name
        self.smarts_string = smarts_string
        self.ast = smarts_grammar.parse(smarts_string)
        if overrides:
            self.overrides = set(overrides)
        else:
            self.overrides = set()

    def __repr__(self):
        return 'Rule({},{},{})'.format(self.name, self.smarts_string, self.overrides)

    def start_atom_expr(self):
        assert self.ast.tail[0].head == 'atom'
        assert self.ast.tail[0].tail
        assert len(self.ast.tail[0].tail) == 1
        return self.ast.tail[0].tail[0]

    def matches(self, atom):
        return self._atom_expr_matches(self.start_atom_expr(), atom)

    def _atom_expr_matches(self, atom_expr, atom):
        if atom_expr.head == 'and_expression':
            return self._atom_expr_matches(atom_expr.tail[0], atom) and self._atom_expr_matches(atom_expr.tail[1], atom)
        elif atom_expr.head == 'or_expression':
            return self._atom_expr_matches(atom_expr.tail[0], atom) or self._atom_expr_matches(atom_expr.tail[1], atom)
        elif atom_expr.head == 'atom_id':
            return self._atom_id_matches(atom_expr.tail[0], atom)
        else:
            raise TypeError('Expected and_expression, or_expression, or atom_id, got {}'.format(atom_expr.head))

    def _atom_id_matches(self, atom_id, atom):
        if atom_id.head == 'any_atom':
            return True
        elif atom_id.head == 'atomic_num':
            return atom.atomic_number == int(atom_id.tail[0])
        elif atom_id.head == 'atom_symbol':
            return atom.atomic_number == pt.AtomicNum[str(atom_id.tail[0])]
        elif atom_id.head == 'has_label':
            return atom_id.tail[0] in atom.whitelist
        elif atom_id.head == 'neighbor_count':
            return len(atom.bond_partners) == int(atom_id.tail[0])
        elif atom_id.head == 'ring_size':
            raise NotImplementedError('ring_size feature is not yet implemented')


if __name__ == '__main__':

    # TODO: set up test cases to catch precedence issues
    # x = '[#1&#2,#3]'
    # show_result(x)
    # y = '[#1,#2&#3]'
    # show_result(y)

    from topos import Topos
    topos = Topos()

    ethanol = topos.load_topo("64-17-5")

    from smarts_opls import opls_rules

    find_atomtypes(ethanol.atoms, opls_rules)

    for atom in ethanol.atoms:
        print("{}: {}".format(atom, atom.whitelist - atom.blacklist))