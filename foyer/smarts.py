import plyplus

_grammar = ("""
    start: string;

    // Rules
    @string: chain nonlastbranch* lastbranch?;
    @chain: atom chain | atom;
    @nonlastbranch: LPAR branch RPAR;
    @lastbranch: branch;
    branch: string;
    atom: (LBRACKET or_expression RBRACKET | atom_symbol) atom_label?;
    atom_label: NUM;
    ?or_expression: (or_expression or_symbol)? and_expression;
    ?and_expression: (and_expression and_symbol)? atom_id;
    @and_symbol: SEMI | AMP;
    @or_symbol: COMMA;
    atom_id: atom_symbol | HASH atomic_num | DOLLAR LPAR matches_string RPAR | has_label | 'X' neighbor_count | 'R' ring_size;
    atom_symbol: SYMBOL | STAR;
    atomic_num: NUM;
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
    LABEL: '\%[A-Za-z_0-9]+';
    // Tokens for chemical elements
    // Optional, custom, non-element underscore-prefixed symbols are pre-pended
    SYMBOL: '{optional}C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|H[eofgas]?|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|U|V|W|Xe|Yb?|Z[nr]';

""")

class SMARTS(object):
    def __init__(self, optional_names=''):
        if optional_names:
            self.grammar = _grammar.format(optional='{}|'.format(
                    '|'.join(optional_names)))
        else:
            self.grammar = _grammar.format(optional='')
        self.PARSER = plyplus.Grammar(self.grammar)

    def parse(self, smarts_string):
        return self.PARSER.parse(smarts_string)
