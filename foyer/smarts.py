import plyplus

from foyer.exceptions import FoyerError

GRAMMAR = (r"""
    start: string;

    // Rules
    @string: chain nonlastbranch* lastbranch?;
    @chain: atom chain | atom;
    @nonlastbranch: LPAR branch RPAR;
    @lastbranch: branch;
    branch: string;
    atom: (LBRACKET weak_and_expression RBRACKET | atom_symbol) atom_label?;
    atom_label: NUM;
    ?weak_and_expression: (weak_and_expression weak_and_symbol)? or_expression;
    ?or_expression: (or_expression or_symbol)? and_expression;
    ?and_expression: (and_expression and_symbol)? (atom_id | not_expression);
    not_expression: not_symbol atom_id;
    @and_symbol: AMP;
    @weak_and_symbol: SEMI;
    @or_symbol: COMMA;
    @not_symbol: EXCL;
    atom_id: atom_symbol
             | HASH atomic_num
             | DOLLAR LPAR matches_string RPAR
             | has_label
             | 'X' neighbor_count
             | 'r' ring_size
             | 'R' ring_count;
    atom_symbol: SYMBOL | STAR;
    atomic_num: NUM;
    matches_string: string ;
    has_label: LABEL ;
    neighbor_count: NUM;
    ring_size: NUM;
    ring_count: NUM;

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
    EXCL: '\!';
    // Tokens for chemical elements
    // Optional, custom, non-element underscore-prefixed symbols are pre-pended
    SYMBOL: '{optional}C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|H[eofgas]?|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|U|V|W|Xe|Yb?|Z[nr]';

""")


class SMARTS(object):
    """A wrapper class for parsing SMARTS grammar using plyplus.

    Provides functionality for injecting optional, custom, non-element symbols
    denoted by an underscore-prefix as additional tokens that the parser can
    recognize.

    """
    def __init__(self, optional_names=''):
        if optional_names:
            for n in optional_names:
                if not n.startswith('_'):
                    raise FoyerError('Non-element types must start with an underscore, you passed {}'.format(', '.join(optional_names)))

            optional_names = sorted(optional_names, reverse=True)
            self.grammar = GRAMMAR.format(optional='{}|'.format(
                '|'.join(optional_names)))

        else:
            self.grammar = GRAMMAR.format(optional='')
        self.PARSER = plyplus.Grammar(self.grammar)

    def parse(self, smarts_string):
        return self.PARSER.parse(smarts_string)
