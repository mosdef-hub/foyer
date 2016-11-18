import itertools
import plyplus

SMARTS_GRAMMAR = plyplus.Grammar("""
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
    ?atom_id: atom_symbol | HASH atomic_num | any_atom | DOLLAR LPAR matches_string RPAR | has_label | 'D' neighbor_count | 'R' ring_size;
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
    LABEL: '\%[a-z_0-9]+';
    // Tokens for chemical elements
    SYMBOL: 'C[laroudsemf]?|Os?|N[eaibdpos]?|S[icernbmg]?|P[drmtboau]?|H[eofgas]?|A[lrsgutcm]|B[eraik]?|Dy|E[urs]|F[erm]?|G[aed]|I[nr]?|Kr?|L[iaur]|M[gnodt]|R[buhenaf]|T[icebmalh]|U|V|W|Xe|Yb?|Z[nr]';

""")

def parse(smarts_string):
    return SMARTS_GRAMMAR.parse(smarts_string)

if __name__ == '__main__':
    ast = parse('O([H&D1])(H)')
    print(ast)
    assert ast.head == "start"
    assert ast.tail[0].head == "atom"
    assert ast.tail[0].tail[0].head == "atom_symbol"
    assert ast.tail[0].tail[0].head == "atom_symbol"
    assert str(ast.tail[0].tail[0].tail[0]) == "O"

    ast = parse('[#6&D2]1CCCCC1')
    # g = ast2graph(ast)
    # assert g.nodes[0].degree == 2


    print(ast)
