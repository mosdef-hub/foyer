
import re, collections
import plyplus

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
    atom_id: atom_symbol | HASH atomic_num | any_atom | DOLLAR LPAR matches_string RPAR | has_label | 'D' neighbor_count | 'R' ring_size;
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
    return tree

def show_result(s):
    t = parse(s)
    print("parsing: {}\nresult :{}\n".format(s,t))



if __name__ == '__main__':
    # """O TIP3P Water """
    opls_111 = '[#12&D2]([#1])[#1]'
    show_result(opls_111)

    # """H TIP3P Water """
    opls_112 = '[#1][%opls_111]'
    show_result(opls_112)

    # """alkane CH3 """
    opls_135 = '[#12&D4]([#12])([#1])([#1])([#1])'
    show_result(opls_135)

    # """alkane CH2 """
    opls_136 = '[#12&D4]([#12])([#12])([#1])([#1])'
    show_result(opls_136)

    # """alkane CH """
    opls_137 = '[#12&D4]([#12])([#12])([#12])([#1])'
    show_result(opls_137)

    # """alkane CH4 """
    opls_138 = '[#12&D4]([#1])([#1])([#1])([#1])'
    show_result(opls_138)

    """alkane C """
    opls_139 = '[#12&D4]([#12])([#12])([#12])([#12])'
    show_result(opls_139)

    # TODO: be specific about alkane C (not just simply C)
    # """alkane H """
    opls_140 = '[#1][#12]'
    show_result(opls_140)

    # NOTE: this catches benzene C as well
    # """alkene C (R2-C=) """
    opls_141 = '[#12&D3]([#12])([#12])([#12])'
    show_result(opls_141)

    # NOTE: this catches benzene C as well
    # """alkene C (RH-C=) """
    opls_142 = '[#12&D3]([#12])([#12])([#1])'
    show_result(opls_142)

    # NOTE: this catches benzene C as well
    # """alkene C (H2-C=) """
    opls_143 = '[#12&D3]([#12])([#1])([#1])'
    show_result(opls_143)

    # """alkene H (H-C=) """
    # NOTE: We make sure that the carbon is an alkene carbon.
    opls_144 = '[#1]([#12&%opls_141,#12&%opls_142,#12&%opls_143])'
    show_result(opls_144)

    # TODO: set up test cases to catch precedence issues
    # x = '[#1&#2,#3]'
    # show_result(x)
    # y = '[#1,#2&#3]'
    # show_result(y)

    # TODO: @Blacklist([141, 142])
    # """Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl """
    opls_145 = '[#12&D3&R6]([#12&R6])([#12&R6])'
    show_result(opls_145)

    # TODO: @Blacklist([145])
    # """Biphenyl C1 """
    opls_145B = '[#12&D3]([#12&%opls_145])([#12&%opls_145])([#12&%opls_145])'
    show_result(opls_145B)


    # TODO: @Blacklist([140, 144])
    # """Benzene H - 12 site. """
    opls_146 = '[#1]([#12&%opls_145])'
    show_result(opls_146)

    # TODO: @Blacklist(135)
    # """C: CH3, toluene """
    opls_148 = '[#12&D4]([#12&%opls_145])([#1])([#1])([#1])'
    show_result(opls_148)

    # TODO: @Blacklist(136)
    # """C: CH2, ethyl benzene """
    opls_149 = '[#12&D4]([#12&%opls_145])([#12])([#1])([#1])'
    show_result(opls_149)


    # TODO: make sure it doesn't catch water
    # """all-atom O: mono alcohols """
    opls_154 = '[#8&D2]([#1])'
    show_result(opls_154)

    # """all-atom H(O): mono alcohols, OP(=O)2 """
    opls_155 = '[H][#8&%opls_154]'
    show_result(opls_155)

    # """alkane CH3 """
    opls_157 = '[C&D4]([H])([H])([O])([O&%opls_154])'
    show_result(opls_154)

    # TODO: @Blacklist(180)
    # """O: anisole """
    opls_179 = '[O&D2]([C&%opls_145])([C])'
    show_result(opls_179)

    # """O: dialkyl ether """
    opls_180 = '[O&D2]([C])([C])'
    show_result(opls_180)

    all_topos()
