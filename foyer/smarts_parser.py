
import re, collections
import plyplus

smarts_grammar = plyplus.Grammar("""
    start: string;

    // Rules
    @string: chain branch | chain;
    branch: LPAR string RPAR;
    @chain: atom chain | atom;
    atom: LBRACKET or_expression RBRACKET;
    ?or_expression: and_expression or_symbol or_expression | and_expression;
    ?and_expression: atom_id and_symbol or_expression | atom_id;
    and_symbol: SEMI | AMP;
    or_symbol: COMMA;
    atom_id: HASH NUM | STAR |DOLLAR LPAR string RPAR;

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
""")


def parse(expr):
    tree = smarts_grammar.parse(expr)
    return tree

if __name__ == '__main__':
    # print(parse('[#1]', 'atom'))
    # print(parse('[#1;#2]', 'atom'))
    # print(parse('[#1]', 'chain'))
    # print(parse('[#1]', 'string'))
    # print(parse('[#1][#2][#3]', 'chain'))
    # print(parse('[#1][#2][#3]', 'string'))
    # print(parse('([#1])', 'branch'))

    t = parse('[#1][#2][#3]')

    # t.to_png_with_pydot('tree.png')
    print t.pretty()
