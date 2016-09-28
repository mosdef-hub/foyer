from smarts_parser import Rule

opls_rules = list()

# """O TIP3P Water """
opls_111 = '[#12&D2]([#1])[#1]'
r111 = Rule('opls_111', opls_111)
opls_rules.append(r111)

# """H TIP3P Water """
opls_112 = '[#1][%opls_111]'
r112 = Rule('opls_112', opls_112)
opls_rules.append(r112)

# """alkane CH3 """
opls_135 = '[#12&D4]([#12])([#1])([#1])([#1])'
r135 = Rule('opls_135', opls_135)
opls_rules.append(r135)

# """alkane CH2 """
opls_136 = '[#12&D4]([#12])([#12])([#1])([#1])'
r136 = Rule('opls_136', opls_136)
opls_rules.append(r136)

# """alkane CH """
opls_137 = '[#12&D4]([#12])([#12])([#12])([#1])'
r137 = Rule('opls_137', opls_137)
opls_rules.append(r137)

# """alkane CH4 """
opls_138 = '[#12&D4]([#1])([#1])([#1])([#1])'

"""alkane C """
opls_139 = '[#12&D4]([#12])([#12])([#12])([#12])'

# TODO: be specific about alkane C (not just simply C)
# """alkane H """
opls_140 = '[#1][#12]'
r140 = Rule('opls_140', opls_140)
opls_rules.append(r140)



# NOTE: this catches benzene C as well
# """alkene C (R2-C=) """
opls_141 = '[#12&D3]([#12])([#12])([#12])'

# NOTE: this catches benzene C as well
# """alkene C (RH-C=) """
opls_142 = '[#12&D3]([#12])([#12])([#1])'

# NOTE: this catches benzene C as well
# """alkene C (H2-C=) """
opls_143 = '[#12&D3]([#12])([#1])([#1])'

# """alkene H (H-C=) """
# NOTE: We make sure that the carbon is an alkene carbon.
opls_144 = '[#1]([#12&%opls_141,#12&%opls_142,#12&%opls_143])'

# TODO: @Blacklist([141, 142])
# """Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl """
opls_145 = '[#12&D3&R6]([#12&R6])([#12&R6])'

# TODO: @Blacklist([145])
# """Biphenyl C1 """
opls_145B = '[#12&D3]([#12&%opls_145])([#12&%opls_145])([#12&%opls_145])'

# TODO: @Blacklist([140, 144])
# """Benzene H - 12 site. """
opls_146 = '[#1]([#12&%opls_145])'

# TODO: @Blacklist(135)
# """C: CH3, toluene """
opls_148 = '[#12&D4]([#12&%opls_145])([#1])([#1])([#1])'

# TODO: @Blacklist(136)
# """C: CH2, ethyl benzene """
opls_149 = '[#12&D4]([#12&%opls_145])([#12])([#1])([#1])'

# TODO: make sure it doesn't catch water
# """all-atom O: mono alcohols """
opls_154 = '[#8&D2]([#1])'
r154 = Rule('opls_154', opls_154)
opls_rules.append(r154)

# """all-atom H(O): mono alcohols, OP(=O)2 """
opls_155 = '[H][#8&%opls_154]'
r155 = Rule('opls_155', opls_155)
opls_rules.append(r155)

# """alkane CH3 """
opls_157 = '[C&D4]([H])([H])([O])([O&%opls_154])'
r157 = Rule('opls_157', opls_157)
opls_rules.append(r157)

# TODO: @Blacklist(180)
# """O: anisole """
opls_179 = '[O&D2]([C&%opls_145])([C])'

# """O: dialkyl ether """
opls_180 = '[O&D2]([C])([C])'
