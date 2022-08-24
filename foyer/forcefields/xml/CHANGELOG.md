# Forcefield Changelog
----
## OPLS-AA

v0.0.1 - April 15, 2021
 - started versioning of forcefield xmls
v0.0.2 - June 24, 2021
 - update SMARTS string for:
    -   opls_182 (from `[C;X4]([O;%opls_180])(H)(H)` to `[C;X4]([O;%opls_180])(H)(H)C`)
    -   opls_282(from `HC[C;%opls_277,%opls_280]` to `HC[C;%opls_277,%opls_280,%opls_465;!%opls_267]`)
    -   opls_468 (from `[C;X4]([O;%opls_467])(H)(H)` to `[C;X4]([O;%opls_467])(H)(H)H`)
    -   opls_469 (from `H[C;%opls_468]` to `H[C;%opls_468,%opls_490]`)
    -   opls_490 (`[C;X4]([O;%opls_467])(H)(H)C`)

 - update overrides for:
    -   opls_279 (from `"opls_144"` to `"opls_185, opls_144"`)
    -   opls_465 (from `""` to `"opls_277"`)
    -   opls_465 (from `""` to `"opls_278"`)

  - update references
    - opls_465 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_466 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_467 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_468 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_469 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)

v0.0.3 - August 7, 2021
 - update SMARTS string for opls_154 (from `[O;X2]H` to `[O;X2](H)([!H])`)

v0.0.4 - August 24, 2022
 - update atomlass for all atomtype entries in the XML (source from https://github.com/gromacs/gromacs/blob/main/share/top/oplsaa.ff/ffnonbonded.itp)
 - Code to transfer the atomclass:
 ```python
 import lxml as etree

 with open("ffnonbonded.itp", "r") as f:
    ref = f.readlines()

skimmed_ref = dict()
for line in ref:
    if "opls_" in line:
        tmp = line.strip().split()[0:2]
        skimmed_ref[tmp[0]] = tmp[1]

oplsaa_xml = etree.parse("oplsaa.xml")
root = xml.getroot()

for child in root.iterchildren():
    if child.tag == "AtomTypes":
        for gchild in child.iterchildren():
            if (type(gchild)==etree.Element and
            gchild.attrib["name"] in skimmed_ref):
            gchild.attrib["class"] = skimmed_ref[gchild.attrib["name"]]

xml.write("oplsaa.xml")
 ```

----
## Trappe-UA

v0.0.1 - April 15, 2021
 - started versioning of forcefield xmls
v0.0.2 - August 9, 2021
 - Updated `combining_rule` to `lorentz`
