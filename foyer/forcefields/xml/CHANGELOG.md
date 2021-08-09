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

----
## Trappe-UA

v0.0.1 - April 15, 2021
 - started versioning of forcefield xmls
v0.0.2 - August 9, 2021
 - Updated `combining_rule` to `lorentz`
