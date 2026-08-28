# Resummation additions
This is H1jet with the ARES resummation coefficients added.
`--obs` sets the observable: `onejet` (default), `Tminor`, `Tthrust`.
`--resum` sets which coefficient is computed: `none` (default, resummation
off), then `h12`, `h11`, `h10`, `h24`, `h23`.

Some of the observable dependence is not covered by `--obs` and has to be
switched by hand. The alternatives are commented out in place, each marked
with the observable it belongs to:
- `initial_state_rad.f90`
- `vboson.f90`
- `cross_sections.f90`
