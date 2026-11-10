# Mathematica decay-chain notebooks (secular-equilibrium ratios)

These notebooks are the origin of the typed ratio column in the `SE (1ppb)` and `Cu SE (1ppb)` tabs of
`../decay_chains_data.xlsx` (the column headed e.g. "232-Th (for 3e-14 per atom)", now column D of `SE (1ppb)`).

They were copied unchanged from the public repository `villano-lab/L3He_paper2023`, folder
`neutron-bknd-spectra/decay_chains/` (commit `126f262`, copied 2026-11-10), where they were written in 2023.

| Notebook | Chain | Method |
|---|---|---|
| `Th232_SE.nb` | Th-232 | symbolic solution; its six outputs equal the typed column to 6 digits |
| `U238_SE_ODE.nb`, `U235_SE_ODE.nb` | U-238, U-235 | NDSolve of the full chain ODEs; the typed U values appear in these |
| `*_SE.nb` (U), `*_ExpMat.nb`, `decay_chains_solver.nb`, others | | earlier or alternative solvers and tests |

## Known problems in the inputs (found 2026-11-10, against ENSDF)

- `U235_SE_ODE.nb`: the Fr-223 alpha branch is entered as 0.006 (a fraction) instead of 6e-5, so At-219 is 100 times too high.
- `U235_SE_ODE.nb`: the At-219 alpha and beta branches are swapped (93.6 % is sent to Rn-219 instead of Bi-215).
- `U235_SE_ODE.nb`: At-215 and Po-211 come out about 4e5 and 360 times too high in the typed column; the cause is not yet traced.
- Decay constants are hard-coded and rounded to about three significant figures.

The corrected ratios are computed live from the half-lives and branching on the `Th232`, `U235` and `U238` tabs
(column C of `SE (1ppb)`); the typed column is kept for comparison.
