# Mathematica decay-chain notebooks (secular-equilibrium ratios)

These notebooks are the origin of the typed ratio column in the `SE (1ppb)` and `Cu SE (1ppb)` tabs of
`../decay_chains_data.xlsx` (column D there, "N/N_seed from Mathematica ...").

They were copied from the public repository `villano-lab/L3He_paper2023`, folder
`neutron-bknd-spectra/decay_chains/` (repository commit `126f262`, copied 2026-11-10), where they were written in 2023.
The first commit of this folder in this repository (`5cbfe89`) holds them unchanged; later commits correct some inputs (below).

| Notebook | Chain | Method |
|---|---|---|
| `Th232_SE.nb` | Th-232 | symbolic solution; its six outputs equal the typed column to 6 digits |
| `U238_SE_ODE.nb`, `U235_SE_ODE.nb` | U-238, U-235 | NDSolve of the full chain ODEs; the typed values come from the last cell of each |
| `U235_SE_ODE_short.nb`, `U235_SE_ExpMat.nb`, `*_SE.nb` (U), `decay_chains_solver.nb`, others | | earlier or alternative solvers and tests |

## Input errors found (2026-11-10, against ENSDF) and what was changed

The original typed values were reproduced from the original inputs with an independent calculation (all 11 U-235-chain rows to 0.1 %), which shows where they come from.

`U235_SE_ODE.nb` (same edits in `U235_SE_ODE_short.nb` and `U235_SE_ExpMat.nb`):

- `b6a` (Fr-223 alpha branch) was 6e-3; ENSDF has 0.006 % = 6e-5. This made At-219 100 times too high. Now `6*^-5`.
- `b8a` / `b8b` (At-219) were swapped (0.064 and 0.936); alpha, to Bi-215, is 0.936. Now swapped back.
- `b11a` (Po-215 beta branch) was 2.3e-4; the sheet and ENSDF value is 0.00023 % = 2.3e-6. Now `23*^-7`.
- The decay constants of the two Po-215 daughters were assigned to the wrong nodes: n12 (fed by the alpha branch `b11b`, i.e. Pb-211) carried At-215's decay constant and n13 carried Pb-211's. `\[Lambda]12` and `\[Lambda]13` are swapped, so n12 is now Pb-211 and n13 is At-215. This had made the typed "At-215" about 4e5 times too high.
- The same swap for the Bi-211 daughters: `\[Lambda]15` and `\[Lambda]16` are swapped, so n15 is Tl-207 and n16 is Po-211. This had made the typed "Po-211" about 360 times too high.
- The plot and result cells at the end now use n13 (At-215) and n16 (Po-211) instead of n12 and n15.
- `U235_SE_ExpMat.nb` only: `b11b` was 0.1, now 1 (its decay constants use a different time unit, which was left alone).

`U238_SE_ODE.nb`: the At-218 decay constant `\[Lambda]9` corresponded to a 1.5 s half-life; now `110*^15` (1.28 s, ENSDF and NuDat 3).

Other: decay constants are hard-coded and rounded to about three significant figures.

## Caveats

- The notebooks were **not** re-evaluated (no Mathematica was available). Saved outputs, equations and plots in the corrected notebooks come from the old inputs and are stale; each corrected notebook starts with a note saying so. The internal cache section of these files was removed so that Mathematica rebuilds it when the notebook is opened.
- Checked instead by an independent calculation: with the corrected inputs the U-235 chain agrees with column C of `SE (1ppb)` (ratios computed from the decay tabs) to within 0.5 %.
- Columns E-J of the SE tabs use column C. Column D keeps the original (uncorrected) notebook values for comparison.
