# Decay-chain cross-checks (archived 2026-11-19)

These are independent checks of the secular-equilibrium ratios N_i/N_seed. They are **not inputs** to the analysis: the ratios
used for the SOURCES inputs are computed from the half-lives and branching ratios on the `Th232`, `U235` and `U238` tabs of
[`decay_chains_data.xlsx`](../../SOURCES/analysis/decay_chains_data.xlsx), and are explained and recomputed in
[`A-Contamination-Estimates.ipynb`](../../A-Contamination-Estimates.ipynb).

| Item | What it is |
|---|---|
| `mathematica_SE/` | The Mathematica ODE notebooks (from `villano-lab/L3He_paper2023`), with the input corrections made in November 2026. See its [README](mathematica_SE/README.md). The corrected U-235 ODE notebook agrees with the workbook to 0.55 %. `U238_SE_ODE.nb` (20 nodes) was edited but not re-run. |
| `A-Decay-Chain-Modeling.ipynb`, `decay_chains.py` | Layne Hellman's Python Monte Carlo / ODE model of the Th-232 and U-235 chains. It uses rounded half-lives, so it is **not** a reference for the numbers. Two bugs were fixed on 2026-11-19: the Ra-226 half-life was 1.602 years instead of 1600 years, and the Fr-223 branches were swapped (beta to Ra-223 is 99.994 %, alpha to At-219 is 0.006 %). Other entries were not audited. Run it from this folder (it imports `decay_chains` from here). |

Nothing in the repository imports these files. They are kept as a cross-check and for the record; the git history before
this move is reachable through `git log --follow`.
