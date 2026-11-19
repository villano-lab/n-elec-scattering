# Analysis history, November 2026 (Flux / SOURCES)

Context note extracted from the Claude Code session "Code repository location" (session `3ec1633e-9d7f-4359-881c-a07585c0d6b6`, 2026-11-09 to 2026-11-19). Paper-writing items (Section III text, LaTeX/VS Code setup, bib edits) are left out. For the workbook's own change list see [changelog.md](changelog.md) and [decay_chains_data_README.md](decay_chains_data_README.md). The daily logs (`villaa/daily-work-log`, `2026/11-November/`) hold the day-by-day detail.

## Repository and process

- Repo: `villano-lab/n-elec-scattering` (public). Default branch `develop`; PRs target `develop`; `main` (renamed from `master`) is release-only and reviewed by Anthony. Branch protection and CODEOWNERS set up; fast-tier CI and a minimal `environment.yml` (`nescatt_2025`) added (PR 1).
- SOURCES4A code (ORNL/LANL, Sheffield-modified libraries) was extracted to the private repo `sources4a-local` and pruned from this repo's public history. Old PR refs and forks may still carry the old history.
- `2-Flux/SOURCES/` reorganised (PR 8): `docs/`, `inputs/`, `results/`, `analysis/`, `notebooks/`, `norite/`. README there has a "Where did files move?" table.
- Cu, Pb, HDPE and Shotcrete Hang material removed from the workbook and archived in the private repo `villano-lab/radioactive-contamination`. The full 14-tab workbook is kept in `L3He_paper2023` (master).

## Workbook `analysis/decay_chains_data.xlsx` (now 9 tabs)

Tabs: Th232, U235, U238 (decay chains); Shotcrete Density, Shotcrete Target Fractions, SE (1ppb); Norite Density, Norite Target Fractions, SE Norite (1ppb).

1. **Decay tabs audited against ENSDF/NuDat 3.** Half-lives (unit conversions, scientific notation) and U-235 branching rows corrected; isotope-name columns added. At-218 (99.95 % alpha, 0.05 % beta to Rn-218) and a Rn-218 row inserted in U238 and propagated; effect is about -1e-7 on Bi-214.
2. **ppb definition.** The sheet originally used per-atom ppb (1e-9 atoms/atom). Adopted the standard mass fraction, 1 ppb = 1 ng/g. Seed density = 1e-9 * rho * N_A / A (shotcrete: Th-232 5.77e12, U-238 5.62e12 atoms/cc before the later density fix). Per-atom and per-mass differ by about 11x for shotcrete.
3. **SE tabs rebuilt.** Mass (amu) column and live secular-equilibrium ratios (column C, from the decay tabs); E-G use the per-atom convention, H-J the ng/g convention; I2/I11 are inputs and I26 is locked to U-238 by natural abundance. The old typed ratios (from Mathematica) had errors in the U-235 block (At-219, At-215, Po-211) and At-218.
4. **Mathematica notebooks** (from `L3He_paper2023` commit 126f262) fixed (percent-vs-fraction slips, swapped branches, swapped lambda assignments), re-run by Anthony, and now in `mathematica_SE/`. The ODE notebook agrees with column C within 0.55 %. `U238_SE_ODE.nb` was extended to 20 nodes but not re-run; short/ExpMat/U235_SE variants not checked.
5. **Shotcrete Density.** Unit error in E2 (1/100000 vs 1/1000000, 10x) fixed; C-13 upper bound corrected to 1.16 %. Composition source: SNOLAB-STR-2007-003 Table 3, shotcrete column (report in `docs/`). Density 2222.22 kg/m3 = 1000 kg / 0.45 m3 from the Sika King MS-D3 datasheet coverage, with caveats. Retained-water fraction added (wet columns K-N, input Q2, default 0 reproduces the old results to 4e-16); SE and Target Fractions tabs link to the wet columns (PR 16).
6. **Norite tabs added** from Layne Hellman's numbers (PR 19): GEOROC mean oxides for orthopyroxene (3.3 g/cm3) and plagioclase (2.7 g/cm3), 50:50 by volume (ratio in B21). Iron counted once as FeOT (the original notebook double-counts, about 104 % total). Implied density 3.049 g/cm3 vs 3000 kg/m3 used. Checked against `norite/inputs/composition.csv`. 1 ppb Th-232 = 7.786E12 and U-238 = 7.589E12 atoms/cm3 at 3.0 g/cm3. Tab order changed in PR 20 (SE after Target Fractions).

## SOURCES4A runs and collaborator exchange

- Built and ran SOURCES4A on this PC with conda-forge gfortran (needs `-O0`, `-fcheck=bounds`, `maxELevPts` overflow noted). A tagged `sources4a-local` version was recorded with its use.
- Correspondence with Vitaly Kudryavtsev (Sheffield) on library updates: his tape3 matches GitHub, his tape4 differs, and he supplied a new `sources_v20261109.f` (excited-state points 200 to 500; Ca-42, Cu, Mo). A reply listing the elements that matter and the tape4 difference was sent.
- Norite SOURCES input built (`norite/norite_sources.py`, iron double-count corrected, global average accepted).

## Open items (as of 2026-11-19)

- Analyze Vitaly's materials (`tape4`, `sources_v20261109.f`) and decide what to adopt.
- Regenerate shotcrete and norite SOURCES inputs/results from column C and the ng/g heads. The norite inputs still take their chain ratios from the shotcrete inputs built from the old values. Notebook B ppm scaling and `python/shielding.py` defaults also depend on the old numbers.
- Run all volumetric-source SOURCES inputs; then shielding model, fluxes, and the plot for Section III; check Table I against the inputs actually used.
- Workbook housekeeping: open the new tabs in Excel (not yet opened after the norite edit); CIAAW bounds for P/Ti/Cr/Ni (single IUPAC value used as both bounds); reconcile the three conflicting norite composition sources (Handbook 3.3 ppm Th / 1.2 U vs the 2007 report averages 5.56 / 1.11 vs Layne's); shotcrete water/H content (core measurement or SNOLAB); the `SE Results` wording in `changelog.md`; the tab lists in the changelog and README may not show the new tab order.
- Notebooks under `notebooks/` and `B-SOURCES4A-Running.ipynb` had only their path strings checked after the reorganisation, not run.
- ppb definition references: mass fraction (ng/g) per the radiopurity/assay literature (e.g. radiopurity.org); no single standards document found. Tracked in the private "Literature Log" artifact.
