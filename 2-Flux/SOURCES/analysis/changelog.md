# Change log of `decay_chains_data.xlsx`

Everything done to the workbook since the version of **January 2026**, and why. The baseline is the file as last changed on
2026-01-14 (commit `173fbb2`, unchanged at `98263c6`, 2026-01-29). That file had 15 tabs (`Th232`, `U235`, `U238`,
`SE Results`, `SE (1ppb)`, `Shotcrete Density`, `Shotcrete Hang`, `Shotcrete Target Fractions`, `HDPE Density`,
`HDPE Target Fractions`, `Cu SE (1ppb)`, `Cu Density`, `Cu Target Fractions`, `Pb Density`, `Pb Target Fractions`).
The current file has nine (`Th232`, `U235`, `U238`, `SE (1ppb)`, `Shotcrete Density`, `Shotcrete Target Fractions`, and the three norite tabs added 2026-11-18: `Norite Density`, `Norite Target Fractions`, `SE Norite (1ppb)`).

The explainer of the current state is [`decay_chains_data_README.md`](decay_chains_data_README.md). The Mathematica notebooks
that produced the original equilibrium ratios are in [`mathematica_SE/`](../../archive/decay-chain-crosschecks/mathematica_SE/README.md). Annotated pictures of every

**Moved 2026-11-19:** the Mathematica notebooks, `A-Decay-Chain-Modeling.ipynb` and `decay_chains.py` are cross-checks, not inputs, and now live in [`2-Flux/archive/decay-chain-crosschecks/`](../../archive/decay-chain-crosschecks/README.md). The equilibrium ratios come from the half-lives and branching ratios on the decay tabs; `2-Flux/A-Contamination-Estimates.ipynb` explains and recomputes them.
changed tab (before against now) are on Trello, notebooks board, list N-MISC-26-002: pg 31 (`Th232`, `U235`, `U238`),
pg 32 (`Shotcrete Density`, `Shotcrete Target Fractions`) and pg 33 (`SE (1ppb)`); the 2026-11-16 daily log has the same pictures.

## Summary: what changed the numbers

| # | Change | Effect on the results |
|---|---|---|
| 1 | Unit factor in `Shotcrete Density` E2 was `1/100000`, now `1/1000000` | Every density and number density 10x too high (22.2 g/cm^3 instead of 2.22). The matrix atom density typed into `SE (1ppb)` F2 and F11 (6.34619E23) was 10x too high (now a live link). |
| 2 | "1 ppb" in `SE (1ppb)` was 1e-9 atoms per atom of matrix | 11.0x (Th-232) and 11.3x (U-238) more seed atoms than the literature ppb (1 ng/g); see the section on ppb below. Together with 1: about 110x too many atoms per "ppb" against a mass-fraction ppb. |
| 3 | U-235 chain equilibrium ratios (typed from `U235_SE_ODE.nb`) had wrong inputs | At-219 100x, At-215 about 4.3e5x, Po-211 about 360x too high; Rn-219 and Bi-215 about 0.1 % off. |
| 4 | At-218 half-life 1.5 s in the U-238 notebook and tab, ENSDF 1.28 s | At-218 17 % too high. |
| 5 | C-13 upper abundance bound 1.6 %, CIAAW range 0.96-1.16 % | C-13 middle target fraction 8.99E-06 to 7.44E-06 (-17 %), upper 1.12E-05 to 8.15E-06 (-27 %). |
| 6 | Half-lives refreshed to ENSDF / NuDat 3 | Th-232, Pa-231, Po-216, Po-212 and Bi-214 changed by 0.6-1.6 % (and with them the ratios in column C). |
| 7 | At-218 branching (99.95 % alpha, 0.05 % beta) and a Rn-218 side branch | Bi-214 -1e-7 relative, Po-214 +2e-11; Rn-218 at 2.4e-26 of the seed. |

## Changes by date

### 2026-11-09 and 2026-11-10 (commits `f1f6b06`, `b2f8128`, `5cbfe89`, `b62a9bf`, `3c32df2`, `dd915ba`, `65e8d6e`)

- **Density tabs (Anthony, 2026-11-09, `f1f6b06`; Claude, 2026-11-10, `b2f8128`).** Divide by 1e6 for g/cm^3 (E2 of
  `Shotcrete Density`, and the same cell on the HDPE, Cu and Pb tabs). The other cells of column F, H follow (F2 1.342 to 0.1342 g/cm^3,
  H2 3.0E22 to 3.0E21 atoms/cm^3). The fraction of atoms (column I) was never affected: the factor 10 cancels.
- **Decay tabs `Th232`, `U235`, `U238`** (`b2f8128`):
  - Half-lives updated to ENSDF / NuDat 3: Th-232 1.41E10 to 1.40E10 y, Po-216 0.145 to 0.144 s, Po-212 299 to 294.3 ns (and the
    Th-232 lambda cell G2 4.92E-11 to 4.95E-11), Pa-231 32,760 to 32,570 y, Bi-214 19.9 to 19.71 min, At-218 1.5 to 1.28 s.
    All other half-lives were checked against ENSDF (within 0.25 %).
  - `U235` branching table corrected: Fr-223 alpha 0.6 % to 0.006 %; the At-219, Po-215 and Bi-211 beta and alpha values were in
    the wrong columns (the table columns are "beta (%)" then "alpha (%)"); Po-215 beta 0.023 % to 0.00023 %.
  - An `Isotope` name column (H, or D beside the branching tables) was added to each decay tab.
- **`SE Results` removed** (`b2f8128`). It held the same chain at the unscaled Mathematica seeds (Th-232 3E-14, U-238 9E-15 atoms per atom;
  E2 = 1.9E10 atoms/cm^3). The commit message does not give a reason; the tab is in the history and in `L3He_paper2023` (commit before `20b4f01`).
- **`SE (1ppb)` rebuilt** (`5cbfe89`, `b62a9bf`, `3c32df2`, `dfceacc`; `Cu SE (1ppb)` the same, since removed). New layout, three blocks (Th-232, U-238,
  U-235 chains):
  - B: isotope mass (amu, IAEA Live Chart of Nuclides).
  - C: N/N_seed at secular equilibrium as a **live formula** from the decay tabs (it was typed Mathematica output).
  - D: the old Mathematica values divided by their typed seed (3E-14, 9E-15, 6.4836E-17) so the seed is 1, marked "uncorrected; not used in E-J".
  - E-G: per-atom ppb (C x F; G = K x E).
  - H-J (new): ng/g ppb. I2 and I11 are the input (mass fraction 1E-9); I26 (U-235) is locked to U-238 by natural abundance
    (0.007204 by atoms, `=0.007204*I11*B26/B11`); H = C x I x rho x N_A / (A_seed x K); J = K x H.
  - K (new): the matrix atom density, a **live link** `='Shotcrete Density'!H13` (it was typed in the January file, F2 = 6.34619E23, 10x too high), L: Q-alpha.
  - E-J use column C, not column D.
  - Notes in rows 40-46 name the exact files and repositories used.
- **Mathematica notebooks** (`mathematica_SE/`): copied unchanged first (`5cbfe89`, from `villano-lab/L3He_paper2023`), then the
  inputs of `U235_SE_ODE.nb`, `U235_SE_ODE_short.nb`, `U235_SE_ExpMat.nb` and `U238_SE_ODE.nb` were corrected (`dd915ba`): Fr-223 alpha branch
  100x too high, At-219 branches swapped, Po-215 beta branch 100x too high, the decay constants of the Po-215 daughters (Pb-211, At-215) and
  the Bi-211 daughters (Tl-207, Po-211) on the wrong nodes, At-218 half-life 1.5 s. `U235_SE_ODE.nb` was re-run in Mathematica on 2026-11-11
  (`e5ae84a`, `60d48c1`): agrees with column C within 0.55 %.
- **`Shotcrete Target Fractions`:** the C-13 upper bound E4 1.6 to 1.16 (`dd915ba`). The same fix and the removal of column I on the HDPE tab (`65e8d6e`);
  `Cu Density` B2 atomic number 6 to 29 (since removed).

### 2026-11-11 (`9be144d`, `b8f326b`, `52a90a0`)

- Explainer `decay_chains_data_README.md` added; the workbook moved from `sources_running/data/raw_sources_IO/` to
  `2-Flux/SOURCES/analysis/` together with the other analysis files (9 exact duplicates dropped, path references updated).

### 2026-11-16 (`41974c9`, `fb3681f`, `31b115d`, `14ddd42`, `65a91f3`, `ac60848`, `5984296`, `5347da7`, `88bfb5c`)

- **At-218 branching and the Rn-218 side branch** (`41974c9`, pull request 9). `U238` has a new branching row 29 (At-218 beta 0.05 %,
  alpha 99.95 %) and a new decay row 20 (decay number "9b", Rn-218, 33.75 ms, 100 % alpha to Po-214). `SE (1ppb)` has a Rn-218 row (23), and new
  formulas for Bi-214, Po-214, Po-210, Bi-210 and Pb-210. `mathematica_SE/U238_SE_ODE.nb` has the 20-node system (edited, **not re-run**; its steady
  state agrees with column C to 0.2 %) and the Bi-214 decay constant is now the 19.71 min value.
- **Tabs removed as not needed for the n-electron scattering paper:** `Cu SE (1ppb)`, `Cu Density`, `Cu Target Fractions`, `Pb Density`,
  `Pb Target Fractions` (`fb3681f`); `HDPE Density`, `HDPE Target Fractions`, with their inputs, results and notebook, and `Cu.ipynb`
  (`31b115d`); `Shotcrete Hang` and column I of `Shotcrete Target Fractions` (which multiplied by Hang's total 0.0474) (`14ddd42`).
  Everything removed is in the git history (`develop` up to `34792c4`), in the private repository `villano-lab/radioactive-contamination`
  and in `villano-lab/L3He_paper2023` (full workbook, commit `20b4f01`).
- **Shotcrete composition and density sourced** (`65a91f3`, `ac60848`, `5984296`):
  - Composition (column C): SNOLAB-STR-2007-003 (Lawson, *Analysis of Rock Samples from the New Laboratory*), Table 3, "Shotcrete" column
    (sample 15, MS-D3 accelerated shotcrete, Creighton Mix 25017). The report is in `2-Flux/SOURCES/docs/`.
  - Density (D2 = 2222.22 kg/m^3): 1000 kg / 0.45 m^3 from the "Coverage" line of the Sika King MS-D3 datasheet (April 2024, v06.02).
    Approximate (the same table gives 2143 and 2141 kg/m^3); a dry-mix yield, not a measured density.
  - References rewritten in A16:A18 of the tab.
- **Retained-water fraction** (`5347da7`). `Shotcrete Density` has the water mass fraction f_w (Q2, default 0), wet columns K-N (mass fraction,
  elemental mass density, number density, fraction of atoms), the water H and O fractions (Q4, Q5), wet atom fractions by Z in C22:C32, and
  three reference lines in A34:A36. At f_w = 0 everything equals the dry columns. The ppm are taken per dry mass: rho stays E2.
- **Links to the wet columns** (`88bfb5c`, pull request 16). `SE (1ppb)` K2, K11, K26 = `'Shotcrete Density'!M13` (were `=...!H13`, the dry value);
  `Shotcrete Target Fractions` column B (atom fractions) = INDEX/MATCH on Z into C22:C32, isotope rows copy their element row (were typed).
  Identical to the previous values at f_w = 0 (to 4E-16).

### 2026-11-18 (`b3eec84`)

- Notes on the ppb/ppm convention in `SE (1ppb)` rows 47-49 and in the README (see the next section).

### 2026-11-18 (norite tabs)

- **Norite added (`Norite Density`, `Norite Target Fractions`, `SE Norite (1ppb)`).** Built from Layne Hellman's GEOROC oxide table and mineral densities
  (`C-Slab-Source-and-Shielding.ipynb`, the numbers already used by `2-Flux/SOURCES/norite/norite_sources.py`), in the layout of the shotcrete tabs. Iron is counted once (FeOT; the
  notebook as written double-counts it). The composition reproduces `norite/inputs/composition.csv` for all 13 elements; the head densities for 1 ng/g (J2 = 7.786E12 Th-232, J11 = 7.589E12 U-238 atoms/cm^3) match the SOURCES inputs.
  P, Ti, Cr and Ni have no CIAAW abundance bounds yet (single values). No existing tab was changed. See the README section "Norite tabs".

## ppb and ppm: the definition used by the literature (mass fraction, not atom density)

The `SE (1ppb)` tab used **1 ppb = 1e-9 atoms of the seed per atom of the matrix** (columns E-G; D2 = 1E-9 in the January file). In the
literature on radiopurity and on U and Th in rock, ppm and ppb are **mass fractions**: 1 ppm = 1 ug/g, 1 ppb = 1 ng/g = 1e-9 g/g. The
two differ by A_matrix / A_seed (mean shotcrete atom 21.1 g/mol; Th-232 232.04, U-238 238.05): per-atom ppb gives 11.0x (Th-232) and 11.3x (U-238)
more seed atoms/cm^3 than 1 ng/g. Columns H-J now follow the literature. Seed atoms/cm^3 per 1 ng/g in shotcrete: 5.77E12 (Th-232), 5.62E12 (U-238);
per 1 per-atom ppb: 6.35E13.

References (checked 2026-11-18):

1. **radiopurity.org** (the PNNL / SNOLAB radiopurity database), conversion chart on its front page,
   <https://www.radiopurity.org/static/logos/conversion_chart.png>: "1 Bq U-238/kg = 81 ppb U (81 x 10^-9 gU/g)" and
   "1 Bq Th-232/kg = 246 ppb Th (246 x 10^-9 gTh/g)". These give 12.3 mBq/kg per ppb U and 4.07 mBq/kg per ppb Th, equal to the
   values calculated here (12.4 and 4.06). (The K-40 row of the same chart gives "32300 ppb K (32300 x 10^-6 gK/g)": the two
   forms disagree by 1000; the chart is consistent only for U and Th.) The site is <https://www.radiopurity.org/>; code at
   <https://github.com/pnnl/Radiopurity-database-assistant>.
2. **J. W. Grate et al. (PNNL)**, "Mass spectrometric analyses of high performance polymers to assess radiopurity as ultra low background materials
   for rare event physics detectors", <https://arxiv.org/abs/2007.03104> (NIM A). Reports "mass concentrations" of Th-232 and U-238 in pg/g and
   converts them to activity: "1 pg/g of 232Th and 238U are equivalent to 4.1 and 12.4 [u]Bq/kg", i.e. 4.1 and 12.4 mBq/kg per ppb, and
   notes the convention of quoting parent activities assuming secular equilibrium.
3. **UK Dark Matter Collaboration, "Useful numbers and relationships"** (radioactivity page of the UKDMC site, Rutherford Appleton Laboratory),
   <https://hepwww.pp.rl.ac.uk/groups/ukdmc/Radioactivity/useful.html>: concentration C in g/g, ppb = parts in 10^9; 1 ppb Th = 0.00406 Bq/kg (parent
   only), 1 ppb natural U = 0.0123 Bq/kg (U-238) + 5.7e-4 Bq/kg (U-235); 0.0406 and 0.179 Bq/kg for the whole chain in equilibrium.
4. **M. Tzortzis and H. Tsertos**, "Determination of thorium, uranium, and potassium elemental concentrations in surface soils in Cyprus",
   <https://arxiv.org/abs/physics/0309008>: elemental concentrations of Th and U in ug/g and K in %, from activity concentrations by
   F = M lambda C A / (N_A f_A) (mass fraction of the element, mass M, decay constant lambda, natural abundance f_A).
5. **Ambiguity of the term itself.** IUPAC defines the amount (mole, number) fraction separately,
   <https://goldbook.iupac.org/terms/view/A00296>, and recommends against using ppm and ppb as units; the background paper
   "Units and quantities for analytical chemistry", <https://pubs.rsc.org/en/content/articlehtml/2020/ay/d0ay90126a>, notes that ppm and ppb are used
   for mass, volume or amount fractions alike and that the quantity should be named. Atom or mole fractions are conventional for gases
   (ppbv, ppt Kr in Xe), not for U and Th in a solid.

No single standards document defining "ppb of U or Th" was found: the mass-fraction reading is the convention of the assay and radiopurity
literature above. The SNOLAB report behind Table I of the paper (SNOLAB-STR-2007-003) quotes the U and Th in ppm from its chemical assay;
its definition of ppm was not checked in the text and should be read before the paper states "ppm by mass".

## Still open

- Re-run `U238_SE_ODE.nb` (and `U235_SE_ODE_short.nb`) in Mathematica; refresh column D of `SE (1ppb)` if wanted (the Rn-218 row is empty there).
- Shotcrete water and hydrogen content (a core measurement, or SNOLAB): choose the f_w cases for the sensitivity run.
- Norite: in the workbook since 2026-11-18, but the composition is the GEOROC average, and three sources of the radioactivity and composition disagree (SNOLAB User's Handbook 3.3 ppm Th / 1.2 ppm U;
  SNOLAB-STR-2007-003 averages 5.56 / 1.11; Table I of the paper 5.44 / 1.06). Enter the CIAAW abundance intervals for P, Ti, Cr and Ni on `Norite Target Fractions`.
- The SOURCES inputs and results were built from the old values (changes 1-3 above) and must be regenerated from column C and the ng/g head densities.
