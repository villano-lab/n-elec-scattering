# `decay_chains_data.xlsx` explained

This workbook holds the inputs that turn "so many ng/g of U or Th in a material" into "so many atoms/cm^3 of each
daughter isotope", and the target-material fractions needed to turn those into (alpha,n) neutron source terms for the
SOURCES-4 inputs. It lives in `2-Flux/SOURCES/analysis/`; the Mathematica notebooks that
produced the original equilibrium ratios are in `mathematica_SE/` (see its [README](mathematica_SE/README.md)).

Last audited 2026-11-10/11 against ENSDF (IAEA Live Chart API, NuDat 3) and CIAAW isotopic abundances; see
"Status of each tab" below.

## Overview of the chain of calculation

```
Th232 / U235 / U238 tabs   half-lives, branching ratios (ENSDF)
          |
          v
SE (1ppb) column C   N_i / N_seed at secular equilibrium (live formula)
          |  x  seed contamination (ppb)   <- Density tab gives rho, atom density
          v
atoms/cm^3 of every chain member in the material (columns E-G per-atom ppb, H-J ng/g ppb)

Density tab -> Target Fractions tab -> (a,n) target isotopes, with abundance bounds   (feeds the SOURCES inputs)
```

## Tabs

| # | Tab | Purpose |
|---|---|---|
| 1 | `Th232` | Th-232 chain: decay number, half-life (yr), mean life, lambda, scaled versions, isotope; branching of Bi-212 (rows 16-17). |
| 2 | `U235` | Same for the U-235 chain; branching table (rows 24-29) for Ac-227, Fr-223, At-219, Po-215, Bi-211. |
| 3 | `U238` | Same for the U-238 chain (Pa-234m isomer included). |
| 4 | `SE (1ppb)` | Shotcrete (norite) contamination: ratios, masses and atoms/cm^3, three blocks (see below). |
| 5 | `Shotcrete Density` | Composition (mass %), density 2.2222 g/cm^3, molar masses, number densities, atom fractions. |
| 6 | `Shotcrete Target Fractions` | For each element, the (a,n) target isotopes with CIAAW abundance bounds and the middle value (fraction of atoms). |

The Cu, Pb and HDPE tabs (`Cu SE (1ppb)`, `Cu Density`, `Cu Target Fractions`, `Pb Density`, `Pb Target Fractions`, `HDPE Density`,
`HDPE Target Fractions`) were removed on 2026-11-16 because they are not needed for the n-electron scattering paper. Nothing else
referenced them. The full 14-tab workbook is in the history of this repository (`develop` up to `34792c4`), in the private repository
`villano-lab/radioactive-contamination` and in `villano-lab/L3He_paper2023` (`neutron-bknd-spectra/decay_chains/`, commit `20b4f01`).

The `Shotcrete Hang` tab (Hang's published atom fractions, renormalised, compared with the atom fractions of `Shotcrete Density`) was also
removed on 2026-11-16, together with column I of `Shotcrete Target Fractions`, which multiplied by that tab's total (0.0474, Hang's
normalisation). Both are in the same three places (the full workbook in the archive repository still has them).

### Decay tabs (1-3)

- Column B is the half-life in years. Where the source gives seconds, minutes or days the cell holds the conversion
  (for example `=3.6319/365`), using a 365-day year. This differs from 365.2422 days by a uniform 0.07 %.
- C = B/ln 2 (mean life), D = 1/C (lambda), E and F are the same scaled to the seed (G2 is the rounded seed lambda).
- H gives the isotope name. All half-lives were checked against ENSDF (within 0.25 %), and the branching tables
  against ENSDF.
- Branching table columns are "beta (%)" and "alpha (%)" (in that order). The U238 table has a row for At-218
  (row 29, listed last so no other row moved): 99.95 % alpha to Bi-214 and 0.05 % beta to Rn-218 (ENSDF; the beta
  branch is as large as its own uncertainty). Rn-218 (33.75 ms, 100 % alpha to Po-214) is a side branch of the chain: it
  has its own row 20 on the U238 tab (decay number "9b") and row 23 on the SE tab. The SE formulas for Bi-214,
  Po-214, Po-210, Bi-210 and Pb-210 use these branches. Effect: Bi-214 -1e-7 relative, Po-214 +2e-11, the rest
  unchanged; Rn-218 is at 2.4e-26 of the seed. The original Mathematica notebooks never had this branch
  (column D of the Rn-218 row is empty); `mathematica_SE/U238_SE_ODE.nb` now has it (edited 2026-11-16, not yet re-run).

### SE tab (4)

Three blocks, one per seed: **Th-232** (header row 1, data rows 2-8), **U-238** (header row 10, data rows
11-22) and **U-235** (header row 25, data rows 26-37). Row labels are `<Z A>-<element>`, e.g. `922350-U`.

| Column | Content |
|---|---|
| A | Isotope. |
| B | Atomic mass (amu), IAEA Live Chart ground-state data (AME-based; edition not checked). |
| C | **N_i/N_seed at secular equilibrium**, a live formula from the decay tabs: (T_1/2 of isotope / T_1/2 of seed) x fraction of seed decays passing through the isotope. The seed is 1. This is the column the others use. |
| D | The original Mathematica values (from `mathematica_SE/`), divided by the typed seed value (3E-14, 9E-15, 6.4836E-17) so that the seed is 1. Kept for comparison; **not used** in E-J. |
| E | `C x F` (per-atom ppb scaling). |
| F | Seed contamination in **atoms per atom of matrix** (head row only; 1e-9 = 1 "per-atom ppb"). |
| G | `K x E`, atoms/cm^3. |
| H | `C x I x rho x N_A / (A_seed x K)`: the ng/g scaling. |
| I | Seed contamination as a **mass fraction** (head row only; 1e-9 = 1 ng/g). The input cell. For U-235 it is locked to U-238: `=0.007204*I11*B26/B11`. |
| J | `K x H`, atoms/cm^3. |
| K | Matrix atom number density (atoms/cm^3), linked to the Density tab (`H13` for shotcrete). |
| L | Q-alpha (MeV). |
| M | Notes. |

### ppb and ppm: the literature convention against the per-atom ppb used so far

In the literature, and in the SNOLAB report behind Table I of the paper, ppm and ppb of U, Th (and K) in a material are
**mass fractions**: 1 ppm = 1 ug/g, 1 ppb = 1 ng/g = 1e-9 g/g. (1 ppb U-238 = 12.4 mBq/kg and 1 ppb Th-232 = 4.06 mBq/kg,
computed here.) Atom (mole) fractions are used for gases and dissolved impurities, for example ppt Kr in Xe, not for U/Th in a
solid.

The workbook carries both definitions side by side:

- **ng/g ppb** (H-J) is the literature definition: 1 ppb = 1e-9 g of the seed per g of matrix. Seed atoms/cm^3 =
  `1e-9 x rho x N_A / A_seed`; for 1 ng/g in shotcrete that is 5.77e12 (Th-232) and 5.62e12 (U-238) atoms/cm^3. **Use these
  columns with a measured ppm x 1000.**
- **per-atom ppb** (E-G) is *not* the literature ppb: 1 ppb = 1e-9 atoms of the seed per atom of matrix, which puts
  A_matrix/A_seed too many seed atoms in. For shotcrete (mean matrix atom 21.1 g/mol) 1 per-atom ppb is 6.35e13 atoms/cm^3,
  11.0 times (Th-232) and 11.3 times (U-238) the ng/g value. E-G are kept only to compare with the earlier results.

Example: 15.24 ppm Th (Table I, shotcrete) is 8.8e16 Th-232 atoms/cm^3 by mass (15240 x J2) but 9.7e17 if read as per atom
(15240 x G2). With the earlier matrix density (typed 10x too high) the absolute source was about 110x too large.

U-235 is tied to U-238 by the natural abundance 0.007204 (by atoms), in both definitions.

Rows 40-46 of the sheet hold the notes naming the exact files and repositories used.

### Density tabs

Mass fraction (C) and overall density (D, kg/m^3; E in g/cm^3) give the elemental mass density (F = E x C/100),
the number density (H = F/G x N_A, with N_A = 6.0221408E+23) and the atom fractions (I = H / total). The atom
fractions are copied by value into the lower table (rows 21+ for shotcrete) that the Target Fractions tab uses.
The shotcrete composition (column C) is Table 3, "Shotcrete" column, of SNOLAB-STR-2007-003 (I. T. Lawson, *Analysis of Rock Samples
from the New Laboratory*, 2007): sample 15, the MS-D3 accelerated shotcrete (Creighton Mix 25017); all 11 elements match. It
sums to 99.66 %; the remaining 0.34 % is the trace elements (Table 2 of the report gives 3582.7 ppm = 0.36 % for that sample).
The same report is the source of the shotcrete Th and U (15.24 and 2.46 ppm) used in the paper. The density (D2 = 2222.22 kg/m^3) is not in that report or in the SNOLAB User's Handbook; it is 1000 kg / 0.45 m^3 from the
"Coverage" line of the Sika King MS-D3 datasheet (April 2024, v06.02: "0.45 m^3 per 1000 kg FIBC"). It is approximate: the same
table gives 0.014 m^3 per 30 kg bag (2143 kg/m^3) and 16.5 ft^3 per 2205 lb (2141 kg/m^3), about 4 % lower, and it is the yield of the
dry mix (no mixing water, no rebound), not a measured density of the cured shotcrete. The 2025 notebook (N-MISC-25-001 pg 26) has the
same calculation. A copy of the report is in `../docs/SNOLAB-STR-2007-003-NoriteRock-Radioactivity.pdf`.

**Water (added 2026-11-17).** Columns K-N, the input `Q2` (f_w, the water mass fraction of the wet material; 0 gives back columns C-I), `Q3` (wet density
= E2/(1-f_w), same dry mass per volume) and `C22:C32` (wet atom fractions by Z) give the composition and number densities with retained water: the dry
composition is scaled by (1-f_w) and water (11.19 % H, 88.81 % O) is added. For example f_w = 5 % gives 2.339 g/cm^3, H 0.60 %, O 50.0 %, and a
total number density 18 % higher (7.52e22 against 6.35e22 atoms/cm^3), most of it hydrogen. `SE (1ppb)` (K2, K11 and K26 = `M13`) and `Shotcrete Target
Fractions` (column B, from `C22:C32`; F-H follow) now use the wet columns, and give the old results exactly at f_w = 0 (checked: every value on
the Target Fractions tab agrees to 4e-16). The ppm are taken per dry mass: rho in the SE tab stays `E2`, so the ng/g results (H-J) are unchanged by
water, and only the per-atom ppb (E-G) and the Target Fractions change. As shot, dry-mix shotcrete holds 8.6-10.8 % water (*Shotcrete*, Winter 2009, pp. 22-27); bound
water is about 0.23 g per g of cement (about 3.5-6 % of the dry mix, 0.4-0.6 % H), against H = 0.04 % in Table 3 of the SNOLAB report. A measurement
on a core (oven-dry at 105 C, ignite at 1000 C) would settle it. The references are also at the bottom of the tab (A34:A36).

### Target Fractions tab (6)

For each element (Z, fraction of atoms), the isotope (A) that acts as an (a,n) target, its lower and upper
abundance bounds (%), and `Middle`, `Lower`, `Upper` as fractions of atoms: `Middle = B x (lower+upper)/200`,
`Lower = B x lower/100`, `Upper = B x upper/100`. Shotcrete bounds follow CIAAW (e.g. C-13 0.96-1.16 %).

## Recomputing by hand

- Equilibrium ratio (for example Pa-231 in the U-235 chain): `T(Pa-231) / T(U-235)` x 1 (every decay of U-235 goes
  through Th-231 to Pa-231) = 3.257e4 / 7.04e8 = 4.626e-5, matching column C.
- ng/g head: `1e-9 x 2.2222 g/cm^3 x 6.0221408e23 / 232.038 g/mol = 5.77e12` Th-232 atoms/cm^3 in shotcrete.

## Status of each tab (audit of commit 65e8d6e, still valid)

OK: Th232, U235, U238, SE (1ppb). Open (works, with a question): Shotcrete Density (the density is
only the datasheet's approximate yield, see above). Shotcrete Target Fractions is OK. No tab has a definite error.

## Known issues and caveats

- Column D of the SE tab is the old Mathematica output; its At-219, At-215 and Po-211 rows disagree with column C
  because of input errors in the original notebooks (listed in `mathematica_SE/README.md`). The corrected notebooks
  were re-run in Mathematica; `U235_SE_ODE.nb` now agrees with column C to within 0.55 %. Column D was kept as it is.
- Downstream SOURCES inputs and results (shotcrete, norite) that used the old U-235-chain values need regenerating
  from column C and the ng/g head densities.
- The workbook is set to recalculate on open.
