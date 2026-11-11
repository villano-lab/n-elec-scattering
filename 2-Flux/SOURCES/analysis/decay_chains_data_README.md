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
SE (1ppb)  /  Cu SE (1ppb) column C   N_i / N_seed at secular equilibrium (live formula)
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
| 6 | `Shotcrete Hang` | Hang's published atom fractions, renormalised, compared with the atom fractions of tab 5 (difference in %). |
| 7 | `Shotcrete Target Fractions` | For each element, the (a,n) target isotopes with CIAAW abundance bounds and the middle value (fraction of atoms). |
| 8 | `HDPE Density` | C2H4 composition from the formula, 961 kg/m^3, number densities. |
| 9 | `HDPE Target Fractions` | Same as tab 7 for HDPE: C-13 (0.96-1.16 %) is the only (a,n) target; H has none. |
| 10 | `Cu SE (1ppb)` | Same layout as tab 4, using the `Cu Density` tab. |
| 11 | `Cu Density` | Cu, Z = 29, 8.96 g/cm^3, 8.49e22 atoms/cm^3. |
| 12 | `Cu Target Fractions` | Cu-63 (69.15 %) and Cu-65 (30.85 %), bounds +/-15 % (a choice, not a measurement). |
| 13 | `Pb Density` | Pb, Z = 82, 11.34 g/cm^3, 207.2 g/mol, 3.296e22 atoms/cm^3. |
| 14 | `Pb Target Fractions` | Pb-204/206/207/208 with bounds (a choice, wider than CIAAW ranges). |

### Decay tabs (1-3)

- Column B is the half-life in years. Where the source gives seconds, minutes or days the cell holds the conversion
  (for example `=3.6319/365`), using a 365-day year. This differs from 365.2422 days by a uniform 0.07 %.
- C = B/ln 2 (mean life), D = 1/C (lambda), E and F are the same scaled to the seed (G2 is the rounded seed lambda).
- H gives the isotope name. All half-lives were checked against ENSDF (within 0.25 %), and the branching tables
  against ENSDF.
- Branching table columns are "beta (%)" and "alpha (%)" (in that order). At-218 in the U-238 chain has no row; the
  SE tabs assume 100 % alpha for it (ENSDF: 99.95 %).

### SE tabs (4 and 10)

Three blocks per tab, one per seed: **Th-232** (header row 1, data rows 2-8), **U-238** (header row 10, data rows
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

Two definitions of "ppb" exist side by side:

- **per-atom ppb** (E-G): 1 ppb = 1e-9 atoms of the seed per atom of matrix.
- **ng/g ppb** (H-J): 1 ppb = 1e-9 g of the seed per g of matrix. Seed atoms/cm^3 = `1e-9 x rho x N_A / A_seed`;
  for 1 ng/g in shotcrete that is 5.77e12 (Th-232) and 5.62e12 (U-238) atoms/cm^3. The two definitions differ by
  about 11x for shotcrete (mean matrix mass about 21 g/mol).

U-235 is tied to U-238 by the natural abundance 0.007204 (by atoms), in both definitions.

Rows 40-46 of the sheet hold the notes naming the exact files and repositories used.

### Density tabs

Mass fraction (C) and overall density (D, kg/m^3; E in g/cm^3) give the elemental mass density (F = E x C/100),
the number density (H = F/G x N_A, with N_A = 6.0221408E+23) and the atom fractions (I = H / total). The atom
fractions are copied by value into the lower table (rows 21+ for shotcrete) that the Target Fractions tabs use.
The shotcrete composition sums to 99.66 % (the Sika datasheet cited does not give composition or density; the other
reference is an inaccessible Confluence page, so this tab is an open question).

### Target Fractions tabs (7, 9, 12, 14)

For each element (Z, fraction of atoms), the isotope (A) that acts as an (a,n) target, its lower and upper
abundance bounds (%), and `Middle`, `Lower`, `Upper` as fractions of atoms: `Middle = B x (lower+upper)/200`,
`Lower = B x lower/100`, `Upper = B x upper/100`. Shotcrete bounds follow CIAAW (e.g. C-13 0.96-1.16 %); the Cu
and Pb bounds are chosen ranges. `Shotcrete Target Fractions` column I (Hang's normalisation) multiplies by the
total of the Shotcrete Hang tab (0.0474); that column only makes sense for shotcrete and was removed from the HDPE
tab.

## Recomputing by hand

- Equilibrium ratio (for example Pa-231 in the U-235 chain): `T(Pa-231) / T(U-235)` x 1 (every decay of U-235 goes
  through Th-231 to Pa-231) = 3.257e4 / 7.04e8 = 4.626e-5, matching column C.
- ng/g head: `1e-9 x 2.2222 g/cm^3 x 6.0221408e23 / 232.038 g/mol = 5.77e12` Th-232 atoms/cm^3 in shotcrete.

## Status of each tab (audit of commit 65e8d6e, still valid)

OK: Th232, U235, U238, SE (1ppb), Cu SE (1ppb), HDPE Density, HDPE Target Fractions, Cu Density, Cu Target
Fractions, Pb Density, Pb Target Fractions. Open (works, with a question): Shotcrete Density (source of
composition/density), Shotcrete Hang (meaning of 0.0474), Shotcrete Target Fractions (column I). No tab has a
definite error.

## Known issues and caveats

- Column D of the SE tabs is the old Mathematica output; its At-219, At-215 and Po-211 rows disagree with column C
  because of input errors in the original notebooks (listed in `mathematica_SE/README.md`). The corrected notebooks
  were re-run in Mathematica; `U235_SE_ODE.nb` now agrees with column C to within 0.55 %. Column D was kept as it is.
- Downstream SOURCES inputs and results (shotcrete, norite) that used the old U-235-chain values need regenerating
  from column C and the ng/g head densities.
- The workbook is set to recalculate on open.
