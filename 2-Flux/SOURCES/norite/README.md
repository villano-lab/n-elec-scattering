# SOURCES for norite

SOURCES 4A inputs and results for the radioactivity of the norite rock around SNOLAB: the Th-232, U-238 and U-235
chains, each for 1 part per billion (by mass) of the head nuclide.

| Path | Contents |
|---|---|
| `inputs/tape1-{Th232,U238,U235}-SNOLAB-Norite-1ppb` | SOURCES inputs (homogeneous problem, neutron spectra, 100 groups from 0 to 10 MeV) |
| `inputs/composition.csv` | element number densities, atom fractions and mass percent used in the inputs |
| `outputs/tape1-*-SNOLAB-Norite-1ppb_tape7` | SOURCES neutron spectra |
| `outputs/*_Norite_summary_1ppb.txt` | the (α,n), spontaneous-fission and total spectra per energy group, all 100 groups |
| `norite_sources.py` | regenerates the inputs and the summaries |

```bash
python norite_sources.py inputs --out inputs
python norite_sources.py summary outputs/tape1-Th232-SNOLAB-Norite-1ppb_tape7 outputs/Th232_Norite_summary_1ppb.txt
```

## Composition

The rock is a 50:50 mix by volume of norite orthopyroxene (3.3 g/cm³) and plagioclase (2.7 g/cm³), density
3.0 g/cm³, with the weight-percent oxides used in `../../C-Slab-Source-and-Shielding.ipynb` (GEOROC means over norite
samples; they reproduce from `saved_data/2024-12-SGFTFN_*.csv`). Compared with that notebook, **iron is counted once**,
as FeOT (total iron expressed as FeO). The notebook also added the FeO and Fe2O3 entries, which double-counts iron and
made its mass fractions sum to 104%. The fractions here are normalized to 1.

| Element | Atom fraction | Mass % |
|---|---|---|
| O | 0.6058 | 44.71 |
| Si | 0.1845 | 23.91 |
| Mg | 0.0731 | 8.19 |
| Al | 0.0635 | 7.91 |
| Fe | 0.0343 | 8.84 |
| Ca | 0.0269 | 4.97 |
| Na | 0.0102 | 1.08 |
| Mn, Ti, Cr, K, Ni, P | 0.0006 to 0.00002 each | 0.15 to 0.002 |

The total atom density is 8.47×10²² cm⁻³ (3.05 g/cm³ from the composition; the inputs use 3.0).

## Source terms

SOURCES takes the atom density of each alpha-emitting nuclide (atoms/cm³, input card 7); the matrix density cancels in
the neutron yield per alpha. Each input is for **1 ppb by mass of the head nuclide in rock of 3.0 g/cm³**, in secular
equilibrium: 7.79×10¹² atoms/cm³ of Th-232, 7.59×10¹² of U-238 or 7.69×10¹² of U-235, with the other chain members
scaled from the ratios in the existing shotcrete inputs. Multiply a result by the concentration in ppb.

The 37 (α,n) targets are every isotope of the norite elements that has data in the modified `tape3` and `tape4`.
Natural isotopic abundances are used, and Ni-58 and Ca-40 have no (α,n) data.

## Results

Neutron production per 1 ppb of the head nuclide (neutrons per second per cm³ of rock):

| Chain | (α,n) | Spontaneous fission |
|---|---|---|
| Th-232 | 5.54×10⁻¹¹ | 3.7×10⁻¹⁶ |
| U-238 | 1.13×10⁻¹⁰ | 4.1×10⁻¹¹ |
| U-235 | 2.01×10⁻⁹ | 9.0×10⁻¹³ |

The mean neutron energy is 2.08 MeV for the Th-232 chain and 1.96 MeV for natural uranium (group midpoints of the
total spectrum).

Natural uranium is 99.29% U-238 and 0.711% U-235 by mass, so 1 ppb of natural U gives 1.67×10⁻¹⁰ n s⁻¹ cm⁻³.
With the rock averages of SNOLAB-STR-2007-003 (Table 1: U 1.11 ± 0.13 ppm from ²²⁶Ra, Th 5.56 ± 0.52 ppm from ²²⁸Th):

| | n s⁻¹ cm⁻³ |
|---|---|
| U chains (1110 ppb) | 1.86×10⁻⁷ |
| Th chain (5560 ppb) | 3.08×10⁻⁷ |
| Total | 4.94×10⁻⁷, or 1.6×10⁻⁴ n s⁻¹ kg⁻¹ |

## How these were run

- SOURCES 4A, Sheffield version (`sources_v20240312.f`, libraries of 5 Feb 2025 with the corrected `tape3` of
  9 Jun 2025), built with `gfortran -O0` (see `../README.md`).
- **`maxELevPts` must be raised from 200 to 500.** `tape4` holds branching tables with up to 434 energy points (Ca-42
  has 434, Cu-65 383, Mo-95 420), and the stock code reads them into a 200-element array. Built with run-time bounds
  checking, the stock code stops with "Index '201' of dimension 1 of array 'ep' above upper bound of 200" when it reaches
  Ca-42; without the check it writes past the array silently. Every run containing Ca-42 is affected, including the
  shotcrete runs; for the shotcrete Th-232 problem the fix changes only the Ca-42 blocks (0.4% and 9.5% in two of
  them) and the total rate by 0.01%. The patch is one line: `sed 's/maxELevPts = 200/maxELevPts = 500/'`.
- With the patch, the 37-target norite problem and the 24-target shotcrete problem run without bounds errors, and the
  bounds-checked and production builds agree exactly.

## Checks and limits

- The ratios between the chains match the existing shotcrete runs (U-235/U-238 is 17.7 here and 17.5 for shotcrete),
  as expected since the chain structure is taken from them.
- 1 ppm of natural U gives about 5.6×10⁻⁸ n s⁻¹ g⁻¹. A rough check: 1 ppm of U emits about 0.1 alpha s⁻¹ g⁻¹ across the
  chain, and an (α,n) yield of about 5×10⁻⁷ neutron per alpha in silicate rock gives the same order.
- The 2013 CDMS norite spectra (`SnoLab_Norite_fast_neutron_{U,Th}.dat` in `HVeV-n-Bknd`) are in arbitrary units and
  much softer: their mean energy is about 1.0 MeV against 2.0 MeV here. The 2013 shotcrete spectra differ from the
  shotcrete results in this repository in the same way (1.09 against 1.80 MeV), so the difference is not specific to
  norite. It is consistent with the older code and libraries, which cap the alpha energy at 6.5 MeV.
- The composition is a global GEOROC norite average, not Sudbury rock. SNOLAB-STR-2007-003 also reports elemental
  analyses of two SNOLAB rock samples that could replace it.
- Densities (3.3, 2.7 g/cm³) and the 50:50 mineral mix are rough, as in the notebook.
