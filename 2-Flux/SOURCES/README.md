# SOURCES: neutron yields from decay chains

[SOURCES 4A](https://www.osti.gov/biblio/15215) computes the neutron production rates and spectra from (α,n)
reactions and spontaneous fission for a given concentration of U or Th chain contamination in a material. We use it
to turn the radioactivity of SNOLAB's shotcrete into neutron source terms for the flux calculations in `2-Flux/`, and
the same chains in the Cu, HDPE and Pb shielding materials. The norite rock has its own inputs and results in
[`norite/`](norite/README.md); `C-Slab-Source-and-Shielding` still applies the shotcrete source spectrum with norite's
composition for the attenuation.

## What is in this repository

| Path | Contents |
|---|---|
| `inputs/` | SOURCES input files (`tape1-*`) for 1 ppb of each chain: `shotcrete/`, `cu/`, `hdpe/`, `pb/`, `hang/`, `pube/`, plus small `examples/`. `archive/` folders hold older versions |
| `results/` | Saved SOURCES outputs (`_tape7` spectra and `*_summary_*.txt`) per material: `shotcrete/`, `cu/`, `hdpe/`, `hang/`, `pube/`, `supersim/`. Dated and `_old` copies are in `archive/` |
| `analysis/` | `decay_chains_data.xlsx` ([explained](analysis/decay_chains_data_README.md)), `shielding.xlsx`, `mathematica_SE/` and the `stanzas2columns*.sh` helpers |
| `notebooks/` | Notebooks that read and plot SOURCES results (`nOutput-SOURCES4A-Sheffield*.ipynb`, `nShielding.ipynb`) |
| `docs/` | SOURCES 4A and 4C manuals, `SOURCES-inputs.pdf`/`.pptx`, `shotcrete_comp.png` |
| `norite/` | Inputs, results and script for the norite rock (self-contained, see its README) |

Everything needed to reproduce the figures from the saved outputs is here. Notebook `../B-SOURCES4A-Running.ipynb`
explains the inputs and reads the saved results; it does not need SOURCES installed.

## What is not in this repository

The SOURCES code and its data libraries (`tape2` to `tape5`) are **not distributed here**. SOURCES is a Los Alamos
code distributed through the Oak Ridge RSICC package, and the libraries we use were extended by V. A. Kudryavtsev's
group at the University of Sheffield and are supplied to people who already hold the RSICC package. To rerun the
calculations you need to obtain them yourself:

1. Request SOURCES from the Oak Ridge RSICC code package collection.
2. Ask the Sheffield group (v.kudryavtsev@sheffield.ac.uk) for their modified code and libraries.

## Checking you have the same version

The results here were produced with the files whose SHA-256 checksums are listed in
[`scripts/get_sources.py`](../../scripts/get_sources.py):

```
src/sources_v20240312.f     SOURCES 4A source, modified (Sheffield, 12 Mar 2024 version)
libs/tape2 ... libs/tape5   libraries (tape3 as corrected on 9 Jun 2025)
```

Put them in `2-Flux/SOURCES/code/` (ignored by git), or point `SOURCES_DIR` at them, and run:

```bash
python scripts/get_sources.py
```

It reports each file as `ok`, `missing` or `different`. A different version can give different yields, particularly
for carbon-13 and above 6.5 MeV, so note any difference when comparing with our results.

## Running a problem

Build the Fortran source with `gfortran -O0 -std=legacy -fallow-argument-mismatch`. **Do not optimize**: at `-O2` the code gives wrong results on problems with many target nuclides (on the shotcrete problem the summed neutron rate is wrong by a factor of about 5e8), while simple problems still look fine. Check a new build by rerunning one of the inputs here and comparing `_tape7` with the saved result in `results/shotcrete/` (the `-250609` files, in `results/shotcrete/archive/` and as `*_summary_1ppb_250609.txt`, are the June 2025 results). Also raise `maxELevPts` from 200 to 500 in the source (`sed 's/maxELevPts = 200/maxELevPts = 500/'`): `tape4` holds branching tables of up to 434 points, which overflow the stock 200-element array (details in [`norite/README.md`](norite/README.md)). SOURCES reads
`tape2` to `tape5` and its inputs from the current directory. List the inputs in a file named `filelist` (first line
the number of files, then one name per line, no spaces), then run the executable. Each input `tape1-X` produces
`tape1-X_outp` (a summary of what was used), `_tape6`, `_tape7` (neutron spectra), `_tape8` and `_tape9`.

Input notes:

- Material fractions are **atom** fractions and should sum to 1.
- Alphas above 6.5 MeV need the extended libraries.
- Cross-section and branching-ratio arrays can be chosen per target in `tape1`; the `_outp` file records which were used.

## Citing

If you use results from SOURCES with these libraries, cite:

- M. Parvu, P. Krawczun and V. A. Kudryavtsev, Appl. Radiat. Isot. **225**, 112035 (2025), arXiv:2408.10910
- V. A. Kudryavtsev, P. Krawczun and R. Bocheva, AIP Conf. Proc. **2908**, 100003 (2023), arXiv:2211.02080
- V. A. Kudryavtsev, P. Zakhary and B. Easeman, Nucl. Instrum. Meth. A **972**, 164095 (2020), arXiv:2005.02499
- W. B. Wilson et al., the SOURCES 4A manual (Los Alamos), linked above


## Where did files move? (reorganisation of 2026-11-11)

The git tag `pre-reorg` marks the layout before the move; `git log --follow <file>` follows a file across it.

| Old location (under `2-Flux/`) | New location (under `2-Flux/SOURCES/`) |
|---|---|
| `SOURCES/SNOLAB_inputs/*`, `SNOLAB_Flux/SOURCES_inputs/tape1-*-SNOLAB-1ppb` | `inputs/shotcrete/` |
| `SOURCES/sources_running/data/raw_sources_IO/tape1-*-SNOLAB-1ppb` (older versions, different content) | `inputs/shotcrete/archive/` |
| `SNOLAB_Flux/SOURCES_inputs/tape1-*-SNOLAB-{Cu,HDPE,Pb}-1ppb` | `inputs/{cu,hdpe,pb}/` |
| `SOURCES/tape1-Th-C*`, `SOURCES/tape1-U-C*`, `SOURCES/filelist` | `inputs/examples/` |
| `SOURCES/.../raw_sources_IO/tape1-Th232-Hang-1ppb`, `tape1-PuBe` | `inputs/hang/`, `inputs/pube/` |
| `SNOLAB_Flux/SOURCES_outputs/*_tape7`, `SOURCES/.../raw_sources_IO/*_summary_1ppb*.txt` | `results/shotcrete/`, `results/cu/`, `results/hdpe/` |
| `SNOLAB_Flux/data/*_summary_1ppb_*.txt` | `results/shotcrete/` (`_250609`), `results/cu/` (`_260130`), `results/hdpe/` (`_260115`); `_260116` in `results/cu/archive/` |
| `.../raw_sources_IO/-250609_tape7`, `SOURCES/.../data/tape1-Th232-SNOLAB_tape7` | `results/shotcrete/archive/` |
| `.../raw_sources_IO/Hang_*`, `*Hang*`, `PuBe*`, `SuperSim_*` | `results/hang/`, `results/pube/`, `results/supersim/` (`_old` in `archive/`) |
| `.../raw_sources_IO/decay_chains_data.xlsx`, `shielding.xlsx`, `mathematica_SE/`, `stanzas2columns*.sh`, `test.txt` | `analysis/` |
| `SOURCES/sources_running/*.ipynb` | `notebooks/` |
| `SNOLAB_Flux/SOURCES_manuals/`, `SNOLAB_Flux/SOURCES_images/`, `SOURCES/SOURCES-inputs.*` | `docs/` |

Nine files that were byte-identical to another copy were deleted (3 `_tape7`, 3 `summary_1ppb_250609`, 3
`stanzas2columns*.sh`); the kept copy is the one in the table above.
