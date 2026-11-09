#!/usr/bin/env python3
"""SOURCES inputs for SNOLAB norite, and a summary of SOURCES output.

  python norite_sources.py inputs  --out DIR     write tape1-{Th232,U238,U235}-SNOLAB-Norite-1ppb and composition.csv
  python norite_sources.py summary TAPE7 OUT     write the (a,n), S.F. and total spectra of a tape7 as a table

Composition. Weight-percent oxides of norite orthopyroxene and plagioclase are the GEOROC means used in
2-Flux/C-Slab-Source-and-Shielding.ipynb (they reproduce from saved_data/2024-12-SGFTFN_*.csv by averaging the
norite rows, gabbronorite excluded). Iron is counted ONCE, as FeOT (total iron expressed as FeO); the notebook also
summed FeO and Fe2O3, which double-counts it. The rock is a 50:50 mix of the two minerals by volume (densities 3.3 and
2.7 g/cm3, giving 3.0 g/cm3), as in the notebook.

Source terms. SOURCES takes the atom density (atoms/cm3) of each alpha-emitting nuclide (manual, card 7). The chains
are written for 1 ppb BY MASS of the head nuclide (Th-232, U-238, U-235) in rock of density 3.0 g/cm3, in secular
equilibrium, by scaling the chain ratios in the existing shotcrete inputs. Multiply the results by the measured
concentration in ppb. The matrix density cancels in the neutron yield per alpha, so only the emitter density matters.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

NA = 6.02214076e23
RHO_NORITE = 3.0  # g/cm3, 50:50 mix of orthopyroxene (3.3) and plagioclase (2.7)
SHOTCRETE_INPUTS = Path(__file__).resolve().parents[2] / "SNOLAB_Flux" / "SOURCES_inputs"

# oxide: (wt%, molar mass g/mol, O atoms per formula unit, metal atoms per formula unit, element)
ORTHOPYROXENE = {
    "SiO2": (52.84, 60.084, 2, 1, "Si"), "TiO2": (0.24, 79.866, 2, 1, "Ti"), "Al2O3": (2.03, 101.961, 3, 2, "Al"),
    "Cr2O3": (0.20, 151.99, 3, 2, "Cr"), "FeOT": (20.87, 71.84, 1, 1, "Fe"), "CaO": (1.3, 56.08, 1, 1, "Ca"),
    "MgO": (25.0, 40.305, 1, 1, "Mg"), "MnO": (0.34, 70.937, 1, 1, "Mn"), "NiO": (0.04, 74.693, 1, 1, "Ni"),
    "K2O": (0.015, 94.196, 1, 2, "K"), "Na2O": (0.033, 61.979, 1, 2, "Na"), "P2O5": (0.01, 141.945, 5, 2, "P"),
}
PLAGIOCLASE = {
    "SiO2": (50.94, 60.084, 2, 1, "Si"), "TiO2": (0.01, 79.866, 2, 1, "Ti"), "Al2O3": (31.27, 101.961, 3, 2, "Al"),
    "FeOT": (0.192, 71.84, 1, 1, "Fe"), "CaO": (14.11, 56.08, 1, 1, "Ca"), "MgO": (0.13, 40.305, 1, 1, "Mg"),
    "MnO": (0.01, 70.937, 1, 1, "Mn"), "NiO": (0.05, 74.693, 1, 1, "Ni"), "K2O": (0.12, 94.196, 1, 2, "K"),
    "Na2O": (3.26, 61.979, 1, 2, "Na"),
}
MINERALS = ((ORTHOPYROXENE, 3.3), (PLAGIOCLASE, 2.7))
ATOMIC_MASS = {"O": 15.999, "Si": 28.086, "Ti": 47.88, "Al": 26.982, "Cr": 51.996, "Fe": 55.845, "Ca": 40.078,
               "Mg": 24.305, "Mn": 54.938, "Ni": 58.693, "K": 39.096, "Na": 22.990, "P": 30.974}
Z = {"O": 8, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "K": 19, "Ca": 20, "Ti": 22, "Cr": 24, "Mn": 25,
     "Fe": 26, "Ni": 28}

# (alpha,n) targets with data in the modified tape3: ZAID -> natural isotopic abundance (IUPAC)
TARGETS = {
    80170: ("O", 0.00038), 80180: ("O", 0.00205), 110230: ("Na", 1.0),
    120240: ("Mg", 0.7899), 120250: ("Mg", 0.1000), 120260: ("Mg", 0.1101), 130270: ("Al", 1.0),
    140280: ("Si", 0.92223), 140290: ("Si", 0.04685), 140300: ("Si", 0.03092), 150310: ("P", 1.0),
    190390: ("K", 0.932581), 190400: ("K", 0.000117), 190410: ("K", 0.067302),
    200420: ("Ca", 0.00647), 200430: ("Ca", 0.00135), 200440: ("Ca", 0.02086), 200460: ("Ca", 0.00004),
    200480: ("Ca", 0.00187),
    220460: ("Ti", 0.0825), 220470: ("Ti", 0.0744), 220480: ("Ti", 0.7372), 220490: ("Ti", 0.0541),
    220500: ("Ti", 0.0518),
    240500: ("Cr", 0.04345), 240520: ("Cr", 0.83789), 240530: ("Cr", 0.09501), 240540: ("Cr", 0.02365),
    250550: ("Mn", 1.0),
    260540: ("Fe", 0.05845), 260560: ("Fe", 0.91754), 260570: ("Fe", 0.02119), 260580: ("Fe", 0.00282),
    280600: ("Ni", 0.26223), 280610: ("Ni", 0.01140), 280620: ("Ni", 0.03635), 280640: ("Ni", 0.00926),
}
# chain name -> (head ZAID, head molar mass, name used in titles)
CHAINS = {"Th232": (902320, 232.0381, "Th-232"), "U238": (922380, 238.0508, "U-238"),
          "U235": (922350, 235.0439, "U-235")}


def composition():
    """Element atom densities (atoms/cm3) of the 50:50 mineral mix, with iron counted once."""
    nd = {}
    for minerals, dens in MINERALS:
        for wt, mm, n_o, n_m, el in minerals.values():
            per_oxide = (wt / 100.0) / mm * dens * NA
            nd["O"] = nd.get("O", 0.0) + 0.5 * per_oxide * n_o
            nd[el] = nd.get(el, 0.0) + 0.5 * per_oxide * n_m
    total = sum(nd.values())
    frac = {el: n / total for el, n in nd.items()}
    mass = sum(nd[el] * ATOMIC_MASS[el] for el in nd)
    mass_frac = {el: 100.0 * nd[el] * ATOMIC_MASS[el] / mass for el in nd}
    return nd, frac, mass_frac, total, mass / NA


def read_chain(path):
    """Return (list of (ZAID, density), head ZAID, head density) from an existing shotcrete input."""
    lines = [l.strip() for l in Path(path).read_text().replace("\r", "").splitlines()]
    i = next(k for k, l in enumerate(lines) if re.fullmatch(r"100 +10\.0 +0\.0", l))
    nq = int(lines[i + 1].split()[0])
    sources = [(int(l.split()[0]), float(l.split()[1])) for l in lines[i + 2:i + 2 + nq]]
    return sources


def make_input(chain, frac):
    head, amass, label = CHAINS[chain]
    template = read_chain(SHOTCRETE_INPUTS / f"tape1-{chain}-SNOLAB-1ppb")
    head_shotcrete = dict(template)[head]
    head_norite = 1e-9 * RHO_NORITE * NA / amass  # atoms/cm3 for 1 ppb by mass
    scale = head_norite / head_shotcrete
    elements = sorted(frac, key=lambda el: Z[el])
    out = [f"Bulk Norite problem 1ppb {label} in secular equilibrium",
           "1 2 1",
           f"{len(elements)} 0"]
    out += [f" {Z[el]}  {frac[el]:.9f}" for el in elements]
    out += ["100 10.0 0.0", f"{len(template)}"]
    out += [f" {zaid}  {n * scale:.6E}" for zaid, n in template]
    out += [f"{len(TARGETS)} 4000"]
    out += [f" {zaid}  {frac[el] * ab:.5E}" for zaid, (el, ab) in sorted(TARGETS.items())]
    return "\n".join(out) + "\n", head_norite


def cmd_inputs(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    nd, frac, mass_frac, total, rho = composition()
    for chain in CHAINS:
        text, head = make_input(chain, frac)
        (outdir / f"tape1-{chain}-SNOLAB-Norite-1ppb").write_text(text, newline="\n")
        print(f"wrote tape1-{chain}-SNOLAB-Norite-1ppb  (1 ppb {chain} head = {head:.4e} atoms/cm3, {len(TARGETS)} targets)")
    with open(outdir / "composition.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["element", "Z", "number_density_per_cm3", "atom_fraction", "mass_percent"])
        for el in sorted(nd, key=lambda e: Z[e]):
            w.writerow([el, Z[el], f"{nd[el]:.5e}", f"{frac[el]:.6f}", f"{mass_frac[el]:.3f}"])
    print(f"total atom density {total:.4e} /cm3 (implied mass density {rho:.3f} g/cm3; nominal {RHO_NORITE})")


def _numbers(lines, start):
    vals = []
    for l in lines[start + 1:]:
        if re.fullmatch(r"\s*([-+0-9.EeD]+\s*)+", l) and l.strip():
            vals += [float(x) for x in l.split()]
        elif vals:
            break
    return vals


def read_spectra(tape7):
    L = Path(tape7).read_text(errors="replace").replace("\r", "").splitlines()
    find = lambda pat: [i for i, l in enumerate(L) if re.search(pat, l)]
    energy = _numbers(L, find(r"Neutron Multigroup Structure")[0])
    an = _numbers(L, find(r"Grand total \(alpha,n\)")[-1])
    total = _numbers(L, find(r"Total Neutron Spectrum")[-1])
    sf_idx = find(r"^\s*S\.F\. neutrons/cc from")
    sf = [0.0] * len(energy)
    for i in sf_idx:
        sf = [a + b for a, b in zip(sf, _numbers(L, i))]
    return energy, an, sf, total


def cmd_summary(tape7, out):
    energy, an, sf, total = read_spectra(tape7)
    with open(out, "w", newline="") as f:
        f.write("Energy\t   (a,n)\t   S.F.\t   tot\n")
        for row in zip(energy, an, sf, total):
            f.write(" ".join(f"{v:.4e}" for v in row) + "\n")
    n = min(len(energy), len(total))
    print(f"wrote {out}: {n} groups; totals (a,n) {sum(an):.4e}, S.F. {sum(sf):.4e}, all {sum(total):.4e} n/s/cm3")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("inputs"); p.add_argument("--out", required=True)
    p = sub.add_parser("summary"); p.add_argument("tape7"); p.add_argument("out")
    a = ap.parse_args()
    cmd_inputs(a.out) if a.cmd == "inputs" else cmd_summary(a.tape7, a.out)


if __name__ == "__main__":
    sys.exit(main())
