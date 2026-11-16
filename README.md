# Mandacaru PAW-LCAO pseudopotentials, version 3

This repository holds the **pseudopotential datasets** used by
[Mandacaru](https://github.com/seixas-research/mandacaru), a framework for
simulating molecules and crystals with variational quantum algorithms
(VQE, ADAPT-VQE) and with density functional theory. It has one file per
element for every element from hydrogen to uranium (Z ≤ 92), in four sets.

Version 3 replaces the version-2 library (`mandacaru-paw-v2`). Every file
was regenerated from scratch in November 2026 by Mandacaru's own generator,
which is written entirely in Rust. Each file passed a fixed list of checks
before it was written. Nothing is converted from another code.

## What changed from version 2

- **Core states in every file (format 3).** Each record stores the orbitals
  and occupations of its frozen core. Hybrid functionals (HSE06) need them for
  the core-valence exchange, and so does Hartree–Fock. Version-2 files had no
  core states.
- **One exchange-correlation scheme.** The reference atom, the unscreening
  of the local potential and the one-center energies all use the same
  functional derivative. Version 2 mixed two schemes across the library.
- **Spin-orbit coupling calibrated against the Dirac atom.** In `lda-fr/` and
  `pbe-fr/`, the spin-orbit term of each file is scaled so that Mandacaru's
  own isolated atom reproduces the all-electron Dirac j-splitting.
  - Version 2 had no such calibration. Without it, the heavy p shells come
    out 8–14 % short: Bi 6p gives 1.607 eV against 1.836 eV.
  - In version 3, `pbe-fr/Bi` gives 1.910 eV against 1.907 eV for its
    Dirac atom (+0.2 %).
  - See [Spin-orbit calibration](#spin-orbit-calibration).
- **Exact Coulomb tail of the local potential.** Beyond the local radius,
  the ionic potential is now exactly −Z_v/r.
  - A new check, G14, refuses any file whose tail is off.
  - This closes a defect that the first fr files of the regeneration had.
    In a crystal, the defect acts as a constant shift per atom: on Bi2Se3 with
    `pbe-fr/`, the stress was −90 GPa off the `pbe-sr/` result and the
    energy +105 eV off.
  - See [The fr tail defect](#the-fr-tail-defect).
- **Recipe updates.** All four sets follow the same recipe tables; a single
  set of tables, keyed by Z, functional and relativistic treatment, covers
  the whole library.
  - Cs, Ba, Fr, Ra and the lanthanides carry their semicore s and p in the
    valence. The lanthanides Ce–Yb use a 4f(n−1) 5d1 6s2 reference.
  - Br, Kr, Te, I and Xe keep a frozen semicore d.
  - The Se d sphere is 2.3 Bohr, with the 3d in the valence.
- **Checks in Rust, applied on release.** Every check runs on the file as
  it is read back, before the file is moved into place. A file that fails
  is never written. See [Checks on every file](#checks-on-every-file).

## Why Mandacaru needs them

Mandacaru computes its Hamiltonian on a real-space grid. A grid fine enough
for valence electrons (0.15–0.30 Å) cannot resolve the core electrons, whose
orbitals vary on a scale of hundredths of an ångström for oxygen and heavier
atoms. Removing the core and replacing the nuclear potential by a smooth one
solves this, and it also shrinks the problem: only valence electrons are left
to put on qubits.

**PAW-LCAO** is the projector augmented-wave (PAW) method of Blöchl in a form
built for a localized basis. The smooth pseudo-atomic orbitals of each
dataset are also the **basis functions** of the calculation. PAW keeps the
information needed to reconstruct the true all-electron orbitals near each
nucleus. It also adds to the Hamiltonian the correction terms that this
reconstruction implies.

## Contents

| Folder | Elements | What it is | Size |
|---|---|---|---|
| `lda-sr/` | 92 (H–U) | scalar-relativistic, LDA — **the default** | ~375 MB |
| `lda-fr/` | 92 (H–U) | the same, plus spin-orbit coupling | ~455 MB |
| `pbe-sr/` | 92 (H–U) | scalar-relativistic, PBE | ~375 MB |
| `pbe-fr/` | 92 (H–U) | the same, plus spin-orbit coupling | ~465 MB |

- **LDA** here is Slater exchange with Perdew–Zunger correlation, with the
  relativistic correction to exchange of MacDonald and Vosko.
- **PBE** is the functional of Perdew, Burke and Ernzerhof, with
  Perdew–Wang correlation and the same relativistic correction to exchange.
- **Scalar-relativistic** (`-sr`) means the reference atoms include
  mass-velocity and Darwin effects, but not spin-orbit coupling.
- **Fully relativistic** (`-fr`) adds spin-orbit coupling from the Dirac
  equation.
- **HSE06** and the other hybrids run on `pbe-sr/` (or `pbe-fr/`). No set is
  generated with a hybrid: the dataset functional moves an HSE06 result by
  at most 0.19 mHa and 16 meV (H2O, SiH4), less than going from LDA to PBE
  does.

## Quick start

The files are format 3, so they need a Mandacaru build that reads format 3.
Older builds refuse them. Install Mandacaru, then clone this repository
and tell Mandacaru where it is:

```bash
git clone https://github.com/seixas-research/mandacaru-paw-v3
mandacaru-paw --set mandacaru-paw-v3   # writes MANDACARU_PAW_PATH to ~/.zshrc or ~/.bashrc
# open a new terminal (or `source` that file), then check:
mandacaru-paw --status
```

PAW-LCAO is chosen as the **basis** of a calculation:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)          # the cell is the real-space box
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis={"name": "PAW-LCAO", "size": "DZP"},
                       grid_spacing=0.25)    # Angstrom
energy = atoms.get_potential_energy()
```

A crystal with spin-orbit coupling and a van der Waals functional:

```python
from ase.io import read

atoms = read("Bi2Se3.cif")
atoms.calc = Mandacaru(method="dft", xc="pbe+rvv10",
                       basis={"name": "PAW-LCAO", "size": "TZP"},
                       paw_dataset="pbe-fr",
                       kpts={"size": (3, 3, 3), "gamma": True},
                       grid_spacing=0.15)
```

If `MANDACARU_PAW_PATH` is not set, a PAW-LCAO calculation stops before it
starts with a `LibraryPathError`, which names the command above.

## Options

**Basis size.** `size` picks how many atomic orbitals each atom contributes.

| `size` | Orbitals per valence shell |
|---|---|
| `"SZ"` (default) | single zeta: one per occupied valence orbital |
| `"DZ"` | double zeta: two |
| `"DZP"` | double zeta plus polarization functions |
| `"TZP"` | triple zeta plus polarization functions |

`SZP`, `DZ2P`, `TZ`, `TZ2P`, `QZ`, `QZP` and `QZ2P` follow the same pattern.
A larger basis is more accurate and needs more qubits. For forces and
stresses in crystals, use DZP at least. In Bi2Se3, going from DZP on a
0.18 Å grid to TZP on a 0.15 Å grid moved the relaxed c by 3.5 %.

**Confinement.** By default each orbital is confined by an energy shift of
0.1 eV. This contracts the diffuse free-atom orbitals toward their size in a
molecule or solid. It is a basis option:
`basis={"name": "PAW-LCAO", "energy_shift": None}` switches it off.

**Grid spacing.** Most files are converged at `grid_spacing` 0.15 Å. A few
files with a hard semicore d (Cu, for example) need 0.12 Å. Check G7 records
each file's working spacing.

**Which set.** The calculator's `paw_dataset` option names a folder of this
repository. It defaults to the set that matches the functional (`"lda-sr"`
for LDA). A Kohn–Sham calculation whose functional differs from its
datasets' warns and names the folder of the matching set. On the command line
the same option is `--paw-dataset pbe-sr`.

**Spin-orbit sets.** Spin-orbit coupling breaks the conservation of the spin
projection S_z.
- **`method="dft"`** on `lda-fr/` or `pbe-fr/` runs spinor (noncollinear)
  Kohn–Sham with the spin-orbit term. It works for `lda`, `pbe`,
  `pbe+rvv10`, `r2scan` and `hse06`, with analytic forces, and also gives
  the stress in crystals.
- **Other methods** on a molecule with a spin-orbit set raise
  `NotImplementedError`, and DFT+U refuses a spin-orbit set.
- **Core states:** those stored in a `-fr` file are j-averaged, so the
  hybrid core-valence exchange does not include the core spin-orbit
  splitting.

**Mixed bases.** A per-element basis may give each atom its own size, as long
as all atoms use PAW-LCAO:

```python
basis={"O": {"name": "PAW-LCAO", "size": "DZP"}, "H": {"name": "PAW-LCAO"}}
```

The full list of options is in Mandacaru's
[pseudopotential guide](https://mandacaru.readthedocs.io/en/latest/guide/pseudopotentials.html).

## How a dataset is built

The method is Blöchl's projector augmented-wave method, in its frozen-core
form. For each element:

1. **Reference atom.** The all-electron atom is solved self-consistently in
   the set's functional. The `-sr` sets use the scalar-relativistic
   (Koelling–Harmon) equation. The `-fr` sets use the radial Dirac equation:
   - each (n, l, κ) level is solved separately;
   - an open subshell is shared between its two j in proportion to 2j+1.
2. **Partial waves.** For each angular momentum l of the valence, the
   all-electron partial waves are taken at two energies: the bound level
   and one scattering energy above it.
3. **Smooth partial waves.** Inside a cutoff radius r_c, each partial wave
   is replaced by a smooth one, expanded in spherical Bessel functions.
4. **Local potential.** A smooth local potential replaces the singular
   nuclear attraction inside the local radius. Beyond it, the stored ionic
   potential is exactly −Z_v/r, joined to the inside by a smooth cubic step.
5. **Projectors.** Projectors are made dual to the smooth partial waves.
6. **One-center terms.** These are the overlap correction q, the kinetic
   and potential differences, the coupling matrix D, the compensation
   charge, and a partial core density for the nonlinear core correction.
7. **Core states.** The orbitals and occupations of every frozen core
   subshell are stored with the core density they reproduce. In the `-fr`
   sets they are the j-averaged large components.
8. **Spin-orbit term (`-fr` only).** The L·S coupling of each channel is
   built from the two j branches, then calibrated against the Dirac atom
   (below).
9. **Release.** The file is written to a temporary path and read back. It
   is moved into place only if every check below passes on the file as
   read.

In a calculation, the smooth partial waves (and their multiple-zeta and
polarization companions) are the basis. The projectors and D add a nonlocal
term to the Hamiltonian, and q makes the basis overlap `S + C q Cᵀ`.

**Two simplifications relative to the full PAW method:**

- **Linearized one-center terms.** The one-center energies are linearized
  around the reference atom, so D is a fixed matrix per element. The error is
  second order in how far the atom's density differs from the reference
  atom's.
- **Frozen core.** The core electrons stay frozen.

**Semicore valence.** These elements carry a filled semicore s and p shell
in the valence (sodium: 9 electrons), so their p channel has projectors. The
list is the same in all four sets:

| Shell | Elements |
|---|---|
| 2s 2p | Na, Mg |
| 3s 3p | K, Ca, Sc, Ti, V, Mn, Fe, Co, Ni, Cu, Zn |
| 4s 4p | Rb, Sr, Y, Zr, Nb, Mo, Tc, Ag, Cd |
| 5s 5p | Cs, Ba, La–Lu, Hf, Ta, W, Re, Ir, Pt, Au, Hg, Pb, Po |
| 6s 6p | Fr, Ra, Ac, Th, Pa, U |

**Frozen shells.** These filled shells are kept in the frozen core instead of
the valence:

| Shell | Elements |
|---|---|
| 3d | Br, Kr |
| 4d | Te, I, Xe |
| 4f | Lu, Hf, Ta, W, Re, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn |

Se carries its 3d in the valence, with a d sphere of 2.3 Bohr.

### Spin-orbit calibration

**Why it is needed.** The basis functions are scalar. Their inner shape is
fixed by the scalar-relativistic partial wave, so they cannot follow the
contraction of the j = l − 1/2 branch. Mandacaru's isolated atom therefore
underestimates the spin-orbit splitting, by more for heavier shells. These
were measured on uncalibrated files during the regeneration, against the
Dirac atom of that stage:

| Shell | Uncalibrated engine | Dirac atom |
|---|---|---|
| Bi 6p | 1.607 eV | 1.836 eV |
| Pb 6p | 1.168 eV | 1.358 eV |
| Tl 6p | 0.867 eV | 0.980 eV |
| Au 5d | 1.567 eV | 1.611 eV |

A model of the same operator in a complete radial basis gives Bi 6p 1.794 eV.
So this is a limit of the basis shape, and the generator cannot remove it
through how it builds the dataset.

**How it is calibrated.**
1. **Measure in the engine's own atom.** The generator solves the atom with
   Mandacaru's own spinor SCF: one atom in a 10 Å box, Γ only, DZP, grid
   0.2 Å, Fermi–Dirac smearing of 0.02 to 0.4 eV. It reads the j-multiplet
   pairs of each bound occupied shell.
2. **Scale the coupling.** It scales the stored coupling of each shell until
   the engine's splitting equals the Dirac splitting to 0.3 % or 1 meV. The
   steps are Newton steps on a measured slope, and the smearing width is held
   fixed after the first round.
3. **Typical factors** are 1.04–1.11 for the lanthanides. The heavy 6p
   shells need the largest ones.

**Shells that cannot be measured** keep the coupling as built and are
marked `reported` by check G13; they are not counted as passed. This covers:
- an open d or f shell, which is fractionally filled, so its splitting moves
  with the SCF occupations;
- an atom whose spinor SCF does not converge at any width.

These are the generator's release rows for the installed files:

| Set | Shells reported (kept as built) | Atoms whose spinor SCF does not converge (every shell kept as built) |
|---|---|---|
| `lda-fr/` | Co d, Er d, Gd d and f, Ho d, La d, Pa d and f, Pt d, Ru d, Sm d, Tb d, Th f, Tm d | none |
| `pbe-fr/` | Ce d and f, Er p and d, La d, Pa d, Pr f, Ru d, Th f, Tm d, Yb d | Cr, Dy, Eu, Gd, Nd, Pm, Sm, Tb |

- **Largest error on a reported shell:** the `pbe-fr/` Er 5p, which is
  −4.6 % from the Dirac atom at every width that converges. The open 4f moves
  it, so it could not be calibrated.
- **What to do:** for spin-orbit properties of these elements, compare the
  result with the scalar set plus a known atomic splitting.

## Generating the datasets

Generating needs a Mandacaru development install whose generator carries the
version-3 changes, with `MANDACARU_PAW_PATH` set. The generator's command is
`mandacaru-build`, run as a module:

```bash
python -m mandacaru.pseudopotentials.build_cli --pp PAW --relativistic --xc LDA --all --workers 2 --ghosts repair --check --install
python -m mandacaru.pseudopotentials.build_cli --pp PAW --dirac        --xc LDA --all --workers 2 --ghosts repair --check --install
python -m mandacaru.pseudopotentials.build_cli --pp PAW --relativistic --xc PBE --all --workers 2 --ghosts repair --check --install
python -m mandacaru.pseudopotentials.build_cli --pp PAW --dirac        --xc PBE --all --workers 2 --ghosts repair --check --install
```

- **Where files go:** `--install` writes into `lda-sr/`, into `lda-fr/` with
  `--dirac`, into `pbe-sr/` with `--xc PBE`, or into `pbe-fr/` with both.
  `--output DIR` writes to a staging folder instead.
- **Single elements:** `--element Fe Cu` rebuilds only those elements.
- **`--ghosts repair`** is the default. When the default construction fails a
  check, the generator searches for a repair: zero norm deficit, a shorter
  local radius, a raised local potential, rebalanced cutoffs or more Bessel
  functions.
- **`--check`** adds the scattering-phase comparison (G3) on the file read
  back.
- **What the command prints:** the status of every check for each element,
  and every `reported` or `waived` row with its reason.

**What a refusal means.** A file that fails a refusing check is not written.
Nothing is left beside the target path, and the previous file stays in
place. The refusal says which check failed and by how much. Fix it through
the recipe (radii, references, channels), never by loosening a threshold.

**Cost.** One element runs per process, and each takes under 2 GB of memory.
Measured per-file times on a loaded 16 GB machine:

| Elements | Time per file |
|---|---|
| O (`pbe-sr/`) | 5 s |
| 5d metals, fr | 5–12 min |
| 6p elements and actinides, fr | 10–30 min |
| Lanthanides, fr | 15–35 min |
| `pbe-fr/` Pa | 58 min |

Most of the time goes to the all-electron atom and the spin-orbit
calibration. Build on an otherwise idle machine: a worker killed under
memory pressure leaves the process pool waiting forever.

**Repairing an existing file** without regenerating it uses two native
entries. Both release the file through every check and leave it untouched if
one fails:
- `_native.paw_reionize(context, path)` enforces the −Z_v/r tail;
- `_native.paw_recalibrate_spin_orbit(context, path)` redoes the spin-orbit
  calibration.

`paw_reionize` is not idempotent near the local radius: apply it once, to
the original file.

## Checks on every file

The generator runs these checks on every file as it is read back, through
`paw::validate`. They are written in Rust and use fixed thresholds. Each row
ends in one of these statuses:

| Status | Meaning |
|---|---|
| `pass` | within the threshold |
| `FAIL` | refuses the file |
| `reported` | measured, no threshold |
| `waived` | recorded with its reason |
| `n/a` | the check does not apply to this file |
| `UNEVALUATED` | the check was skipped; never counted as a pass |

| Check | What it measures | Threshold |
|---|---|---|
| G1 reload | the file read back reproduces every channel's bound levels | 1e-5 Ha |
| G2 ghosts | no spurious bound state in any channel, j branch or angular momentum without projectors | none deeper than 1e-4 Ha below the atom |
| G3 scattering | phase of the logarithmic derivative against the all-electron atom | 0.05 rad within ε ± 0.5 Ha, 0.3 rad within ε ± 1 Ha; channels without projectors 0.2 rad |
| G4 states | missing bound states; an intruding hydrogen 1s | missing levels < 1e-2 Ha (5e-2 without projectors); 1s miss < 1 |
| G5 core | format 3; every core subshell stored; Σ occupations = Z − Z_v; the stored orbitals reproduce the stored core density | 1e-8 e; 1e-8 of the peak density |
| G6 conditioning | projector duality, symmetry of D, smallest eigenvalue of the overlap operator | 1e-8; 1e-4 Ha; ≥ 0.1 |
| G7 hardness and grid | kinetic energy above the cutoff q_cut; the working grid spacing | ≤ 1 Ha; 0.15 Å, a finer one recorded |
| G8 sphere and partition | widest augmentation sphere against the shortest bond; the core/valence split equals the twin set's | below the bond; equal |
| G9 smooth core | integrated partial core charge, against the stored value | reported; 1e-4 e |
| G10 transferability | double ionization energy against an all-electron X2C atom | measured on a pilot set, not a per-file check |
| G11 grid | radial points at r_i = (i+1) h | 1e-4 of h |
| G13 spin-orbit (`-fr`) | the engine's isolated-atom j-splitting of every bound valence shell against the Dirac atom | 2 % or 10 meV; reported for open d/f shells |
| G14 ionic tail | max \|r v_local + Z_v\| from 1.6 r_c to the end of the table; ∫(v_local + Z_v/r) d³r over the same range | 2e-3 e; 0.1 Ha·bohr³ |

G13 and G14 are new in version 3. There is no G12 row in a file's report:
G12 is the library check on Bi2Se3 described below.

**Re-checking a file.** Pass `twin=` to also compare the core/valence split
with the twin file:

```python
from mandacaru import _native
from mandacaru.pseudopotentials.native_context import context

report = _native.paw_validate_file(context(), "pbe-fr/Bi.parquet",
                                   twin="pbe-sr/Bi.parquet")
print(report["passed"])
for row in report["rows"]:
    print(row["gate"], row["check"], row["value"], row["threshold"], row["status"])
```

A file read alone gives `n/a` for G1, because it has no generated record to
compare with.

`paw_validate_file`, `paw_reionize` and `paw_recalibrate_spin_orbit` come
with the version-3 generator. A Mandacaru build without it does not have
them.

## Checks that a single file cannot pass alone

An isolated atom in a box cannot see some defects. The library is therefore
also checked in crystals, against the twin set.

### The fr tail defect

**Symptom.** With the first fr files of the regeneration, `pbe-fr/` Bi2Se3
gave a stress of −90 GPa and an energy 105 eV above `pbe-sr/`. Spin-orbit
coupling can change the stress of Bi2Se3 by about 1 GPa at most.

**Cause.**
- The Dirac reference atom filled an open subshell one j level at a time,
  while the partial waves are the spherical average. Version 3 shares an
  open subshell between its two j levels in proportion to 2j+1, and builds
  the density from the averaged amplitudes.
- The stored local potential therefore carried the Hartree potential of the
  density difference out to 20 Bohr, instead of −Z_v/r.
- That adds a constant α = ∫(v_local + Z_v/r) d³r to every atom, and the
  crystal energy shifts by N_val·Δα/V. Against their sr twins, Δα was
  30–110 Ha·bohr³ for Bi, Pb and Sn.
- An isolated atom cannot see it: the fr − sr energy of Bi in a box of 8 Å
  and of 12 Å agrees to 0.6 meV, with the defect and without it.

**Fix.** The 2j+1 sharing of open subshells, plus the exact −Z_v/r tail
beyond the local radius. Check G14 now refuses any file whose tail or
constant is off.

**Check on Bi2Se3.** At the same cell (DZP, k 3×3×3, grid 0.18 Å), the
stress difference against `pbe-sr/`:

| `pbe-fr/` files | σ_xx = σ_yy | σ_zz | E − E(`pbe-sr/`) |
|---|---|---|---|
| first regeneration | −89 GPa | −103 GPa | +105 eV |
| 2j+1 sharing; the tail checked only to 15 Bohr | −8.65 GPa | −8.67 GPa | +13.5 eV |
| version 3: exact tail over the whole table | −0.93 GPa | −0.81 GPa | +9.2 eV |

- **The remaining stress:** −0.9 GPa is within the spin-orbit effect
  expected.
- **The remaining energy:** +9.2 eV is a constant per atom between the two
  sets. Compare energies only within one set.
- **The next term to check:** the ionic constant inside the local radius
  (+0.475 Ha·bohr³ for Bi), if a tighter twin agreement is needed.

### Spin-orbit topology (G12)

Bulk Bi2Se3 is a strong topological insulator only with spin-orbit
coupling. The Z2 index comes from the parities at the eight
time-reversal-invariant momenta, and from the flow of the Wannier centers on
the six time-reversal-invariant planes.

**Settings:** PBE+rVV10, DZP, k 3×3×3, grid 0.15 Å.
**Cells:** relaxed with PBE+rVV10, TZP, k 3×3×3 and grid 0.151 Å (64 points
along each axis, held fixed during the relaxation), to a maximum force below
0.01 eV/Å and a maximum stress below 0.02 GPa.

| Cell | a (Å) | c (Å) | V (Å³) | u(Bi) | u(Se2) |
|---|---|---|---|---|---|
| `pbe-sr/`, scalar | 4.190 | 28.147 | 142.67 | 0.4016 | 0.2081 |
| `pbe-fr/`, spin-orbit | 4.218 | 28.145 | 144.52 | 0.4015 | 0.2080 |
| Experiment (Nakajima) | 4.143 | 28.636 | 141.89 | 0.4008 | 0.2117 |

These cells belong to the library of 2026-11-14/15, not to this release.
The scalar cell was relaxed with `pbe-sr/` Bi and Se of 2026-11-13. The
spin-orbit cell was relaxed with `pbe-fr/Bi` of 2026-11-14 and the first
`pbe-fr/Se` (content digest `ecc33f30`). That Se file missed the ionic tail
check and was replaced on 2026-11-15.

A single point with this release, at the same cells and settings (TZP, 64
grid points per axis), gives:
- **`pbe-sr/`, scalar cell:** σxx = −0.0034 GPa and σzz = +0.0027 GPa, so the
  cell is still relaxed to below 0.02 GPa.
- **`pbe-fr/`, spin-orbit cell:** σxx = +0.113 GPa and σzz = +0.118 GPa,
  above the 0.02 GPa it was relaxed to. The files installed just before this
  release, which already had the replaced Se, give +0.114 and +0.119 GPa on
  the same build, so this regeneration moves it by 0.001 GPa only. The shift
  matches the Se replacement, which moved this stress by about 0.13 GPa on a
  DZP check.

The spin-orbit cell would need a new relaxation to stand for this release;
the table keeps the cell as relaxed.

The spin-orbit cell is 1.8 % longer in a and 1.7 % shorter in c than
experiment. On a 0.12 Å grid the stress along c is about 1 GPa larger and
points to c ≈ 27.4 Å, so c carries an uncertainty of about 0.7 Å at this grid.

| Set | Cell | Z2 indices | Band gap at Γ | Method |
|---|---|---|---|---|
| `pbe-sr/` | scalar | (0;000), trivial | +0.13 eV, normal order | parity and flow |
| `pbe-fr/` | scalar | (1;000), strong topological insulator | −0.47 eV, inverted | parity and flow |
| `pbe-fr/` | spin-orbit | (1;000), strong topological insulator | −0.46 eV, inverted | parity and flow |

With spin-orbit coupling the flow gives Z2 = 1 on the three planes k = 0 and
0 on the three planes k = π.

Rerun with this release (parity method, same settings and cells): `pbe-fr/`
at the spin-orbit cell gives (1;000) with a band gap at Γ of −0.466 eV,
inverted, and `pbe-sr/` at the scalar cell gives (0;000) with +0.131 eV,
normal order.

**Relaxing a crystal.** With a fixed grid spacing, the number of grid
points jumps as the cell changes. In Bi2Se3 each jump shifts the energy by
about 0.2 eV. The analytic stress is the derivative at a fixed number of
points, so it does not see these jumps. Across a jump, the optimizer's
curvature estimates become inconsistent and the relaxation stalls. Hold the
number of points fixed while the cell moves, then re-mesh at the final cell.

## Validation of this release

All 368 files are format 3, carry their core states, and record no defect.
- **Regenerated (2026-11-16):** 295 files, all 184 scalar files and the 111
  fr files that an earlier repair step (`paw_reionize`) had released. They
  were built by `mandacaru-build`'s route from one build of the generator.
  Each file passed the Rust release checks G1–G14 before it was written.
- **Kept:** the other 73 fr files were already made by this generator. Five
  of them were regenerated as a check, and all five came out byte-identical:
  `lda-fr/Au`, `pbe-fr/Fr` and `pbe-fr/Gd` with the build that made the 295, and
  `lda-fr/Pt` and `pbe-fr/Se` on the day the generator was ported. The 73
  installed files are therefore kept.
- **Native validation:** every one of the 368 installed files is then
  validated natively, G1–G14 against its twin. The result is listed below.

The digest of every file before and after the regeneration is kept with the
generation logs.

<!-- generated: validation -->
**368 of 368 files pass** the native validation, G1–G14 against their twins.

| Set | Pass | G13 `reported` (spin-orbit shell not measurable) |
|---|---|---|
| `lda-sr/` | 92 of 92 | — (scalar) |
| `pbe-sr/` | 92 of 92 | — (scalar) |
| `lda-fr/` | 92 of 92 | Co, Er, Gd, Ho, La, Pa, Pt, Ru, Sm, Tb, Th, Tm |
| `pbe-fr/` | 92 of 92 | Ce, Cr, Dy, Er, Eu, Gd, La, Nd, Pa, Pm, Pr, Ru, Sm, Tb, Th, Ti, Tm, Yb |

**Ionic tail (G14).** In every file the local potential beyond 1.6 r_c is
−Z_v/r to within 1.7e-8 e. The scalar files of the previous install reached
1.06e-3 e (Lu). Their tail constant, up to 6.9e-2 Ha bohr³ (Pt), is now
at most 6.7e-7 Ha bohr³ in any file.

**Crystal check of the fr ionic potential (G14).** In Bi2Se3 (PBE+rVV10,
DZP, k 3×3×3, grid 0.18 Å), switching both elements from `pbe-sr/` to
`pbe-fr/` changes the stress by −0.775 GPa (xx) and −0.665 GPa (zz). Both are
inside the 1 GPa target. The previous install gave −0.777 and −0.668 GPa.

**Single points against the previous install.** Same program build, same
settings, energy and stress changes from the old files to the new ones:

| Crystal | Set | ΔE (eV) | Δσ (GPa) |
|---|---|---|---|
| Si, diamond, a = 5.43 Å (PBE, DZP, k 4×4×4, grid 0.15 Å) | `pbe-sr/` | −0.007 | +0.019 |
| Cu, fcc, a = 3.61 Å (PBE, DZP, k 8×8×8, grid 0.15 Å) | `pbe-sr/` | +0.074 | −0.847 |
| Bi2Se3, scalar cell (PBE+rVV10, TZP) | `pbe-sr/` | +0.017 | −0.004 (xx), −0.004 (zz) |
| Bi2Se3, spin-orbit cell (PBE+rVV10, TZP) | `pbe-fr/` | −0.003 | −0.001 (xx), −0.001 (zz) |

The scalar changes come from the exact ionic tail. The old Cu file carried a
tail constant of −4.8e-3 Ha bohr³, which a crystal adds as a constant energy
per volume. Its change has the same sign as the measured one, and the same
order of magnitude.
<!-- end generated: validation -->

A record with no defects means that the checks above passed when the file
was generated. It is not a certificate for every chemical environment. For
molecules containing elements with a `reported` row, check bond lengths
against an energy scan rather than trusting a relaxed geometry alone.

## File format

Each `<Symbol>.parquet` is a self-describing Mandacaru record:

- format `mandacaru-pseudopotential`, version 3, family `paw-lcao`;
- a uniform radial grid starting one step from the origin and ending at
  30 Bohr, with every quantity in atomic units. The grid sizes are:
  - smallest: 1500 points 0.020 Bohr apart (H, Be);
  - largest in the scalar sets: 34500 points 0.00087 Bohr apart (U);
  - finer in a few fr files: `lda-fr/Te` 39000 points, `pbe-fr/Ge` 48000 and
    `pbe-fr/Se` 51000.

A record holds:

- the all-electron and smooth partial waves;
- the projectors and the one-center matrices;
- the local potential;
- the core and partial core densities, **and the core states**;
- the compensation charge and the frozen one-center energy;
- the generation record: functional, relativity, relativistic exchange,
  core correction, frozen and semicore subshells;
- any recorded defects;
- in `lda-fr/` and `pbe-fr/`, the calibrated spin-orbit term.

**Digests.** A file has two fingerprints:
- **File SHA-256:** the SHA-256 of the file itself.
- **Content digest:** `_native.paw_content_sha256(context, path)` hashes
  the physics content: the payload and tables, without the format version
  and the core states.

Compare a file with a recorded digest using the same kind of digest. Use a
content digest from the same Mandacaru build: its definition is versioned
(`mandacaru-paw-content/2`).

## References

- P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994) — the projector
  augmented-wave method.
- D. D. Koelling and B. N. Harmon, *J. Phys. C* **10**, 3107 (1977) — the
  scalar-relativistic equation.
- J. P. Perdew and A. Zunger, *Phys. Rev. B* **23**, 5048 (1981) — LDA
  correlation.
- J. P. Perdew, K. Burke and M. Ernzerhof, *Phys. Rev. Lett.* **77**, 3865
  (1996) — the PBE functional.
- A. H. MacDonald and S. H. Vosko, *J. Phys. C* **12**, 2977 (1979) — the
  relativistic correction to exchange.
- S. G. Louie, S. Froyen and M. L. Cohen, *Phys. Rev. B* **26**, 1738 (1982)
  — the nonlinear core correction.
- J. Heyd, G. E. Scuseria and M. Ernzerhof, *J. Chem. Phys.* **118**, 8207
  (2003); erratum **124**, 219906 (2006) — the HSE screened hybrid.
- L. Fu and C. L. Kane, *Phys. Rev. B* **76**, 045302 (2007) — the Z2
  index from parities.
- A. A. Soluyanov and D. Vanderbilt, *Phys. Rev. B* **83**, 235401 (2011) —
  the Z2 index from the flow of Wannier centers.
- H. Zhang *et al.*, *Nat. Phys.* **5**, 438 (2009) — Bi2Se3 as a
  topological insulator.
- S. Nakajima, *J. Phys. Chem. Solids* **24**, 479 (1963) — the Bi2Se3
  crystal structure.

## License

MIT, see [LICENSE](LICENSE).
