<div align="center">

```
        █████ █████ ██████   █████ █████ ██████████
       ░░███ ░░███ ░░██████ ░░███ ░░███ ░░███░░░░░█
        ░███  ░███  ░███░███ ░███  ░███  ░███  █ ░
        ░███  ░███  ░███░░███░███  ░███  ░██████
        ░███  ░███  ░███ ░░██████  ░███  ░███░░█
  ███   ░███  ░███  ░███  ░░█████  ░███  ░███ ░   █
 ░░████████   █████ █████  ░░█████ █████ ██████████
  ░░░░░░░░   ░░░░░ ░░░░░    ░░░░░ ░░░░░ ░░░░░░░░░░
```

**Fortran toolkit for the analysis of Li₃PS₄ molecular dynamics trajectories**

*Program for the Analysis of LPS — Like a Genie !*

<img width="600" height="450" alt="JINIE Screenshot" src="https://github.com/user-attachments/assets/25e1d820-c0cb-4b06-967a-c023d6fdd8ed" />

<br>

![Fortran](https://img.shields.io/badge/Fortran-734F96?style=flat-square&logo=fortran&logoColor=white)
![Intel](https://img.shields.io/badge/Intel%20ifx-0071C5?style=flat-square&logo=intel&logoColor=white)
![License](https://img.shields.io/github/license/wlsdk9803/JINIE?style=flat-square)

</div>

---

## Overview

JINIE is a Fortran-based analysis toolkit designed for Li₃PS₄ (LPS) glass and crystal MD simulations. It covers the full analysis pipeline — from structural properties and ion dynamics to Li-ion hopping kinetics and PS₄ rotational motion — through an interactive terminal UI powered by `fzf`.

```
@..@             ▄▖    ▜     ▘     ▐▘   ▖ ▄▖▄▖
(----)           ▌▌▛▌▀▌▐ ▌▌▛▘▌▛▘ ▛▌▜▘  ▌ ▙▌▚
( >__< )         ▛▌▌▌█▌▐▖▙▌▄▌▌▄▌ ▙▌▐   ▙▖▌ ▄▌
^^ ~~ ^^                 ▄▌
```

---

## Requirements

| Item | Requirement |
|------|-------------|
| Compiler | Intel Fortran Compiler (`ifx`) |
| Interactive UI | [`fzf`](https://github.com/junegunn/fzf) |
| Trajectory format | `.dcd` · `.lammpstraj` |
| System | Li₃PS₄ (3 Li : 1 P : 4 S atom ratio) |

> `READER.f90` hard-codes the LPS atom ratio (3/8 Li, 1/8 P, 1/2 S). Ensure this matches the trajectory being analyzed.

---

## Installation

```bash
git clone https://github.com/wlsdk9803/JINIE.git
cd JINIE
make              # build all → bin/
make bin/MSD      # build single program
make clean
```

Add the wrapper to your PATH:

```bash
export PATH="$PATH:/path/to/JINIE"
```

---

## Usage

Run from the directory containing your trajectory files:

```bash
jinie traj.dcd
jinie traj1.dcd traj2.dcd     # concatenated trajectories
```

The interactive menu launches automatically. Navigate with `← →` (section), `↑ ↓` (property), and `Enter` to run.

```
[type] filter  ·  [← →] section  ·  [↑↓] navigate  ·  [Enter] select  ·  [Esc] cancel
```

Jobs are launched via `nohup` in the background. Logs are written to `./log/<PROP>.log`.

---

## Configuration Files

Place these in `inp/` under your working directory before running.

### `inp/bin.inp` — used by structure & dynamics analyses

```
0.01    # bin size (Å)
1000    # bin number
1000    # time window
0.05    # bin size for q  [R2X / R2N]
500     # bin number for q [R2X / R2N]
2.5     # P–S cutoff (Å)  [CLS / CLF]
4.5     # P–P cutoff (Å)  [CLS / CLF]
3.7     # S–S cutoff (Å)  [CLS / CLF]
500     # time window      [CLS / CLF]
```

### `inp/hop.inp` — used by hopping analyses

```
4       # deltaT  (frames)
1.5     # h_c     (hop function threshold)
20      # nsurr   (neighbor shell count)
```

---

## Program Reference

### Structure

| Program | Description | Output |
|---------|-------------|--------|
| `AGD` | Angle Distribution | `AGD.out` |
| `CLS` | Polyanion Clustering | `CLS.lammpstraj` |
| `CLF` | Polyanion Cluster Fraction | `CLF.out` |
| `RDF` | Radial Distribution Function $g(r)$ | `RDF.out` |
| `R2X` | RDF → X-ray Structure Factor $S(q)$ | `XSF.out` |
| `R2N` | RDF → Neutron Structure Factor $S(q)$ | `NSF.out` |

### Dynamics

| Program | Description | Output |
|---------|-------------|--------|
| `MSD` | Mean Squared Displacement $\langle r^2(t) \rangle$ | `MSD.out` |
| `NGP` | Non-Gaussian Parameter $\alpha_2(t) = \dfrac{3\langle r^4\rangle}{5\langle r^2\rangle^2} - 1$ | `NGP.out` |
| `ICD` | Ionic Conductivity $\sigma$ (from MSD) | `ICD.out` |
| `GSRT` | Self-part Van Hove Correlation $G_s(r,t)$ | `GSRT_Li/P/S.out` |

### Hopping

| Program | Requires | Description | Output |
|---------|----------|-------------|--------|
| `HOP` | — | Hop Function $h_i(t)$ | `HOP.out` |
| `EVT` | `HOP` | Binary event classification (0 = cage, 1 = hop) | `EVT.out` |
| `HSH_G` | `HOP` | Hop Shell Analysis — Glass phase | `HSH_G.out` |
| `HSH_C` | `HOP` | Hop Shell Analysis — Crystal phase | `HSH_C.out` |
| `CDDF_Li` | `HOP` | Conditional Distance Distribution — Li | `CDDF_Li.out` |
| `CDDF_P` | `HOP` | Conditional Distance Distribution — P | `CDDF_P.out` |
| `CDDF_S` | `HOP` | Conditional Distance Distribution — S | `CDDF_S.out` |
| `CRDF` | `HOP` | Conditional Radial Distribution Function | `CRDF_Li/P/S.out` |
| `EFH` | `CLR` | Effective (Forward/Backward) Hopping Ratio | `EFH.out` |
| `CST` | `EVT` | Hop Lifetime & Survival Probability $C_s(t)$ | `CST_Lifetime/Survival/Randomness.out` |

### Rotation (PS₄ Tetrahedra)

| Program | Description | Output |
|---------|-------------|--------|
| `RGRT` | Rotational Van Hove Correlation $G_s(\theta, t)$ | `RGRT.out` |
| `RNGP` | Rotational Non-Gaussian Parameter $\alpha_2^{\rm rot}(t) = \dfrac{\langle\theta^4\rangle}{3\langle\theta^2\rangle^2} - 1$ | `RNGP.out` |

### Visualization & Utilities

| Program | Requires | Description | Output |
|---------|----------|-------------|--------|
| `D2L` | — | DCD → LAMMPS trajectory | `out.lammpstraj` |
| `HVIS` | — | Hop value / event / sieve overlay on trajectory | `hop.lammpstraj` |
| `CLR` | `EVT` | Cage-averaged trajectory (vibrations sieved) | `clear.lammpstraj` |
| `NVIS` | — | Nearest-neighbor assignment at hop events | `neighbor.lammpstraj` |

---

## Hopping Analysis Workflow

```
HOP  ──→  EVT  ──┬──→  CLR  ──→  EFH
                 ├──→  HSH_G / HSH_C
                 ├──→  CDDF_Li / CDDF_P / CDDF_S
                 ├──→  CRDF
                 └──→  CST
                 
HVIS  (independent, uses raw trajectory + hop.inp)
NVIS  (independent, uses raw trajectory + HOP.out)
```

---

## Auxiliary Tools

| Tool | Description |
|------|-------------|
| `jqueue [user]` | Real-time SLURM job monitor with progress bars |
| `wrapping` | Moves all `*.out` files to `OUTPUT/` |

---

## Module Architecture

```
READER.f90            ← single-trajectory module (most programs)
ENSEMBLE_READER.f90   ← multi-ensemble module (CRDF, CDDF_*, CST)
```

All analysis programs `use` one of these two modules. The module is compiled first and linked at build time via the Makefile.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
