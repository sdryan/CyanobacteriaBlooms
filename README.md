# CyanobacteriaBlooms
Files and Data associated with Brown, Munther, Ryan Article

# Modeling Multi-Phase Growth in Cyanobacterial Harmful Algal Blooms

Simulation and analysis code for the manuscript by A. M. Brown, D. S. Munther, and S. D. Ryan.

The model is a one-dimensional, four-field reaction–diffusion system on a no-flux domain representing a surface transect of the western basin of Lake Erie:

- `u` — *Microcystis aeruginosa* ("Mike"), non–nitrogen-fixing, nitrogen-limited
- `w` — *Anabaena circinalis* / *Dolichospermum circinale* ("Ana"), nitrogen-fixing, phosphorus-limited
- `N`, `P` — dissolved nitrogen and phosphorus

Growth follows Michaelis–Menten kinetics; nutrients are supplied by a localized source term `f` that switches off at `T_c` (springtime agricultural runoff). The system is solved in non-dimensional form by explicit finite differences.

## Requirements

- **MATLAB** R2025a — base simulations and all manuscript figures.
- **Python** 3.9+ with `numpy` and `matplotlib` — analytical-bound verification, sensitivity analysis, and the Lake Erie phenology comparison.
  ```bash
  pip install numpy matplotlib
  ```

## Repository layout

```
current-code/                MATLAB source (model + figure generation)
Decay_Spreading_Analysis/    Python source (verification, sensitivity, data comparison) + analysis notes
JMBSubmission/               Manuscript (LaTeX, svjour3) + figures
PNAS_Submission_2025/        Manuscript (LaTeX, PNAS) + figures
```

## MATLAB code (`current-code/`)

| File | Purpose | Manuscript figure |
|------|---------|-------------------|
| `algae_1d_ext_nondim_paper_2025.m` | Base non-dimensional simulation; bloom onset and multi-phase (Mike→Ana) dynamics | Figs 1–4 |
| `new3dmovies-Nov2025/algae_1d_code_2025.m` | 3-D space–time evolution of each field; frame export for movies | Figs 1, 3 |
| `new3dmovies-Nov2025/threedplotgen_algae.m` | Renders the 3-D surface plots from saved data | Figs 1, 3 |
| `competition/algae_1d_comp_2025.m` | Sweeps competition strength `beta` (31 values); critical times vs `beta` | Fig 5A |
| `competition/algae_1d_qmax_2025.m` | Sweeps efficiency ratio `qmax/kmax` (61 values) | Fig 5B |
| `nutrientmanagement/algae_1d_f0_source_2025.m` | Sweeps nutrient load `f0` (21 values); onset/decay/overtake times and the bifurcation threshold | Fig 6 |
| `Phase Diagrams/algae_1d_phasediagrams_2025.m` | Generates phase-diagram data (`phase_time_form.txt`, `phase_time_decay.txt`, `phase_time_overtake.txt`) over `(f0, beta, qmax/kmax)` | Fig 7 (data) |
| `Phase Diagrams/algae_1d_phasediagrams_largef0_2025.m` | Phase-diagram data for the high-`f0` regime | Fig 7 (data) |
| `Phase Diagrams/bloom_phasediagrams.m` | Reads the `phase_time_*.txt` files and renders the phase diagrams | Fig 7 |
| `Bounds On Depletion Time/algae_1d_ext_nondim_paper_2025.m` | Decay / depletion-time runs after source shut-off | context for Fig (decay), Theorem/Corollary |


**Running:** open a script in MATLAB and press *Run*. Simulation scripts write image frames (`Algae_mov_###.jpg`) and, for the sweep scripts, the timing data used by the figure scripts. For the phase diagrams, run the two `algae_1d_phasediagrams_*` generators first (they write the `phase_time_*.txt` files), then `bloom_phasediagrams.m`.

**Key parameters** are set at the top of each script: maximal uptake `kmax`, `qmax`; half-saturation `kn`, `qn`; mortality `delta`; nutrient degradation `gamma`; uptake rates `eta1`, `zeta1`, `zeta2` (with `eta2 = 0`, Mike does not consume phosphorus); competition `beta1`, `beta2`; source strength `f0`; domain length `L`, grid points `M`, step `dt`, total time `Time`, and shut-off time `T_c`. The non-dimensional groups are computed in-script from these.

## Python code (`Decay_Spreading_Analysis/`)

| File | Purpose | Output |
|------|---------|--------|
| `algae_sim.py` | Vectorized port of the non-dimensional model; verifies the bloom-extinction theorem (post-shut-off decay rate → mortality `delta`) and the decay-time bound, and evaluates the Fisher–KPP spreading speed | `decay_verification.png` + printed summary |
| `sensitivity.py` | Local normalized sensitivities (elasticities, central differences at ±5%) of peak bloom, decay time, and overtake time to the model parameters | `sensitivity_analysis.png` + printed table |
| `model_vs_lakeerie.py` | Aligns the base run to a calendar (runoff onset = 1 July) and compares bloom onset/peak/decay and the Mike→Ana succession against documented western-basin phenology | `model_vs_lakeerie.png` + printed comparison |

**Running:**
```bash
cd Decay_Spreading_Analysis
python3 algae_sim.py
python3 sensitivity.py
python3 model_vs_lakeerie.py
```
These scripts are self-contained (parameters mirror the MATLAB base case) and use a coarser grid than the MATLAB runs, which is justified because the discrete diffusion number is `O(10^-9)` on the 400 km domain (diffusion is dynamically negligible relative to reaction at lake scale).

## Data

The phenology comparison is anchored to the NOAA GLERL / CIGLR western Lake Erie HAB monitoring program (weekly, May–October, since 2012): NCEI accession doi:10.25921/11da-3x54, described in Boegehold et al. (2023), *Earth Syst. Sci. Data* 15, 3853–3868. The raw data are openly available from NCEI and are not redistributed in this repository.

## Notes

- Variables are non-dimensionalized with length `l0`, time `t0 = 1/kmax`, cyanobacteria scaled by carrying capacity, and each nutrient by its half-saturation constant (see the manuscript appendix for the full derivation).
- When forming the non-dimensional diffusivities, keep the diffusivity and rate units consistent (diffusivities are in µm²·s⁻¹ while rate constants are in hr⁻¹); the simulation output is insensitive to this because diffusion is negligible at lake scale, but the reported `alpha*` values depend on it.

## Citation

A. M. Brown, D. S. Munther, S. D. Ryan, *Modeling Multi-Phase Growth in Cyanobacterial Harmful Algal Blooms*.
