# GLS tutorial

This is a MATLAB tutorial that walks through the implementation of generalized least squares (GLS) inference. 
It assumes no linear algebra background. 
It walks through the fitting routine and ends with the student building their own least-squares fabric estimator and running it on four real quad-pol ApRES sites.

It is intended to run on your personal computer with MATLAB preinstalled.

## Setup

1. Download `GHOST24_Polarimetric_pRES_OZ.zip` (the GHOST 2023/24 polarimetric
   pRES data, about 130 MB) and unzip it into `data/`, so you have
   `data/GHOST24_Polarimetric_pRES_OZ/PpRES_20240551_001/...`.
2. In MATLAB, `cd` to this folder. Every lesson calls `tutorial_setup`
   itself.

Each lesson is a plain `.m` script split into `%%` sections. Run one section
at a time (Ctrl+Enter, or Cmd+Enter on a Mac) and read what it prints.

## Lessons

Each lesson will generate it's own data, identifies a porblem that uses these data and then solves the problem, and checks the error estimates using Monte Carlo methods.

| # | file | description |
|---|---|---|---|
| 0 | `lesson00_matlab_matrices.m` | vectors, `G*m`, `'`, `.*`, `\`, variance, covariance |
| 1 | `lesson01_ols.m` | least squares by grid search, then normal equations; `C_M`; Monte Carlo |
| 2 | `lesson02_wls.m` | 1/variance weights; the coherence-phase CRB |
| 3 | `lesson03_gls_correlated.m` | correlated noise, `C_d`, whitening with `chol`, variance inflation |
| 4 | `lesson04_chi2.m` | reduced chi-square |
| 5 | `lesson05_prior_traveltime.m` | ill-posed inversion, priors, resolution, abstention | 
| 6 | `lesson06_nonlinear.m` | Jacobian, Gauss-Newton, local minima vs. global minimum |
| 7 | `lesson07_fabricGLS.m` | reading `fabricGLS` block by block, then testing its error bars |
| 8 | `lesson08_apres_radar.m` | how ApRES forms a range profile, the four-shot polarimetric measurement, and how these data differ from the Ridge A survey |
| 9 | `lesson09_apres_data.m` | one GHOST site from raw files to the coherence field |

## The project: build your own estimator

`project/project_guide.m` walks the student through writing six functions in `project/student/`. 
After each one, `check_step(k)` compares it with a hidden reference on data with a known answer, prints PASS or FAIL, and gives some suggestions about what might be wrong.

| step | function | description |
|---|---|---|
| 1 | `coh_model` | the single-column HH-VV coherence model |
| 2 | `window_cost` | the weighted misfit, with the coherence scale solved in closed form |
| 3 | `fit_window` | grid search, then Gauss-Newton with a line search |
| 4 | `window_errors` | linearized error bars and the chi-square check |
| 5 | `window_jackknife` | leave-one-row-out error bars |
| 6 | `fit_column` | the whole column, converted to Δλ, abstaining where it cannot know |

Along the way, the guide makes the student test their own error bars. Their
step-4 σ fails the Monte Carlo test (the synthesized azimuths have correlated
errors): θ's σ is about 4 times too large and the phase gradient's σ about 1.3
times too small, which motivates step 5. The guide ends by running their
estimator on all four sites next to `ptt.ershadiFabric`, followed by open
questions.

**For instructors:**
- `check_step(k, 'ref')` runs a test on the reference solution in
  `project/solutions/+ref/`.
- To grade a student's folder, put it first on the path and run
  `check_step(1)` through `check_step(6)`.
- The reference agrees with `ershadiFabric` on the axis at all four sites.

### What the reference estimator can and can't do

Measured on synthetic sites (`apres.syntheticSite`), with 60 m windows:

- **Δλ:** within about 0.005 of the truth, and its jackknife σ is about 1.4
  times cautious.
- **θ:** its jackknife σ is about 2 times too small. It is unbiased when the
  top 40 m is isotropic. The phase calibration assumes that firn; with fabric
  right up to the surface it biases θ by about −1.1° (Δλ = 0.08, 12 seeds),
  which is where the "speckle bias" of earlier versions came from.
- **Phase nodes:** where the birefringent phase passes π, |C| drops sharply
  with azimuth, the model misfits (χ²/dof ≫ 5) and those windows abstain.
- **Rotating axes:** it assumes one axis for the whole column above each
  window. Where the axis turns with depth most windows abstain, and the few
  that pass get θ wrong; that is the case `fabricGLS` exists for.
- **Weak fabric:** windows abstain when Δλ is below about 0.07 (too little
  phase turn in 60 m at 300 MHz).
- **Choice of branch:** the phase calibration can't tell θ from −θ (a mirror
  about north). It picks the branch with smaller antenna offsets and records
  that it did (`apres.calibratePhase`).

## Code

| folder | contents |
|---|---|
| `+apres/` | ApRES reading (`loadBurst`, `rangeProcess`, checked against an independent Python reader to 1e-7), `loadQuadpolSite` (polarization from file names, shared attenuator setting, sub-bin co-registration), `calibratePhase` (antenna phase offsets from reciprocity plus isotropic firn), `observables` and `coherenceField`, and `syntheticSite` |
| `vendor/+ptt/` | ten unmodified `fabric_anisotropy` files at a pinned commit; see `vendor/README.md` |
| `project/` | the guide, `check_step`, the student skeletons, and the reference solution |

## Running everything

```sh
matlab -batch "run_all"          # lessons 0-9, figures -> figs/  (about 2 minutes)
matlab -batch "run_all([6 8])"   # just some
```

## Data

The GHOST 2023/24 polarimetric pRES measurements on Thwaites Glacier were
made by Ole Zeising. The field log `GHOST24_PpRES_Log_OZ.pdf` ships with the
data and gives the antenna layout used throughout.
