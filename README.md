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
| 2 | `window_cost` | the GLS misfit: each depth row whitened with its full covariance, coherence scale solved in closed form |
| 3 | `fit_window` | grid search, then Gauss-Newton with a line search |
| 4 | `window_errors` | formula error bars and the chi-square check |
| 5 | `fit_column` | the whole column, converted to Δλ, abstaining where it cannot know |
| 6 | `bootstrap_column` | error bars by repeating each window's measurement in simulation |

Along the way the guide has the student test their error bars. Step 4
compares the formula with the scatter over 20 synthetic sites; step 6
repeats each window in simulation and compares with the true scatter. The
guide ends by running the estimator on all four sites next to
`ptt.ershadiFabric`, followed by open questions.

**For instructors:**
- `check_step(k, 'ref')` runs a test on the reference solution in
  `project/solutions/+ref/`.
- To grade a student's folder, put it first on the path and run
  `check_step(1)` through `check_step(6)`.
- The reference agrees with `ershadiFabric` on the axis at all four sites.

### What the reference estimator can and can't do

The four radar channels are the measurements; the 18 synthesized azimuths
in each depth row are made from them and carry only about 7 independent
numbers. `apres.coherenceField` therefore gives each row its full covariance
(`apres.coherenceCov`, checked against simulation to within 6% at 10 looks)
and a whitening matrix, and the estimator is a proper GLS fit. Measured on
synthetic sites (`apres.syntheticSite`) with 60 m windows:

- **Accuracy:** θ is unbiased (mean error 0.002°) and Δλ within 0.001-0.002
  of the truth. Modelling each row at its brightness-weighted depth
  (`F.z_eff`) instead of its centre cut the Δλ scatter tenfold.
- **Formula error bars (step 4):** cautious, never optimistic: about 1.7× the
  true scatter for θ and 5× for Δλ. The caution comes from `cond_floor`,
  which raises every variance in a row to at least 5% of the row's largest,
  so no direction of the data is trusted as near-perfect; without it, small
  model errors near phase nodes wrecked fits.
- **Simulated repeats (step 6):** within about 1-2× of the true scatter,
  window by window, for both θ and Δλ. These are the error bars to quote.
- **Phase nodes:** where the birefringent phase passes π the model misfits
  (χ²/dof ≫ 5) and those windows abstain.
- **Rotating axes:** the model assumes one axis for the whole column above
  each window. Where the axis turns with depth most windows abstain and the
  rest get θ wrong; that is the case `fabricGLS` exists for.
- **Weak fabric:** windows abstain when Δλ is below about 0.07 (too little
  phase turn in 60 m at 300 MHz).
- **Real sites:** 9-11 of 16-18 windows reported, median χ²/dof 1.7-3.7: the
  single-column model fits real ice less well than synthetics, and the
  error bars are widened by √(χ²/dof).
- **Calibration:** the phase calibration needs cross-polarized signal between
  40 and 800 m and assumes the top 40 m is isotropic firn. It cannot tell θ
  from −θ (a mirror about north); it picks the branch with smaller antenna
  offsets and records that it did (`apres.calibratePhase`).

## Code

| folder | contents |
|---|---|
| `+apres/` | ApRES reading (`loadBurst`, `rangeProcess`, checked against an independent Python reader to 1e-7), `loadQuadpolSite` (polarization from file names, shared attenuator setting, sub-bin co-registration), `calibratePhase` (antenna phase offsets from reciprocity plus isotropic firn), `observables`, `coherenceField` and `coherenceCov` (the coherence field, its per-row covariance and whitening), and `syntheticSite` (known-answer sites, with adjustable speckle, signal and VV coherence) |
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
