# Audit scripts

Each script reproduces one finding of the 2026-10-13 review of GEMB's physics. The numbering
follows that review's findings table (kept with the GlobalGlacierAnalysis.jl project as
`notes/gemb_pipeline_review_2026-10-13.md`); the script header states the finding in full, so
the scripts stand on their own. Every script runs unchanged against `main` (0cda823) and
against `fix/physics-audit`, prints which one it found, and states in its header what each
should print.

## Running

The scripts need only GEMB and the standard library. From a checkout of this branch:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/audit/01_longwave_absorption.jl
```

Against `main`, point `--project` at a worktree of 0cda823 and keep running the scripts from
this branch, since `main` does not have them:

```sh
git worktree add /tmp/gemb-main 0cda823
julia --project=/tmp/gemb-main -e 'using Pkg; Pkg.instantiate()'
julia --project=/tmp/gemb-main scripts/audit/01_longwave_absorption.jl
```

Pass `--startup-file=no` if your `startup.jl` activates the project in the working directory.
`common.jl` holds the version detection and the wrappers for the few internal functions whose
signatures differ between the two versions.

## Index

| Script | Finding | `main` | `fix/physics-audit` |
|---|---|---|---|
| `01_longwave_absorption.jl` | Absorbed LW↓ instead of ε·LW↓ | +0.45 K/day at radiative equilibrium | 0 K |
| `02_bare_ice_albedo.jl` | 1: bare-ice albedo held near 0.58 | albedo 0.58 for both `albedo_ice`, SW ratio 1.00 | 0.30 / 0.48, ratio 1.346 |
| `03_stability.jl` | 4, 16: stable-branch jump, effective height, Exner factor | SHF/ΔT jumps ×85 at ΔT ≈ 1.6 K; pressure ratio 0.94 | < 1% per 0.01 K; ratio 1.0000 |
| `04_rain.jl` | 3 and rain accounting: rain compacts snow, refreeze uncounted | cell 1 at 450 kg m⁻³; cold-cell refreeze 0 | 350 kg m⁻³; refreeze 1.0745 kg m⁻² |
| `05_vapor_phase.jl` | 9: vapor at a melting surface booked to ice | energy residual ≈ −4.6e4 J m⁻² | 0 |
| `06_explicit_step.jl` | Explicit sub-step vs surface relaxation; low-density divergence | 30% flux change with sub-step; surface 138 K | 0.2%; within 2 K of implicit |
| `07_wet_cell.jl` | 17: wet cell not a latent buffer | LW↑ 302 W m⁻² (cools below 0 °C) | 312.54, the 0 °C value |
| `08_shortwave_budget.jl` | Subsurface SW ignores the final and diffuse albedos | absorbed off by 22–27 W m⁻² | matches for every method |
| `09_radiation_geometry.jl` | 5: forcing SZA, τ and diffuse replaced by `mp` zeros | albedo identical with and without geometry | higher with geometry |
| `10_grain_growth.jl` | 15: `:Arthern` default barely coarsens surface snow | r_e ≤ 0.21 mm after 60 d | 0.30–0.46 mm at 30 K m⁻¹ |
