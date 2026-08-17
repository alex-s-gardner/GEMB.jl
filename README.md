<img src="docs/src/assets/logo.png" align="right" width="130" alt="GEMB.jl logo: a firn column of grains coarsening with depth, in the Julia colors">

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alex-s-gardner.github.io/GEMB.jl/dev/)
[![Build Status](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alex-s-gardner/GEMB.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# GEMB.jl: Glacier Energy and Mass Balance Model

GEMB.jl is a Julia implementation of the Glacier Energy and Mass Balance (GEMB, the "B" is
silent) model — a one-dimensional model of the surface energy balance and vertical firn
evolution of glaciers and ice sheets. It couples atmospheric forcing to subsurface
thermodynamics and densification physics to resolve temperature, density, water content,
grain properties, and layer age through the column.

GEMB is a column model of intermediate complexity with no horizontal communication,
prioritizing computational efficiency to accommodate the multi-millennial spin-ups required
to initialize deep firn columns. It is used for interpreting satellite altimetry, firn and
ice-core studies, surface mass balance inversion, and uncertainty quantification. A complete
description of version 1.0 of the model is given in
[*Gardner et al*., 2023](https://doi.org/10.5194/gmd-16-2277-2023).

GEMB.jl is the reference implementation and is developed independently. It began as a
translation of an earlier MATLAB version; since v2.0.0 its physics, defaults, and numerics
are set on their own merits — against the published literature and against the
[Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) and
[IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM). The
[MATLAB version](https://github.com/alex-s-gardner/GEMB) is maintained separately and the two
are no longer expected to agree numerically; departures are recorded under
[Physics notes](https://alex-s-gardner.github.io/GEMB.jl/dev/physics_notes).

## Installation

GEMB.jl is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/alex-s-gardner/GEMB.jl")
```

Climate forcing lives in the companion
[GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl), also
unregistered. It is optional — GEMB consumes any conforming `DimStack` — but it is what the
examples use, and it is required to run the test suite:

```julia
Pkg.add(url="https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl")
```

Julia ≥ 1.11 is required. A Makie backend (CairoMakie, GLMakie, or WGLMakie) enables
`gemb_plot_output`.

## Quick Start

```julia
using GEMB
using GEMB_ClimateForcing

# Climate forcing: synthetic 3-hourly series (a DimStack) → ClimateForcing
cf = GEMB.initialize_forcing(simulate_climate_forcing("test_1", 3))

mp = initialize_parameters(output_frequency=:daily)   # parameters, validated at construction
profile = initialize_profile(mp, cf)                  # the firn/ice column
output = gemb(profile, cf, mp)                        # returns a DimStack

T_surface = surface_timeseries(output[:temperature])
```

For research applications, spin the column up to a quasi-steady state first with
[`forcing_climatology`](https://alex-s-gardner.github.io/GEMB.jl/dev/#Spinup) and `gemb_spinup`.

Runnable scripts are in [`examples/`](examples): `synthetic_example.jl` (spinup and run on
synthetic forcing) and `era5_example.jl` (the same against ERA5 reanalysis).

## Documentation

Full documentation is at **[alex-s-gardner.github.io/GEMB.jl](https://alex-s-gardner.github.io/GEMB.jl/dev/)**:

| Page | Contents |
|---|---|
| [Home](https://alex-s-gardner.github.io/GEMB.jl/dev/) | Installation, quick start, output structure, spinup |
| [Model Architecture](https://alex-s-gardner.github.io/GEMB.jl/dev/architecture) | Data flow, physics modules, the vertical grid, design principles, surface energy balance numerics, meltwater percolation |
| [Model Parameters](https://alex-s-gardner.github.io/GEMB.jl/dev/parameters) | Every option, its default, and the intercomparison or model the default follows |
| [Physics Notes](https://alex-s-gardner.github.io/GEMB.jl/dev/physics_notes) | Where GEMB.jl's physics departs from an earlier form or from a published law, and why |
| [Variable Reference](https://alex-s-gardner.github.io/GEMB.jl/dev/variables) | Every input and output variable with units |
| [API Reference](https://alex-s-gardner.github.io/GEMB.jl/dev/api) | Exported functions and types |

## Testing

```julia
using Pkg
Pkg.test("GEMB")
```

Validation comes from four sources: closed-form checks where a scheme has an analytic
solution, published values from the paper a scheme is cited to, agreement with the
[CFM](https://alex-s-gardner.github.io/GEMB.jl/dev/cfm_comparison) and
[IMAU-FDM](https://alex-s-gardner.github.io/GEMB.jl/dev/imau_fdm_comparison) where the algebra
is the same, and invariants — conservation of mass and energy per timestep and over a whole
run, monotonicity, and bounds.

The test target needs GEMB_ClimateForcing.jl checked out as a sibling directory, resolved via
`[sources]` in `test/Project.toml` (Julia ≥ 1.11).

## Performance

`bench/opt_bench.jl` runs the representative hot path — a 75-year climatological spinup
followed by a 32-year transient run, **~107 model-years at 3-hourly resolution** on a single
column — in **5.7 s** (Apple M2 Max, minimum of 3 runs after warmup). It also emits a
per-field numerical snapshot, so an optimization can be checked for bit-level equivalence
against a baseline rather than only for speed.

## Citation

If you use GEMB.jl in your research, please cite:

> Gardner, A. S., Schlegel, N.-J., and Larour, E.: Glacier Energy and Mass Balance (GEMB): a
> model of firn processes for cryosphere research, *Geosci. Model Dev.*, 16, 2277–2302,
> [https://doi.org/10.5194/gmd-16-2277-2023](https://doi.org/10.5194/gmd-16-2277-2023), 2023.

## Related repositories

- [GEMB_ClimateForcing.jl](https://github.com/alex-s-gardner/GEMB_ClimateForcing.jl) — climate
  forcing for GEMB: real (ERA5-Land, via CDS) and synthetic, produced as a `DimStack`
- [GEMB_GlacierSims.jl](https://github.com/alex-s-gardner/GEMB_GlacierSims.jl) — multi-site and
  regional simulation workflows
- [GEMB (MATLAB)](https://github.com/alex-s-gardner/GEMB) — the separately maintained MATLAB
  implementation
- [Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) and
  [IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM) — other open firn models, used
  as independent cross-checks on GEMB's physics

## Contributing

Contributions are welcome. Please:

1. Justify any change to the physics against the literature or against an independent
   implementation, and record it under
   [Physics notes](https://alex-s-gardner.github.io/GEMB.jl/dev/physics_notes) if it moves
   model output
2. Add tests that pin the new behaviour — against a closed form, a published value, or an
   independent implementation, whichever applies
3. Follow Julia style guidelines
4. Update documentation for new features

## License

MIT — see [LICENSE](LICENSE).

## Authors

GEMB was created by Alex Gardner, with contributions from Nicole-Jeanne Schlegel and Chad
Greene. The Julia implementation was developed by Alex Gardner.

For questions or issues, please open an issue on GitHub.
