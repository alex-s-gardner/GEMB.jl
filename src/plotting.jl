"""
    gemb_plot_output(output::DimStack; datelims=nothing, depthlims=(-10, 0),
                     variables=nothing, title="")

Create a diagnostic dashboard for a GEMB model run (`output = gemb(profile, cf, mp)`).

The figure has two columns sharing a common (linked) time axis:

- **Left — profile fields.** Each 2-D variable (depth × time, e.g. `temperature`,
  `density`) is interpolated from the evolving Lagrangian grid onto a fixed vertical
  grid with [`gemb_interp`](@ref) and drawn as a heatmap. Colormaps are chosen per
  field and color limits are clipped to robust (2–98%) percentiles for maximum
  contrast. Long runs are strided in time so the raster stays a sensible width.
- **Right — scalar fields.** 1-D variables are grouped by physical theme and shared
  unit (energy fluxes, mass fluxes, densification, …) and overplotted as time series
  on a common axis, with a legend when a panel holds more than one series.

Redundant decorations are removed: x ticks and the "year" label appear only on the
bottom panel of each column, and a header banner records provenance (GEMB version,
time span, sampling cadence, max depth, striding, generation time).

# Arguments
- `output::DimStack`: The `DimStack` returned by [`gemb`](@ref).

# Keywords
- `datelims`: Optional `(start, stop)` limits for the time axis. Accepts a tuple of
  `DateTime`s or of decimal years.
- `depthlims`: `(low, high)` depth limits (metres, negative down). Sets both the
  regridding range and the depth axis of the heatmaps. Defaults to `(-10, 0)` — the
  near-surface zone where most of the dynamics live. Pass a wider range (or derive it
  from the data) to see the full column.
- `variables`: Iterable of variable name `Symbol`s to include. Defaults to every
  variable in `output` except the profile fields `dz`, `albedo`, and `albedo_diffuse`
  (a grid diagnostic and two fields already represented by the surface-albedo series);
  pass any of them explicitly here to include them.
- `title`: Optional figure title. Defaults to `"GEMB diagnostic output"`.

# Returns
A Makie `Figure`.

!!! note "Requires Makie"
    This function lives in a package extension. Load a Makie backend
    (e.g. `using CairoMakie` or `using GLMakie`) before calling it.

# Example
```julia
using GEMB, CairoMakie
output = gemb(profile_spunup, cf, mp)
fig = gemb_plot_output(output; depthlims=(-5, 0))
save("gemb_diagnostics.png", fig)
```
"""
function gemb_plot_output end
