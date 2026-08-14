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
  unit (air temperature, energy fluxes, mass fluxes, …) and overplotted as time series
  on a common axis, with a legend when a panel holds more than one series. Air
  temperature is placed first (top), aligned with the temperature profile heatmap.

Axis, colorbar, and legend labels take their units from each layer's CF metadata (see
[`GEMB_CF_ATTRIBUTES`](@ref)) rather than from a hard-coded table, so a panel cannot
disagree with the data it draws. Fields stored in kelvin are converted to °C for
display; everything else is shown in its stored units.

Redundant decorations are removed: x ticks and the "year" label appear only on the
bottom panel of each column, and a header banner records provenance (GEMB version,
forcing source, time span, sampling cadence, spinup climatology window + convergence
— or "no spinup" — mass budget, generation time).

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
  variable in `output` except the profile fields `dz`, `albedo`, `albedo_diffuse`,
  `grain_dendricity`, and `grain_sphericity`, and the densification scalar group (a
  grid diagnostic, fields already represented by the surface-albedo/grain-radius
  panels, and a lower-priority group); pass any of them explicitly here to include them.
- `title`: Optional figure title. Defaults to `"GEMB.jl glacier firn model output"`.

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
function gemb_plot_output(args...; kwargs...)
    error("""
        `gemb_plot_output` requires a Makie backend. Load one before calling this function:

            using CairoMakie   # for file output (PNG, SVG, PDF)
            using GLMakie      # for interactive windows
            using WGLMakie     # for Jupyter / Pluto notebooks
        """)
end
