"""
    plot_output(output::DimStack; datelims=nothing, depthlims=(-10, 0),
                     variables=nothing, title="", reference_frame=:surface,
                     detrend=:spinup)

Create a diagnostic dashboard for a GEMB model run (`output = gemb(profile, cf, mp)`).

The figure has two columns sharing a common (linked) time axis:

- **Left — profile fields.** Each 2-D variable (depth × time, e.g. `temperature`,
  `density`) is interpolated from the evolving Lagrangian grid onto a fixed vertical
  grid with [`gemb_interp`](@ref) and drawn as a heatmap. Colormaps are chosen per
  field and color limits are clipped to robust (2–98%) percentiles for maximum
  contrast. Long runs are strided in time so the raster stays a sensible width.
  Panels are ordered `temperature`, `age`, `density`, `grain_radius`, `water` —
  physical order, not the order the fields happen to be stored in.
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
  from the data) to see the full column. In the `:datum` frame the upper bound is raised
  to the highest cell centre present, so the surface excursion the frame exists to show is
  never cropped; the lower bound is used as given.
- `variables`: Iterable of variable name `Symbol`s to include. Defaults to every
  variable in `output` except the profile fields `dz`, `grain_dendricity`, and
  `grain_sphericity`, and the densification scalar group (a grid diagnostic, fields
  already represented by the grain-radius panel, and a lower-priority group); pass any
  of them explicitly here to include them.
- `title`: Optional figure title. Defaults to `"GEMB.jl glacier firn model output"`.
- `reference_frame`: Vertical coordinate the profile heatmaps (and the isochrone overlay)
  are drawn in — `:surface` (default) or `:datum`. See
  [Vertical reference frames](#Vertical-reference-frames) below.
- `detrend`: Rate [m per year] subtracted from the `:datum` frame's surface height, making it
  an anomaly about the mean state. Ignored for `:surface`. Needed because GEMB has no ice
  dynamics — see [Vertical reference frames](#Vertical-reference-frames).
  - `:spinup` (default) — the `spinup_smb_rate` provenance the run inherits from
    [`gemb_spinup`](@ref), measured over the climatological cycle the column was spun up on.
    Errors if the run carries none. Preferred because detrending by the *climatological* rate
    leaves the plotted record's own departure from that climatology visible.
  - `:transient` — the same quantity measured over the run being plotted. Agrees closely with
    `:spinup` unless the record's climate departs from the climatology, in which case this
    removes that departure along with the mean.
  - a number — used verbatim. `0.0` gives the raw datum frame, which climbs at the mean
    accumulation rate; worth asking for on a short record.

  Both symbolic options give `(precipitation + evaporation_condensation - runoff) /
  density_ice` over their record divided by its length — SMB by definition. The rate and its
  source are recorded in the header banner, so a figure cannot leave the reader unable to tell
  an anomaly from an absolute height.

# Vertical reference frames

The model's own vertical coordinate is depth below the **instantaneous** surface. Because
the column has a fixed total depth ([`trim_bottom!`](@ref)), accumulation is expressed as
mass leaving through the base rather than as a rising surface, so in that coordinate a
parcel of firn appears to sink even though it is not moving, and annual layers slope
downward across the panel. `reference_frame` shifts each timestep's cell centres by a
per-timestep offset before regridding, to undo that:

- `:surface` — no offset; the model coordinate, and the historical behaviour.
- `:datum` — add back the cumulative basal ice flux, `-cumsum(ice_flux)`, then subtract
  `detrend·(t-t₀)`. Layers hold their elevation as they should, and the surface rises and
  falls the way a glacier surface does against an earth datum such as the ellipsoid, as an
  anomaly about the mean state.

The first term is not an approximation of elevation; it is elevation, exactly:

```
-cumsum(ice_flux) = SMB / density_ice + Δ(firn_air_content)
```

the standard altimetry decomposition of surface elevation change into a mass term and a
compaction term. The identity holds because the basal flux is whatever it must be to hold
`Σdz` fixed while both terms act on the column (see [`trim_bottom!`](@ref)), and it is
checked numerically in `test/test_ablation_regime.jl`.

The second term is what makes the result an *anomaly*, and it is not cosmetic. A column in
equilibrium with its climate — ice dynamics included — loses as much mass to ice flux over long
averages as it gains from SMB, or the reverse. GEMB has no ice dynamics, and the basal flux it
does perform is a numerical device holding `Σdz` fixed rather than that export, so nothing here
balances the mean SMB and the raw height walks off at the mean accumulation rate. Subtracting
the mean SMB rate supplies the missing equilibrium flux, leaving the variability about the mean
state — which is what an altimeter measures against an earth datum once the long-term mean is
removed.

Note the asymmetry, which follows from the identity: the mean SMB *rate* is subtracted, but
cumulative SMB must **not** be *added*, since it is already one of the two terms of
`-cumsum(ice_flux)` and adding it would count it twice. `ice_flux` is a model construct — what
crosses the domain boundary to hold the column depth fixed — while SMB is a climate-forcing
quantity; the identity relates the two but does not make them interchangeable.

# Returns
A Makie `Figure`.

!!! note "Requires Makie"
    This function lives in a package extension. Load a Makie backend
    (e.g. `using CairoMakie` or `using GLMakie`) before calling it.

# Example
```julia
using GEMB, CairoMakie
output = gemb(profile_spunup, cf, mp)
fig = plot_output(output; depthlims=(-5, 0))
save("gemb_diagnostics.png", fig)
```
"""
function plot_output(args...; kwargs...)
    error("""
        `plot_output` requires a Makie backend. Load one before calling this function:

            using CairoMakie   # for file output (PNG, SVG, PDF)
            using GLMakie      # for interactive windows
            using WGLMakie     # for Jupyter / Pluto notebooks
        """)
end
