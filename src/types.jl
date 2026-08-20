using DimensionalData
using DimensionalData: AbstractDimStack, data_eltype
import DimensionalData as DD
using Dates

"""
    RESTART_LAYERS

The layers that constitute a complete column state — exactly what [`gemb`](@ref) reads off the
`profile` it is handed, and exactly what [`gemb_profile`](@ref) writes back out.

This is the restart contract. Every layer is required: `gemb` indexes them by name to build its
initial state, so a profile missing one fails on the first timestep rather than silently
defaulting. `age` is part of the state even though no physics reads it — it is the column's only
clock, and dropping it would restart every continuation's age at whatever the saved column held.

Exported so a caller that persists a column between runs (writing it to disk and reading it back
later) can validate a restored profile against the set `gemb` actually requires, instead of
against its own copy of this list. A copy cannot detect the case that matters: a layer added to
the state here would leave the caller's list stale in precisely the same way as the old file it is
checking, and the validation would pass.

`z_center` is deliberately absent — `gemb_profile` omits it so an extracted profile matches the
output layout; recompute it with `dz2z(profile[:dz])` where needed.
"""
const RESTART_LAYERS = (:dz, :temperature, :density, :water, :grain_radius,
                        :grain_dendricity, :grain_sphericity, :age)

"""
    DERIVED_PARAMETERS

The [`ModelParameters`](@ref) fields [`gemb`](@ref) computes for itself rather than reading from
the caller, and overwrites on entry with values derived from the forcing timestep.

These are outputs of a run, not settings of it. Exported so a caller comparing two runs'
configurations — deciding whether a saved record may be extended, say — can exclude them and not
report a difference neither run asked for. Recording one would also pin a value the caller never
chose and cannot reproduce without knowing the forcing's timestep.
"""
const DERIVED_PARAMETERS = (:dt_divisors,)

"""
    AbstractThermalSolver

Supertype for the schemes that advance the subsurface temperature profile over one forcing
timestep. The concrete subtype held by `ModelParameters.thermal_solver` selects the scheme by
dispatch, so a new scheme is a new subtype plus one `GEMB._thermal_solve!` method — there is no
central branch to extend.

The interface a subtype must implement is

    _thermal_solve!(solver, temperature, dz, density, K, shortwave_flux,
                    water_surface, grain_radius, sfc, cfs, mp, verbose, workspace)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent,
            heat_flux_basal, evaporation_condensation)

updating `temperature` in place and returning forcing-step averages (the first four in W m-2,
the last in kg m-2 accumulated over the step). Three invariants bind every implementation: the
bottom cell is a Dirichlet reservoir whose temperature must be returned bit-unchanged; the
returned fluxes must be the ones actually applied to the column, since `gemb_core` closes its
energy and mass budgets against them; and every buffer entry drawn from `workspace` must be
written before it is read, since the workspace arrives holding a previous timestep's values.

Subtypes are singletons: `isbits`, immutable, and carrying no state. A scheme's scratch buffers
live in a separate per-run [`ThermalWorkspace`](@ref) rather than in the solver, so a
`ModelParameters` stays a pure, freely shareable description of *what* to compute — one `mp` can
be read by any number of threads stepping columns concurrently. The workspace is the mutable
half, and each concurrent run needs its own.

See [`ExplicitThermal`](@ref) and [`ImplicitThermal`](@ref).
"""
abstract type AbstractThermalSolver end

"""
    _ExplicitWorkspace()

Scratch buffers for `_thermal_solve!(::ExplicitThermal, ...)`. Named for the quantities in that
scheme's algebra; see the solver for what each holds.

Empty on construction and grown to the column by [`_resize_workspace!`](@ref) on first use.
"""
struct _ExplicitWorkspace
    M_inv::Vector{Float64}   # 1 / (cell mass * enthalpy scale)
    H::Vector{Float64}       # cell enthalpy [J m-2]
    Q_sw::Vector{Float64}    # SW absorbed per sub-step [J m-2]
    A_face::Vector{Float64}  # face conductance * dt [J m-2 K-1]
end

_ExplicitWorkspace() = _ExplicitWorkspace(Float64[], Float64[], Float64[], Float64[])

"""
    _ImplicitWorkspace()

Scratch buffers for `_thermal_solve!(::ImplicitThermal, ...)`. All seven are sized to the `m-1`
unknowns (the bottom cell is a Dirichlet reservoir, outside the solve).

Empty on construction and grown to the column by [`_resize_workspace!`](@ref) on first use.
"""
struct _ImplicitWorkspace
    M::Vector{Float64}       # cell mass [kg m-2]
    H::Vector{Float64}       # cell enthalpy [J m-2] — the prognostic
    Q_sw::Vector{Float64}    # SW absorbed per sub-step [J m-2]
    A_face::Vector{Float64}  # face conductance * dt [J m-2 K-1]
    rhs::Vector{Float64}     # right-hand side, then the Thomas solution
    sup::Vector{Float64}     # eliminated superdiagonal
    T_it::Vector{Float64}    # current Newton iterate
end

_ImplicitWorkspace() = _ImplicitWorkspace(ntuple(_ -> Float64[], 7)...)

"""
    _resize_workspace!(ws, n) -> ws

Grow every buffer in `ws` to length `n`, or leave it alone if it is already that long.

`resize!` keeps the underlying capacity, so the column length is reached on the first timestep
and every later call is a no-op. The buffers hold no state between timesteps — each is fully
written before it is read — so the contents after a grow are deliberately undefined.

Callers index by `m`/`n` computed from the column, never by `length(buffer)`, so a buffer that is
longer than needed (a column that shrank) is still correct.
"""
function _resize_workspace!(ws::Union{_ExplicitWorkspace,_ImplicitWorkspace}, n::Int)
    for buffer in _workspace_buffers(ws)
        length(buffer) < n && resize!(buffer, n)
    end
    return ws
end

# Split out so `_resize_workspace!` stays a single method over both workspace types. `fieldnames`
# would need a generated function to stay type-stable; an explicit tuple is clearer and inlines.
_workspace_buffers(ws::_ExplicitWorkspace) = (ws.M_inv, ws.H, ws.Q_sw, ws.A_face)
_workspace_buffers(ws::_ImplicitWorkspace) =
    (ws.M, ws.H, ws.Q_sw, ws.A_face, ws.rhs, ws.sup, ws.T_it)

"""
    ExplicitThermal()

Explicit finite-volume thermal scheme, sub-stepped to the Von Neumann stability limit
(see `GEMB._max_safe_dt`). The default. Conserves energy to the last bit — diffusion is applied
as one flux per face with opposite signs to the two adjacent cells — at the cost of a sub-step
count set by the single stiffest cell in the column.

A stateless singleton; the scheme's scratch buffers live in a [`ThermalWorkspace`](@ref).

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing,
  Ch. 3-4.
"""
struct ExplicitThermal <: AbstractThermalSolver end

"""
    ImplicitThermal()

Backward-Euler thermal scheme on a tridiagonal system, solved by the Thomas algorithm.
Unconditionally stable, so it needs no stability sub-stepping and is insensitive to thin
layers; `mp.dt_divisors` and `GEMB._max_safe_dt` are unused on this path.

A stateless singleton; the scheme's scratch buffers live in a [`ThermalWorkspace`](@ref).

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing,
  Ch. 4.
- Thomas, L. H. (1949). *Elliptic Problems in Linear Difference Equations over a Network*.
  Watson Scientific Computing Laboratory, Columbia University.
"""
struct ImplicitThermal <: AbstractThermalSolver end

"""
    ThermalWorkspace()

Per-run scratch space for the thermal solver, holding the buffers each scheme's algebra needs so
that a timestep allocates nothing after the first.

This is the mutable half of the split that keeps [`ModelParameters`](@ref) shareable: `mp`
describes *what* to compute and can be read concurrently by any number of threads, while a
`ThermalWorkspace` is *where one run scratches* and must not be shared between concurrent runs.

Both schemes' buffer sets are carried, so a workspace is valid for either
[`ExplicitThermal`](@ref) or [`ImplicitThermal`](@ref) and costs nothing for the scheme not in
use (both start empty; only the one actually run is ever grown).

`gemb` creates one per call automatically, so single-threaded use never mentions this type.
Callers stepping columns concurrently should give each thread its own:

```julia
Threads.@threads for i in eachindex(columns)
    ws = ThermalWorkspace()          # one per thread, not one per package
    outputs[i] = gemb(columns[i], forcings[i], mp; thermal_workspace=ws)
end
```

Reusing one across *sequential* runs is fine and saves the first-timestep growth.
"""
struct ThermalWorkspace
    explicit::_ExplicitWorkspace
    implicit::_ImplicitWorkspace
end

ThermalWorkspace() = ThermalWorkspace(_ExplicitWorkspace(), _ImplicitWorkspace())

# Which buffer set a scheme draws from. One method per solver, so a new scheme adds its own
# workspace type and one accessor rather than editing a branch.
_solver_workspace(::ExplicitThermal, ws::ThermalWorkspace) = ws.explicit
_solver_workspace(::ImplicitThermal, ws::ThermalWorkspace) = ws.implicit

"""
    ModelParameters

All GEMB model configuration parameters with validation.
Construct with keyword arguments; unspecified fields use defaults.

Defaults follow the recommendations of the two firn-model intercomparisons (Vandecrux et
al., 2020, RetMIP; Lundin et al., 2017, FirnMICE) and, where those are silent, the shipped
defaults of the Community Firn Model and IMAU-FDM. See the parameter table in `README.md`.

Parameterized on the thermal solver type so that `thermal_solver` is concretely typed and the
scheme is resolved at compile time. `::ModelParameters` in a signature still matches every
instance.

# References

Works cited in the per-field notes below. Each option's own physics is documented, with its
equation numbers, on the function named in that field's comment — this list resolves the
author-year mentions to full entries in one place.

- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010).
  In situ measurements of Antarctic snow compaction compared with predictions of models.
  *J. Geophys. Res.* 115, F03011. (`densification_method`, `grain_growth_method`)
- Barnola, J.-M., Pimienta, P., Raynaud, D., and Korotkevich, Y. S. (1991). CO2-climate
  relationship as deduced from the Vostok ice core: a re-examination based on new measurements
  and on a re-evaluation of the air dating. *Tellus B* 43, 83-90. (`density_ice`)
- Calonne, N., Flin, F., Morin, S., Lesaffre, B., Rolland du Roscoat, S., and Geindreau, C.
  (2011). Numerical and experimental investigations of the effective thermal conductivity of
  snow. *Geophys. Res. Lett.* 38, L23501. (`thermal_conductivity_method = :Calonne`)
- Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F., and
  Geindreau, C. (2019). Thermal conductivity of snow, firn, and porous ice from 3-D
  image-based computations. *Geophys. Res. Lett.* 46, 13079-13089.
  (`thermal_conductivity_method` = `:Calonne2019`/`:Calonne2019Air`, and `density_ice`)
- Colbeck, S. C. (1974). Water flow through snow overlying an impermeable boundary.
  *Water Resour. Res.* 10, 119-123. (`water_irreducible_method = :constant`)
- Coléou, C. and Lesaffre, B. (1998). Irreducible water saturation in snow: experimental
  results in a cold laboratory. *Ann. Glaciol.* 26, 64-68.
  (`water_irreducible_method = :ColeouLesaffre`)
- Cuffey, K. M. and Paterson, W. S. B. (2010). *The Physics of Glaciers*, 4th ed.
  Butterworth-Heinemann. eq. 9.1. (`heat_capacity_method = :CuffeyPaterson`)
- Fausto, R. S., et al. (2018). A snow density dataset for improving surface boundary
  conditions in Greenland ice sheet firn modeling. *Front. Earth Sci.* 6, 51.
  (`new_snow_method` = `:Fausto`/`:FaustoFit`)
- Fourteau, K., Brondex, J., Cancès, C., and Dumont, M. (2026). Numerical strategies for
  representing Richards' equation and its couplings in snowpack models. *Geosci. Model Dev.*
  19, 3193-3212. Sect. 2.3. (`melt_geometry = :density`)
- Gordon, M., Simon, K., and Taylor, P. A. (2006). On snow depth predictions with the Canadian
  land surface scheme including a parametrization of blowing snow sublimation.
  *Atmos.-Ocean* 44, 239-255. (`blowing_snow_sublimation`)
- Gregory, S. A., Albert, M. R., and Baker, I. (2014). Impact of physical properties and
  accumulation rate on pore close-off in layered firn. *The Cryosphere* 8, 91-105.
  (`impermeable_density`)
- Lafaysse, M., et al. (2026). Version 3.0.2 of the Crocus snowpack model.
  *Geosci. Model Dev.* 19, 6273-6334. (`blowing_snow_method`, `blowing_snow_sublimation`,
  `new_snow_method = :Pahaut`)
- Lundin, J. M. D., et al. (2017). Firn Model Intercomparison Experiment (FirnMICE).
  *J. Glaciol.* 63, 401-422.
- Marchenko, S., et al. (2019). Thermal conductivity of firn at Lomonosovfonna, Svalbard,
  derived from subsurface temperature measurements. *The Cryosphere* 13, 1843-1859.
  (`thermal_conductivity_method = :Marchenko2019`)
- Pahaut, E. (1975). *La métamorphose des cristaux de neige*. Monographies de la Météorologie
  Nationale 96, Météo-France. (`new_snow_method = :Pahaut`)
- Sturm, M., Holmgren, J., König, M., and Morris, K. (1997). The thermal conductivity of
  seasonal snow. *J. Glaciol.* 43, 26-41. (`thermal_conductivity_method = :Sturm`)
- Vandecrux, B., et al. (2020). The firn meltwater Retention Model Intercomparison Project
  (RetMIP): evaluation of nine firn models at four weather station sites on the Greenland ice
  sheet. *The Cryosphere* 14, 3785-3810. (`thermal_conductivity_method`,
  `water_irreducible_method`, `impermeable_density`, `runoff_method`)
- Vionnet, V., et al. (2012). The detailed snowpack scheme Crocus and its implementation in
  SURFEX v7.2. *Geosci. Model Dev.* 5, 773-791. (`densification_method`,
  `blowing_snow_method`)
- Zuo, Z. and Oerlemans, J. (1996). Modelling albedo and specific balance of the Greenland ice
  sheet: calculations for the Søndre Strømfjord transect. *J. Glaciol.* 42, 305-317.
  (`runoff_method = :ZuoOerlemans`)
"""
Base.@kwdef struct ModelParameters{S<:AbstractThermalSolver}
    # --- Density & Densification ---
    densification_method::Symbol = :Arthern
    densification_coeffs_M01::Symbol = :Gre_RACMO_GS_SW0
    # Which climatological mass flux the accumulation-driven schemes compact against.
    # `:accumulation` uses `cf.accumulation_mean` (snowfall only) and is the default:
    # `precipitation_mean` includes rain, which does not bury the column, so it overstates the
    # burial rate at any raining site. `:precipitation` compacts against total precipitation.
    # Inert on schemes whose rate does not depend on accumulation (`:ArthernB`,
    # `:Barnola1991`, `:Crocus`, `:CrocusPure`). See `calculate_density`.
    densification_accumulation::Symbol = :accumulation
    # Which mean temperature the Arrhenius rate factors evaluate at. `exp(E/RT)` is convex in
    # 1/T, so `exp(E/R⟨T⟩) ≠ ⟨exp(E/RT)⟩` — at a site with a ±15 K seasonal swing the two
    # differ by tens of percent, biasing toward under-densification. `:arrhenius` uses
    # `cf.temperature_air_effective`, which averages the exponential itself; `:arithmetic`
    # (the default) uses the plain mean `cf.temperature_air_mean`.
    mean_temperature_method::Symbol = :arithmetic
    # Fresh-snow density parameterization. The default `:Constant350` is a constant 350 kg m-3,
    # matching the Community Firn Model's shipped `rhos0: 350.0`. `:Constant150` and
    # `:Constant315` are the other bare constants. `:Fausto` is the Fausto et al.
    # (2018) Greenland fit frozen at a single temperature (315 kg m-3, that fit evaluated at
    # T ≈ 256.2 K); `:FaustoFit` is the same fit carrying its temperature dependence, as
    # IMAU-FDM implements it — though note IMAU-FDM drives it from a one-year running-mean
    # temperature on its Greenland domain, where `:FaustoFit` uses the instantaneous value.
    # Both also select the Crocus wind-dependent fresh-grain properties, as does `:Pahaut`.
    # `:Pahaut` is the alpine-seasonal-snow alternative to those polar fits — Pahaut (1975) as
    # Crocus implements it, temperature- and wind-dependent at the instantaneous timestep, and
    # much lighter over the same range. Prefer it for temperate and mid-latitude glaciers.
    # See `fresh_snow_density`.
    new_snow_method::Symbol = :Constant350
    # Density of pure ice [kg m-3]. 917 is the value used by both the Community Firn Model
    # (`constants.py`) and IMAU-FDM (`constants.toml`), and is the pure-ice density the
    # Calonne (2019) conductivity and Barnola (1991) densification fits were built against.
    density_ice::Float64 = 917.0
    rain_temperature_threshold::Float64 = 273.15

    # --- Longwave Emissivity ---
    emissivity_method::Symbol = :uniform
    emissivity::Float64 = 0.97
    emissivity_grain_radius_large::Float64 = 0.97
    emissivity_grain_radius_threshold::Float64 = 10.0
    surface_roughness_effective_ratio::Float64 = 0.10

    # --- Thermal Conductivity ---
    # `:Sturm` (Sturm et al., 1997) and `:Calonne` (Calonne et al., 2011) are density-only
    # quadratics. `:Calonne2019`, `:Calonne2019Air` and `:Marchenko2019` are the
    # temperature-dependent forms Vandecrux et al. (2020, RetMIP Sects. 5.1 and 7) recommend
    # over the older fits, having attributed part of a multi-model cold bias at Summit to
    # conductivity parameterization. `:Calonne2019` is the default and is what the Community
    # Firn Model ships (`conductivity: "Calonne2019"`); IMAU-FDM hardcodes the same equation.
    # `:Sturm` gives the lowest conductivity of the five — 0.56 against `:Calonne2019`'s 0.91
    # W m-1 K-1 at 243 K and ρ = 550 — so selecting it damps the seasonal wave.
    # The two 2019 forms differ only in whether the snow branch carries the air-conductivity
    # ratio (`:Calonne2019Air` does, matching IMAU-FDM; `:Calonne2019` does not, matching the
    # Community Firn Model), which matters only for cold, low-density snow. See
    # `thermal_conductivity`.
    thermal_conductivity_method::Symbol = :Calonne2019

    # --- Heat Capacity ---
    # Specific heat capacity of the load-bearing ice matrix, applied to snow, firn, and ice
    # alike (pore air and pore water contribute nothing to this term). `:constant` (the default)
    # uses `heat_capacity_ice` at every temperature, matching the Community Firn Model's
    # `enthalpyDiff` melt path, which sets `c_firn = CP_I = 2097`. `:CuffeyPaterson` makes it
    # temperature-dependent (`152.5 + 7.122·T`, as IMAU-FDM uses throughout) and ignores
    # `heat_capacity_ice`. See `heat_capacity`.
    heat_capacity_method::Symbol = :constant
    heat_capacity_ice::Float64 = HEAT_CAPACITY_ICE_DEFAULT
    # Heat capacity carrying the *sensible* heat of above-freezing liquid entering the column,
    # which only rain is. `:water` (the default) uses `HEAT_CAPACITY_WATER`; `:ice` uses the
    # matrix value, understating it roughly twofold. Everything else in the model
    # — pore water, refreezing, runoff — is isothermal at the melting point and unaffected by
    # this field. See `specific_enthalpy_water`.
    rain_heat_capacity::Symbol = :water

    # --- Grain Growth ---
    # Selects the non-dendritic *dry* grain-growth law. `:Arthern` (the default, and what the
    # Community Firn Model ships as `GrGrowPhysics: "Arthern"`) uses Arthern et al. (2010)
    # `dr²/dt = kgr·exp(-Eg/RT)` at every density. `:Marbouty` stops growing grains at
    # `DENSITY_MARBOUTY_MAX`, which leaves the radius frozen in deep firn — where `:ArthernB`
    # and `:Crocus` densification read it as `1/r²`. `:hybrid` switches from the latter to the
    # former at that density. See `calculate_grain_size`.
    grain_growth_method::Symbol = :Arthern

    # --- Melt & Water ---
    # `:ColeouLesaffre` (the default) makes the irreducible water content density-dependent
    # after Coléou & Lesaffre (1998) and ignores `water_irreducible_saturation`. It is what the
    # Community Firn Model ships (`ColeouLesaffre: true`) and what IMAU-FDM implements, and
    # RetMIP Sect. 5.4 found it gave realistic percolation depths at Dye-2. `:constant` holds
    # `water_irreducible_saturation` at every density (Colbeck, 1974). See
    # `irreducible_saturation`.
    water_irreducible_method::Symbol = :ColeouLesaffre
    water_irreducible_saturation::Float64 = 0.07
    # What melting does to a cell's geometry. Refreezing is at constant thickness under both
    # settings; this field only concerns the melt half of the phase change.
    #
    # `:thickness` (the default) holds density fixed and shrinks `dz`, which is Crocus's
    # treatment and what GEMB has always done. `:density` holds `dz` fixed and lowers density,
    # which is SNOWPACK's treatment and the one Fourteau et al. (2026) Sect. 2.3 argues for:
    # the phase change happens *within* the microstructure, so it removes ice from the pore
    # walls rather than collapsing the layer, and the high density of wet snow is better
    # explained by its low viscosity under overburden than by melt-driven thinning. `:density`
    # also removes an asymmetry — under `:thickness` refreezing raises density at constant
    # thickness while melting does not lower it, so a melt-refreeze cycle that returns the
    # cell's mass does not return its geometry. See `calculate_melt`.
    melt_geometry::Symbol = :thickness
    # The density-based impermeability criterion of the bucket scheme: percolating water is
    # routed to runoff at a contiguous run of cells at or above `impermeable_density` that is
    # thicker than `impermeable_thickness`. Below that thickness water passes through, on the
    # grounds that thin lenses are laterally discontinuous at the scale the model represents.
    #
    # The defaults (830, 0.1) are exactly what the Community Firn Model ships as `RhoImp` and
    # `ThickImp`. Both are genuinely uncertain and are the dominant control
    # on bucket-scheme behaviour at ice-slab sites — across the RetMIP models (Vandecrux et
    # al., 2020) the density threshold ranges from 810 kg m-3 (DMIHH, after Gregory et al.,
    # 2014, which gave the lowest firn-temperature RMSE at KAN_U) to 917 (DTU, which
    # over-produced runoff badly enough to be excluded from the paper's multi-model mean).
    # Exposed here so that range can be explored. See `calculate_melt`.
    impermeable_density::Float64 = DENSITY_PORE_CLOSEOFF
    impermeable_thickness::Float64 = 0.1
    # How water that percolation cannot pass leaves the column.
    #
    # `:instantaneous` is the default — a plain bucket scheme, as both the Community Firn Model
    # (`liquid: "bucket"`, `Ponding: false`) and IMAU-FDM ship, and which RetMIP Sect. 5.2 found
    # performs no worse than the deep-percolation schemes: blocked water runs off within the
    # timestep, and no cell can ever hold more than its irreducible water. The other two let it
    # pond above the barrier (see `calculate_melt`) and drain laterally over a finite timescale
    # (see `apply_lateral_drainage!`), which is what makes standing water and firn aquifers
    # representable at all. Both are scaled by `surface_slope`, and both were used by the two
    # models with the lowest firn-temperature error at the RetMIP ice-slab site KAN_U
    # (Vandecrux et al., 2020): `:ZuoOerlemans` by DMIHH, `:Darcy` by GEUS.
    runoff_method::Symbol = :instantaneous
    # Cap on how full a cell may get, as a fraction of its pore space. 1.0 permits full
    # saturation. Only read when `runoff_method !== :instantaneous`. See `pore_capacity`.
    pore_saturation_max::Float64 = 1.0

    # --- Albedo & Radiation ---
    albedo_method::Symbol = :GardnerSharp
    albedo_density_threshold::Float64 = Inf
    shortwave_subsurface_absorption::Bool = false
    albedo_snow::Float64 = 0.85
    albedo_ice::Float64 = 0.48
    albedo_fixed::Float64 = 0.85
    shortwave_downward_diffuse::Float64 = 0.0
    solar_zenith_angle::Float64 = 0.0
    cloud_optical_thickness::Float64 = 0.0
    black_carbon_snow::Float64 = 0.0
    black_carbon_ice::Float64 = 0.0
    cloud_fraction::Float64 = 0.1

    # --- Output Controls ---
    output_frequency::Symbol = :all
    # Age epoch of an initialized column. `:steady_state` inherits the residence time
    # `steady_state_profile` marched; `:zero` starts every cell at 0 (the instant of
    # initialization is the epoch). No physics reads `age`, so this is output-only.
    # See `initialize_profile`.
    initialize_age::Symbol = :steady_state
    # Add a `viscosity` profile layer [Pa s] to the `gemb` output. Off by default: it costs
    # one length-`m` allocation per timestep, and it is populated only under
    # `densification_method = :Crocus`/`:CrocusPure`, the only schemes that form an effective
    # viscosity. Under any other scheme the layer exists but is all `NaN` — see
    # `calculate_density`.
    output_viscosity::Bool = false

    # --- Grid Geometry ---
    column_ztop::Float64 = 10.0
    column_dztop::Float64 = 0.05
    column_dzmin::Float64 = 0.025
    column_dzmax::Float64 = 0.075
    # A *ceiling* on the depth of the constructed grid, not the depth itself:
    # `initialize_profile` sizes the column to the depth the forcing climate needs
    # (`_derive_column_depth`) and clamps it here, so an ablation site is not made to carry
    # hundreds of metres of solid ice. The depth actually chosen is `sum(profile[:dz])`,
    # which is what stays fixed for the run — nothing downstream reads this field.
    column_depth_max::Float64 = 250.0
    column_zy::Float64 = 1.10

    # --- Ice Dynamics ---
    # Trace of the horizontal strain-rate tensor, ε̇_xx + ε̇_yy [yr-1]. Positive is
    # horizontal divergence, which thins every layer at constant density; negative is
    # convergence, which thickens. 0.0 (the default) disables the term entirely and leaves
    # output bit-identical to a run without it. See `apply_horizontal_strain!`.
    horizontal_strain_rate::Float64 = 0.0
    # Surface slope [m m-1], the lateral hydraulic gradient driving runoff under
    # `runoff_method` `:ZuoOerlemans` or `:Darcy`; unread under `:instantaneous`. The RetMIP
    # sites span 0 (Summit) to 0.6 degrees (FA), i.e. 0 to ~0.010 m m-1 — note the units are
    # a gradient, not degrees.
    surface_slope::Float64 = 0.0

    # --- Blowing Snow ---
    # Which scheme, if any, lets wind rework snow that has already landed. `:none` (the
    # default) disables the term entirely and leaves output bit-identical to a run without it.
    # `:Crocus` is the SURFEX/Crocus `SNOWDRIFT` scheme (Vionnet et al., 2012; Lafaysse et al.,
    # 2026 eqs. 59-66): near-surface layers are compacted toward `DRIFT_DENSITY_MAX` and their
    # grains fragmented at a rate set by a wind-driven driftability index that decays with
    # depth. Mass-conserving on its own — `dz` shrinks in step with the density rise.
    # Independent of `drift_rate`, which is a prescribed mass sink rather than a scheme, and of
    # `blowing_snow_sublimation`, which is the only part that removes mass. This mirrors
    # SURFEX's two independent namelist flags (`OSNOWDRIFT`, `OSNOWDRIFT_SUBLIM`). See
    # `apply_blowing_snow!`.
    blowing_snow_method::Symbol = :none
    # Whether the blowing-snow scheme also sublimates suspended snow after Gordon et al.
    # (2006), as Crocus implements it (Lafaysse et al., 2026 eqs. 67-68). Unlike the compaction
    # path this removes mass from the surface cell and consumes latent heat, so it appears in
    # the mass and energy budgets. `false` (the default) matches Crocus's own default. Read
    # only when `blowing_snow_method !== :none`. See `blowing_snow_sublimation_rate`.
    blowing_snow_sublimation::Bool = false
    # Prescribed net drift divergence [kg m-2 yr-1], positive for erosion (mass leaves) and
    # negative for deposition. The IMAU-FDM path: that model reads `SnowDrif` from RACMO and
    # applies it as a surface mass sink, eroding at the surface cell's own density and
    # depositing at the fresh-snow density. Use it when a drift product is available; it is
    # independent of `blowing_snow_method` and can be combined with it, though doing so risks
    # double-counting if the product already includes the compaction the scheme computes.
    # Superseded per-timestep by the optional `snow_drift` forcing layer when that is present.
    # 0.0 (the default) disables it. See `apply_prescribed_drift!`.
    drift_rate::Float64 = 0.0

    # --- Thermal Solver ---
    # Which scheme advances the temperature profile. Selected by type, not by symbol, so the
    # choice is resolved at compile time and a new scheme needs no branch here — see
    # `AbstractThermalSolver`. `ExplicitThermal` is the default and the only scheme whose
    # output the regression fingerprints are pinned to.
    thermal_solver::S = ExplicitThermal()

    # --- Thermal Time Stepping ---
    # Read by `ExplicitThermal` only. That scheme's sub-step comes from the stability limit of
    # the scheme actually being solved, whose face conductances are harmonic means over two
    # half-cells: `dt ≤ ρᵢcᵢdzᵢ / (Gᵢ + Gᵢ₋₁)` — not the textbook uniform-grid form
    # `0.5·ρ·c·dz²/K`, which is not conservative on a graded grid. See `_max_safe_dt`.
    # `ImplicitThermal` is unconditionally stable and never reads this.
    dt_divisors::Vector{Float64} = Float64[]  # pre-computed divisors for thermo sub-stepping; set by gemb driver
    # Fraction of `_max_safe_dt` the explicit sub-step is allowed to use. `ExplicitThermal`
    # only. 1.0 would run at the diffusive limit exactly; the default is below it because
    # `_max_safe_dt` covers diffusion alone and the surface cell carries a second, unbounded
    # feedback — see `THERMAL_EXPLICIT_SAFETY_FACTOR`, which documents what the margin is for
    # and what was measured before choosing it.
    thermal_explicit_safety_factor::Float64 = THERMAL_EXPLICIT_SAFETY_FACTOR
end

"""
    ClimateForcing <: DimensionalData.AbstractDimStack

Time-series surface meteorological forcing for GEMB, implemented as an
`AbstractDimStack`. It therefore supports the full DimensionalData stack API —
`keys`, `length`, `dims`, `map`, `layers`, and DimStack-style indexing such as
date-range subsetting `cf[Ti(a .. b)]` and single-time selection `cf[Ti(At(t))]`
— all returning a `ClimateForcing` (a sub-stack), never a scalar step.

# Layers (14 time-series `DimArray`s sharing a common `Ti` dimension)
Required forcing: `temperature_air`, `pressure_air`, `precipitation`,
`wind_speed`, `shortwave_downward`, `longwave_downward`, `vapor_pressure`.
Time-varying model parameters (typically `Fill`-backed): `black_carbon_snow`,
`black_carbon_ice`, `cloud_optical_thickness`, `solar_zenith_angle`,
`shortwave_downward_diffuse`, `cloud_fraction`.
Optional forcing: `snow_drift` [kg m-2 yr-1, positive for erosion], the prescribed
drifting-snow divergence read by [`apply_prescribed_drift!`](@ref); `Fill`-backed zeros
when the forcing carries no drift product, and superseded by `mp.drift_rate` in that case.

# Metadata (scalars, carried in the stack `metadata` as a `NamedTuple`)
`time_step::Int` [s], plus `temperature_air_mean`, `wind_speed_mean`,
`precipitation_mean`, `temperature_observation_height`, `wind_observation_height`,
source provenance (`dataset`, `latitude`, `longitude`, `elevation_offset`) and
climatology provenance (`climatology_window_start`, `climatology_window_stop`,
`climatology_n_years`, `climatology_steps_per_year`; set by
[`forcing_climatology`](@ref)).

Both layers and scalar metadata are reachable by name via property access
(`cf.temperature_air` returns the layer `DimArray`; `cf.time_step` returns the
concrete scalar), so existing `cf.<field>` code continues to work unchanged.
Construct via the 19-argument positional constructor (13 `DimArray`s then the 6
scalars, in the order listed above) or, more commonly, via [`initialize_forcing`](@ref).
"""
struct ClimateForcing{K,T,N,L,D<:Tuple,R<:Tuple,LD,M,LM} <: AbstractDimStack{K,T,N,L,D}
    data::L
    dims::D
    refdims::R
    layerdims::NamedTuple{K,LD}
    metadata::M
    layermetadata::NamedTuple{K,LM}
    # Inner constructors mirror `DimStack` so the generic AbstractDimStack
    # `rebuild`/`stacktype` machinery (which reconstructs via
    # `basetypeof(s){K,...}(data, dims, refdims, layerdims, metadata, layermetadata)`)
    # finds a matching signature.
    function ClimateForcing(
        data, dims, refdims, layerdims::LD, metadata, layermetadata::NamedTuple{K}
    ) where LD <: NamedTuple{K} where K
        T = data_eltype(data)
        N = length(dims)
        ClimateForcing{K,T,N}(data, dims, refdims, layerdims, metadata, layermetadata)
    end
    function ClimateForcing{K,T,N}(
        data::L, dims::D, refdims::R, layerdims::NamedTuple{K,LD}, metadata::M, layermetadata::NamedTuple{K,LM}
    ) where {K,T,N,L,D,R,LD,M,LM}
        new{K,T,N,L,D,R,LD,M,LM}(data, dims, refdims, layerdims, metadata, layermetadata)
    end
end

# The 14 forcing layers. The first 13 are in positional-constructor order; `snow_drift` is
# keyword-only and trailing, so every existing 19-argument positional call site stays valid
# and forcing built without a drift product carries a `Fill`-backed zero layer.
const CLIMATE_FORCING_LAYER_KEYS = (
    :temperature_air, :pressure_air, :precipitation, :wind_speed,
    :shortwave_downward, :longwave_downward, :vapor_pressure,
    :black_carbon_snow, :black_carbon_ice, :cloud_optical_thickness,
    :solar_zenith_angle, :shortwave_downward_diffuse, :cloud_fraction,
    :snow_drift,
)
# The scalar metadata carried alongside the layers. The first six are the physics
# scalars (positional in the constructor); the trailing four are provenance
# (passed as keyword args) that trace where the forcing came from and how it was
# adjusted, so they can be surfaced in `gemb` output and diagnostic plots.
const CLIMATE_FORCING_META_KEYS = (
    :time_step, :temperature_air_mean, :wind_speed_mean, :precipitation_mean,
    :temperature_observation_height, :wind_observation_height,
    :accumulation_mean, :temperature_air_effective,
    :dataset, :latitude, :longitude, :elevation, :elevation_native, :elevation_offset,
    :climatology_window_start, :climatology_window_stop,
    :climatology_n_years, :climatology_steps_per_year,
)

"""
    ClimateForcing(temperature_air, pressure_air, precipitation, wind_speed,
                   shortwave_downward, longwave_downward, vapor_pressure,
                   black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
                   solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
                   time_step, temperature_air_mean, wind_speed_mean,
                   precipitation_mean, temperature_observation_height,
                   wind_observation_height)

Positional constructor: 13 forcing `DimArray`s (sharing a `Ti` dimension) followed
by the 6 scalar metadata values. Builds the stack layers and stores the scalars in
the stack `metadata` as a `NamedTuple` (keeping their concrete types).

The 14th layer, `snow_drift` [kg m-2 yr-1, positive for erosion], is keyword-only and
defaults to `Fill`-backed zeros, so a forcing built without a drift product is valid and
inert. Supply it to drive [`apply_prescribed_drift!`](@ref) from a drift product
(e.g. RACMO's `SnowDrif`) rather than from the constant `mp.drift_rate`.

Optional provenance keywords — `dataset` (source name, e.g. `"ERA5-Land"`),
`latitude`, `longitude`, `elevation` (absolute target elevation, m; the surface the
forcing represents after any elevation adjustment), `elevation_native` (the source
dataset's own surface elevation, m, before any adjustment), and `elevation_offset`
(m the forcing was elevation-adjusted by) — are stored alongside the physics scalars
and default to `""`/`NaN`/`NaN`/`NaN`/`NaN`/`0.0`.

Two optional refinements of the physics scalars, both defaulting to `NaN`, in which case
they fall back to the scalar they refine so the forcing behaves exactly as one built
without them:

- `accumulation_mean` [kg m-2 yr-1] — mean *snowfall*, i.e. `precipitation_mean` with rain
  excluded on `mp.rain_temperature_threshold`. Read by densification under
  `mp.densification_accumulation = :accumulation`; falls back to `precipitation_mean`.
- `temperature_air_effective` [K] — the Arrhenius-weighted mean temperature,
  `Eg / (R·log(mean(exp(Eg/(R·T)))))`, which is what the grain-growth Arrhenius factors
  actually need (`exp(E/R⟨T⟩) ≠ ⟨exp(E/RT)⟩`). Read under
  `mp.mean_temperature_method = :arrhenius`; falls back to `temperature_air_mean`.

`elevation_offset` is forced to `0.0` when `elevation` is not finite: with no known
target elevation there is nothing an offset could be measured against, and carrying a
bare offset invites downstream code to derive a native elevation from `NaN - offset`.

Climatology provenance keywords — `climatology_window_start`,
`climatology_window_stop` (the requested averaging window, `DateTime` or
`nothing`), `climatology_n_years` (number of complete years averaged), and
`climatology_steps_per_year` — are set by [`forcing_climatology`](@ref) and
default to `nothing`/`0` for non-climatological forcing.
"""
function ClimateForcing(
    temperature_air, pressure_air, precipitation, wind_speed,
    shortwave_downward, longwave_downward, vapor_pressure,
    black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
    solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
    time_step, temperature_air_mean, wind_speed_mean, precipitation_mean,
    temperature_observation_height, wind_observation_height;
    snow_drift=nothing,
    accumulation_mean::Real=NaN, temperature_air_effective::Real=NaN,
    dataset::AbstractString="", latitude::Real=NaN, longitude::Real=NaN,
    elevation::Real=NaN, elevation_native::Real=NaN, elevation_offset::Real=0.0,
    climatology_window_start=nothing, climatology_window_stop=nothing,
    climatology_n_years::Integer=0, climatology_steps_per_year::Integer=0,
)
    # A `Fill`-backed zero layer when no drift product was supplied, on the `Ti` dimension the
    # required layers already carry — the same treatment the time-varying model parameters get.
    if snow_drift === nothing
        tdim = dims(temperature_air, Ti)
        snow_drift = DimArray(Fill(0.0, length(tdim)), (tdim,))
    end
    das = NamedTuple{CLIMATE_FORCING_LAYER_KEYS}((
        temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure,
        black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
        solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
        snow_drift,
    ))
    meta = NamedTuple{CLIMATE_FORCING_META_KEYS}((
        time_step, temperature_air_mean, wind_speed_mean, precipitation_mean,
        temperature_observation_height, wind_observation_height,
        # Both fall back to the total-precipitation / arithmetic-mean scalars they refine, so
        # a forcing built without them behaves exactly as before and the two
        # `ModelParameters` gates that read them are inert rather than reading a NaN.
        isnan(accumulation_mean) ? Float64(precipitation_mean) : Float64(accumulation_mean),
        isnan(temperature_air_effective) ? Float64(temperature_air_mean) :
            Float64(temperature_air_effective),
        String(dataset), Float64(latitude), Float64(longitude),
        Float64(elevation), Float64(elevation_native),
        # No target elevation ⇒ no meaningful offset (see docstring).
        isfinite(elevation) ? Float64(elevation_offset) : 0.0,
        climatology_window_start, climatology_window_stop,
        Int(climatology_n_years), Int(climatology_steps_per_year),
    ))
    stackdims = DD.combinedims(collect(das))
    data = map(parent, das)
    layerdims = map(DD.basedims, das)
    layermetadata = map(DD.metadata, das)
    ClimateForcing(data, stackdims, (), layerdims, meta, layermetadata)
end

# Preserve the `cf.<field>` interface on top of the stack layout: forcing fields
# resolve to their layer `DimArray`, scalar metadata to the concrete stored value,
# and everything else (`data`, `dims`, `refdims`, ...) to the underlying struct field.
Base.@constprop :aggressive function Base.getproperty(cf::ClimateForcing, k::Symbol)
    if k in CLIMATE_FORCING_LAYER_KEYS
        return cf[k]
    elseif k in CLIMATE_FORCING_META_KEYS
        return getfield(cf, :metadata)[k]
    else
        return getfield(cf, k)
    end
end

"""
    ClimateForcingStep

Single time-step forcing values extracted from ClimateForcing.
Plain struct of scalars for efficient use in the physics loop.
"""
struct ClimateForcingStep
    dt::Float64
    temperature_air::Float64
    pressure_air::Float64
    precipitation::Float64
    wind_speed::Float64
    shortwave_downward::Float64
    longwave_downward::Float64
    vapor_pressure::Float64
    # Metadata
    temperature_air_mean::Float64
    wind_speed_mean::Float64
    precipitation_mean::Float64
    temperature_observation_height::Float64
    wind_observation_height::Float64
    # Time-varying model params (from forcing or ModelParam defaults)
    black_carbon_snow::Float64
    black_carbon_ice::Float64
    cloud_optical_thickness::Float64
    solar_zenith_angle::Float64
    shortwave_downward_diffuse::Float64
    cloud_fraction::Float64
    # Refinements of `precipitation_mean` and `temperature_air_mean`, selected by
    # `mp.densification_accumulation` and `mp.mean_temperature_method`. Last, not beside the
    # scalars they refine, so the 19-argument positional form stays valid; the two-argument
    # tail constructor below defaults each to the scalar it refines, which is what makes an
    # existing positional call site behave exactly as it did.
    accumulation_mean::Float64
    temperature_air_effective::Float64
    # Prescribed drifting-snow divergence [kg m-2 yr-1], positive for erosion. Last for the
    # same reason the two refinements above are: it keeps the 19- and 21-argument positional
    # forms valid. Defaults to 0.0 in both tail constructors, which is inert in
    # `apply_prescribed_drift!`.
    snow_drift::Float64
end

# 19-argument positional form: the two refinements default to the scalars they refine, so a
# step built this way reads identically on either setting of the two gates, and `snow_drift`
# defaults to zero (no drift).
ClimateForcingStep(dt, temperature_air, pressure_air, precipitation, wind_speed,
    shortwave_downward, longwave_downward, vapor_pressure,
    temperature_air_mean, wind_speed_mean, precipitation_mean,
    temperature_observation_height, wind_observation_height,
    black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
    solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction) =
    ClimateForcingStep(dt, temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure,
        temperature_air_mean, wind_speed_mean, precipitation_mean,
        temperature_observation_height, wind_observation_height,
        black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
        solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
        precipitation_mean, temperature_air_mean, 0.0)

# 21-argument positional form: both refinements given explicitly, `snow_drift` still zero.
ClimateForcingStep(dt, temperature_air, pressure_air, precipitation, wind_speed,
    shortwave_downward, longwave_downward, vapor_pressure,
    temperature_air_mean, wind_speed_mean, precipitation_mean,
    temperature_observation_height, wind_observation_height,
    black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
    solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
    accumulation_mean, temperature_air_effective) =
    ClimateForcingStep(dt, temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure,
        temperature_air_mean, wind_speed_mean, precipitation_mean,
        temperature_observation_height, wind_observation_height,
        black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
        solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
        accumulation_mean, temperature_air_effective, 0.0)

"""
    ClimateForcingStep(; dt, temperature_air, ...)

Keyword constructor, defaulting every field to zero.

The physics loop builds these positionally from the forcing vectors, but the
initialization helpers each need only the handful of fields their callee reads
(grain growth reads `dt`; the surface energy balance reads the air state; albedo
reads the optical properties). Naming those fields keeps a field insertion or
reorder in the struct above from silently mis-assigning values — every field is a
`Float64`, so a positional mistake is not a type error.
"""
function ClimateForcingStep(;
    dt::Real=0.0,
    temperature_air::Real=0.0, pressure_air::Real=0.0, precipitation::Real=0.0,
    wind_speed::Real=0.0, shortwave_downward::Real=0.0, longwave_downward::Real=0.0,
    vapor_pressure::Real=0.0,
    temperature_air_mean::Real=0.0, wind_speed_mean::Real=0.0,
    precipitation_mean::Real=0.0, temperature_observation_height::Real=0.0,
    wind_observation_height::Real=0.0,
    black_carbon_snow::Real=0.0, black_carbon_ice::Real=0.0,
    cloud_optical_thickness::Real=0.0, solar_zenith_angle::Real=0.0,
    shortwave_downward_diffuse::Real=0.0, cloud_fraction::Real=0.0,
    accumulation_mean::Real=precipitation_mean,
    temperature_air_effective::Real=temperature_air_mean,
    snow_drift::Real=0.0)

    return ClimateForcingStep(
        dt, temperature_air, pressure_air, precipitation, wind_speed,
        shortwave_downward, longwave_downward, vapor_pressure,
        temperature_air_mean, wind_speed_mean, precipitation_mean,
        temperature_observation_height, wind_observation_height,
        black_carbon_snow, black_carbon_ice, cloud_optical_thickness,
        solar_zenith_angle, shortwave_downward_diffuse, cloud_fraction,
        accumulation_mean, temperature_air_effective, snow_drift,
    )
end
