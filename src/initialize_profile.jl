"""
    initialize_profile(mp::ModelParameters, cf::ClimateForcing;
                       steady_state=true, depth_autoadjust=true,
                       ddf_snow=3.0, melt_accum_ratio=0.6)
        -> (profile::DimStack, mp::ModelParameters)

Initialize a GEMB firn column profile as a DimStack.

Returns a tuple `(profile, mp)`. The returned `mp` equals the input parameters
unless `depth_autoadjust` shrinks the grid depth (see below), in which case it is a
modified copy — pass this returned `mp` to [`gemb`](@ref) / [`gemb_spinup`](@ref) so the
adjusted `column_depth` is honored at runtime as well as in the initial grid.

By default the density profile is initialized to the **Herron–Langway steady
state** ([`herron_langway_steady_state`](@ref)) derived from the mean annual
temperature (`cf.temperature_air_mean`) and accumulation (`cf.precipitation_mean`)
carried on the forcing. This starts a column close to equilibrium so
[`gemb_spinup`](@ref) converges in far fewer cycles. The surface density is taken
from [`fresh_snow_density`](@ref) for the configured `mp.new_snow_method`.

Where melt is expected to dominate — annual potential melt (a positive-degree-day
estimate, [`annual_pdd_melt`](@ref)) reaches `melt_accum_ratio` of the annual
accumulation — or accumulation is non-positive, no firn forms and the column is
initialized to **pure ice** (`mp.density_ice`), matching the historical behavior.
An ice column is also initialized with non-dendritic, faceted, large grains
(`grain_dendricity = 0`, `grain_sphericity = 0`, `grain_radius = 2.5` mm) and ice
albedo (`mp.albedo_ice`); the firn/snow regime keeps fresh-snow grains and
`mp.albedo_snow`.
The Herron–Langway / Arthern firn models are calibrated for the cold, dry-snow
zone (mean-annual T ≈ 247–257 K, accumulation ≈ 130–1040 kg m-2 yr-1; Arthern et
al. 2010) and are not valid under significant melt.

Temperature is initialized uniformly to the mean annual temperature, clamped to
the melt point (`min(T_mean, 273.15)`).

# Keyword arguments
- `steady_state`: use the Herron–Langway steady-state density profile (default
  `true`); set `false` to force the legacy pure-ice initialization.
- `depth_autoadjust`: when `true` (default) and the column is inferred to be ice
  (the ablation-zone / pure-ice regime — same condition that fills the column with
  `mp.density_ice`), shrink the grid to `column_depth = 25 m` (a deep firn grid is
  unnecessary for solid ice). The adjusted value is returned in `mp`. Set `false` to
  keep the configured depth (required for MATLAB-fidelity and regression tests).
- `constant_density`: when `true`, initialize the density profile to pure ice
  (`mp.density_ice`) everywhere, bypassing the steady-state/ablation regime
  logic (default `false`).
- `constant_temperature`: when `true`, initialize the temperature profile to the
  mean annual air temperature (`cf.temperature_air_mean`) everywhere, without
  clamping to the melt point (default `false`).
- `ddf_snow`: snow degree-day factor [mm w.e. °C-1 d-1] for the melt estimate
  (default 3.0).
- `melt_accum_ratio`: initialize pure ice when `melt ≥ melt_accum_ratio × accum`
  (default 0.6).

Setting both `constant_density=true` and `constant_temperature=true` reproduces
the MATLAB `model_initialize_profile` initialization (pure ice, uniform
mean-annual temperature).

Returns a DimStack with Z dimension containing:
- dz, temperature, density, water, grain_radius,
  grain_dendricity, grain_sphericity, albedo, albedo_diffuse
"""
function initialize_profile(mp::ModelParameters, cf::ClimateForcing;
    steady_state::Bool=true, depth_autoadjust::Bool=true, constant_density::Bool=false,
    constant_temperature::Bool=false, ddf_snow::Real=3.0, melt_accum_ratio::Real=0.6)

    T_mean = cf.temperature_air_mean

    @assert T_mean > 0 "temperature_air_mean must exceed 0 K."
    if T_mean < 100
        @warn "temperature_air_mean should be in kelvin, but is below 100, suggesting an error."
    end

    # Decide firn (steady-state density) vs. ablation/ice (pure ice) regime. This
    # is computed before the grid is built so it can also drive depth_autoadjust.
    accum = cf.precipitation_mean                        # [kg m-2 yr-1]
    melt = annual_pdd_melt(cf; ddf_snow=ddf_snow)        # [kg m-2 yr-1]
    is_ice = constant_density || !steady_state || accum <= 0.0 || melt >= melt_accum_ratio * accum

    # For an ice column there is no deep firn to resolve, so shrink the grid depth. The
    # rebuilt mp is returned because column_depth fixes the run's column depth.
    if depth_autoadjust && is_ice
        mp = ModelParameters(;
            (f => getfield(mp, f) for f in fieldnames(ModelParameters)
                 if f != :column_depth)...,
            column_depth=25.0)
    end

    # Initialize grid
    dz = initialize_grid(mp)
    z_center = dz2z(dz)
    m = length(dz)

    # Temperature: mean annual, clamped to the melt point (matters in the warm /
    # ablation regime where T_mean approaches or exceeds 273.15 K). When
    # constant_temperature=true, use the unclamped mean annual temperature
    # everywhere (matches the MATLAB initialization).
    T_init = constant_temperature ? T_mean : min(T_mean, CtoK)

    if is_ice
        density = fill(mp.density_ice, m)
    else
        # Warn if outside the Arthern et al. (2010) dry-firn calibration envelope.
        if T_mean > 258.0 || accum < 100.0 || accum > 1100.0
            @warn "initialize_profile: mean T ($(round(T_mean, digits=1)) K) / accumulation " *
                  "($(round(accum, digits=1)) kg m-2 yr-1) is outside the Herron–Langway / " *
                  "Arthern (2010) dry-firn calibration envelope (T≈247–257 K, b≈130–1040); " *
                  "steady-state density profile may be inaccurate."
        end
        ρ0 = fresh_snow_density(mp, T_mean, accum, cf.wind_speed_mean)
        density = herron_langway_steady_state(z_center, T_init, accum, ρ0, mp.density_ice)
    end

    # Create Z dimension
    zdim = Z(1:m; metadata=cf_layer_index_attributes())

    # Grain and albedo initialization is regime-dependent. An ice column is
    # non-dendritic (dendricity 0), faceted (sphericity 0), large-grained
    # (grain_radius 2.5 mm = 5 mm diameter, the non-spherical cap in
    # calculate_grain_size), and low-albedo (mp.albedo_ice). Firn/snow starts as
    # fresh dendritic snow with snow albedo. grain_radius is 2.5 mm in both cases.
    if is_ice
        gdn0, gsp0, re0, alb0 = 0.0, 0.0, 2.5, mp.albedo_ice
    else
        gdn0, gsp0, re0, alb0 = 1.0, 0.5, 2.5, mp.albedo_snow
    end

    # Note: `z_center` is intentionally not stored — it is a pure function of `dz`
    # (via `dz2z`) and is recomputed on demand, so the profile stays a clean slice
    # of the `gemb` output layout (which likewise carries `dz`, not `z_center`).
    layers = (
        dz=DimArray(dz, (zdim,)),
        temperature=DimArray(fill(T_init, m), (zdim,)),
        density=DimArray(density, (zdim,)),
        water=DimArray(zeros(m), (zdim,)),
        grain_radius=DimArray(fill(re0, m), (zdim,)),
        grain_dendricity=DimArray(fill(gdn0, m), (zdim,)),
        grain_sphericity=DimArray(fill(gsp0, m), (zdim,)),
        albedo=DimArray(fill(alb0, m), (zdim,)),
        albedo_diffuse=DimArray(fill(alb0, m), (zdim,)),
    )

    # Same CF attributes a profile extracted from `gemb` output carries (`gemb_profile`),
    # so a profile describes itself the same way whichever end of the pipeline it came
    # from. No `cell_methods`: there is no time dimension here.
    profile = DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))

    return profile, mp
end

"""
    initialize_grid(mp::ModelParameters)

Generate the initial vertical grid layer thicknesses.
Matches MATLAB's `model_initialize_grid` (local function in model_initialize_profile.m).

Returns a Vector{Float64} of layer thicknesses from surface to depth.
"""
function initialize_grid(mp::ModelParameters)

    # Calculate number of top grid points
    n_top = mp.column_ztop / mp.column_dztop
    @assert mod(n_top, 1) == 0 "Top grid cell structure length does not go evenly into specified top structure depth."

    n_top = Int(n_top)

    if mp.column_dztop < 0.05 - D_TOLERANCE
        @warn "Initial top grid cell length (column_dztop) is < 0.05 m."
    end

    # Initialize top grid (constant spacing)
    dzT = fill(mp.column_dztop, n_top)

    # Build bottom grid (geometrically stretched)
    dzB = Float64[]
    gp0 = mp.column_dztop
    z0 = mp.column_ztop

    while mp.column_depth > (z0 + D_TOLERANCE)
        dz_new = gp0 * mp.column_zy
        push!(dzB, dz_new)
        gp0 = dz_new
        z0 += gp0
    end

    # Combine top and bottom
    return vcat(dzT, dzB)
end
