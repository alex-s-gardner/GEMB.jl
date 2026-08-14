"""
    initialize_profile(mp::ModelParameters, cf::ClimateForcing;
                       constant_density=false, constant_temperature=false)
        -> profile::DimStack

Initialize a GEMB firn column profile as a DimStack, as a steady-state guess
derived entirely from the climate forcing.

The grid is sized to the depth this climate needs (see
[`_derive_column_depth`](@ref)) rather than to the configured `mp.column_depth_max`,
which acts as a ceiling. Nothing else has to be told about that: the derived depth
lives in the returned `dz`, and both [`gemb`](@ref) and [`gemb_spinup`](@ref) take
the column depth they hold fixed from `sum(profile[:dz])`, not from `mp`.

# The scheme

[`initialize_climate_summary`](@ref) reduces the forcing to a handful of scalars
in one pass: snowfall and rainfall (partitioned on
`mp.rain_temperature_threshold`, so rain is not counted as accumulating mass), a
surface-energy-balance melt estimate, the cold-content-capped refreeze, the fitted
annual temperature harmonic, and the net annual mass balance
`b = snowfall + refreeze − melt`.

[`steady_state_profile`](@ref) then marches one parcel of snow forward in age with
`b` as its burial rate, recording every state variable as it is buried: density
(relaxing toward `mp.density_ice` under the run's own `mp.densification_method`),
temperature (a damped annual wave about the latent-heat-warmed mean), grain size
(evolved by [`calculate_grain_size`](@ref) itself) and irreducible water. Surface
albedo comes from [`calculate_albedo`](@ref) applied to the resulting surface
state.

**There is no regime threshold.** The sign of `b` is the only discriminator, and
it needs none: where `b ≤ 0` nothing is ever buried, the march terminates
immediately, and the column is the ice it exposes. Sites near `b = 0` get
intermediate profiles rather than a cliff between a deep firn column and a block
of solid ice.

# Keyword arguments

Both flags exist **only to reproduce the MATLAB `model_initialize_profile`
initialization** for the fidelity and regression tests. They are not physically
meaningful choices for a real run — the climate-derived guess is strictly better —
and either one builds the grid on the configured `mp.column_depth_max`, bypassing the
derived depth.

- `constant_density`: fill density with `mp.density_ice`, with the matching
  bare-ice grain state (`grain_dendricity = 0`, `grain_sphericity = 0`,
  `grain_radius = 2.5` mm) and no pore water.
- `constant_temperature`: fill temperature with `cf.temperature_air_mean`,
  unclamped.

Setting **both** takes a separate early-return path
([`_uniform_ice_profile`](@ref)) that reproduces the MATLAB column exactly: pure
ice, uniform mean-annual temperature, ice albedo, no march and no climate summary.

Returns a DimStack with Z dimension containing:
- dz, temperature, density, water, grain_radius,
  grain_dendricity, grain_sphericity, albedo, albedo_diffuse
"""
function initialize_profile(mp::ModelParameters, cf::ClimateForcing;
    constant_density::Bool=false, constant_temperature::Bool=false)

    T_mean = Float64(cf.temperature_air_mean)

    @assert T_mean > 0 "temperature_air_mean must exceed 0 K."
    if T_mean < 100
        @warn "temperature_air_mean should be in kelvin, but is below 100, suggesting an error."
    end

    # MATLAB-fidelity path: nothing here is climate-derived, so short-circuit
    # before the summary and the march so neither can perturb it.
    if constant_density && constant_temperature
        return _uniform_ice_profile(mp, T_mean)
    end

    cs = initialize_climate_summary(cf, mp)

    # Size the grid to the depth this climate needs, unless a fidelity flag pins it
    # to the configured value. The depth is not returned: it is realized in `dz`,
    # which is what fixes the column depth for the run.
    depth = (constant_density || constant_temperature) ? mp.column_depth_max :
            _derive_column_depth(mp, cs)

    dz = initialize_grid(mp, depth)
    m = length(dz)

    ss = steady_state_profile(dz, cs, mp)

    density = constant_density ? fill(mp.density_ice, m) : ss.density
    temperature = constant_temperature ? fill(T_mean, m) : ss.temperature
    water = constant_density ? zeros(m) : ss.water
    grain_radius = constant_density ? fill(GRAIN_RADIUS_ICE, m) : ss.grain_radius
    grain_dendricity = constant_density ? zeros(m) : ss.grain_dendricity
    grain_sphericity = constant_density ? zeros(m) : ss.grain_sphericity

    # Surface albedo from the initialized surface state rather than a binary
    # snow/ice pick: `calculate_albedo` already transitions continuously with
    # density and grain size, so a bare-ice column lands near `mp.albedo_ice` and a
    # snow column near `mp.albedo_snow` without being told which it is.
    albedo = _initial_albedo(temperature, dz, density, water, grain_radius, cs, mp)

    # Create Z dimension
    zdim = Z(1:m; metadata=cf_layer_index_attributes())

    # Note: `z_center` is intentionally not stored — it is a pure function of `dz`
    # (via `dz2z`) and is recomputed on demand, so the profile stays a clean slice
    # of the `gemb` output layout (which likewise carries `dz`, not `z_center`).
    layers = (
        dz=DimArray(dz, (zdim,)),
        temperature=DimArray(temperature, (zdim,)),
        density=DimArray(density, (zdim,)),
        water=DimArray(water, (zdim,)),
        grain_radius=DimArray(grain_radius, (zdim,)),
        grain_dendricity=DimArray(grain_dendricity, (zdim,)),
        grain_sphericity=DimArray(grain_sphericity, (zdim,)),
        albedo=DimArray(fill(albedo, m), (zdim,)),
        albedo_diffuse=DimArray(fill(albedo, m), (zdim,)),
    )

    # Same CF attributes a profile extracted from `gemb` output carries (`gemb_profile`),
    # so a profile describes itself the same way whichever end of the pipeline it came
    # from. No `cell_methods`: there is no time dimension here.
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end

"""
    _uniform_ice_profile(mp::ModelParameters, T_mean) -> DimStack

The MATLAB `model_initialize_profile` column: pure ice at a uniform temperature,
with non-dendritic faceted ice grains (`grain_radius = 2.5` mm, the non-spherical
cap in [`calculate_grain_size`](@ref)), ice albedo and no pore water.

Reached only when both fidelity flags of [`initialize_profile`](@ref) are set.
Kept as a separate function taking no climate-derived input so no future change to
the steady-state scheme can perturb it.
"""
function _uniform_ice_profile(mp::ModelParameters, T_mean::Real)
    dz = initialize_grid(mp)
    m = length(dz)
    zdim = Z(1:m; metadata=cf_layer_index_attributes())
    layers = (
        dz=DimArray(dz, (zdim,)),
        temperature=DimArray(fill(Float64(T_mean), m), (zdim,)),
        density=DimArray(fill(mp.density_ice, m), (zdim,)),
        water=DimArray(zeros(m), (zdim,)),
        grain_radius=DimArray(fill(GRAIN_RADIUS_ICE, m), (zdim,)),
        grain_dendricity=DimArray(zeros(m), (zdim,)),
        grain_sphericity=DimArray(zeros(m), (zdim,)),
        albedo=DimArray(fill(mp.albedo_ice, m), (zdim,)),
        albedo_diffuse=DimArray(fill(mp.albedo_ice, m), (zdim,)),
    )
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end

"""
    _derive_column_depth(mp::ModelParameters, cs::ClimateSummary;
                         margin=1.2, thermal_depths=4)
        -> depth [m]

Total column depth this climate needs, replacing the old binary 250-vs-25 m switch
with a continuous derivation. Used to build the grid in
[`initialize_profile`](@ref); it goes no further than that, since the grid is what
carries the depth forward.

The configured `mp.column_depth_max` is treated as a **ceiling** rather than a target,
so a shallow-firn or ablation site is not forced to carry hundreds of metres of
solid ice, while a deep cold site is unaffected. Two requirements set the depth:

1. **Resolve the firn.** Below the depth where the steady-state profile reaches
   `mp.density_ice` the column carries no further firn information, so
   `margin × ss.ice_depth` is deep enough — the margin leaves room for the real
   column to densify more slowly than the guess.
2. **Resolve the annual thermal wave.** Even a column of pure ice must be deep
   enough that its lower boundary does not clamp the seasonal temperature cycle.
   `thermal_depths` damping depths (`d = sqrt(2K/(ρ·c_p·ω))` at ice density, ≈3.3 m,
   so ≈13 m) is the floor. This is what keeps an ablation site — where
   `ice_depth = 0` because nothing is ever buried — from collapsing to a
   physically useless grid, and it recovers the old 25 m intent from the thermal
   physics rather than a constant.

The result is `clamp(max(firn requirement, thermal floor), mp.column_ztop,
mp.column_depth_max)`. Depths are found by marching on a coarse 1 m probe grid, which
is all the resolution two scalars need and keeps this independent of the final grid
it is used to build.
"""
function _derive_column_depth(mp::ModelParameters, cs::ClimateSummary;
    margin::Real=1.2, thermal_depths::Real=4)

    ceiling = mp.column_depth_max

    # Coarse uniform probe grid spanning the ceiling.
    probe_dz = fill(1.0, max(ceil(Int, ceiling), 1))
    ss = steady_state_profile(probe_dz, cs, mp)

    # Annual thermal-wave damping depth at ice density.
    d_thermal = thermal_damping_depth(cs.temperature_air_mean, mp.density_ice, mp)

    depth = max(margin * ss.ice_depth, thermal_depths * d_thermal)
    return Float64(clamp(depth, mp.column_ztop, ceiling))
end

"""
    _initial_albedo(temperature, dz, density, water, grain_radius, cs, mp) -> α

Surface albedo of the initialized column, from [`calculate_albedo`](@ref) applied
to the initialized surface state.

This replaces a binary `mp.albedo_snow` / `mp.albedo_ice` pick with the same
function the model uses at runtime, so the first timestep does not begin from an
albedo the physics immediately contradicts.

The seed is [`_snow_cover_albedo`](@ref) — the same continuous snow/ice blend the
melt estimate was solved against, so it is consistent with the climate summary and
carries no threshold. This matters for the methods that *decay from* the seed
rather than recomputing it (`:Bougamont2005`): a seed picked on a density
threshold would jump by `albedo_snow − albedo_ice` as a site crossed pore
close-off, which is exactly the cliff this scheme removes. For `:GardnerSharp`
(the default) `calculate_albedo` overwrites element `[1]` outright, so the seed is
immaterial there.
"""
function _initial_albedo(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64}, grain_radius::Vector{Float64},
    cs::ClimateSummary, mp::ModelParameters)

    seed = _snow_cover_albedo(cs.accumulation_effective, cs.melt, mp)
    albedo = fill(seed, length(density))
    albedo_diffuse = fill(seed, length(density))

    cfs = _albedo_forcing_step(cs, mp)
    albedo, _ = calculate_albedo(temperature, dz, density, water, grain_radius,
        albedo, albedo_diffuse, 0.0, 0.0, cfs, mp)

    return albedo[1]
end

"""
    _albedo_forcing_step(cs::ClimateSummary, mp::ModelParameters) -> ClimateForcingStep

Minimal [`ClimateForcingStep`](@ref) for the initial [`calculate_albedo`](@ref)
call. The optical properties come from the `mp` defaults (the same values the
forcing's `Fill`-backed layers carry when the forcing does not supply them), and
`precipitation` is zero: there is no fresh snow event at initialization.
"""
function _albedo_forcing_step(cs::ClimateSummary, mp::ModelParameters)
    return ClimateForcingStep(; dt=SECONDS_PER_YEAR,
        temperature_air=cs.temperature_air_mean,
        pressure_air=cs.pressure_air_mean,
        wind_speed=cs.wind_speed_mean,
        temperature_air_mean=cs.temperature_air_mean,
        wind_speed_mean=cs.wind_speed_mean,
        precipitation_mean=cs.accumulation_effective,
        temperature_observation_height=cs.temperature_observation_height,
        wind_observation_height=cs.wind_observation_height,
        black_carbon_snow=mp.black_carbon_snow,
        black_carbon_ice=mp.black_carbon_ice,
        cloud_optical_thickness=mp.cloud_optical_thickness,
        solar_zenith_angle=mp.solar_zenith_angle,
        shortwave_downward_diffuse=mp.shortwave_downward_diffuse,
        cloud_fraction=mp.cloud_fraction)
end

"""
    initialize_grid(mp::ModelParameters, depth=mp.column_depth_max)

Generate the initial vertical grid layer thicknesses, spanning `depth`.
Matches MATLAB's `model_initialize_grid` (local function in model_initialize_profile.m).

`depth` is an argument rather than read from `mp` because [`initialize_profile`](@ref)
builds the grid on the climate-derived depth ([`_derive_column_depth`](@ref)), for
which `mp.column_depth_max` is only the ceiling. It is the depth the grid actually
spans, hence `depth` and not `column_depth_max`.

Returns a Vector{Float64} of layer thicknesses from surface to depth.
"""
initialize_grid(mp::ModelParameters) = initialize_grid(mp, mp.column_depth_max)

function initialize_grid(mp::ModelParameters, depth::Real)

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

    while depth > (z0 + D_TOLERANCE)
        dz_new = gp0 * mp.column_zy
        push!(dzB, dz_new)
        gp0 = dz_new
        z0 += gp0
    end

    # Combine top and bottom
    return vcat(dzT, dzB)
end
