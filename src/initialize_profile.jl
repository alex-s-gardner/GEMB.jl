"""
    initialize_profile(mp::ModelParameters, cf::ClimateForcing;
                       constant_density=false, constant_temperature=false,
                       climate_summary=nothing)
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
temperature (a damped annual wave about the mean *surface* temperature, warmed near
the surface by refreezing latent heat), grain size
(evolved by [`calculate_grain_size`](@ref) itself) and irreducible water. No albedo
is stored: [`calculate_albedo`](@ref) diagnoses it from the column at the start of
every timestep, including the first.

**There is no regime threshold.** The sign of `b` is the only discriminator, and
it needs none: where `b ≤ 0` nothing is ever buried, the march terminates
immediately, and the column is the ice it exposes. Sites near `b = 0` get
intermediate profiles rather than a cliff between a deep firn column and a block
of solid ice.

# Keyword arguments

Both flags exist to build a **simple, fully specified starting column** — the form the
tests use, where a climate-derived initial guess would make an expected value depend on the
whole initialization chain. They are not physically meaningful choices for a real run (the
climate-derived guess is strictly better), and either one builds the grid on the configured
`mp.column_depth_max`, bypassing the derived depth.

- `constant_density`: fill density with `mp.density_ice`, with the matching
  bare-ice grain state (`grain_dendricity = 0`, `grain_sphericity = 0`,
  `grain_radius = 2.5` mm) and no pore water.
- `constant_temperature`: fill temperature with `cf.temperature_air_mean`, clamped
  to the melt point as every other path is.

Setting **both** takes a separate early-return path
([`_uniform_ice_profile`](@ref)): pure ice at a uniform mean-annual temperature, with no
march and no climate summary.

- `climate_summary`: a [`ClimateSummary`](@ref) for `cf`, when the caller already has one.
  Purely an optimization — the summary is a pure function of `(cf, mp)`, so passing one cannot
  change the result, and `nothing` (the default) builds it here as before. Exists because
  `forcing_climatology(method=:representative)` needs the record's summary for its accumulation
  guard and would otherwise pay for an identical second one: measured at 346 ms of a 3020 ms call
  on a 32-year 3-hourly record.

Every path here returns temperature at or below the melt point (273.15 K): the
climate-derived march clamps in [`_steady_state_temperature`](@ref), and the
two flag paths clamp `cf.temperature_air_mean` (with a warning) before filling.

Returns a DimStack with Z dimension containing:
- dz, temperature, density, water, grain_radius,
  grain_dendricity, grain_sphericity, age

`age` [d] follows `mp.initialize_age`. Under `:steady_state` (the default) it is the
residence time the march itself integrated — the accumulated burial time at each depth,
increasing downward — which makes age and residence-time diagnostics meaningful from the
first timestep instead of only after a multi-century spinup has flushed the column. Under
`:zero` the instant of initialization is the epoch and every cell starts at 0, which is
what the flag paths do regardless (`constant_density` discards the march's density, so
its age describes a column that was replaced). Nothing in the physics reads `age`, so this
choice is output-only.
"""
function initialize_profile(mp::ModelParameters, cf::ClimateForcing;
    constant_density::Bool=false, constant_temperature::Bool=false,
    climate_summary=nothing)

    T_mean = Float64(cf.temperature_air_mean)

    @assert T_mean > 0 "temperature_air_mean must exceed 0 K."
    if T_mean < 100
        @warn "temperature_air_mean should be in kelvin, but is below 100, suggesting an error."
    end

    # No initialized cell may start above the melt point: the column carries no
    # enthalpy above `CtoK` (that energy is melt, and there is no water here to hold
    # it), so a warmer cell would be an unphysical heat reservoir that the first
    # timesteps discharge downward. The climate-derived path is already clamped in
    # `_steady_state_temperature`; this clamp covers the two flag paths, which fill
    # `T_mean` verbatim.
    #
    # Applied only where `T_mean` is actually used — under `constant_temperature` — so a
    # temperate site does not get told its temperature was clamped on a path that never
    # reads `T_mean` at all.
    if constant_temperature && T_mean > CtoK
        @warn "temperature_air_mean is above the melt point; initialized temperature clamped to 0 °C." T_mean
        T_mean = CtoK
    end

    # Fully specified column: nothing here is climate-derived, so short-circuit before the
    # summary and the march so neither can perturb it.
    if constant_density && constant_temperature
        return _uniform_ice_profile(mp, T_mean)
    end

    # Reuse a caller's summary when given one. `forcing_climatology(method=:representative)`
    # already needs the record's summary for its accumulation guard, and building it is not cheap
    # (346 ms on a 32-year 3-hourly record), so recomputing an identical one here was ~11% of that
    # call. `nothing` — every other caller — behaves exactly as before.
    cs = climate_summary === nothing ? initialize_climate_summary(cf, mp) : climate_summary

    # Size the grid to the depth this climate needs, unless a flag pins it
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
    # `mp.initialize_age` picks the epoch convention; `constant_density` forces `:zero`
    # because that path discards the march's density, and the marched age describes the
    # firn column the march built, not the block of ice this flag substitutes for it.
    age = (constant_density || mp.initialize_age === :zero) ? zeros(m) : ss.age

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
        age=DimArray(age, (zdim,)),
    )

    # Same CF attributes a profile extracted from `gemb` output carries (`gemb_profile`),
    # so a profile describes itself the same way whichever end of the pipeline it came
    # from. No `cell_methods`: there is no time dimension here.
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end

"""
    _uniform_ice_profile(mp::ModelParameters, T_mean) -> DimStack

A pure-ice column at a uniform temperature, with non-dendritic faceted ice grains
(`grain_radius = 2.5` mm, the non-spherical cap in [`calculate_grain_size`](@ref)) and no
pore water.

Reached only when both keyword flags of [`initialize_profile`](@ref) are set.
Kept as a separate function taking no climate-derived input so no future change to
the steady-state scheme can perturb it.

`T_mean` is expected already clamped to `CtoK` by the caller — the clamp lives
there so it covers both flag paths in one place.
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
        age=DimArray(zeros(m), (zdim,)),
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
    initialize_grid(mp::ModelParameters, depth=mp.column_depth_max)

Generate the initial vertical grid layer thicknesses, spanning `depth`.

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
