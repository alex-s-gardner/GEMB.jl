"""
    fresh_snow_density(mp::ModelParameters, T_air_mean, precip_mean, wind_speed_mean,
                       T_air=T_air_mean, wind_speed=wind_speed_mean)

Return the density of fresh snow [kg m-3] for the configured `mp.new_snow_method`.

The first three forcing arguments are climatological means (temperature [K], accumulation
[kg m-2 yr-1], wind speed [m s-1]), so this is shared by `calculate_accumulation` (per
timestep) and `initialize_profile` (as the surface density ρ₀ of the steady-state firn
profile).

`T_air` and `wind_speed` are the *instantaneous* air temperature [K] and wind speed
[m s-1], read only by the methods whose fits are meant to apply at the moment of snowfall
rather than to a climatology: `T_air` by `:FaustoFit` and `:Pahaut`, `wind_speed` by
`:Pahaut`. They default to their climatological counterparts so that the steady-state
march — which has no instantaneous forcing to offer — gets a well-defined climatological
ρ₀ from the same call, and so that callers of the older methods need not pass them.

Methods:
- `:Constant150` → 150, `:Constant315` → 315, `:Constant350` → 350: constants.
- `:Fausto` → 315: the Fausto et al. (2018) Greenland fit frozen at one temperature (it is
  `:FaustoFit` evaluated at T ≈ 256.2 K). Inherited from MATLAB and kept as the default for
  reference fidelity. Numerically equal to `:Constant315`, but `:Fausto` additionally selects
  the Crocus wind-dependent fresh-grain properties below, where `:Constant315` does not — use
  `:Constant315` for a bare 315 kg m-3 with the default grain properties.
- `:FaustoFit` → `362.1 + 2.78·(T_air - CtoK)`: that same fit carrying its temperature
  dependence, as implemented in IMAU-FDM (`initialise_model.f90`) for its Greenland domain.
  Unbounded below as published, so callers clamp — see `calculate_accumulation` and
  `steady_state_profile`.
- `:Pahaut` → `max(50, 109 + 6·(T_air - CtoK) + 26·√wind_speed)`: Pahaut (1975) as
  implemented in Crocus, via Lafaysse et al. (2026) eq. 35 (SURFEX/Crocus v3.0.2, their
  `:V12` fresh-snow option). The only method here fitted to *alpine seasonal snow* rather
  than to a polar ice sheet, and the only one carrying an explicit wind dependence at the
  instantaneous timestep. Both dependencies are physical: warmer snowfall gives denser
  crystals, and wind fragments them before and during deposition. It is much lighter than
  the polar fits over their common range — at T = -10 °C and 5 m s-1 it gives 107 kg m-3
  against `:FaustoFit`'s 334 — because alpine snowfall is warmer, wetter, and much less
  wind-packed than the katabatic-scoured surfaces the Greenland and Antarctic fits were
  regressed on. Prefer it for temperate and mid-latitude glaciers; the polar fits will
  overestimate fresh-snow density there, which suppresses the albedo of new snow and speeds
  its burial. The published floor of 50 kg m-3 is retained (it binds below about -10 °C in
  calm air).
- `:Kaspers` → `(7.36e-2 + 1.06e-3·T + 6.69e-2·P/1000 + 4.77e-3·U)·1000`: the Antarctic
  fit of Kaspers et al. (2004), the only method here depending on all three climatological
  means. Temperature is capped at the melting point, as published.
- `:KuipersMunneke` → `481.0 + 4.834·(T - CtoK)`: the Antarctic fit used by Kuipers Munneke
  et al. (2015). Both take *climatological* means rather than instantaneous forcing.

# References
- Pahaut, E. (1975). *La métamorphose des cristaux de neige (Snow crystal metamorphosis)*.
  Monographies de la Météorologie Nationale 96, Météo-France.
- Lafaysse, M., et al. (2026). Version 3.0.2 of the Crocus snowpack model.
  *Geosci. Model Dev.* 19, 6273-6334. https://doi.org/10.5194/gmd-19-6273-2026
- Fausto, R. S., et al. (2018). A snow density dataset for improving surface boundary
  conditions in Greenland ice sheet firn modeling. *Front. Earth Sci.* 6, 51.
- Vionnet, V., et al. (2012). The detailed snowpack scheme Crocus and its implementation in
  SURFEX v7.2. *Geosci. Model Dev.* 5, 773-791.
- Kaspers, K. A., et al. (2004). Model calculations of the age of firn air across the
  Antarctic continent. *Atmos. Chem. Phys.* 4, 1365-1380. (`:Kaspers`)
- Kuipers Munneke, P., et al. (2015). Elevation change of the Greenland Ice Sheet due to
  surface mass balance and firn processes, 1960-2014. *The Cryosphere* 9, 2009-2025.
  (`:KuipersMunneke`)
"""
function fresh_snow_density(mp::ModelParameters, T_air_mean::Real,
    precip_mean::Real, wind_speed_mean::Real, T_air::Real=T_air_mean,
    wind_speed::Real=wind_speed_mean)
    if mp.new_snow_method == :Constant150
        return 150.0
    elseif mp.new_snow_method == :Constant315
        return 315.0
    elseif mp.new_snow_method == :Constant350
        return 350.0
    elseif mp.new_snow_method == :Fausto
        return DENSITY_NEW_SNOW_FAUSTO_CONSTANT
    elseif mp.new_snow_method == :FaustoFit
        return DENSITY_NEW_SNOW_FAUSTO_A + DENSITY_NEW_SNOW_FAUSTO_B * (T_air - CtoK)
    elseif mp.new_snow_method == :Pahaut
        return max(DENSITY_NEW_SNOW_PAHAUT_MIN,
            DENSITY_NEW_SNOW_PAHAUT_A + DENSITY_NEW_SNOW_PAHAUT_B * (T_air - CtoK) +
            DENSITY_NEW_SNOW_PAHAUT_C * sqrt(max(wind_speed, 0.0)))
    elseif mp.new_snow_method == :Kaspers
        return (7.36e-2 + 1.06e-3 * min(T_air_mean, CtoK - T_TOLERANCE) +
                6.69e-2 * precip_mean / 1000.0 + 4.77e-3 * wind_speed_mean) * 1000.0
    elseif mp.new_snow_method == :KuipersMunneke
        return 481.0 + 4.834 * (T_air_mean - CtoK)
    end
    return 0.0
end

"""
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Add precipitation and deposition to the model column.

Precipitation is classified as snow or rain based on `mp.rain_temperature_threshold`.
Snow is added as a new layer (if depth > dzmin) or merged into the top cell.
Rain is added by increasing the mass and temperature of the top grid cell,
with temperature adjusted to account for latent heat of fusion.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
age, rain)`. Arrays grow by one cell (via [`open_slot!`](@ref)) when snow depth > dzmin; the
extra cell is reclaimed later in the timestep by [`manage_layer_thickness`](@ref)'s count
controller.

This is the model's principal age source. All three paths add mass at age 0: a fresh cell is
set to 0 outright, and the two merge-into-cell-1 paths dilute `age[1]` by the arriving mass
via [`dilute_age`](@ref).

Albedo is absent: it is diagnosed from the column at the top of the *next* timestep (see
[`calculate_albedo`](@ref)), which reads the fresh-snow grain size and density this
function sets, so new snow brightens the surface without an albedo being stored here.
"""
function calculate_accumulation(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, age::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

    # Note: arrays are modified in place, and grow by one cell via `open_slot!` when new
    # snow is deep enough to warrant its own layer.

    # Specify constants (shared with the steady-state initial guess, so a freshly
    # initialized column and freshly fallen snow agree by construction)
    re_new_snow = RE_NEW_SNOW    # new snow grain size [mm]
    gdn_new_snow = GDN_NEW_SNOW  # new snow dendricity
    gsp_new_snow = GSP_NEW_SNOW  # new snow sphericity
    rain = 0.0               # rainfall [mm w.e. or kg m^-2]

    if verbose
        M = dz .* density
        M_total_initial = sum(M)
        E_total_initial = column_enthalpy(mp, M, temperature, water)
    end

    # Density of fresh snow [kg m-3] (shared with initialize_profile). The instantaneous air
    # temperature and wind speed are passed as the fifth and sixth arguments, for `:FaustoFit`
    # and `:Pahaut`; the other methods ignore them. Clamped for the same reason
    # `steady_state_profile` clamps its own call: the Fausto fit is a linear regression with no
    # bounds, so a sufficiently cold forcing step would otherwise hand a non-positive density
    # to the `dz = mass/density` below. (`:Pahaut` carries its own published floor.)
    density_new_snow = clamp(fresh_snow_density(mp, cfs.temperature_air_mean,
            cfs.precipitation_mean, cfs.wind_speed_mean, cfs.temperature_air,
            cfs.wind_speed),
        1.0, mp.density_ice)
    if mp.new_snow_method == :Fausto || mp.new_snow_method == :FaustoFit ||
       mp.new_snow_method == :Pahaut
        # From Vionnet et al., 2012 (Crocus): wind-dependent grain properties. `:Pahaut` takes
        # these too, being the fresh-snow density Crocus itself pairs them with.
        gdn_new_snow = min(max(1.29 - 0.17 * cfs.wind_speed, 0.20), 1.0)
        gsp_new_snow = min(max(0.08 * cfs.wind_speed + 0.38, 0.5), 0.9)
        re_new_snow = max(1e-1 * (gdn_new_snow / 0.99 + (1.0 - 1.0 * gdn_new_snow / 0.99) * (gsp_new_snow / 0.99 * 3.0 + (1.0 - gsp_new_snow / 0.99) * 4.0)) / 2.0, GDN_TOLERANCE)
    end

    M_surface = dz[1] * density[1]
    # Total mass of the surface cell, the weight for the age dilutions below. Captured here,
    # before any branch mutates `dz[1]` or `density[1]`.
    S_surface = M_surface + water[1]

    if cfs.precipitation > 0
        # if snow
        if cfs.temperature_air <= (mp.rain_temperature_threshold + T_TOLERANCE)
            z_snow = cfs.precipitation / density_new_snow          # depth of snow
            dfall = gdn_new_snow
            sfall = gsp_new_snow
            refall = re_new_snow

            # if snow depth is greater than specified min dz, new cell created
            if z_snow > mp.column_dzmin + D_TOLERANCE
                cols = column_state(temperature, dz, density, water, grain_radius,
                    grain_dendricity, grain_sphericity, age)
                open_slot!(cols, 1)
                @inbounds begin
                    temperature[1] = cfs.temperature_air
                    dz[1] = z_snow
                    density[1] = density_new_snow
                    water[1] = 0.0
                    grain_radius[1] = refall
                    grain_dendricity[1] = dfall
                    grain_sphericity[1] = sfall
                    # `open_slot!` duplicated old cell 1 here; fresh snow must not inherit
                    # its age.
                    age[1] = 0.0
                end
            else
                # if snow depth is less than specified minimum dz
                M_surface_new = M_surface + cfs.precipitation
                dz[1] = dz[1] + cfs.precipitation / density_new_snow
                density[1] = M_surface_new / dz[1]

                # adjust temperature (assume precipitation is same temp as air)
                temperature[1] = mix_temperature(mp, cfs.temperature_air, cfs.precipitation,
                    temperature[1], M_surface)

                age[1] = dilute_age(age[1], S_surface, cfs.precipitation)

                grain_dendricity[1] = dfall
                grain_sphericity[1] = sfall
                grain_radius[1] = max(0.1 * (grain_dendricity[1] / 0.99 + (1.0 - 1.0 * grain_dendricity[1] / 0.99) * (grain_sphericity[1] / 0.99 * 3.0 + (1.0 - grain_sphericity[1] / 0.99) * 4.0)) / 2, GDN_TOLERANCE)
            end

        else
            # rain
            # grid cell adjusted mass
            M_surface_new = M_surface + cfs.precipitation

            # adjust temperature (liquid: must account for latent heat of fusion)
            temperature[1] = mix_temperature_liquid(mp, cfs.temperature_air, cfs.precipitation,
                temperature[1], M_surface)

            # Rain mass joins the ice/firn matrix here (not `water`), so it is new age-zero
            # mass in the cell on the same footing as snowfall.
            age[1] = dilute_age(age[1], S_surface, cfs.precipitation)

            # adjust grid cell density
            density[1] = M_surface_new / dz[1]

            # if density > the density of ice
            if density[1] > mp.density_ice - D_TOLERANCE
                density[1] = mp.density_ice
                dz[1] = M_surface_new / density[1]
            end

            rain = cfs.precipitation
        end

        if verbose
            # Check for conservation of mass
            M = dz .* density
            M_total_final = sum(M)
            M_delta = M_total_final - M_total_initial - cfs.precipitation

            E_total_final = column_enthalpy(mp, M, temperature, water)

            E_snow = (cfs.precipitation - rain) * specific_enthalpy(mp, cfs.temperature_air)
            E_rain = rain * specific_enthalpy_water(mp, cfs.temperature_air)

            E_delta = E_total_final - E_total_initial - E_snow - E_rain

            E_tol = energy_tolerance(E_total_initial)

            if (abs(M_delta) > M_TOLERANCE) || (abs(E_delta) > E_tol)
                error("Mass and/or energy are not conserved:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
            end
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
    age, rain
end
