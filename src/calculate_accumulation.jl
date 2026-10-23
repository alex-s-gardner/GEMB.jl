"""
    fresh_snow_density(mp::ModelParameters, T_air_mean, precip_mean, wind_speed_mean,
                       T_air=T_air_mean, wind_speed=wind_speed_mean;
                       wind_observation_height=10.0, z0=Z0_SNOW_DRY)

Return the density of fresh snow [kg m-3] for the configured `mp.new_snow_method`.

The first three forcing arguments are climatological means (temperature [K], accumulation
[kg m-2 yr-1], wind speed [m s-1]), so this is shared by `calculate_accumulation` (per
timestep) and `initialize_profile` (as the surface density ρ₀ of the steady-state firn
profile).

`T_air` and `wind_speed` are the *instantaneous* air temperature [K] and wind speed
[m s-1], read only by the methods whose fits are meant to apply at the moment of snowfall
rather than to a climatology: `T_air` by `:FaustoFit`, `:Pahaut`, `:vanKampenhout2017` and
`:vanKampenhout2017SNOTEL`, `wind_speed` by `:Pahaut` and the two `:vanKampenhout2017`
methods. They default to their climatological counterparts so that the steady-state march — which has no instantaneous forcing to offer —
gets a well-defined climatological ρ₀ from the same call, and so that callers of the older
methods need not pass them.

Methods:
- `:Constant150` → 150, `:Constant315` → 315, `:Constant350` → 350: constants.
- `:Constant` → `mp.new_snow_density`: a constant of the caller's choosing, validated to
  `[NEW_SNOW_DENSITY_MIN, NEW_SNOW_DENSITY_MAX]`, with the same fresh-grain properties as the
  three fixed constants.
- `:Fausto` → 315: the Fausto et al. (2018) Greenland fit frozen at one temperature (it is
  `:FaustoFit` evaluated at T ≈ 256.2 K). Numerically equal to `:Constant315`, but `:Fausto`
  additionally selects the Crocus wind-dependent fresh-grain properties below, where
  `:Constant315` does not — use `:Constant315` for a bare 315 kg m-3 with the default grain
  properties.
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
- `:vanKampenhout2017` → `ρ_T(T_air) + ρ_W(U₁₀)`: van Kampenhout et al. (2017) eqs. 2-4,
  p. 2587, the CLM5 deposition density; see [`vankampenhout_density_temperature`](@ref) and
  [`vankampenhout_density_wind`](@ref). It reads the instantaneous `T_air` and the 10 m wind
  speed, which is `wind_speed` carried from `wind_observation_height` to 10 m by
  [`wind_speed_10m`](@ref) over roughness `z0`. It keeps GEMB's default fresh-grain
  properties. It is a deposition density (50-436 kg m-3), so its near-surface density
  depends on the settling scheme: pair it with `densification_method = :Crocus`
  (destructive metamorphism is on by default there). Under an accumulation-driven
  densification scheme the snow stays near its deposition density for years.
- `:vanKampenhout2017SNOTEL` → `ρ_T,SNOTEL(T_air) + ρ_W(U₁₀)`: `:vanKampenhout2017` with the
  SNOTEL-calibrated temperature term [`snotel_density_temperature`](@ref); the wind term and
  fresh-grain properties are unchanged.

`wind_observation_height` [m] and `z0` [m] are read only by `:vanKampenhout2017` and
`:vanKampenhout2017SNOTEL`. Their defaults treat `wind_speed` as a 10 m wind over dry snow.

Throws an `ErrorException` for an unrecognized `mp.new_snow_method`.

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
- van Kampenhout, L., Lenaerts, J. T. M., Lipscomb, W. H., Sacks, W. J., Lawrence, D. M.,
  Slater, A. G., and van den Broeke, M. R. (2017). Improving the representation of polar snow
  and firn in the Community Earth System Model. *J. Adv. Model. Earth Syst.* 9, 2583-2600.
  (`:vanKampenhout2017`)
"""
function fresh_snow_density(mp::ModelParameters, T_air_mean::Real,
    precip_mean::Real, wind_speed_mean::Real, T_air::Real=T_air_mean,
    wind_speed::Real=wind_speed_mean;
    wind_observation_height::Real=WIND_HEIGHT_10M,
    z0::Real=Z0_SNOW_DRY)
    if mp.new_snow_method == :Constant150
        return 150.0
    elseif mp.new_snow_method == :Constant315
        return 315.0
    elseif mp.new_snow_method == :Constant350
        return 350.0
    elseif mp.new_snow_method == :Constant
        return mp.new_snow_density
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
    elseif mp.new_snow_method in (:vanKampenhout2017, :vanKampenhout2017SNOTEL)
        ρ_T = mp.new_snow_method === :vanKampenhout2017SNOTEL ?
              snotel_density_temperature(T_air) : vankampenhout_density_temperature(T_air)
        return ρ_T + vankampenhout_density_wind(
            wind_speed_10m(wind_speed, wind_observation_height, z0))
    end
    error("unrecognized new_snow_method: $(mp.new_snow_method)")
end

# Throws an `ArgumentError` unless `T_air` [K] is finite and positive.
function _air_temperature_check(T_air::Real)
    (isfinite(T_air) && T_air > 0) ||
        throw(ArgumentError("air temperature must be finite and positive [K], got $T_air"))
    return nothing
end

"""
    vankampenhout_density_temperature(T_air) -> ρ_T [kg m-3]

Temperature term of the van Kampenhout et al. (2017) fresh-snow density, eq. 3, p. 2587, with
`T_air` the near-surface air temperature [K] and `ΔT = T_air - 273.15`:

    ρ_T = 50 + 1.7·17^1.5                  ΔT > 2
          50 + 1.7·(ΔT + 15)^1.5           -15 < ΔT ≤ 2
          -3.8328·ΔT - 0.0333·ΔT²          ΔT ≤ -15

The middle branch is Anderson (1976). The warm branch holds it at its +2 °C value
(169.16 kg m-3). The cold branch is a fit to unpublished data by A. G. Slater, and gives a
density that rises again as the temperature falls (100.03 kg m-3 at -40 °C). The middle branch
is continuous at +2 °C. At -15 °C it gives 50 and the cold branch gives 49.9995, so ρ_T steps
down by 5e-4 kg m-3 there, as published.

Throws an `ArgumentError` unless `T_air` is finite and positive.
"""
function vankampenhout_density_temperature(T_air::Real)
    _air_temperature_check(T_air)
    ΔT = T_air - CtoK
    if ΔT > DENSITY_NEW_SNOW_VK17_DT_WARM
        return DENSITY_NEW_SNOW_VK17_A + DENSITY_NEW_SNOW_VK17_B *
               (DENSITY_NEW_SNOW_VK17_DT_WARM - DENSITY_NEW_SNOW_VK17_DT_COLD)^1.5
    elseif ΔT > DENSITY_NEW_SNOW_VK17_DT_COLD
        return DENSITY_NEW_SNOW_VK17_A + DENSITY_NEW_SNOW_VK17_B *
               (ΔT - DENSITY_NEW_SNOW_VK17_DT_COLD)^1.5
    else
        return DENSITY_NEW_SNOW_VK17_C1 * ΔT + DENSITY_NEW_SNOW_VK17_C2 * ΔT^2
    end
end

"""
    snotel_density_temperature(T_air) -> ρ_T [kg m-3]

Temperature term of `new_snow_method = :vanKampenhout2017SNOTEL`, with `T_air` the
near-surface air temperature [K] and `ΔT = T_air - 273.15`:

    ρ_T = 56 + 0.85·(clamp(ΔT, -8, 2) + 15)^1.5

The form of the van Kampenhout et al. (2017) middle branch
([`vankampenhout_density_temperature`](@ref)) with refitted coefficients and a floor in place
of its cold branch: ρ_T is 71.74 kg m-3 at and below -8 °C and is held at its +2 °C value
above. It is calibrated with `densification_method = :Crocus`, destructive metamorphism and
`destructive_metamorphism_density = 100`; see the Physics Notes. The fit has no data below
about -25 °C, where the published cold branch exceeds the floor.

Throws an `ArgumentError` unless `T_air` is finite and positive.
"""
function snotel_density_temperature(T_air::Real)
    _air_temperature_check(T_air)
    ΔT = clamp(T_air - CtoK, DENSITY_NEW_SNOW_SNOTEL_DT_FLOOR, DENSITY_NEW_SNOW_SNOTEL_DT_WARM)
    return DENSITY_NEW_SNOW_SNOTEL_A +
           DENSITY_NEW_SNOW_SNOTEL_B * (ΔT - DENSITY_NEW_SNOW_SNOTEL_DT_REF)^1.5
end

"""
    vankampenhout_density_wind(wind_speed_10m) -> ρ_W [kg m-3]

Wind term of the van Kampenhout et al. (2017) fresh-snow density, eq. 4, p. 2587, with `U` the
10 m wind speed [m s-1]:

    ρ_W = 266.861·(½(1 + tanh(U/5)))^8.8

It rises from 0.599 kg m-3 at U = 0 through 87.34 at 5 m s-1 toward 266.861, and is half its
limit at 6.25 m s-1. The paper describes it as a smooth transition at 5 m s-1, built on the
experiments of Liston et al. (2007).

Throws an `ArgumentError` unless `wind_speed_10m` is finite and non-negative.
"""
function vankampenhout_density_wind(wind_speed_10m::Real)
    (isfinite(wind_speed_10m) && wind_speed_10m >= 0) ||
        throw(ArgumentError("wind speed must be finite and non-negative [m s-1], got $wind_speed_10m"))
    return DENSITY_NEW_SNOW_VK17_W *
           (0.5 * (1.0 + tanh(wind_speed_10m / DENSITY_NEW_SNOW_VK17_U_SCALE)))^
           DENSITY_NEW_SNOW_VK17_W_EXPONENT
end

"""
    fresh_snow_properties(cfs, mp, z0) -> (density, grain_radius, grain_dendricity, grain_sphericity)

Density [kg m-3] and grain state of the snow falling in step `cfs`, onto a surface of roughness
`z0` [m] (read by the two `new_snow_method = :vanKampenhout2017` methods only):
[`fresh_snow_density`](@ref) clamped to `[1, mp.density_ice]`, and GEMB's fresh-grain constants, or under `:Fausto`,
`:FaustoFit` and `:Pahaut` the wind-dependent grains of Crocus (Vionnet et al., 2012). Shared
by [`calculate_accumulation`](@ref) and the snow-free ground of `lower_boundary = :ground`.
"""
function fresh_snow_properties(cfs::ClimateForcingStep, mp::ModelParameters, z0::Real)
    # Specify constants (shared with the steady-state initial guess, so a freshly
    # initialized column and freshly fallen snow agree by construction)
    re_new_snow = RE_NEW_SNOW    # new snow grain size [mm]
    gdn_new_snow = GDN_NEW_SNOW  # new snow dendricity
    gsp_new_snow = GSP_NEW_SNOW  # new snow sphericity

    # Density of fresh snow [kg m-3] (shared with initialize_profile). The instantaneous air
    # temperature and wind speed are passed as the fifth and sixth arguments, for `:FaustoFit`,
    # `:Pahaut` and the two `:vanKampenhout2017` methods; the other methods ignore them. The wind
    # height and the roughness of the surface the snow lands on are read by those two only.
    # Clamped for the same reason `steady_state_profile` clamps its own call: the Fausto fit
    # is a linear regression with no bounds, so a sufficiently cold forcing step would
    # otherwise hand a non-positive density to the `dz = mass/density` below. (`:Pahaut`
    # carries its own published floor.)
    density_new_snow = clamp(fresh_snow_density(mp, cfs.temperature_air_mean,
            cfs.precipitation_mean, cfs.wind_speed_mean, cfs.temperature_air,
            cfs.wind_speed; wind_observation_height=cfs.wind_observation_height,
            z0=z0),
        1.0, mp.density_ice)
    if mp.new_snow_method == :Fausto || mp.new_snow_method == :FaustoFit ||
       mp.new_snow_method == :Pahaut
        # From Vionnet et al., 2012 (Crocus): wind-dependent grain properties. `:Pahaut` takes
        # these too, being the fresh-snow density Crocus itself pairs them with.
        gdn_new_snow = min(max(1.29 - 0.17 * cfs.wind_speed, 0.20), 1.0)
        gsp_new_snow = min(max(0.08 * cfs.wind_speed + 0.38, 0.5), 0.9)
        re_new_snow = max(1e-1 * (gdn_new_snow / 0.99 + (1.0 - 1.0 * gdn_new_snow / 0.99) * (gsp_new_snow / 0.99 * 3.0 + (1.0 - gsp_new_snow / 0.99) * 4.0)) / 2.0, GDN_TOLERANCE)
    end

    return density_new_snow, re_new_snow, gdn_new_snow, gsp_new_snow
end

"""
    calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)

Add precipitation and deposition to the model column.

Precipitation is classified as snow or rain based on `mp.rain_temperature_threshold`.

Snowfall follows the layering of Crocus (Vionnet et al., 2012, p. 778; Lafaysse et al., 2026,
eq. 11, p. 6278): it joins the top cell only when that cell is similar snow, and otherwise
forms a layer of its own. Snow of mass `P` and depth `z = P/ρ_new` becomes a new top cell with
fresh-snow properties when `z > mp.column_dzmin`, or when `P` is at least
[`snow_layer_min_mass`](@ref)`(cfs.dt)` and it is distinct from the top cell
([`snow_layers_distinct`](@ref) at the [`snow_layer_distance`](@ref) `D(n, 1)` of the top cell
and the fresh snow). Over firn or ice the distance is infinite, so such
snow always forms its own cell. The new cell may be thinner than `mp.column_dzmin`;
[`manage_layer_thickness`](@ref) keeps it while it stays distinct from the cell below and holds
at least that mass.

Otherwise the snow merges into the top cell, of matrix mass `M = dz[1]·density[1]`, whether
that is snow, firn or ice: thickness grows by `P/ρ_new` and density becomes `(M + P)/dz[1]`
(mass conserved); temperature is the enthalpy-weighted mix of the cell and the snow at the air
temperature ([`mix_temperature`](@ref); energy conserved); age is diluted by the age-0 snow
([`dilute_age`](@ref)); and the grain state is mixed by mass ([`mix_grains`](@ref), weight
`f = P/(M + P)` on the fresh snow). Snowfall below the minimum mass onto firn or ice
therefore lowers the cell's density, so repeated light snowfall builds up snow there.

`P = 0` changes nothing. A trace snowfall merges, being below the minimum mass. It moves
the grain state by `O(f)`, and not at all where the mixed values round back to the cell's own
(5e-324 kg m-2 on any cell), and the density by `O(f)` (at most an ulp for 5e-324 kg m-2). The
surface albedo, a function of the top cell's grain radius and density and of the snow mass
above firn or ice, therefore moves by `O(f)` at most.

Rain is added to the pore water of the top cell, `water[1]`, at the melt-point water enthalpy,
and its sensible heat above the melting point, `P·(h_w(T_air) − h_w(0 °C))`, is added to the
cell's ice matrix. The cell's density and thickness are unchanged. In a cold cell the rain
refreezes in [`calculate_melt`](@ref), where it is counted as refreeze; in a cell at the
melting point its sensible heat melts ice there.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
age, rain)`. Arrays grow by one cell (via [`open_slot!`](@ref)) when the snow forms its own
cell; the extra cell is reclaimed later in the timestep by [`manage_layer_thickness`](@ref)'s
count controller.

This is the model's principal age source. Every path adds mass at age 0: a fresh cell is
set to 0 outright, and the merge-into-cell-1 paths dilute `age[1]` (the mean age of the
cell's matrix and pore water together) by the arriving mass via [`dilute_age`](@ref).

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

    rain = 0.0               # rainfall [mm w.e. or kg m^-2]

    if verbose
        M = dz .* density
        M_total_initial = sum(M) + sum(water)
        E_total_initial = column_enthalpy(mp, M, temperature, water)
    end

    density_new_snow, re_new_snow, gdn_new_snow, gsp_new_snow =
        fresh_snow_properties(cfs, mp, surface_roughness(mp, density[1], water[1]))

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

            # Fresh snow forms its own top cell when it is deeper than `column_dzmin`, or when
            # it holds at least `snow_layer_min_mass` and is distinct from the top cell
            # (always so over firn or ice); see the docstring. D(n, 1) takes the top cell as
            # the upper layer and the new snow as the lower.
            D = snow_layer_distance(grain_dendricity[1], grain_sphericity[1], grain_radius[1],
                water[1] > WATER_TOLERANCE, age[1], density[1],
                dfall, sfall, refall, false, 0.0, density_new_snow)
            if (z_snow > mp.column_dzmin + D_TOLERANCE) ||
               ((cfs.precipitation >= snow_layer_min_mass(cfs.dt)) &&
                snow_layers_distinct(z_snow, D, mp))
                cols = column_state(temperature, dz, density, water, grain_radius,
                    grain_dendricity, grain_sphericity, age)
                open_slot!(cols, 1)
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
            else
                # Merged into the top cell (adjust temperature: precipitation is at the air
                # temperature).
                temperature[1] = mix_temperature(mp, cfs.temperature_air, cfs.precipitation,
                    temperature[1], M_surface)
                age[1] = dilute_age(age[1], S_surface, cfs.precipitation)
                dz[1] = dz[1] + cfs.precipitation / density_new_snow
                density[1] = (M_surface + cfs.precipitation) / dz[1]
                grain_dendricity[1], grain_sphericity[1], grain_radius[1] = mix_grains(
                    grain_dendricity[1], grain_sphericity[1], grain_radius[1], M_surface,
                    dfall, sfall, refall, cfs.precipitation, independent_grain_radius(mp))
            end

        else
            # Rain joins the pore water at the melt-point enthalpy; its sensible heat above
            # the melting point warms the matrix, so the total enthalpy added is exactly
            # P·h_w(T_air). Density and thickness are untouched.
            rain = cfs.precipitation
            water[1] += rain
            Q_sensible = rain * (specific_enthalpy_water(mp, cfs.temperature_air) -
                                 specific_enthalpy_water(mp))
            temperature[1] = add_energy_temperature(mp, temperature[1], M_surface, Q_sensible)

            # New age-zero mass in the cell on the same footing as snowfall.
            age[1] = dilute_age(age[1], S_surface, rain)
        end

        if verbose
            # Check for conservation of mass
            M = dz .* density
            M_total_final = sum(M) + sum(water)
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
