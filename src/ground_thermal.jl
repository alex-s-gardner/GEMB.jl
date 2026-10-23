# Thermal solve of a `lower_boundary = :ground` column: `ns` snow cells (possibly none) over the
# fixed soil column, solved together as one column whose bottom (the deepest ground cell) is the
# Dirichlet reservoir.

"""
    _ground_turbulent_heat_flux(T1, sfc, cfs, β)
        -> (heat_flux_sensible, heat_flux_latent, latent_heat, dshf, dlhf)

Turbulent fluxes [W m-2, positive into the surface] of a snow-free ground surface at `T1` [K] and
their slopes at fixed wind, from the stability table of `sfc`, with the latent flux scaled by
the bare-soil evaporation factor `β`.

The wind in the transfer law is `U_e = √(U² + w*²)`, with `U` the forcing wind (floored as in
`sfc`) and `w* = (z_i·g·Q_0v/T_a)^{1/3}` the free-convection velocity scale of a heated surface,
`Q_0v = −H/(ρ_a c_p) − 0.61·T_a·LE/(ρ_a L)` the upward kinematic virtual-temperature flux and
`z_i = FREE_CONVECTION_HEIGHT` (Beljaars, 1995; ECMWF, 2018, eqs. 3.17-3.18). Without it the
transfer coefficient of a calm hour vanishes with `U`, and sunlit dry soil has no turbulent
loss at all. `w*` depends on the fluxes it sets, so it is found by fixed-point iteration from
`w* = 1 m s-1`; on a stable or neutral surface it is zero and `U_e = U`. A zero transfer factor
`sfc.C` is a surface sealed from the air: every flux and slope is zero.

# References
- Beljaars, A. C. M. (1995). The parametrization of surface fluxes in large-scale models under
  free convection. *Q. J. R. Meteorol. Soc.* 121, 255-270.
- ECMWF (2018). IFS Documentation CY45R1, Part IV: Physical Processes, §3.2.1.
"""
@inline function _ground_turbulent_heat_flux(T1::Float64, sfc::_ThermalSurface,
    cfs::ClimateForcingStep, β::Float64)
    # A zero transfer factor is a surface that exchanges nothing with the air.
    sfc.C == 0.0 && return 0.0, 0.0, LS, 0.0, 0.0
    U2 = sfc.C / VON_KARMAN^2
    U2 *= U2
    k_q = GRAVITY * FREE_CONVECTION_HEIGHT / cfs.temperature_air
    ρcp = sfc.density_air * HEAT_CAPACITY_AIR
    w = 1.0
    for _ in 1:FREE_CONVECTION_ITERATIONS
        Ue = sqrt(U2 + w * w)
        H, LE, Lh, _, _ = _table_turbulent_heat_flux(sfc.table, T1, sfc.density_air,
            VON_KARMAN^2 * Ue, GRAVITY * cfs.wind_observation_height / (Ue * Ue), cfs)
        Q0v = -H / ρcp - 0.61 * cfs.temperature_air * β * LE / (sfc.density_air * Lh)
        w_new = Q0v > 0.0 ? cbrt(k_q * Q0v) : 0.0
        converged = abs(w_new - w) <= FREE_CONVECTION_TOLERANCE
        w = w_new
        converged && break
    end
    Ue = sqrt(U2 + w * w)
    heat_flux_sensible, heat_flux_latent, latent_heat, dshf, dlhf =
        _table_turbulent_heat_flux(sfc.table, T1, sfc.density_air, VON_KARMAN^2 * Ue,
            GRAVITY * cfs.wind_observation_height / (Ue * Ue), cfs)
    return heat_flux_sensible, β * heat_flux_latent, latent_heat, dshf, β * dlhf
end

"""
    _ground_surface_fluxes_and_slopes(T1, emissivity, sfc, cfs, β)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, dlw, dshf, dlhf)

Surface fluxes [W m-2, positive into the column] of a ground surface at surface-cell temperature
`T1` [K] and their slopes with respect to `T1`, as [`_surface_fluxes_and_slopes`](@ref) gives them
for snow, with three differences: the surface temperature is `T1` itself, not `min(T1, 273.15)`,
since ground warms above the melting point; the latent heat flux is scaled by the bare-soil
evaporation factor `β`; and the turbulent transfer includes free convection
([`_ground_turbulent_heat_flux`](@ref)). The latent heat is `LV` with saturation over water at
or above the melting point and `LS` over ice below it. The turbulent slope pair is set to zero
when its sum would be positive, so the total slope is `≤ 0`.
"""
@inline function _ground_surface_fluxes_and_slopes(T1::Float64, emissivity::Float64,
    sfc::_ThermalSurface, cfs::ClimateForcingStep, β::Float64)
    heat_flux_sensible, heat_flux_latent, latent_heat, dshf, dlhf =
        _ground_turbulent_heat_flux(T1, sfc, cfs, β)
    T2 = T1 * T1
    longwave_upward = -(SB * T2 * T2 * emissivity + (1.0 - emissivity) * cfs.longwave_downward)
    dlw = -4 * SB * emissivity * T2 * T1
    if dshf + dlhf > 0.0
        dshf = dlhf = 0.0
    end
    return longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, dlw, dshf, dlhf
end

"""
    _GroundNodes

The ground half of the thermal solve: per node the solids heat capacity `Cs` [J m-2 K-1], the
freezable water `W` [kg m-2] (the pore water, plus the snow lumped into the top node, see
[`gemb_core_ground`](@ref)), the enthalpy `H` [J m-2] (advanced in place), the conductivity `K`
[W m-1 K-1] and the thickness `dz` [m]. The last node is the Dirichlet reservoir at `T_bottom`.
"""
struct _GroundNodes
    Cs::Vector{Float64}
    W::Vector{Float64}
    H::Vector{Float64}
    K::Vector{Float64}
    dz::Vector{Float64}
    T_bottom::Float64
end

# Temperature of cell `i` of the combined column from its enthalpy `H`: snow cells `1:ns` through
# GEMB's dry and wet-cell maps, ground nodes through `ground_state`.
@inline function _combined_cell_temperature(mp::ModelParameters, gn::_GroundNodes, ns::Int,
    i::Int, H::Float64, M::Vector{Float64}, W::Vector{Float64})
    if i <= ns
        return W[i] > 0.0 ? _wet_cell_temperature(mp, H, M[i], W[i]) :
               temperature_from_specific_enthalpy(mp, H / M[i])
    end
    j = i - ns
    return ground_state(mp, gn.Cs[j], gn.W[j], H)[1]
end

"""
    _thermal_solve_ground!(temperature, dz, density, water, gn, shortwave_flux, sfc,
                           emissivity, emissivity_melt_switch, β, cfs, mp)
        -> (longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal,
            evaporation_condensation, evaporation_condensation_liquid, heat_flux_ground, refrozen)

Advance the snow-and-ground column over one forcing step with the explicit finite-volume
enthalpy scheme of `_thermal_solve!(::ExplicitThermal, ...)`: one flux per face, so the column
total is conserved exactly; harmonic-mean face conductances, which join the snow and the soil at
their interface as Crocus joins its snow base to ISBA's first soil layer (Lafaysse et al., 2026,
eqs. 95-96); a surface flux linearized in the surface cell's temperature change; and the deepest
ground node a Dirichlet reservoir whose enthalpy is left unchanged.

`temperature`, `dz`, `density` and `water` are the `ns` snow cells (possibly none), updated in
place. `gn` holds the ground nodes, whose enthalpies `gn.H` are advanced in place.
`shortwave_flux` [W m-2] covers all `ns + ng` cells. The surface cell is the top snow cell, whose
fluxes are the snow surface balance of [`_surface_fluxes_and_slopes`](@ref), or with no snow cell
the top ground node, whose fluxes are [`_ground_surface_fluxes_and_slopes`](@ref) with
evaporation factor `β`.

Snow cells use GEMB's enthalpy maps, including the freezing plateau of a wet cell
([`_wet_cell_temperature`](@ref)) and [`_wet_cell_handoff!`](@ref) at the end. Ground nodes use
[`ground_state`](@ref), whose plateau freezes and thaws their water at 273.15 K.

The flux across the snow–ground face takes the snow side at `min(T, 273.15)`. A snow cell's
temperature above the melting point is energy that [`calculate_melt`](@ref) turns into melt
after the solve; conducting it into the soil instead would heat the ground from snow at 0 °C.
Crocus evaluates its basal flux the same way (`snowcro.F90`, "Lower boundary flux": `PGBAS` with
`MIN(XTT, ZSNOWTEMP)`).

`heat_flux_ground` [W m-2, positive upward] is the step mean of the conductive flux from the top
ground node into the snow cell above it; with no snow cell it is minus the energy the surface
balance and shortwave deliver to the ground surface. Shortwave transmitted through the snow and
absorbed by the ground is not part of it. `refrozen` [kg m-2] is the snow pore water frozen by
`_wet_cell_handoff!`.

The sub-step is the explicit stability limit of the combined column, a ground node's sensible
heat capacity taken with its water frozen (the smaller of the two), and the surface relaxation
bound of [`_surface_relaxation_dt`](@ref).

# References
- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing,
  Ch. 3-4.
- Lafaysse, M., et al. (2026). Version 3.0.2 of the Crocus snowpack model. *Geosci. Model Dev.*
  19, 6273-6334.
"""
function _thermal_solve_ground!(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64}, gn::_GroundNodes,
    shortwave_flux::Vector{Float64}, sfc::_ThermalSurface, emissivity::Float64,
    emissivity_melt_switch::Bool, β::Float64, cfs::ClimateForcingStep, mp::ModelParameters)

    ns = length(dz)
    ng = length(gn.H)
    m = ns + ng
    length(shortwave_flux) == m || throw(DimensionMismatch(
        "shortwave_flux has $(length(shortwave_flux)) cells, the column $(m)"))

    M = zeros(m)                        # snow matrix mass [kg m-2]
    W = zeros(m)                        # snow pore water [kg m-2]
    H = Vector{Float64}(undef, m)       # cell enthalpy [J m-2]
    K = Vector{Float64}(undef, m)
    dzc = Vector{Float64}(undef, m)
    Cmin = Vector{Float64}(undef, m)    # sensible heat capacity for the stability bound [J m-2 K-1]
    T = Vector{Float64}(undef, m)       # surface-cell temperature lives in T[1]
    h_water = specific_enthalpy_water(mp)

    if ns > 0
        K_snow = thermal_conductivity(temperature, density, mp)
        for i in 1:ns
            M[i] = density[i] * dz[i]
            W[i] = water[i]
            H[i] = M[i] * specific_enthalpy(mp, temperature[i]) + water[i] * h_water
            K[i] = K_snow[i]
            dzc[i] = dz[i]
            Cmin[i] = density[i] * heat_capacity(mp, temperature[i]) * dz[i]
            T[i] = temperature[i]
        end
    end
    for j in 1:ng
        i = ns + j
        H[i] = gn.H[j]
        K[i] = gn.K[j]
        dzc[i] = gn.dz[j]
        Tj = ground_state(mp, gn.Cs[j], gn.W[j], gn.H[j])[1]
        Cmin[i] = gn.Cs[j] + gn.W[j] * heat_capacity(mp, min(Tj, CtoK))
        T[i] = Tj
    end
    T_bottom = gn.T_bottom
    # A wet surface snow cell is held at its solve temperature from the first sub-step.
    if ns > 0 && W[1] > 0.0
        T[1] = _wet_cell_temperature(mp, H[1], M[1], W[1])
    end

    G = Vector{Float64}(undef, m - 1)   # face conductance [W m-2 K-1]
    for i in 1:m-1
        G[i] = 1.0 / (dzc[i+1] / (2 * K[i+1]) + dzc[i] / (2 * K[i]))
    end

    # Stability: diffusion over cells 1 … m-1, and the surface relaxation bound.
    max_safe_dt = Inf
    G_prev = 0.0
    for i in 1:m-1
        max_safe_dt = min(max_safe_dt, Cmin[i] / (G[i] + G_prev))
        G_prev = G[i]
    end
    Λ0 = if ns > 0
        _surface_energy_balance_and_slope(T[1], emissivity, sfc, cfs)[6]
    else
        r = _ground_surface_fluxes_and_slopes(T[1], emissivity, sfc, cfs, β)
        r[5] + r[6] + r[7]
    end
    C1_0 = ns > 0 ? (M[1] + W[1]) * heat_capacity(mp, T[1]) : Cmin[ns+1]
    max_safe_dt = min(max_safe_dt, THERMAL_SURFACE_STEP_FRACTION * C1_0 / (G[1] + abs(Λ0)))
    dt = _find_dt_divisor(max_safe_dt * mp.thermal_explicit_safety_factor, mp.dt_divisors)
    n_steps = round(Int, cfs.dt / dt)

    Q_sw = [shortwave_flux[i] * dt for i in 1:m]
    Q_sw[m-1] += Q_sw[m]
    Q_sw[m] = 0.0
    A = G .* dt
    longwave_downward = cfs.longwave_downward * dt

    lw_cum = 0.0; shf_cum = 0.0; lhf_cum = 0.0; basal_cum = 0.0; ground_cum = 0.0
    EC_cum = 0.0; EC_liquid_cum = 0.0
    T_pre = Vector{Float64}(undef, m)

    for _ in 1:n_steps
        T1_0 = T[1]
        if ns > 0
            liquid = _surface_liquid(mp, H[1], M[1], W[1]) + EC_liquid_cum > 0.0
            lw_up, shf, lhf, latent_heat, _, dlw, dshf, dlhf =
                _surface_fluxes_and_slopes(T1_0, emissivity, sfc, cfs, liquid)
            C1 = (M[1] + W[1]) * heat_capacity(mp, T1_0)
        else
            lw_up, shf, lhf, latent_heat, dlw, dshf, dlhf =
                _ground_surface_fluxes_and_slopes(T1_0, emissivity, sfc, cfs, β)
            C1 = ground_heat_capacity(mp, gn.Cs[1], gn.W[1], T1_0)
        end
        Λ = dlw + dshf + dlhf
        if Λ < 0.0
            F1 = A[1] * (_combined_cell_temperature(mp, gn, ns, 2, H[2], M, W) - T1_0)
            Q0 = longwave_downward + (lw_up + shf + lhf) * dt
            ΔT1 = (Q_sw[1] + Q0 + F1) / (C1 - Λ * dt)
            lw_up += dlw * ΔT1
            shf += dshf * ΔT1
            lhf += dlhf * ΔT1
        end
        lhf_cum += lhf * dt
        shf_cum += shf * dt
        lw_cum -= lw_up * dt
        EC = lhf / latent_heat * dt
        EC_cum += EC
        latent_heat == LV && (EC_liquid_cum += EC)
        Q_surface = longwave_downward + lw_up * dt + (shf + lhf) * dt

        for i in 1:m-1
            T_pre[i] = _combined_cell_temperature(mp, gn, ns, i, H[i], M, W)
        end
        H[1] += Q_sw[1] + Q_surface
        ns == 0 && (ground_cum -= Q_sw[1] + Q_surface)
        F_prev = 0.0
        for i in 1:m-2
            # The snow side of the snow–ground face is at most the melting point.
            F = A[i] * (T_pre[i+1] - (i == ns ? min(T_pre[i], CtoK) : T_pre[i]))
            H[i] += F - F_prev
            H[i+1] += Q_sw[i+1]
            i == ns && (ground_cum += F)
            F_prev = F
        end
        F_basal = A[m-1] * (T_bottom - (m - 1 == ns ? min(T_pre[m-1], CtoK) : T_pre[m-1]))
        m - 1 == ns && (ground_cum += F_basal)
        basal_cum += F_basal
        H[m-1] += F_basal - F_prev

        T[1] = _combined_cell_temperature(mp, gn, ns, 1, H[1], M, W)
        if emissivity_melt_switch
            emissivity = T[1] < (CtoK - T_MELT_SWITCH_TOLERANCE) ? mp.emissivity :
                         mp.emissivity_grain_radius_large
        end
    end

    refrozen = zeros(ns)
    for i in 1:ns
        if water[i] > 0.0
            _wet_cell_handoff!(mp, temperature, dz, density, water, refrozen, i, H[i], M[i])
        else
            temperature[i] = temperature_from_specific_enthalpy(mp, H[i] / M[i])
        end
    end
    for j in 1:ng-1
        gn.H[j] = H[ns+j]
    end

    return lw_cum / cfs.dt, shf_cum / cfs.dt, lhf_cum / cfs.dt, basal_cum / cfs.dt, EC_cum,
        EC_liquid_cum, ground_cum / cfs.dt, sum(refrozen)
end
