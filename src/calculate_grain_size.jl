"""
    calculate_grain_size(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, cfs::ClimateForcingStep, mp::ModelParameters)

Model the evolution of effective snow grain size, dendricity, and sphericity.

Accounts for different physical processes depending on snow state:
- Dendritic Snow (fresh): Evolves based on temperature gradients and liquid water content.
- Nondendritic Dry Snow: Temperature gradient metamorphism using Marbouty (1980), or
  normal grain growth after Arthern et al. (2010), selected by `mp.grain_growth_method`.
- Wet Snow: Rapid grain growth due to liquid water using Brun (1989).

`mp.grain_growth_method` governs the non-dendritic *dry* branch only; the dendritic and wet
branches, and the sphericity caps, are the same on every setting.

- `:Marbouty` (default, MATLAB's behaviour) — Marbouty (1980) everywhere. Its density factor
  `H` vanishes at `DENSITY_MARBOUTY_MAX`, so grain radius is frozen throughout the firn
  column below a few metres.
- `:Arthern` — `dr²/dt = kgr·exp(-Eg/RT)` (`GRAIN_GROWTH_KGR`, `GRAIN_GROWTH_EG`) at every
  density, dropping the temperature-gradient dependence entirely.
- `:hybrid` — Marbouty below `DENSITY_MARBOUTY_MAX`, Arthern at or above it. Seasonal snow
  keeps the temperature-gradient physics; firn grains keep growing, which matters because
  `:ArthernB` densification goes as `1/r²`.

Runs unconditionally. Metamorphism is not contingent on model configuration — grains coarsen
in a real snowpack whichever albedo scheme is selected — and `grain_radius` is read by four
schemes across three modules ([`calculate_albedo`](@ref),
[`calculate_shortwave_radiation`](@ref), `:ArthernB`/`:Crocus`/`:CrocusPure` densification in
[`calculate_density`](@ref), and the grain-radius emissivity methods in
[`calculate_temperature`](@ref)). MATLAB skips this function unless `albedo_method` is
`:GardnerSharp` or `:BrunLefebre`, which covers only the first two, so the other
configurations silently ran with grain size frozen at the initial profile. Enumerating the
consumers to skip the work costs 11% on the runs that can skip it and re-breaks silently
whenever a new consumer is added, so the work is simply always done.

Returns `(grain_radius, grain_dendricity, grain_sphericity)`. **All three are updated in
place**, and the returned `grain_radius` is the same array that was passed in — pass a copy
if the incoming values are still needed. `grain_radius` used to be returned as a fresh array;
it now shares the mutation convention of the other two, which removed the largest allocation
in the function (see the `gsz` comment below).

This is a scalar-loop implementation that is numerically identical, element by
element, to the reference vectorized MATLAB translation, but avoids the ~30 mask
/ gather / broadcast temporaries the vectorized form allocated per call.

# References
- Brun, E., et al. (1992). J. Glaciol., 38, 13-22.
- Marbouty, D. (1980). J. Glaciol., 26, 303-312.
- Brun, E. (1989). Ann. Glaciol., 13, 22-26.
"""
function calculate_grain_size(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    m = length(temperature)

    # Convert dt from seconds to days
    dt_days = cfs.dt / 86400.0

    # Work in grain size (diameter) [mm], converting `grain_radius` in place and halving it
    # back at the end. Nothing reads the incoming radius after this doubling, and both
    # callers rebind the returned vector to the same array they passed in, so the previous
    # `similar(grain_radius)` scratch array was a copy no one needed — it was the largest
    # allocation in the function (8 bytes/cell, ~73% of the total).
    gsz = grain_radius
    @inbounds @simd for i in 1:m
        gsz[i] *= 2
    end

    # Classify cells once (fixed for the whole call, matching the vectorized
    # G/J masks that were computed from the *initial* dendricity). This cannot be
    # recomputed on the fly: the dendritic branch clamps `grain_dendricity` to 0 below, so
    # by the time the non-dendritic branches test membership the distinction is gone.
    # Deliberately a `Vector{Bool}` (1 byte/cell) rather than a `BitVector` (1 bit/cell): the
    # BitVector was measured 15x slower to fill at this column size (556 ns vs 38 ns) and
    # allocates one object *more*, since its wrapper and chunk array are separate. Masked
    # bit `setindex!` does not pay for itself against 7 bytes/cell.
    isG = Vector{Bool}(undef, m)
    anyG = false
    anyJ = false
    @inbounds for i in 1:m
        g = grain_dendricity[i] > 0 + GDN_TOLERANCE
        isG[i] = g
        anyG |= g
        anyJ |= !g
    end

    ## DENDRITIC SNOW METAMORPHISM (dendricity > 0)
    if anyG
        @inbounds for i in 1:m
            isG[i] || continue
            dTi = _grain_gradient(temperature, dz, i, m)
            wi = water[i]
            if wi <= 0 + WATER_TOLERANCE
                ex = exp(-6e3 / temperature[i])
                if dTi <= 5 + T_TOLERANCE
                    # dT < 5 degC m-1
                    A = -2e8 * ex * dt_days
                    B = 1e9 * ex * dt_days
                    grain_dendricity[i] += A
                    grain_sphericity[i] += B
                else
                    # dT >= 5 degC m-1
                    C = (-2e8 * ex * dt_days) * (dTi^0.4)
                    grain_dendricity[i] += C
                    grain_sphericity[i] += C
                end
            else
                # dendritic wet snow
                lwci = _grain_lwc(water[i], density[i], dz[i])
                D = (1 / 16) * (lwci^3) * dt_days
                grain_dendricity[i] -= D
                grain_sphericity[i] += D
            end
        end

        # dendricity and sphericity can not be > 1 or < 0 (applies to all cells)
        @inbounds for i in 1:m
            gd = grain_dendricity[i]
            if gd <= 0 + GDN_TOLERANCE
                grain_dendricity[i] = 0.0
            elseif gd >= 1 - GDN_TOLERANCE
                grain_dendricity[i] = 1.0
            end
            gs = grain_sphericity[i]
            if gs <= 0 + GDN_TOLERANCE
                grain_sphericity[i] = 0.0
            elseif gs >= 1 - GDN_TOLERANCE
                grain_sphericity[i] = 1.0
            end
        end

        # new grain size (mm) for dendritic cells, using clamped values
        @inbounds for i in 1:m
            isG[i] || continue
            gd = grain_dendricity[i]
            gs = grain_sphericity[i]
            gsz[i] = max(0.1 * (gd / 0.99 + (1.0 - 1.0 * gd / 0.99) *
                        (gs / 0.99 * 3.0 + (1.0 - gs / 0.99) * 4.0)), GDN_TOLERANCE * 2)
        end
    end

    ## NONDENDRITIC SNOW (dendricity == 0)
    if anyJ
        # Wet-snow grains (class 6) sphericity evolution (Brun eq. 4 regimes)
        @inbounds for i in 1:m
            isG[i] && continue
            gs = grain_sphericity[i]
            if gs > GDN_TOLERANCE && gs < 1 - GDN_TOLERANCE
                dTi = _grain_gradient(temperature, dz, i, m)
                if dTi > 5 + T_TOLERANCE
                    F1 = (-2e8 * exp(-6e3 / temperature[i]) * dt_days) * dTi^0.4
                    grain_sphericity[i] += F1
                elseif water[i] > 0 + WATER_TOLERANCE
                    lwci = _grain_lwc(water[i], density[i], dz[i])
                    F2 = (1.0 / 16.0) * lwci^3.0 * dt_days
                    grain_sphericity[i] += F2
                else
                    F3 = 1e9 * exp(-6e3 / temperature[i]) * dt_days
                    grain_sphericity[i] += F3
                end
            end
        end

        # sphericity can not be > 1 or < 0 (applies to all cells)
        @inbounds for i in 1:m
            gs = grain_sphericity[i]
            if gs <= 0 + GDN_TOLERANCE
                grain_sphericity[i] = 0.0
            elseif gs >= 1 - GDN_TOLERANCE
                grain_sphericity[i] = 1.0
            end
        end

        # Dry-snow (Marbouty 1980 / Arthern 2010) and wet-snow (Brun 1989) grain growth
        dt_seconds = cfs.dt
        method = mp.grain_growth_method
        @inbounds for i in 1:m
            isG[i] && continue
            dTi = _grain_gradient(temperature, dz, i, m)
            if (water[i] <= 0 + WATER_TOLERANCE) ||
               ((grain_sphericity[i] <= 0 + GDN_TOLERANCE) && (dTi > 5 + T_TOLERANCE))
                # DRY SNOW METAMORPHISM
                if method === :Marbouty ||
                   (method === :hybrid && density[i] < DENSITY_MARBOUTY_MAX - D_TOLERANCE)
                    # Marbouty (1980) temperature-gradient metamorphism
                    Q = _marbouty_Q(temperature[i], density[i], dTi)
                    gsz[i] += Q * dt_days
                else
                    # Arthern et al. (2010): dr²/dt = kgr·exp(-Eg/RT), in radius, so halve
                    # the working diameter and double back.
                    r = gsz[i] / 2
                    gsz[i] = 2 * sqrt(r * r + GRAIN_GROWTH_KGR *
                                            exp(-GRAIN_GROWTH_EG / (R_GAS * temperature[i])) *
                                            dt_seconds * 1e6)
                end
            else
                # WET SNOW METAMORPHISM (Brun, 1989)
                lwci = _grain_lwc(water[i], density[i], dz[i])
                E = (1.28e-8 + 4.22e-10 * (lwci^3)) * (dt_days * 86400)   # [mm^3 s^-1]
                gsz[i] = 2 * (3 / (pi * 4) * ((4 / 3) * pi * (gsz[i] / 2)^3 + E))^(1 / 3)
            end
        end

        # grain-size caps by sphericity (Brun, 1992), applied to all cells
        @inbounds for i in 1:m
            if abs(grain_sphericity[i] - 1) < WATER_TOLERANCE
                # spherical grains: <= 2 mm
                if gsz[i] > 2 - WATER_TOLERANCE
                    gsz[i] = 2.0
                end
            else
                # non-spherical grains: <= 5 mm
                if gsz[i] > 5 - WATER_TOLERANCE
                    gsz[i] = 5.0
                end
            end
        end
    end

    # Back from diameter to effective radius, in the caller's own array.
    @inbounds @simd for i in 1:m
        gsz[i] /= 2
    end

    return gsz, grain_dendricity, grain_sphericity
end

"""
    _grain_gradient(temperature, dz, i, m)

Absolute temperature gradient [degC m-1] at cell `i` of a column of `m` cells,
using the same finite-difference stencil as the reference implementation. The
grid-point-center separations reduce to local `dz` combinations, so no cumulative
depth array is needed.
"""
@inline function _grain_gradient(temperature::Vector{Float64}, dz::Vector{Float64}, i::Int, m::Int)
    @inbounds begin
        if m <= 1
            return 0.0
        elseif m == 2
            denom = dz[1] / 2 + dz[2] / 2
            return abs((temperature[2] - temperature[1]) / denom)
        else
            if i == 1
                denom = dz[1] / 2 + dz[2] + dz[3] / 2
                return abs((temperature[3] - temperature[1]) / denom)
            elseif i == m
                denom = dz[m-2] / 2 + dz[m-1] + dz[m] / 2
                return abs((temperature[m] - temperature[m-2]) / denom)
            else
                denom = dz[i-1] / 2 + dz[i] + dz[i+1] / 2
                return abs((temperature[i+1] - temperature[i-1]) / denom)
            end
        end
    end
end

"""
    _grain_lwc(water, density, dz)

Liquid-water content as a mass fraction [%] for a single cell, capped at 9%
(Brun, 1980).
"""
@inline function _grain_lwc(water::Float64, density::Float64, dz::Float64)
    lwc = water / (density * dz) * 100
    return lwc > (9 + WATER_TOLERANCE) ? 9.0 : lwc
end

"""
    _marbouty_Q(temperature, density, dT)

Grain-growth rate coefficient Q [mm d-1] for a single cell per Fig. 9 of
Marbouty (1980). No grain growth above `DENSITY_MARBOUTY_MAX` (H = 0).
Scalar equivalent of the reference `_Marbouty`.
"""
@inline function _marbouty_Q(temperature::Float64, density::Float64, dT::Float64)

    E = 0.09       # model time growth constant [mm d-1]
    T = temperature - 273.15   # K to C
    dTc = dT / 100.0           # degC/m to degC/cm

    ## Temperature coefficient F
    F = 0.0
    if T > -6 + T_TOLERANCE
        F = 0.7 + ((T / -6) * 0.3)
    elseif T > -22 + T_TOLERANCE
        F = 1 - ((T + 6) / -16 * 0.8)
    elseif T > -40 + T_TOLERANCE
        F = 0.2 - ((T + 22) / -18 * 0.2)
    end

    ## Density coefficient H
    # Ramps 1 → 0 over [150, DENSITY_MARBOUTY_MAX]; the span is derived from the ceiling
    # rather than restated as 250 so the two cannot drift apart.
    H = 0.0
    if density < 150 - D_TOLERANCE
        H = 1.0
    elseif density < DENSITY_MARBOUTY_MAX - D_TOLERANCE
        H = 1 - ((density - 150) / (DENSITY_MARBOUTY_MAX - 150))
    end

    ## Temperature gradient coefficient G
    G = 0.0
    if dTc >= 0.7 - T_TOLERANCE
        G = 1.0
    elseif dTc >= 0.50 - T_TOLERANCE
        G = 0.90 + (((dTc - 0.50) / 0.20) * 0.1)
    elseif dTc >= 0.40 - T_TOLERANCE
        G = 0.67 + (((dTc - 0.40) / 0.10) * 0.23)
    elseif dTc >= 0.25 - T_TOLERANCE
        G = 0.10 + (((dTc - 0.25) / 0.15) * 0.57)
    elseif dTc >= 0.16 - T_TOLERANCE
        G = ((dTc - 0.16) / 0.09) * 0.1
    end

    return F * H * G * E
end
