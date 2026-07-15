"""
    calculate_density(temperature, dz, density, grain_radius, cfs::ClimateForcingStep, mp::ModelParameters)

Compute the densification of snow/firn using one of several models:
- `:HerronLangway`: Herron and Langway (1980)
- `:Arthern`: semi-empirical model of Arthern et al. (2010) [default]
- `:ArthernB`: physical model from Appendix B of Arthern et al. (2010)
- `:LiZwally`: empirical model of Li and Zwally (2004)
- `:Helsen`: modified empirical model by Helsen et al. (2008)
- `:Ligtenberg`: semi-empirical model of Ligtenberg et al. (2011)

Returns `(dz, density)`. `density` is updated in place; `dz` is returned as a
new array (recomputed from the conserved cell mass).

This is a scalar-loop implementation that is numerically identical, element by
element, to the reference vectorized MATLAB translation, but avoids the
mask / gather / broadcast temporaries (`mass_init`, `idx`, `H`, `c0`, `c1`, ...)
the vectorized form allocated per call.

# References
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
- Li, J. and Zwally, H. (2004). Ann. Glaciol., 38, 309-313.
- Helsen, M. M., et al. (2008). Science, 320, 1626-1629.
- Ligtenberg, S. R. M., et al. (2011). The Cryosphere, 5, 809-819.
"""
function calculate_density(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, grain_radius::Vector{Float64},
    cfs::ClimateForcingStep, mp::ModelParameters)

    d_tolerance = 1e-11

    # specify constants
    dt = cfs.dt / 86400.0   # convert from [s] to [d]
    R = 8.314               # gas constant [mol-1 K-1]

    m = length(density)
    density_ice = mp.density_ice
    pm = cfs.precipitation_mean
    tam = cfs.temperature_air_mean

    # New grid-cell lengths (fresh array; density is updated in place).
    dz_out = similar(dz)

    method = mp.densification_method

    if method == :HerronLangway
        pm_frac = pm / 1000
        sqrt_term = (pm / 1000)^0.5
        @inbounds for i in 1:m
            T = temperature[i]
            if density[i] <= 550.0 + d_tolerance
                c = (11 * exp(-10160 / (T * R))) * pm / 1000
            else
                c = (575 * exp(-21400 / (T * R))) * sqrt_term
            end
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :Arthern
        precip_force = pm * 9.81
        @inbounds for i in 1:m
            T = temperature[i]
            H = exp((-60000.0 / (T * R)) + (42400.0 / (tam * R))) * precip_force
            c = (density[i] <= 550.0 + d_tolerance ? 0.07 : 0.03) * H
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :ArthernB
        # Overburden pressure, replicating the reference exactly:
        #   obp[1] = 0; obp[i] = (cumulative dz through i-1) * density[i-1]
        # i.e. cumulative depth times the density of the immediately overlying cell.
        cumdz = 0.0      # sum of dz over cells above the current one
        prev_d = 0.0     # original density of the immediately overlying cell
        @inbounds for i in 1:m
            d0 = density[i]
            T = temperature[i]
            gr = grain_radius[i] / 1000
            obp = i == 1 ? 0.0 : cumdz * prev_d
            H = exp((-60000.0 / (T * R))) * obp / gr^2
            c = (d0 <= 550.0 + d_tolerance ? 9.2e-9 : 3.7e-9) * H
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
            # advance the running overburden terms using original values
            cumdz += dz[i]
            prev_d = d0
        end

    elseif method == :LiZwally
        base = (pm / density_ice) * max(139.21 - 0.542 * tam, 1.0) * 8.36
        @inbounds for i in 1:m
            c = base * max(CtoK - temperature[i], 1.0)^(-2.061)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :Helsen
        base = (pm / density_ice) * max(76.138 - 0.28965 * tam, 1.0) * 8.36
        @inbounds for i in 1:m
            c = base * max(CtoK - temperature[i], 1.0)^(-2.061)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :Ligtenberg
        precip_force = pm * 9.81
        M01 = densification_lookup_M01(mp.densification_coeffs_M01)
        if length(M01) == 4
            M0 = max(M01[1] - (M01[2] * log(pm)), 0.25)
            M1 = max(M01[3] - (M01[4] * log(pm)), 0.25)
        else
            if abs(density_ice - 820.0) < d_tolerance
                M0 = max(M01[1, 1] - (M01[1, 2] * log(pm)), 0.25)
                M1 = max(M01[1, 3] - (M01[1, 4] * log(pm)), 0.25)
            else
                M0 = max(M01[2, 1] - (M01[2, 2] * log(pm)), 0.25)
                M1 = max(M01[2, 3] - (M01[2, 4] * log(pm)), 0.25)
            end
        end
        @inbounds for i in 1:m
            T = temperature[i]
            H = exp((-60000.0 / (T * R)) + (42400.0 / (tam * R))) * precip_force
            c = density[i] <= 550.0 + d_tolerance ? M0 * (0.07 * H) : M1 * (0.03 * H)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    else
        error("unrecognized densification method")
    end

    return dz_out, density
end

"""
    _densify_cell!(density, dz_in, dz_out, i, c, dt, density_ice, d_tolerance)

Apply the densification increment for cell `i` given rate coefficient `c`:
update `density[i]` in place, clamp it to the density of ice, and write the
mass-conserving new grid-cell length to `dz_out[i]`.
"""
@inline function _densify_cell!(density::Vector{Float64}, dz_in::Vector{Float64},
    dz_out::Vector{Float64}, i::Int, c::Float64, dt::Float64,
    density_ice::Float64, d_tolerance::Float64)
    @inbounds begin
        d0 = density[i]
        mass = d0 * dz_in[i]
        d = d0 + (c * (density_ice - d0) / 365 * dt)
        if d > density_ice - d_tolerance
            d = density_ice
        end
        density[i] = d
        dz_out[i] = mass / d
    end
end
