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
        # Rate coefficients `_hl_c0`/`_hl_c1` are shared with `steady_state_density`.
        @inbounds for i in 1:m
            T = temperature[i]
            c = density[i] <= 550.0 + d_tolerance ? _hl_c0(T, pm) : _hl_c1(T, pm)
            _densify_cell!(density, dz, dz_out, i, c, dt, density_ice, d_tolerance)
        end

    elseif method == :Arthern
        precip_force = pm * 9.81
        @inbounds for i in 1:m
            c = _arthern_c(density[i], temperature[i], tam, precip_force)
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
        # M0/M1 depend only on the climatological accumulation, so hoist them via
        # DensificationCoeffs (shared with `steady_state_density`).
        p = DensificationCoeffs(mp, pm, tam)
        @inbounds for i in 1:m
            c = _ligtenberg_c(density[i], temperature[i], tam, p.precip_force, p.M0, p.M1)
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

# ---------------------------------------------------------------------------
# Densification rate coefficients `c` [yr-1]
#
# These are the single source of truth for the three methods `validate_parameters`
# permits: they are called both by the transient branches of `calculate_density`
# above and by the steady-state marcher (`steady_state_density`), so the two can
# never drift. Each is a plain `@inline` function of the layer state, and the
# loop-invariant terms are hoisted by the caller into `DensificationCoeffs`, so
# routing the hot loop through them costs nothing relative to the inline form.
#
# `T` is layer temperature [K], `ρ` layer density [kg m-3], `pm` mean annual
# accumulation [kg m-2 yr-1], `tam` mean annual air temperature [K].
# `R_GAS = 8.314 J mol-1 K-1`.
# ---------------------------------------------------------------------------

# Herron & Langway (1980) two-stage coefficients, equivalent to Arthern et al.
# (2010) eq. (2).
@inline _hl_c0(T::Float64, pm::Float64) = (11.0 * exp(-10160.0 / (T * R_GAS))) * (pm / 1000.0)          # ρ <= 550
@inline _hl_c1(T::Float64, pm::Float64) = (575.0 * exp(-21400.0 / (T * R_GAS))) * sqrt(pm / 1000.0)     # ρ > 550

# Arthern et al. (2010) semi-empirical. `precip_force = pm * 9.81`, hoisted.
@inline function _arthern_c(ρ::Float64, T::Float64, tam::Float64, precip_force::Float64)
    H = exp((-60000.0 / (T * R_GAS)) + (42400.0 / (tam * R_GAS))) * precip_force
    return (ρ <= 550.0 + D_TOLERANCE ? 0.07 : 0.03) * H
end

# Ligtenberg et al. (2011): Arthern scaled by the accumulation-dependent M0/M1.
@inline function _ligtenberg_c(ρ::Float64, T::Float64, tam::Float64,
    precip_force::Float64, M0::Float64, M1::Float64)
    H = exp((-60000.0 / (T * R_GAS)) + (42400.0 / (tam * R_GAS))) * precip_force
    return ρ <= 550.0 + D_TOLERANCE ? M0 * (0.07 * H) : M1 * (0.03 * H)
end

"""
    DensificationCoeffs(mp::ModelParameters, pm, tam)

Loop-invariant densification terms for `mp.densification_method`, hoisted out of
the per-cell / per-age-step loop: the accumulation forcing `pm * 9.81` and the
Ligtenberg `M0`/`M1` scalings. Shared by [`calculate_density`](@ref) and
[`steady_state_density`](@ref) so both derive `c` from identical inputs.

`pm` is mean annual accumulation [kg m-2 yr-1] and `tam` mean annual air
temperature [K].
"""
struct DensificationCoeffs
    method::Symbol
    pm::Float64
    tam::Float64
    precip_force::Float64
    M0::Float64
    M1::Float64
end

function DensificationCoeffs(mp::ModelParameters, pm::Real, tam::Real)
    pm = Float64(pm)
    tam = Float64(tam)
    M0 = M1 = NaN
    if mp.densification_method == :Ligtenberg
        M01 = densification_lookup_M01(mp.densification_coeffs_M01)
        if length(M01) == 4
            M0 = max(M01[1] - (M01[2] * log(pm)), 0.25)
            M1 = max(M01[3] - (M01[4] * log(pm)), 0.25)
        elseif abs(mp.density_ice - 820.0) < D_TOLERANCE
            M0 = max(M01[1, 1] - (M01[1, 2] * log(pm)), 0.25)
            M1 = max(M01[1, 3] - (M01[1, 4] * log(pm)), 0.25)
        else
            M0 = max(M01[2, 1] - (M01[2, 2] * log(pm)), 0.25)
            M1 = max(M01[2, 3] - (M01[2, 4] * log(pm)), 0.25)
        end
    end
    return DensificationCoeffs(mp.densification_method, pm, tam, pm * 9.81, M0, M1)
end

"""
    _densification_rate(p::DensificationCoeffs, ρ, T) -> c [yr-1]

Densification rate coefficient for the method in `p`, dispatching on
`p.method`. Used by [`steady_state_density`](@ref), where the branch is taken
once per age step rather than once per cell. `calculate_density`'s hot loop calls
the per-method functions directly, so it takes the branch once per call.
"""
@inline function _densification_rate(p::DensificationCoeffs, ρ::Float64, T::Float64)
    if p.method == :HerronLangway
        return ρ <= 550.0 + D_TOLERANCE ? _hl_c0(T, p.pm) : _hl_c1(T, p.pm)
    elseif p.method == :Arthern
        return _arthern_c(ρ, T, p.tam, p.precip_force)
    elseif p.method == :Ligtenberg
        return _ligtenberg_c(ρ, T, p.tam, p.precip_force, p.M0, p.M1)
    end
    error("unrecognized densification method: $(p.method)")
end

"""
    herron_langway_steady_state(z_center, T_mean, accumulation, ρ0, ρ_ice;
                                n_age=2000)

Analytic steady-state Herron & Langway (1980) firn density profile.

Returns a `Vector{Float64}` of density [kg m-3] at each depth in `z_center`
(cell-center heights, negative below the surface as produced by [`dz2z`](@ref)).

This is the fixed point of the same densification law GEMB integrates transiently
(`∂ρ/∂t = c(ρᵢ − ρ)`, Arthern et al. 2010 eq. 1) using the shared
[`_hl_c0`](@ref)/[`_hl_c1`](@ref) rate coefficients. Under steady accumulation the
column is described by the age of a parcel `t` [yr]:

- Density integrates the two-stage exponential relaxation
  `ρ(t) = ρᵢ − (ρᵢ − ρ₀)·exp(−c·t)`, with `c = c₀` while `ρ ≤ 550` and `c = c₁`
  (restarted from `ρ₅₅₀`) above.
- Depth follows from mass balance, `z(t) = A·∫₀ᵗ ds/ρ(s)`, where `A = accumulation`
  [kg m-2 yr-1] is the mass of firn deposited per year.

The `(z, ρ)` curve is marched on a fine age grid and interpolated onto `z_center`.
With negligible accumulation (`accumulation ≤ 0`) no firn forms and the profile is
returned as pure ice (`ρ_ice`).

# Arguments
- `z_center`: cell-center heights [m], negative below surface.
- `T_mean`: mean annual temperature [K].
- `accumulation`: mean annual accumulation [kg m-2 yr-1].
- `ρ0`: surface (fresh-snow) density [kg m-3].
- `ρ_ice`: density of ice [kg m-3], used as the asymptote / clamp.

# References
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
"""
function herron_langway_steady_state(z_center::AbstractVector, T_mean::Real,
    accumulation::Real, ρ0::Real, ρ_ice::Real; n_age::Int=2000)

    depth = -Float64.(z_center)                     # positive downward depths [m]
    m = length(depth)

    # No accumulation → no firn column forms; return pure ice.
    if accumulation <= 0.0
        return fill(Float64(ρ_ice), m)
    end

    T = Float64(T_mean)
    A = Float64(accumulation)
    ρi = Float64(ρ_ice)
    ρ = clamp(Float64(ρ0), 1.0, ρi)                 # surface density, guarded

    c0 = _hl_c0(T, A)
    c1 = _hl_c1(T, A)
    d_tol = 1e-11

    max_depth = maximum(depth)

    # March a parcel forward in age, accumulating (z, ρ). Depth increment over a
    # small age step dt uses the trapezoidal thickness A*dt / mean(ρ).
    z_curve = Float64[0.0]
    ρ_curve = Float64[ρ]
    t = 0.0
    z = 0.0
    dt = 0.02                                        # age step [yr]
    iter = 0
    max_iter = n_age * 50
    @inbounds while z < max_depth && ρ < ρi - d_tol && iter < max_iter
        c = ρ <= 550.0 + d_tol ? c0 : c1
        ρ_next = ρi - (ρi - ρ) * exp(-c * dt)
        if ρ_next > ρi
            ρ_next = ρi
        end
        ρ_mid = 0.5 * (ρ + ρ_next)
        z += A * dt / ρ_mid
        t += dt
        ρ = ρ_next
        push!(z_curve, z)
        push!(ρ_curve, ρ)
        iter += 1
    end
    # Ensure the curve spans past the deepest requested cell.
    if z_curve[end] < max_depth
        push!(z_curve, max_depth + 1.0)
        push!(ρ_curve, ρi)
    end

    # Interpolate ρ(z) onto the model depths; constant extrapolation handles the
    # z=0 surface (z_curve[1]=0) and any tiny overshoot at depth.
    itp = DataInterpolations.LinearInterpolation(ρ_curve, z_curve;
        extrapolation=DataInterpolations.ExtrapolationType.Constant)

    out = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        out[i] = clamp(itp(depth[i]), Float64(ρ0), ρi)
    end
    return out
end
