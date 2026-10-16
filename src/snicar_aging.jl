# Snow aging of SNICAR as implemented in CLM5.0 (SnowSnicarMod.F90, subroutine SnowAge_grain;
# Lawrence et al., 2018, CLM5.0 Technical Note, section 2.3.2.3, eqs. 2.3.72-2.3.75), selected by
# `grain_growth_method = :FlannerZender2006`. The dry-snow law is the radius form of the SSA decay
# of Flanner and Zender (2006, eq. 16), with parameters looked up per cell from the table Flanner
# fitted to their microphysical model over 223-273 K, 0-300 K m-1 and 50-400 kg m-3.

const SNICAR_AGING_TABLE_FILE = joinpath(@__DIR__, "..", "data", "snicar_drdt_bst_fit_60_c070416.txt")
include_dependency(SNICAR_AGING_TABLE_FILE)

# Nodes of the aging table: snow density, temperature gradient and temperature, in the index
# order of CLM's `snowage_tau(rhos_idx, Tgrd_idx, T_idx)`.
const SNICAR_AGING_DENSITY = range(50.0, 400.0; length = 8)      # kg m-3
const SNICAR_AGING_GRADIENT = range(0.0, 300.0; length = 31)     # K m-1
const SNICAR_AGING_TEMPERATURE = range(223.15, 273.15; length = 11)  # K

"""
    _read_snicar_aging_table(path) -> (tau, kappa, drdt0)

Read the SNICAR aging parameters τ [µm], κ [-] and (dr/dt)₀ [µm h-1] from the text copy of
`snicar_drdt_bst_fit_60_c070416.nc` at `path` (one row per node: temperature, gradient,
density, τ, κ, (dr/dt)₀), into arrays indexed `[density, gradient, temperature]` on the nodes
`SNICAR_AGING_DENSITY`, `SNICAR_AGING_GRADIENT` and `SNICAR_AGING_TEMPERATURE`. Throws if a row
is malformed, lies off the grid, repeats a node, or if any node is missing.
"""
function _read_snicar_aging_table(path::AbstractString)
    dims = (length(SNICAR_AGING_DENSITY), length(SNICAR_AGING_GRADIENT),
            length(SNICAR_AGING_TEMPERATURE))
    tau = fill(NaN, dims)
    kappa = fill(NaN, dims)
    drdt0 = fill(NaN, dims)
    node(x, grid) = (k = round(Int, (x - first(grid)) / step(grid)) + 1;
                     (k in eachindex(grid) && abs(grid[k] - x) <= 1e-9) ? k :
                     error("$path: $x is not a node of $grid"))
    n = 0
    for line in eachline(path)
        (isempty(line) || startswith(line, '#')) && continue
        v = parse.(Float64, split(line))
        length(v) == 6 || error("$path: expected 6 columns, got $(length(v)): $line")
        k = node(v[1], SNICAR_AGING_TEMPERATURE)
        j = node(v[2], SNICAR_AGING_GRADIENT)
        i = node(v[3], SNICAR_AGING_DENSITY)
        isnan(tau[i, j, k]) || error("$path: node (T, dT/dz, ρ) = $(v[1:3]) repeats")
        tau[i, j, k], kappa[i, j, k], drdt0[i, j, k] = v[4], v[5], v[6]
        n += 1
    end
    n == prod(dims) || error("$path: $n rows for $(prod(dims)) nodes")
    return tau, kappa, drdt0
end

const SNICAR_AGING_TAU, SNICAR_AGING_KAPPA, SNICAR_AGING_DRDT0 =
    _read_snicar_aging_table(SNICAR_AGING_TABLE_FILE)

# Fortran `nint`: nearest integer, halves away from zero.
@inline _nint(x::Float64) = round(Int, x, RoundNearestTiesAway)

"""
    snicar_aging_parameters(T, gradient, density) -> (τ [µm], κ, (dr/dt)₀ [µm h-1])

Aging parameters at snow temperature `T` [K], temperature gradient `gradient` [K m-1] and
density `density` [kg m-3], taken from the nearest node of the SNICAR table as CLM does
(`T_idx = nint((T − 223)/5) + 1`, `Tgrd_idx = nint(dT/dz/10) + 1`,
`rhos_idx = nint((max(ρ, 50) − 50)/50) + 1`, each clamped to the table). Values outside the
table take its edge nodes. There is no interpolation between nodes.
"""
@inline function snicar_aging_parameters(T::Float64, gradient::Float64, density::Float64)
    nT, nG, nR = length(SNICAR_AGING_TEMPERATURE), length(SNICAR_AGING_GRADIENT),
                 length(SNICAR_AGING_DENSITY)
    # Clamped before rounding so that a large finite argument cannot overflow `Int`; a NaN
    # still throws in `round`.
    k = _nint(clamp((T - 223) / 5, -1.0, Float64(nT))) + 1
    j = _nint(clamp(gradient / 10, -1.0, Float64(nG))) + 1
    i = _nint(clamp((max(density, 50.0) - 50) / 50, -1.0, Float64(nR))) + 1
    k = clamp(k, 1, nT)
    j = clamp(j, 1, nG)
    i = clamp(i, 1, nR)
    return SNICAR_AGING_TAU[i, j, k], SNICAR_AGING_KAPPA[i, j, k], SNICAR_AGING_DRDT0[i, j, k]
end

"""
    snicar_aging_radius(r, T, gradient, density, liquid_fraction, dt) -> r_new [µm]

Effective radius [µm] of snow of radius `r` [µm] after a step of `dt` seconds, by the dry and
wet growth of CLM5.0 (Technical Note eqs. 2.3.72, 2.3.73 and 2.3.75):

    dr_dry = (dr/dt)₀·(τ/(r − r₀ + τ))^(1/κ)·dt/3600
    dr_wet = 10¹⁸·C₂·f³·dt/(4π·r²)

with τ, κ, (dr/dt)₀ from [`snicar_aging_parameters`](@ref) at the cell's temperature `T` [K],
temperature gradient `gradient` [K m-1] and density `density` [kg m-3] (ice and liquid),
`r₀ = SNICAR_RADIUS_FRESH`, `C₂ = SNICAR_LIQUID_GROWTH` [m3 s-1], and `f` the liquid mass
fraction `liquid_fraction` capped at `SNICAR_LIQUID_FRACTION_MAX`. The two are summed: CLM
drops the dry term `C₁` of Brun (1989) because dry aging is counted separately. The result is
capped at `SNICAR_RADIUS_MAX`.

Snow finer than `r₀` (GEMB's fresh snow is `1000·RE_NEW_SNOW` = 50 µm) grows at the fresh-snow
rate `(dr/dt)₀`. CLM raises such snow to `r₀` instead.
"""
@inline function snicar_aging_radius(r::Float64, T::Float64, gradient::Float64,
    density::Float64, liquid_fraction::Float64, dt::Float64)
    τ, κ, drdt0 = snicar_aging_parameters(T, gradient, density)
    dr = drdt0 * (τ / (max(r - SNICAR_RADIUS_FRESH, 0.0) + τ))^(1 / κ) * (dt / 3600)
    f = min(liquid_fraction, SNICAR_LIQUID_FRACTION_MAX)
    dr += 1e18 * dt * SNICAR_LIQUID_GROWTH * f^3 / (4 * pi * r * r)
    return min(r + dr, SNICAR_RADIUS_MAX)
end

"""
    snicar_refrozen_radius(r, refrozen, mass) -> r_new [µm]

Effective radius [µm] of a cell of total mass `mass` [kg m-2] (ice and liquid) and radius `r`
[µm] after `refrozen` [kg m-2] of its liquid water has frozen: the mass-weighted mean of `r`
and `SNICAR_RADIUS_REFROZEN` (CLM5.0 Technical Note eq. 2.3.72), with the refrozen fraction
capped at 1.
"""
@inline function snicar_refrozen_radius(r::Float64, refrozen::Float64, mass::Float64)
    f = min(refrozen / mass, 1.0)
    return r + f * (SNICAR_RADIUS_REFROZEN - r)
end

"""
    snicar_refreeze!(grain_radius, grain_dendricity, freeze, M, water, density)

Apply [`snicar_refrozen_radius`](@ref) to every cell that [`calculate_grain_size`](@ref) ages by
the SNICAR law (dendritic, or below `DENSITY_MARBOUTY_MAX`), for the `freeze[i]` [kg m-2] of
liquid water that froze in cell `i` this step. `grain_radius` is in mm and updated in place;
`M` and `water` [kg m-2] are the cell's ice and liquid mass after freezing.
"""
function snicar_refreeze!(grain_radius::AbstractVector{Float64},
    grain_dendricity::AbstractVector{Float64}, freeze::AbstractVector{Float64},
    M::AbstractVector{Float64}, water::AbstractVector{Float64},
    density::AbstractVector{Float64})
    for i in eachindex(grain_radius, grain_dendricity, freeze, M, water, density)
        freeze[i] > 0.0 || continue
        (grain_dendricity[i] > GDN_TOLERANCE ||
         density[i] < DENSITY_MARBOUTY_MAX - D_TOLERANCE) || continue
        grain_radius[i] = snicar_refrozen_radius(1000 * grain_radius[i], freeze[i],
                                                 M[i] + water[i]) / 1000
    end
    return grain_radius
end
