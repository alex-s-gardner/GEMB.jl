# Thermal properties of the soil column of `lower_boundary = :ground`.
#
# A ground cell holds mineral solids (volume fraction 1 − φ) and a fixed mass of pore water
# `W = θ·ρ_w·dz`, liquid or frozen. Its enthalpy counts the solids at `c_s·T` (0 K reference) and the
# pore water as ice below the melting point (GEMB's `specific_enthalpy`), so the pore water
# freezes on a latent plateau at 273.15 K exactly as a wet snow cell's water does in
# `_wet_cell_temperature`.

# Conductivities of ice and liquid water in the saturated-soil geometric mean [W m-1 K-1]:
# ISBA `XCONDI` (`ini_csts.F90:155`) and `XCONDWTR` (`modd_isba_par.F90:94`).
const GROUND_CONDUCTIVITY_ICE = 2.22
const GROUND_CONDUCTIVITY_WATER = 0.57
# Least degree of saturation in the Kersten number (ISBA `soil.F90:251`), which keeps log10 finite.
const GROUND_SATURATION_MIN = 0.1

# Depths of the lower faces of the ground cells [m]. Below 0.1 m these are the ISBA-DIF 20-layer
# grid (`XSOILGRID` of the Crocus reference namelist, `OPTIONS_V9_reference.nam:210`); its four
# cells above 0.1 m (0.01, 0.03, 0.06, 0.1) are two 5 cm cells here, whose centres bracket the
# 5 cm soil-temperature depth and whose heat capacity bounds the explicit sub-step. Below 10 m the
# grid continues in 2 m cells.
const GROUND_LAYER_BOTTOMS = (0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0,
    4.0, 5.0, 6.5, 8.0, 10.0)
const GROUND_LAYER_DEEP_DZ = 2.0

# Depths below the ground surface of the soil-temperature outputs [m], the SNOTEL soil sensors
# (2, 4, 8, 20 and 40 in.).
const GROUND_OUTPUT_DEPTHS = (0.05, 0.1, 0.2, 0.5, 1.0)

# Material codes of the `material` profile layer of a `lower_boundary = :ground` column.
const MATERIAL_SNOW = 0.0     # snow, firn or ice above the ground
const MATERIAL_GROUND = 1.0   # soil
const MATERIAL_EMPTY = 2.0    # an unused snow slot: zero thickness, mass and enthalpy

"""
    ground_grid(depth) -> dz [m]

Thicknesses of the ground cells of a soil column `depth` metres deep, surface first: the faces of
`GROUND_LAYER_BOTTOMS` above `depth`, then 2 m cells, the last cell ending at `depth`.
"""
function ground_grid(depth::Real)
    bottoms = Float64[]
    for b in GROUND_LAYER_BOTTOMS
        b < depth - D_TOLERANCE && push!(bottoms, b)
    end
    b = isempty(bottoms) ? 0.0 : bottoms[end]
    while b + GROUND_LAYER_DEEP_DZ < depth - D_TOLERANCE
        b = max(b, GROUND_LAYER_BOTTOMS[end]) + GROUND_LAYER_DEEP_DZ
        push!(bottoms, b)
    end
    push!(bottoms, Float64(depth))
    dz = diff(vcat(0.0, bottoms))
    length(dz) >= 2 || error("ground_grid: a ground column of $(depth) m has fewer than 2 cells")
    return dz
end

"""
    GroundColumn

The fixed soil column of a `lower_boundary = :ground` run, surface first: cell thickness `dz`
[m], solids heat capacity per unit area `C_solid = (1 − φ)·ρ_s·c_s·dz` [J m-2 K-1], pore-water
mass `water_total = θ·ρ_w·dz` [kg m-2], and the solids and dry-soil conductivities shared by
every cell [W m-1 K-1]. Built by [`ground_column`](@ref) from `ModelParameters` and `dz`; holds
no state that changes during a run.
"""
struct GroundColumn
    dz::Vector{Float64}
    C_solid::Vector{Float64}
    water_total::Vector{Float64}
    density_bulk::Float64       # solids and all pore water per unit volume [kg m-3]
    conductivity_solid::Float64
    conductivity_dry::Float64
    saturation::Float64         # degree of saturation θ/φ
    porosity::Float64
end

"""
    ground_dry_conductivity(porosity, solid_density) -> k_dry [W m-1 K-1]

Dry-soil conductivity of Johansen (1975) in the form ISBA implements
(`thrmcondz.F90:110-111`, after Peters-Lidard et al., 1998):
`k_dry = (0.135·γ_d + 64.7)/(ρ_s − 0.947·γ_d)`, `γ_d = (1 − φ)·ρ_s` the dry bulk density.

# References
- Johansen, O. (1975). *Thermal conductivity of soils*. PhD thesis, University of Trondheim.
  (CRREL Draft Translation 637, 1977.)
- Peters-Lidard, C. D., Blackburn, E., Liang, X., and Wood, E. F. (1998). The effect of soil
  thermal conductivity parameterization on surface energy fluxes and temperatures.
  *J. Atmos. Sci.* 55, 1209-1224.
"""
@inline function ground_dry_conductivity(porosity::Real, solid_density::Real)
    γd = (1 - porosity) * solid_density
    return (0.135 * γd + 64.7) / (solid_density - 0.947 * γd)
end

"""
    ground_column(mp, dz) -> GroundColumn

The soil column of `dz` [m] under `mp`'s `ground_*` parameters.
"""
function ground_column(mp::ModelParameters, dz::AbstractVector{<:Real})
    φ = mp.ground_porosity
    θ = mp.ground_water
    dzv = Vector{Float64}(dz)
    C_solid = [(1 - φ) * mp.ground_solid_density * mp.ground_solid_heat_capacity * d for d in dzv]
    water_total = [θ * DENSITY_WATER * d for d in dzv]
    density_bulk = (1 - φ) * mp.ground_solid_density + θ * DENSITY_WATER
    return GroundColumn(dzv, C_solid, water_total, density_bulk, mp.ground_solid_conductivity,
        ground_dry_conductivity(φ, mp.ground_solid_density), θ / φ, φ)
end

"""
    ground_conductivity(g, liquid, total) -> k [W m-1 K-1]

Conductivity of a soil cell holding `liquid` of its `total` pore water [kg m-2] unfrozen, after
Johansen (1975) as ISBA implements it (`soil.F90:236-266`):

    k     = Ke·(k_sat − k_dry) + k_dry
    k_sat = k_s^(1−φ) · k_i^(φ−u) · k_w^u,      u = (1 − f)·φ
    Ke    = (1 − f)·(log10(Sr) + 1) + f·Sr,      Sr = max(0.1, θ/φ)

with `f` the frozen fraction of the pore water, `k_s` the solids conductivity,
`k_i = 2.22` and `k_w = 0.57`. The Kersten number `Ke` is the unfrozen form `log10(Sr) + 1` of
Johansen for fine soil, blended linearly in `f` with the frozen form `Sr`, as ISBA does
("a smooth transition from thawed to frozen soils"). A cell without pore water is dry soil.

# References
- Johansen, O. (1975). *Thermal conductivity of soils*. PhD thesis, University of Trondheim.
- Peters-Lidard, C. D., et al. (1998). *J. Atmos. Sci.* 55, 1209-1224.
- ISBA-DIF in SURFEX v8 (`soil.F90`, `thrmcondz.F90`), as distributed with Crocus 3.0.2.
"""
@inline function ground_conductivity(g::GroundColumn, liquid::Float64, total::Float64)
    total > 0.0 || return g.conductivity_dry
    f = clamp(1.0 - liquid / total, 0.0, 1.0)
    φ = g.porosity
    u = (1.0 - f) * φ
    k_sat = g.conductivity_solid^(1.0 - φ) * GROUND_CONDUCTIVITY_ICE^(φ - u) *
            GROUND_CONDUCTIVITY_WATER^u
    Sr = max(GROUND_SATURATION_MIN, g.saturation)
    Ke = (1.0 - f) * (log10(Sr) + 1.0) + f * Sr
    return Ke * (k_sat - g.conductivity_dry) + g.conductivity_dry
end

# Enthalpy bounds of the freezing plateau of a ground node of solids heat capacity `Cs`
# [J m-2 K-1] and pore water `W` [kg m-2]: all of the water frozen at the melting point (`H_lo`)
# and none of it (`H_hi = H_lo + W·LF`) [J m-2].
@inline function _ground_plateau(mp::ModelParameters, Cs::Float64, W::Float64)
    H_lo = Cs * CtoK + W * specific_enthalpy(mp, CtoK)
    return H_lo, H_lo + W * LF
end

"""
    ground_enthalpy(mp, Cs, W, T, liquid) -> H [J m-2]

Enthalpy of a ground node of solids heat capacity `Cs` [J m-2 K-1] and pore water `W` [kg m-2]
at temperature `T` [K] with `liquid` [kg m-2] of that water unfrozen: the solids at `Cs·T`; the
pore water at the ice enthalpy `h(T)` below the melting point, at `h(273.15) + LF·(liquid/W)`
per kg on the plateau, and as liquid water `h(273.15) + LF + c_w·(T − 273.15)` above it, with
`c_w` = `HEAT_CAPACITY_WATER`. `liquid` is read only at `T = 273.15 K`; elsewhere the phase
follows from `T`.

The volumetric heat capacity off the plateau is that of the design,
`(1 − φ)ρ_s c_s + θ_l ρ_w c_w + θ_i ρ_i c_i`, with `c_i` GEMB's ice heat capacity.
"""
@inline function ground_enthalpy(mp::ModelParameters, Cs::Float64, W::Float64, T::Float64,
    liquid::Float64)
    if T < CtoK
        return Cs * T + W * specific_enthalpy(mp, T)
    elseif T > CtoK
        return Cs * T + W * (specific_enthalpy_water(mp) + HEAT_CAPACITY_WATER * (T - CtoK))
    end
    H_lo, _ = _ground_plateau(mp, Cs, W)
    return H_lo + clamp(liquid, 0.0, W) * LF
end

"""
    ground_state(mp, Cs, W, H) -> (T [K], liquid [kg m-2])

Inverse of [`ground_enthalpy`](@ref): the temperature and unfrozen pore water of a ground node
with enthalpy `H`. Above the plateau the map is linear; below it, with `h = aT + (b/2)T²`
(`:CuffeyPaterson`), `Cs·T + W·h(T) = H` is solved by the cancellation-free root
`T = 2H/(B + √(B² + 2bWH))`, `B = Cs + aW`.
"""
@inline function ground_state(mp::ModelParameters, Cs::Float64, W::Float64, H::Float64)
    H_lo, H_hi = _ground_plateau(mp, Cs, W)
    if H >= H_hi
        return CtoK + (H - H_hi) / (Cs + W * HEAT_CAPACITY_WATER), W
    elseif H >= H_lo
        return CtoK, (H - H_lo) / LF
    end
    if mp.heat_capacity_method === :constant
        return H / (Cs + W * mp.heat_capacity_ice), 0.0
    end
    B = Cs + HEAT_CAPACITY_CUFFEY_A * W
    return 2 * H / (B + sqrt(B * B + 2 * HEAT_CAPACITY_CUFFEY_B * W * H)), 0.0
end

"""
    ground_heat_capacity(mp, Cs, W, T) -> C [J m-2 K-1]

Sensible heat capacity of a ground node at `T` off the plateau: the solids plus the pore water
as ice below the melting point, as liquid at or above it.
"""
@inline ground_heat_capacity(mp::ModelParameters, Cs::Float64, W::Float64, T::Float64) =
    Cs + W * (T < CtoK ? heat_capacity(mp, T) : HEAT_CAPACITY_WATER)
