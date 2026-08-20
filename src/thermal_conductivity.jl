"""
    thermal_conductivity(temperature, density, mp::ModelParameters)

Compute thermal conductivity for snow/firn/ice based on density and temperature.
Matches MATLAB's `thermal_conductivity.m` for the two methods MATLAB carries.

For snow/firn (density < density_ice), by `mp.thermal_conductivity_method`:

- `:Sturm` — Sturm et al. (1997): `K = 0.138 - 1.01e-3*ρ + 3.233e-6*ρ^2`
- `:Calonne` — Calonne et al. (2011): `K = 0.024 - 1.23e-4*ρ + 2.5e-6*ρ^2`
- `:Calonne2019` — Calonne et al. (2019) eq. 5, a temperature-dependent sigmoid blend of a
  snow and a firn regime, with `K_air` pinned at its reference value (matching the
  Community Firn Model); see [`_thermal_conductivity_calonne2019`](@ref)
- `:Calonne2019Air` — the same equation carrying the `K_air(T)/K_air(T_ref)` factor on its
  snow branch as well (matching IMAU-FDM). Gives lower conductivity than `:Calonne2019`
  for cold low-density snow — up to ~20% at 220 K — and is identical to it at and above
  ρ ≈ 550 and at the 270.15 K reference temperature
- `:Marchenko2019` — Marchenko et al. (2019) eq. 30, a linear firn fit floored by
  `:Calonne`; see [`_thermal_conductivity_marchenko2019`](@ref)

For ice (density >= density_ice):
- K = 9.828 * exp(-5.7e-3 * T)

The two 2019 forms are the exception: each is continuous into ice by construction (both
return `K_ice(T)` at ρ = 917) and so is evaluated at all densities rather than
short-circuiting to the ice branch.

`:Sturm` and `:Calonne` are the MATLAB-equivalent options and either reproduces the
reference to 1e-12. The `:Calonne2019` pair and `:Marchenko2019` are additions with no
MATLAB counterpart, all recommended over the older fits by Vandecrux et al. (2020, RetMIP
Sect. 5.1), who attribute part of a multi-model cold bias at Summit and Dye-2 to
conductivity parameterizations.

Returns vector of thermal conductivities [W m-1 K-1].

# References
- Sturm, M., Holmgren, J., König, M., and Morris, K. (1997). The thermal conductivity of
  seasonal snow. *J. Glaciol.* 43, 26-41.
- Calonne, N., Flin, F., Morin, S., Lesaffre, B., Rolland du Roscoat, S., and Geindreau, C.
  (2011). Numerical and experimental investigations of the effective thermal conductivity of
  snow. *Geophys. Res. Lett.* 38, L23501.
- Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F., and
  Geindreau, C. (2019). Thermal conductivity of snow, firn, and porous ice from 3-D
  image-based computations. *Geophys. Res. Lett.* 46, 13079-13089.
- Marchenko, S., Cheng, G., Lötstedt, P., Pohjola, V., Pettersson, R., van Pelt, W., and
  Reijmer, C. (2019). Thermal conductivity of firn at Lomonosovfonna, Svalbard, derived from
  subsurface temperature measurements. *The Cryosphere* 13, 1843-1859.
- Yen, Y.-C. (1981). *Review of thermal properties of snow, ice and sea ice*. CRREL Report
  81-10, U.S. Army Cold Regions Research and Engineering Laboratory. (Ice branch.)
- Reid, R. C., Prausnitz, J. M., and Sherwood, T. K. (1966). *The Properties of Gases and
  Liquids*. McGraw-Hill. (Air conductivity, `:Calonne2019Air` only.)
- Vandecrux, B., et al. (2020). The firn meltwater Retention Model Intercomparison Project
  (RetMIP). *The Cryosphere* 14, 3785-3810. Sects. 5.1 and 7.
"""
function thermal_conductivity(temperature::AbstractVector, density::AbstractVector, mp::ModelParameters)
    K = Vector{Float64}(undef, length(density))
    method = _conductivity_code(mp.thermal_conductivity_method)
    @inbounds for i in eachindex(density)
        K[i] = _thermal_conductivity_scalar(temperature[i], density[i], mp.density_ice, method)
    end
    return K
end

"""
    thermal_conductivity(temperature::Real, density::Real, mp) -> K [W m-1 K-1]

Scalar method, for callers with a single (temperature, density) pair — the
steady-state initial guess marches one parcel and would otherwise allocate a
one-element array per step.

The vector method above loops over this one, so the branch and the calibrated
coefficients have a single definition.
"""
function thermal_conductivity(temperature::Real, density::Real, mp::ModelParameters)
    return _thermal_conductivity_scalar(Float64(temperature), Float64(density),
        mp.density_ice, _conductivity_code(mp.thermal_conductivity_method))
end

# The method is resolved to a small integer once, outside the vector method's loop, rather
# than compared as a `Symbol` per cell. (This was a `Bool` when there were only two methods.)
const _K_STURM = 1
const _K_CALONNE = 2
const _K_CALONNE2019 = 3
const _K_MARCHENKO2019 = 4
const _K_CALONNE2019AIR = 5

@inline function _conductivity_code(method::Symbol)
    method === :Sturm && return _K_STURM
    method === :Calonne && return _K_CALONNE
    method === :Calonne2019 && return _K_CALONNE2019
    method === :Marchenko2019 && return _K_MARCHENKO2019
    method === :Calonne2019Air && return _K_CALONNE2019AIR
    error("unknown thermal_conductivity_method: $(method)")
end

"""
    _thermal_conductivity_ice(T) -> K [W m-1 K-1]

Thermal conductivity of pure ice, `9.828*exp(-5.7e-3*T)` (Yen, 1981). Used both for the
ρ >= `density_ice` branch and as the temperature scaling inside the
[`_thermal_conductivity_calonne2019`](@ref) blend.
"""
@inline _thermal_conductivity_ice(T::Real) = 9.828 * exp(-5.7e-3 * T)

"""
    _thermal_conductivity_air(T) -> K [W m-1 K-1]

Thermal conductivity of air, `a*T^1.5/(b + T)` (Reid, 1966), with coefficients
`K_AIR_REID_A`/`K_AIR_REID_B`. Read only as the ratio `K_air(T)/K_air(T_ref)` inside
[`_thermal_conductivity_calonne2019`](@ref)'s `:Calonne2019Air` form, so its absolute
calibration never enters the result — see the note there.

`T^1.5` is written `T*sqrt(T)`, which the hardware does in one instruction where the
generic `pow` takes ~15x longer — measured at 14.9 ns vs 1.0 ns per element, which was
the whole of a 3.5x slowdown on this function when the air form was added. The two agree
exactly for three quarters of the doubles in 150-320 K and differ by 1 ulp (~2e-16
relative) for the rest, twelve orders of magnitude below the fit's own accuracy.
"""
@inline _thermal_conductivity_air(T::Real) = K_AIR_REID_A * (T * sqrt(T)) / (K_AIR_REID_B + T)

# The reference-temperature denominator of the air ratio, hoisted to a `const` so the
# `sqrt` and two flops behind it are not recomputed per cell per timestep. Kept as a
# denominator rather than a precomputed reciprocal: multiplying by `1/x` rounds differently
# from dividing by `x`, and the published equation is a ratio.
const K_AIR_REF = K_AIR_REID_A * (CONDUCTIVITY_T_REF * sqrt(CONDUCTIVITY_T_REF)) /
                  (K_AIR_REID_B + CONDUCTIVITY_T_REF)

"""
    _thermal_conductivity_calonne2019(T, ρ; air_factor::Bool) -> K [W m-1 K-1]

Calonne et al. (2019) eq. 5. A logistic weight `θ` in density blends two regime-specific
fits, each rescaled by the constituent conductivities at `T` relative to their values at
the temperature the fits were calibrated at (`CONDUCTIVITY_T_REF` = 270.15 K):

    θ      = 1/(1 + exp(-2a(ρ - ρ_t))),      a = 0.02, ρ_t = 450 kg m-3
    k_snow = 0.024 - 1.23e-4ρ + 2.5e-6ρ^2   (the :Calonne 2011 fit, ρ -> 0 gives air)
    k_firn = 2.107 + 0.003618(ρ - 917)       (ρ -> 917 gives ice)
    K      = (1-θ)·(K_ice/k_i)·(K_air/k_a)·k_snow + θ·(K_ice/k_i)·k_firn

The two regimes carry different scalings because they are physically different: the firn
fit describes a connected ice skeleton and so scales with ice alone, while the snow fit
describes grains in air and scales with both constituents. The reference values are
internally consistent with that reading — `k_i = 2.107` is `_thermal_conductivity_ice`
at 270.15 K to within 1.2e-4 relative, and `k_a = 0.024` is `_thermal_conductivity_air`
there to within 6.6e-3 — so the ratios are 1 at the reference temperature by construction
and the fits return their published values.

`air_factor` selects which of the two published readings of the snow term is used:

- `false` (`:Calonne2019`) pins `K_air` at its reference value, so the snow prefactor
  collapses to `K_ice/k_i` and the whole blend takes a single factor. This matches the
  Community Firn Model, whose `diffusion.py` hardcodes `K_air = kref_a` with a standing
  TODO to find the temperature dependence, and it is the form GEMB has always used.
- `true` (`:Calonne2019Air`) carries `K_air(T)/K_air(T_ref)` as well. This matches
  IMAU-FDM (`firn_physics.f90`, `Thermal_Cond`). Since it enters only as a ratio, Reid's
  absolute calibration cancels and only its temperature *shape* matters — which is why the
  0.7% gap between `K_air(T_ref)` and `k_a` is irrelevant here.

Air conducts more poorly as it cools, so the ratio is below 1 below the reference
temperature (0.83 at 220 K, 0.88 at 233 K, 0.94 at 253 K) and `:Calonne2019` therefore
returns a *higher* conductivity than `:Calonne2019Air` for cold snow — by up to ~20% at
220 K for ρ ≲ 300, ~10% at ρ = 450, and under 0.3% by ρ = 550 where `θ → 1` and the
air-free firn branch takes over. The default is `:Calonne2019`; `:Calonne2019Air` is the
fuller form of the published equation.

Both forms are continuous into ice: at ρ = 917 the weight is 1 to within 1e-8 and
`k_firn = k_i`, so `K = K_ice(T)` to the same relative precision and no ice branch is
needed. The air factor cannot disturb this because it is weighted by `(1-θ)`.

Unlike most density constants in GEMB this uses the literal 917, not `mp.density_ice`: it is
a fitted coefficient of the published regression, not a configurable property of the column
(the same reasoning as the `:Barnola1991` note in `initialize_parameters.jl`).

Cross-checked against the Community Firn Model's `Calonne2019` conductivity in
`CFM_main/diffusion.py` and IMAU-FDM's `Thermal_Cond` in `source/firn_physics.f90`.
"""
@inline function _thermal_conductivity_calonne2019(T::Real, ρ::Real; air_factor::Bool)
    a = 0.02
    ρ_transition = 450.0
    k_i = 2.107      # reference ice conductivity of the fit [W m-1 K-1]
    θ = 1.0 / (1.0 + exp(-2.0 * a * (ρ - ρ_transition)))
    k_snow = 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
    k_firn = 2.107 + 0.003618 * (ρ - 917.0)
    # `air_factor` is a compile-time constant at both call sites, so this branch is folded
    # away and neither form pays for the other's arithmetic.
    if air_factor
        k_snow *= _thermal_conductivity_air(T) / K_AIR_REF
    end
    return (_thermal_conductivity_ice(T) / k_i) * ((1.0 - θ) * k_snow + θ * k_firn)
end

"""
    _thermal_conductivity_marchenko2019(ρ) -> K [W m-1 K-1]

Marchenko et al. (2019) eq. 30, `K = 0.301e-2ρ - 0.724`, calibrated by fitting modelled to
observed subsurface temperatures at Lomonosovfonna, Svalbard, over ρ ≈ 350–900 kg m-3. It
gives higher firn conductivity than `:Sturm` or `:Calonne`, which is the direction the
multi-model cold bias in RetMIP implies.

The fit is linear and crosses zero at ρ ≈ 241, so it cannot be used at low density. It is
floored by the `:Calonne` (2011) snow fit, which the two curves cross near ρ ≈ 321 — just
below the calibration range — so the floor takes over continuously and exactly where the
regression stops being supported by data. Temperature-independent, so the ice branch
(`_thermal_conductivity_ice`) still applies at ρ >= density_ice.
"""
@inline function _thermal_conductivity_marchenko2019(ρ::Real)
    k_marchenko = 0.301e-2 * ρ - 0.724
    k_calonne = 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
    return max(k_marchenko, k_calonne)
end

@inline function _thermal_conductivity_scalar(T::Real, ρ::Real, density_ice::Float64,
    method::Int)
    # Both 2019 forms are continuous into ice by construction, so neither is short-circuited
    # to the ice branch below.
    method == _K_CALONNE2019 && return _thermal_conductivity_calonne2019(T, ρ; air_factor=false)
    method == _K_CALONNE2019AIR && return _thermal_conductivity_calonne2019(T, ρ; air_factor=true)

    if ρ < density_ice - D_TOLERANCE
        method == _K_STURM && return 0.138 - 1.01e-3 * ρ + 3.233e-6 * ρ^2
        method == _K_CALONNE && return 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
        return _thermal_conductivity_marchenko2019(ρ)
    end
    return _thermal_conductivity_ice(T)
end
