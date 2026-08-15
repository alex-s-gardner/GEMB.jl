"""
    thermal_conductivity(temperature, density, mp::ModelParameters)

Compute thermal conductivity for snow/firn/ice based on density and temperature.
Matches MATLAB's `thermal_conductivity.m` for the two methods MATLAB carries.

For snow/firn (density < density_ice), by `mp.thermal_conductivity_method`:

- `:Sturm` — Sturm et al. (1997): `K = 0.138 - 1.01e-3*ρ + 3.233e-6*ρ^2`
- `:Calonne` — Calonne et al. (2011): `K = 0.024 - 1.23e-4*ρ + 2.5e-6*ρ^2`
- `:Calonne2019` — Calonne et al. (2019) eq. 5, a temperature-dependent sigmoid blend of a
  snow and a firn regime; see [`_thermal_conductivity_calonne2019`](@ref)
- `:Marchenko2019` — Marchenko et al. (2019) eq. 30, a linear firn fit floored by
  `:Calonne`; see [`_thermal_conductivity_marchenko2019`](@ref)

For ice (density >= density_ice):
- K = 9.828 * exp(-5.7e-3 * T)

`:Calonne2019` is the exception: it is continuous into ice by construction (it returns
`K_ice(T)` exactly at ρ = 917) and so is evaluated at all densities rather than
short-circuiting to the ice branch.

`:Sturm` and `:Calonne` are the MATLAB-equivalent options and either reproduces the
reference to 1e-12. `:Calonne2019` and `:Marchenko2019` are additions with no MATLAB
counterpart, both recommended over the older fits by Vandecrux et al. (2020, RetMIP
Sect. 5.1), who attribute part of a multi-model cold bias at Summit and Dye-2 to
conductivity parameterizations.

Returns vector of thermal conductivities [W m-1 K-1].
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

@inline function _conductivity_code(method::Symbol)
    method === :Sturm && return _K_STURM
    method === :Calonne && return _K_CALONNE
    method === :Calonne2019 && return _K_CALONNE2019
    method === :Marchenko2019 && return _K_MARCHENKO2019
    error("unknown thermal_conductivity_method: $(method)")
end

"""
    _thermal_conductivity_ice(T) -> K [W m-1 K-1]

Thermal conductivity of pure ice, `9.828*exp(-5.7e-3*T)` (Yen, 1981). Used both for the
ρ >= density_ice branch and as the temperature scaling inside
[`_thermal_conductivity_calonne2019`](@ref).
"""
@inline _thermal_conductivity_ice(T::Real) = 9.828 * exp(-5.7e-3 * T)

"""
    _thermal_conductivity_calonne2019(T, ρ) -> K [W m-1 K-1]

Calonne et al. (2019) eq. 5. A logistic weight `θ` in density blends two regime-specific
fits, each rescaled by the temperature-dependent ice conductivity:

    θ      = 1/(1 + exp(-2a(ρ - ρ_t))),      a = 0.02, ρ_t = 450 kg m-3
    k_snow = 0.024 - 1.23e-4ρ + 2.5e-6ρ^2   (the :Calonne 2011 fit, ρ -> 0 gives air)
    k_firn = 2.107 + 0.003618(ρ - 917)       (ρ -> 917 gives ice)
    K      = (1-θ)·(K_ice·K_air)/(k_i·k_a)·k_snow + θ·(K_ice/k_i)·k_firn

with `k_i = 2.107` and `k_a = 0.024` the reference ice and air conductivities the two fits
were built at. Because `K_air` is held at its reference value `k_a`, the snow term's
prefactor reduces to `K_ice/k_i`, the same scaling as the firn term — so this is written as
a single `K_ice/k_i` factor on the blend. At ρ = 917 the weight is 1 to within 1e-8 and
`k_firn = k_i`, so `K = K_ice(T)` to the same relative precision: the parameterization is
continuous into ice and needs no ice branch.

Unlike most density constants in GEMB this uses the literal 917, not `mp.density_ice`: it is
a fitted coefficient of the published regression, not a configurable property of the column
(the same reasoning as the `:Barnola1991` note in `initialize_parameters.jl`).

Cross-checked against the Community Firn Model's `Calonne2019` conductivity in
`CFM_main/diffusion.py`.
"""
@inline function _thermal_conductivity_calonne2019(T::Real, ρ::Real)
    a = 0.02
    ρ_transition = 450.0
    k_i = 2.107      # reference ice conductivity of the fit [W m-1 K-1]
    θ = 1.0 / (1.0 + exp(-2.0 * a * (ρ - ρ_transition)))
    k_snow = 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
    k_firn = 2.107 + 0.003618 * (ρ - 917.0)
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
    # Continuous into ice by construction, so it is not short-circuited below.
    method == _K_CALONNE2019 && return _thermal_conductivity_calonne2019(T, ρ)

    if ρ < density_ice - D_TOLERANCE
        method == _K_STURM && return 0.138 - 1.01e-3 * ρ + 3.233e-6 * ρ^2
        method == _K_CALONNE && return 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2
        return _thermal_conductivity_marchenko2019(ρ)
    end
    return _thermal_conductivity_ice(T)
end
