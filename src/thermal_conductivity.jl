"""
    thermal_conductivity(temperature, density, mp::ModelParameters)

Compute thermal conductivity for snow/firn/ice based on density and temperature.
Matches MATLAB's `thermal_conductivity.m`.

For snow/firn (density < density_ice):
- "Sturm": K = 0.138 - 1.01e-3*d + 3.233e-6*d^2
- "Calonne": K = 0.024 - 1.23e-4*d + 2.5e-6*d^2

For ice (density >= density_ice):
- K = 9.828 * exp(-5.7e-3 * T)

Returns vector of thermal conductivities [W m-1 K-1].
"""
function thermal_conductivity(temperature::AbstractVector, density::AbstractVector, mp::ModelParameters)
    K = Vector{Float64}(undef, length(density))
    calonne = mp.thermal_conductivity_method == :Calonne
    @inbounds for i in eachindex(density)
        K[i] = _thermal_conductivity_scalar(temperature[i], density[i], mp.density_ice, calonne)
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
        mp.density_ice, mp.thermal_conductivity_method == :Calonne)
end

# `calonne` is passed as a `Bool` rather than re-read from `mp` so the vector
# method hoists the symbol comparison out of its loop.
@inline function _thermal_conductivity_scalar(T::Real, ρ::Real, density_ice::Float64,
    calonne::Bool)
    d_tolerance = 1e-11
    if ρ < density_ice - d_tolerance
        return calonne ? 0.024 - 1.23e-4 * ρ + 2.5e-6 * ρ^2 :
                         0.138 - 1.01e-3 * ρ + 3.233e-6 * ρ^2
    end
    return 9.828 * exp(-5.7e-3 * T)
end
