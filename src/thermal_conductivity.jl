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
    d_tolerance = 1e-11

    # Vectorized version using broadcasting and ifelse
    snow_mask = density .< (mp.density_ice - d_tolerance)

    if mp.thermal_conductivity_method == :Calonne
        K = @. ifelse(snow_mask,
                      0.024 - 1.23e-4 * density + 2.5e-6 * density^2,
                      9.828 * exp(-5.7e-3 * temperature))
    else  # "Sturm"
        K = @. ifelse(snow_mask,
                      0.138 - 1.01e-3 * density + 3.233e-6 * density^2,
                      9.828 * exp(-5.7e-3 * temperature))
    end

    return K
end
