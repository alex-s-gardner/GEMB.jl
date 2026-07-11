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
    m = length(density)
    K = zeros(m)

    for i in 1:m
        if density[i] < mp.density_ice - d_tolerance
            # Snow/firn
            if mp.thermal_conductivity_method == "Calonne"
                K[i] = 0.024 - 1.23e-4 * density[i] + 2.5e-6 * density[i]^2
            else  # "Sturm"
                K[i] = 0.138 - 1.01e-3 * density[i] + 3.233e-6 * density[i]^2
            end
        else
            # Ice
            K[i] = 9.828 * exp(-5.7e-3 * temperature[i])
        end
    end

    return K
end
