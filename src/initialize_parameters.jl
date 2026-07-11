"""
    initialize_parameters(; kwargs...)

Create and validate a `ModelParameters` struct.
Matches MATLAB's `model_initialize_parameters.m`.

All parameters have defaults matching the MATLAB version.
Validation checks are performed after construction.
"""
function initialize_parameters(; kwargs...)
    mp = ModelParameters(; kwargs...)
    validate_parameters(mp)
    return mp
end

function validate_parameters(mp::ModelParameters)
    # Densification method
    @assert mp.densification_method in ("HerronLangway", "Arthern", "Ligtenberg") "densification_method must be one of: HerronLangway, Arthern, Ligtenberg"

    # Densification coefficients
    valid_coeffs = ("Ant_ERA5_GS_SW0", "Ant_ERA5v4_Paolo23", "Ant_ERA5_BF_SW1",
        "Ant_RACMO_GS_SW0", "Ant_Ligtenberg", "Gre_ERA5_GS_SW0",
        "Gre_RACMO_GS_SW0", "Gre_RACMO_GB_SW1", "Gre_KuipersMunneke")
    @assert mp.densification_coeffs_M01 in valid_coeffs "Invalid densification_coeffs_M01"

    # New snow method
    @assert mp.new_snow_method in ("150kgm2", "350kgm2", "Fausto", "Kaspers", "KuipersMunneke") "Invalid new_snow_method"

    # Density of ice
    @assert 800 <= mp.density_ice <= 950 "density_ice must be in [800, 950]"

    # Rain temperature threshold
    @assert 270.15 <= mp.rain_temperature_threshold <= 276.15 "rain_temperature_threshold must be in [270.15, 276.15]"

    # Emissivity method
    @assert mp.emissivity_method in ("uniform", "grain_radius_threshold", "grain_radius_w_threshold") "Invalid emissivity_method"

    # Emissivity values
    @assert 0 <= mp.emissivity <= 1 "emissivity must be in [0, 1]"
    @assert 0 <= mp.emissivity_grain_radius_large <= 1 "emissivity_grain_radius_large must be in [0, 1]"
    @assert 0 <= mp.emissivity_grain_radius_threshold <= 100 "emissivity_grain_radius_threshold must be in [0, 100]"

    # Surface roughness
    @assert 0 <= mp.surface_roughness_effective_ratio <= 3 "surface_roughness_effective_ratio must be in [0, 3]"

    # Thermal conductivity
    @assert mp.thermal_conductivity_method in ("Sturm", "Calonne") "thermal_conductivity_method must be Sturm or Calonne"

    # Water
    @assert 0 <= mp.water_irreducible_saturation <= 0.2 "water_irreducible_saturation must be in [0, 0.2]"

    # Albedo method
    @assert mp.albedo_method in ("None", "GardnerSharp", "BrunLefebre", "GreuellKonzelmann", "Bougamont2005") "Invalid albedo_method"

    # Albedo values
    @assert mp.albedo_density_threshold >= 0 "albedo_density_threshold must be >= 0"
    @assert 0.5 <= mp.albedo_snow <= 0.95 "albedo_snow must be in [0.5, 0.95]"
    @assert 0.2 <= mp.albedo_ice <= 0.6 "albedo_ice must be in [0.2, 0.6]"
    @assert 0.2 <= mp.albedo_fixed <= 0.95 "albedo_fixed must be in [0.2, 0.95]"

    # Radiation parameters
    @assert 0 <= mp.shortwave_downward_diffuse <= 1000
    @assert 0 <= mp.solar_zenith_angle <= 90
    @assert 0 <= mp.cloud_optical_thickness <= 30
    @assert 0 <= mp.black_carbon_snow <= 2
    @assert 0 <= mp.black_carbon_ice <= 2
    @assert 0 <= mp.cloud_fraction <= 1

    # Bougamont2005 parameters
    @assert 5 <= mp.albedo_wet_snow_t0 <= 25
    @assert 20 <= mp.albedo_dry_snow_t0 <= 40
    @assert 2 <= mp.albedo_K <= 12

    # Output
    @assert mp.output_frequency in ("all", "monthly", "daily", "last") "Invalid output_frequency"
    @assert 0 <= mp.output_padding <= 10000

    # Grid geometry
    @assert 0 <= mp.column_ztop <= 100
    @assert 0 <= mp.column_dztop <= 0.2
    @assert 0 <= mp.column_dzmin <= 0.2
    @assert 0 <= mp.column_dzmax <= 0.2
    @assert 0 <= mp.column_zmax <= 1000
    @assert 0 <= mp.column_zmin <= 1000
    @assert 1 <= mp.column_zy <= 2
end
