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
    # `:LiZwally` and `:Helsen` stay gated: their accumulation units are ambiguous in the
    # literature (m i.e. vs m w.e., a factor of ~1.099), so the coefficients cannot be
    # confirmed against an independent implementation.
    @assert mp.densification_method in (:HerronLangway, :Arthern, :ArthernB, :Barnola1991, :Crocus, :CrocusPure, :GSFC2020, :Simonsen2013, :Ligtenberg) "densification_method must be one of: HerronLangway, Arthern, ArthernB, Barnola1991, Crocus, CrocusPure, GSFC2020, Simonsen2013, Ligtenberg"
    @assert mp.densification_accumulation in (:accumulation, :precipitation) "densification_accumulation must be accumulation or precipitation"
    @assert mp.mean_temperature_method in (:arithmetic, :arrhenius) "mean_temperature_method must be arithmetic or arrhenius"

    # Densification coefficients
    valid_coeffs = (:Ant_ERA5_GS_SW0, :Ant_ERA5v4_Paolo23, :Ant_ERA5_BF_SW1,
        :Ant_RACMO_GS_SW0, :Ant_Ligtenberg, :Gre_ERA5_GS_SW0,
        :Gre_RACMO_GS_SW0, :Gre_RACMO_GB_SW1, :Gre_KuipersMunneke)
    @assert mp.densification_coeffs_M01 in valid_coeffs "Invalid densification_coeffs_M01"

    # New snow method
    @assert mp.new_snow_method in (Symbol("150kgm2"), Symbol("350kgm2"), :Fausto, :FaustoFit, :Kaspers, :KuipersMunneke) "Invalid new_snow_method"

    # Density of ice
    @assert 800 <= mp.density_ice <= 950 "density_ice must be in [800, 950]"

    # `:Barnola1991` carries a narrower `density_ice` range than the rest of the model. Its
    # `f(ρ)` switches from a fitted polynomial to a closed-pore form at `BARNOLA_CLOSEOFF`
    # (800), and only the closed-pore branch scales with `density_ice`: the polynomial's
    # coefficients were fitted against the fit's own ice density of ~920 and cannot be
    # rescaled. Below ~900 the two branches diverge badly at the handover (-27% at 900, -94%
    # at 820) and the polynomial is applied to cells whose true porosity is nearly zero,
    # overstating the rate by up to ~16x; at 800 the polynomial covers the whole column and
    # `f` never vanishes, so the scheme loses the self-limiting property it exists to provide.
    if mp.densification_method === :Barnola1991
        @assert mp.density_ice >= BARNOLA_MIN_DENSITY_ICE "densification_method=:Barnola1991 requires density_ice >= $(BARNOLA_MIN_DENSITY_ICE) (its f(ρ) polynomial was fitted against an ice density of ~920 kg m-3 and cannot be rescaled); got $(mp.density_ice)"
    end

    # Rain temperature threshold
    @assert 270.15 <= mp.rain_temperature_threshold <= 276.15 "rain_temperature_threshold must be in [270.15, 276.15]"

    # Emissivity method
    @assert mp.emissivity_method in (:uniform, :grain_radius_threshold, :grain_radius_w_threshold) "Invalid emissivity_method"

    # Emissivity values
    @assert 0 <= mp.emissivity <= 1 "emissivity must be in [0, 1]"
    @assert 0 <= mp.emissivity_grain_radius_large <= 1 "emissivity_grain_radius_large must be in [0, 1]"
    @assert 0 <= mp.emissivity_grain_radius_threshold <= 100 "emissivity_grain_radius_threshold must be in [0, 100]"

    # Surface roughness
    @assert 0 <= mp.surface_roughness_effective_ratio <= 3 "surface_roughness_effective_ratio must be in [0, 3]"

    # Thermal conductivity
    @assert mp.thermal_conductivity_method in (:Sturm, :Calonne, :Calonne2019, :Calonne2019Air, :Marchenko2019) "thermal_conductivity_method must be one of: Sturm, Calonne, Calonne2019, Calonne2019Air, Marchenko2019"

    # Heat capacity
    @assert mp.heat_capacity_method in (:constant, :CuffeyPaterson) "heat_capacity_method must be constant or CuffeyPaterson"
    @assert 1500 <= mp.heat_capacity_ice <= 2500 "heat_capacity_ice must be in [1500, 2500]"
    @assert mp.rain_heat_capacity in (:water, :ice) "rain_heat_capacity must be water or ice"

    # Grain growth
    @assert mp.grain_growth_method in (:Marbouty, :Arthern, :hybrid) "grain_growth_method must be one of: Marbouty, Arthern, hybrid"

    # Water
    @assert mp.water_irreducible_method in (:constant, :ColeouLesaffre) "water_irreducible_method must be constant or ColeouLesaffre"
    @assert 0 <= mp.water_irreducible_saturation <= 0.2 "water_irreducible_saturation must be in [0, 0.2]"
    # Bounded at 917 (pure ice, and the largest value used by any RetMIP model) rather than at
    # `mp.density_ice`: a threshold above the column's ice density is inert, not invalid — a
    # cell at `density_ice` is already unconditionally impermeable in `calculate_melt`
    # regardless of run thickness. Tying the bound to `mp.density_ice` would reject the default
    # criterion whenever a caller lowers `density_ice`, a legitimate sensitivity test.
    @assert 700 <= mp.impermeable_density <= 917 "impermeable_density must be in [700, 917] (got $(mp.impermeable_density))"
    @assert mp.impermeable_thickness >= 0 "impermeable_thickness must be >= 0"

    # Runoff
    @assert mp.runoff_method in (:instantaneous, :ZuoOerlemans, :Darcy) "runoff_method must be one of: instantaneous, ZuoOerlemans, Darcy"
    @assert 0 < mp.pore_saturation_max <= 1 "pore_saturation_max must be in (0, 1]"
    # `:Darcy`'s flux is directly proportional to the slope, so a zero slope makes it drain
    # nothing at all while still allowing water to pond — an aquifer that only ever fills.
    # `:ZuoOerlemans` is well defined at zero slope (its timescale saturates at c1 + c2).
    if mp.runoff_method === :Darcy
        @assert mp.surface_slope > 0 "runoff_method=:Darcy requires surface_slope > 0; its flux is proportional to slope, so a zero slope drains nothing while still permitting ponding"
    end

    # Albedo method
    @assert mp.albedo_method in (:None, :GardnerSharp, :BrunLefebre, :GreuellKonzelmann) "Invalid albedo_method"

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

    # Output
    @assert mp.output_frequency in (:all, :monthly, :weekly, :daily, :last) "Invalid output_frequency"
    @assert mp.initialize_age in (:zero, :steady_state) "initialize_age must be zero or steady_state"

    # Grid geometry
    @assert 0 <= mp.column_ztop <= 100
    @assert 0 <= mp.column_dztop <= 0.2
    @assert 0 <= mp.column_dzmin <= 0.2
    @assert 0 <= mp.column_dzmax <= 0.2
    @assert 0 <= mp.column_depth_max <= 1000
    @assert 1 <= mp.column_zy <= 2

    # Ice dynamics. The bound is orders of magnitude above any real ice-sheet strain rate
    # (fast shear margins reach ~1e-2 yr-1) and exists to catch a units mistake — an s-1
    # value, or a percent — rather than to constrain the physics.
    @assert -1 <= mp.horizontal_strain_rate <= 1 "horizontal_strain_rate must be in [-1, 1] yr-1"

    # Slope as a gradient, not degrees. The bound is far above any ice-sheet surface slope
    # (the RetMIP sites reach ~0.010 m m-1) and is there to catch degrees passed by mistake.
    @assert 0 <= mp.surface_slope <= 1 "surface_slope must be in [0, 1] m m-1 (a gradient, not degrees)"
end
