"""
    initialize_parameters(; kwargs...)

Create and validate a `ModelParameters` struct.

Defaults are documented per field on [`ModelParameters`](@ref); several were changed in
v2.0.0 on the recommendation of the firn-model intercomparisons, so they are not the MATLAB
model's. Validation checks are performed after construction.
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

    # Destructive metamorphism is a strain rate summed with the viscous one, so it exists only
    # under the two schemes that form a strain rate from a viscosity.
    if mp.destructive_metamorphism
        @assert mp.densification_method in (:Crocus, :CrocusPure) "destructive_metamorphism=true requires densification_method :Crocus or :CrocusPure (it is added to their viscous strain rate); got :$(mp.densification_method)"
    end
    @assert 50 <= mp.destructive_metamorphism_density <= 450 "destructive_metamorphism_density must be in [50, 450] kg m-3 (got $(mp.destructive_metamorphism_density))"

    # The Crocus handoff. Its settings are read only under `:Crocus`, so a non-default value
    # with any other scheme would be silently ignored.
    @assert mp.crocus_handoff_law in (:GSFC2020, :Arthern) "crocus_handoff_law must be GSFC2020 or Arthern (got :$(mp.crocus_handoff_law))"
    @assert mp.crocus_handoff_rule in (:density, :rate_match) "crocus_handoff_rule must be density or rate_match (got :$(mp.crocus_handoff_rule))"
    @assert 200 <= mp.crocus_handoff_density <= DENSITY_PORE_CLOSEOFF "crocus_handoff_density must be in [200, $(DENSITY_PORE_CLOSEOFF)] kg m-3 (got $(mp.crocus_handoff_density))"
    @assert 0 < mp.crocus_handoff_k0 <= 10 && 0 < mp.crocus_handoff_k1 <= 10 "crocus_handoff_k0 and crocus_handoff_k1 must be in (0, 10] (got $(mp.crocus_handoff_k0), $(mp.crocus_handoff_k1))"
    @assert abs(mp.crocus_handoff_activation_energy_shift) <= 60000 "crocus_handoff_activation_energy_shift must be in [-60000, 60000] J mol-1 (got $(mp.crocus_handoff_activation_energy_shift))"
    @assert 0 < mp.crocus_handoff_accumulation_exponent_scale <= 3 "crocus_handoff_accumulation_exponent_scale must be in (0, 3] (got $(mp.crocus_handoff_accumulation_exponent_scale))"
    @assert 400 <= mp.crocus_handoff_stage_density < DENSITY_PORE_CLOSEOFF "crocus_handoff_stage_density must be in [400, $(DENSITY_PORE_CLOSEOFF)) kg m-3 (got $(mp.crocus_handoff_stage_density))"
    @assert 0 < mp.crocus_overburden_factor <= 10 "crocus_overburden_factor must be in (0, 10] (got $(mp.crocus_overburden_factor))"
    if mp.densification_method !== :Crocus
        changed = filter(f -> getfield(mp, f) != getfield(CROCUS_FIRN_CALIBRATION, f),
                         keys(CROCUS_FIRN_CALIBRATION))
        @assert isempty(changed) "$(join(changed, ", ")) are read only by densification_method = :Crocus; got :$(mp.densification_method)"
    end

    # Densification coefficients
    valid_coeffs = (:Ant_ERA5_GS_SW0, :Ant_ERA5v4_Paolo23, :Ant_ERA5_BF_SW1,
        :Ant_RACMO_GS_SW0, :Ant_Ligtenberg, :Gre_ERA5_GS_SW0,
        :Gre_RACMO_GS_SW0, :Gre_RACMO_GB_SW1, :Gre_KuipersMunneke)
    @assert mp.densification_coeffs_M01 in valid_coeffs "Invalid densification_coeffs_M01"

    # New snow method
    @assert mp.new_snow_method in (:Constant, :Constant150, :Constant315, :Constant350, :Fausto, :FaustoFit, :Pahaut, :Kaspers, :KuipersMunneke, :vanKampenhout2017, :vanKampenhout2017SNOTEL) "Invalid new_snow_method"
    @assert NEW_SNOW_DENSITY_MIN <= mp.new_snow_density <= NEW_SNOW_DENSITY_MAX "new_snow_density must be in [$(NEW_SNOW_DENSITY_MIN), $(NEW_SNOW_DENSITY_MAX)] kg m-3 (got $(mp.new_snow_density))"

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

    # Explicit thermal sub-step safety factor. Bounded in (0, 1]: it scales a limit that is
    # already an upper bound, so exceeding 1 would run the explicit scheme knowingly unstable,
    # and 0 would collapse `dt` to the smallest divisor available. 1.0 is admitted — it runs at
    # the start-of-step limit exactly, leaving no margin for the surface slope changing within
    # the step. See `THERMAL_EXPLICIT_SAFETY_FACTOR`.
    @assert 0 < mp.thermal_explicit_safety_factor <= 1 "thermal_explicit_safety_factor must be in (0, 1] (got $(mp.thermal_explicit_safety_factor)); it scales an upper bound on the stable sub-step, so a value above 1 runs `ExplicitThermal` knowingly unstable — use `thermal_solver=ImplicitThermal()` for a scheme with no step limit"

    # Heat capacity
    @assert mp.heat_capacity_method in (:constant, :CuffeyPaterson) "heat_capacity_method must be constant or CuffeyPaterson"
    @assert 1500 <= mp.heat_capacity_ice <= 2500 "heat_capacity_ice must be in [1500, 2500]"
    @assert mp.rain_heat_capacity in (:water, :ice) "rain_heat_capacity must be water or ice"

    # Grain growth
    @assert mp.grain_growth_method in (:Marbouty, :Arthern, :hybrid, :Taillandier2007, :FlannerZender2006) "grain_growth_method must be one of: Marbouty, Arthern, hybrid, Taillandier2007, FlannerZender2006"

    # Water
    @assert mp.water_irreducible_method in (:constant, :ColeouLesaffre) "water_irreducible_method must be constant or ColeouLesaffre"
    @assert mp.melt_geometry in (:thickness, :density) "melt_geometry must be thickness or density"
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

    # Blowing snow. The `drift_rate` bound is far above any observed net drift divergence
    # (Lenaerts et al. (2012) report up to ~200 kg m-2 yr-1 of drifting-snow erosion in the
    # windiest Antarctic escarpment cells) and exists to catch a per-timestep value passed
    # where an annual rate belongs. Signed: negative is net deposition.
    @assert mp.blowing_snow_method in (:none, :Crocus, :Keenan2021) "blowing_snow_method must be none, Crocus or Keenan2021"
    # Suspension sublimation belongs to the Crocus scheme; the Keenan et al. (2021) scheme
    # moves no mass across the surface.
    if mp.blowing_snow_method === :Keenan2021
        @assert !mp.blowing_snow_sublimation "blowing_snow_sublimation is part of blowing_snow_method=:Crocus and cannot be combined with :Keenan2021"
    end
    @assert 0 < mp.drift_fetch_length <= 1000 "drift_fetch_length must be in (0, 1000] m (got $(mp.drift_fetch_length))"
    @assert -2000 <= mp.drift_rate <= 2000 "drift_rate must be in [-2000, 2000] kg m-2 yr-1 (an annual rate, positive for erosion)"

    # Slope as a gradient, not degrees. The bound is far above any ice-sheet surface slope
    # (the RetMIP sites reach ~0.010 m m-1) and is there to catch degrees passed by mistake.
    @assert 0 <= mp.surface_slope <= 1 "surface_slope must be in [0, 1] m m-1 (a gradient, not degrees)"

    # Lower boundary
    @assert mp.lower_boundary in (:ice, :ground, :prescribed) "lower_boundary must be ice, ground or prescribed"
    mp.lower_boundary === :ground && validate_ground_parameters(mp)
    mp.lower_boundary === :prescribed && validate_prescribed_parameters(mp)
    return nothing
end

# The `:prescribed` base has no soil column, so only the snow-over-ground settings apply; the
# base temperature comes from the forcing.
function validate_prescribed_parameters(mp::ModelParameters)
    @assert mp.thermal_solver isa ExplicitThermal "lower_boundary=:prescribed requires thermal_solver=ExplicitThermal()"
    @assert mp.horizontal_strain_rate == 0 "lower_boundary=:prescribed requires horizontal_strain_rate = 0"
    @assert 0.05 <= mp.ground_albedo <= 0.5 "ground_albedo must be in [0.05, 0.5]"
    @assert 0.8 <= mp.ground_emissivity <= 1 "ground_emissivity must be in [0.8, 1]"
    @assert 1e-4 <= mp.ground_roughness <= 0.5 "ground_roughness must be in [1e-4, 0.5] m"
    @assert 0 <= mp.ground_evaporation_factor <= 1 "ground_evaporation_factor must be in [0, 1]"
    @assert 3 <= mp.ground_snow_cells_max <= 1000 "ground_snow_cells_max must be in [3, 1000] (the count controller merges pairs below the surface cell)"
    @assert mp.drift_rate == 0 "lower_boundary=:prescribed requires drift_rate = 0 (prescribed erosion of a thin seasonal pack is not supported)"
    return nothing
end

function validate_ground_parameters(mp::ModelParameters)
    @assert mp.thermal_solver isa ExplicitThermal "lower_boundary=:ground requires thermal_solver=ExplicitThermal()"
    @assert mp.horizontal_strain_rate == 0 "lower_boundary=:ground requires horizontal_strain_rate = 0 (the ground does not strain)"
    @assert isfinite(mp.ground_temperature_deep) && 200 <= mp.ground_temperature_deep <= 320 "lower_boundary=:ground requires ground_temperature_deep, the mean annual ground temperature [K], in [200, 320] (got $(mp.ground_temperature_deep))"
    @assert 1 <= mp.ground_depth <= 100 "ground_depth must be in [1, 100] m (got $(mp.ground_depth))"
    @assert 0 < mp.ground_porosity < 1 "ground_porosity must be in (0, 1) (got $(mp.ground_porosity))"
    @assert 0 <= mp.ground_water <= mp.ground_porosity "ground_water must be in [0, ground_porosity] (got $(mp.ground_water), porosity $(mp.ground_porosity))"
    @assert 0.5 <= mp.ground_solid_conductivity <= 8 "ground_solid_conductivity must be in [0.5, 8] W m-1 K-1 (quartz is 7.7)"
    @assert 1000 <= mp.ground_solid_density <= 3500 "ground_solid_density must be in [1000, 3500] kg m-3"
    @assert 500 <= mp.ground_solid_heat_capacity <= 2000 "ground_solid_heat_capacity must be in [500, 2000] J kg-1 K-1"
    @assert 0.05 <= mp.ground_albedo <= 0.5 "ground_albedo must be in [0.05, 0.5]"
    @assert 0.8 <= mp.ground_emissivity <= 1 "ground_emissivity must be in [0.8, 1]"
    @assert 1e-4 <= mp.ground_roughness <= 0.5 "ground_roughness must be in [1e-4, 0.5] m"
    @assert 0 <= mp.ground_evaporation_factor <= 1 "ground_evaporation_factor must be in [0, 1]"
    @assert 3 <= mp.ground_snow_cells_max <= 1000 "ground_snow_cells_max must be in [3, 1000] (the count controller merges pairs below the surface cell)"
    @assert mp.drift_rate == 0 "lower_boundary=:ground requires drift_rate = 0 (prescribed erosion of a thin seasonal pack is not supported)"
    return nothing
end
