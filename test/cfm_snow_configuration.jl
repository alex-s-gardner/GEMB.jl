# The Community Firn Model's snow configuration, which tests of the Arthern, Constant350 and hybrid
# grain-growth schemes run under (the defaults are a different set). The CFM has no wind packing.
const CFM_SNOW = (; new_snow_method = :Constant350, densification_method = :Arthern,
    grain_growth_method = :hybrid, destructive_metamorphism = false,
    destructive_metamorphism_density = 175.0, blowing_snow_method = :none)
