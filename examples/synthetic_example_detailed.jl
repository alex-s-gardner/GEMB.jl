# Detailed synthetic example with output comparable to MATLAB
# Saves results for cross-validation

using GEMB
using GEMB_ClimateForcing
using Statistics
using Dates

println("="^60)
println("GEMB.jl Synthetic Example - Detailed Output")
println("="^60)

## Set up the model and run it:

# Generate 3-hourly synthetic climate forcing data:
println("\n1. Generating synthetic climate forcing...")
time_step_hours = 3
ds = simulate_climate_forcing("test_1", time_step_hours)
cf = initialize_forcing(ds)

println("   Climate forcing generated:")
println("     Time steps: $(length(cf.time))")
println("     Duration: $(cf.time[end] - cf.time[1]) days")
println("     Variables: $(keys(cf))")

# Initialize model parameters:
println("\n2. Initializing model parameters...")
mp = initialize_parameters(output_frequency=:daily)
println("   Output frequency: $(mp.output_frequency)")
println("   Densification model: $(mp.densification_model)")

# Initialize a column. Its grid is sized to the depth this climate needs.
println("\n3. Initializing profile...")
profile = initialize_profile(mp, cf)
println("   Initial layers: $(length(profile.dz))")
println("   Initial column height: $(sum(profile.dz)) m")
println("   Initial mean temperature: $(round(mean(profile.temperature), digits=2)) K")
println("   Initial mean density: $(round(mean(profile.density), digits=1)) kg/m³")

# Create a climatological average time series:
println("\n4. Creating climatology...")
cf_climatology = forcing_climatology(cf)
println("   Climatology length: $(length(cf_climatology.temperature_air)) steps")

# Spin up a profile for up to 75 years of average forcing (exits early on convergence).
# gemb_spinup internally forces output_frequency=:last, so no separate :last params
# object is needed.
println("\n5. Running spinup (up to 75 years)...")
profile_spunup = gemb_spinup(profile, cf_climatology, mp;
                             max_iterations=75, convergence_delta_density=0.01)
# Provenance recorded on the spun-up profile (see DimensionalData.metadata):
prov = metadata(profile_spunup)
println("   Spinup complete! $(prov[:spinup_cycles]) cycles, converged=$(prov[:spinup_converged])")
println("   Climatology years averaged: $(prov[:climatology_n_years])")
println("   Post-spinup layers: $(length(profile_spunup.dz))")
println("   Post-spinup column height: $(round(sum(profile_spunup.dz), digits=2)) m")
println("   Post-spinup mean temperature: $(round(mean(profile_spunup.temperature), digits=2)) K")
println("   Post-spinup mean density: $(round(mean(profile_spunup.density), digits=1)) kg/m³")

# Run GEMB with the spun-up profile:
println("\n6. Running main simulation...")
output = gemb(profile_spunup, cf, mp)
println("   Simulation complete!")

## Examine results:

println("\n" * "="^60)
println("SIMULATION RESULTS")
println("="^60)

# Get a 2D matrix of grid cell centers:
z_center = dz2z(parent(output[:dz]))

# Print summary statistics:
println("\nBasic Statistics:")
println("  Time steps: $(size(output[:melt], 1))")
println("  Profile layers: $(size(output[:temperature], 1))")
println("  Mean surface albedo: $(round(mean(parent(output[:albedo_broadband])), digits=3))")
println("  Total melt: $(round(sum(parent(output[:melt])), digits=2)) kg/m²")
println("  Total runoff: $(round(sum(parent(output[:runoff])), digits=2)) kg/m²")
println("  Total refreezing: $(round(sum(parent(output[:refreezing])), digits=2)) kg/m²")
println("  Total accumulation: $(round(sum(parent(output[:accumulation])), digits=2)) kg/m²")

# Surface properties
println("\nSurface Properties:")
temp_surface = parent(output[:temperature])[1, :]
albedo_broadband = parent(output[:albedo_broadband])
println("  Mean surface temperature: $(round(mean(temp_surface), digits=2)) K")
println("  Min surface temperature: $(round(minimum(temp_surface), digits=2)) K")
println("  Max surface temperature: $(round(maximum(temp_surface), digits=2)) K")
println("  Albedo range: $(round(minimum(albedo_broadband), digits=3)) - $(round(maximum(albedo_broadband), digits=3))")

# Energy balance
println("\nEnergy Balance:")
sw_net = parent(output[:shortwave_net])
lw_net = parent(output[:longwave_net])
sensible = parent(output[:sensible_heat_flux])
latent = parent(output[:latent_heat_flux])
println("  Mean net shortwave: $(round(mean(sw_net), digits=1)) W/m²")
println("  Mean net longwave: $(round(mean(lw_net), digits=1)) W/m²")
println("  Mean sensible heat flux: $(round(mean(sensible), digits=1)) W/m²")
println("  Mean latent heat flux: $(round(mean(latent), digits=1)) W/m²")

# Column evolution
println("\nColumn Evolution:")
# `ice_flux` is the per-interval basal ice flux; its cumulative sum, negated, is surface
# elevation change against a datum fixed in the ice at the start of the run.
cumulative_ice_flux = cumsum(parent(output[:ice_flux]))
println("  Cumulative basal ice flux: $(round(cumulative_ice_flux[end], digits=2)) m of ice")
println("  Surface elevation change: $(round(-cumulative_ice_flux[end], digits=2)) m")

# Profile statistics at final timestep
println("\nFinal Profile (last timestep):")
temp_final = parent(output[:temperature])[:, end]
density_final = parent(output[:density])[:, end]
grain_final = parent(output[:grain_radius])[:, end]
temp_valid = temp_final[isfinite.(temp_final)]
density_valid = density_final[isfinite.(density_final)]
grain_valid = grain_final[isfinite.(grain_final)]
println("  Active layers: $(length(temp_valid))")
println("  Mean temperature: $(round(mean(temp_valid), digits=2)) K")
println("  Mean density: $(round(mean(density_valid), digits=1)) kg/m³")
println("  Mean grain radius: $(round(mean(grain_valid)*1000, digits=3)) mm")

println("\n" * "="^60)
println("Simulation completed successfully!")
println("="^60)

# Save output for cross-validation
using MAT
output_dict = Dict(
    "time" => collect(output.time),
    "temperature" => parent(output[:temperature]),
    "density" => parent(output[:density]),
    "melt" => parent(output[:melt]),
    "runoff" => parent(output[:runoff]),
    "albedo_broadband" => parent(output[:albedo_broadband]),
    "ice_flux" => parent(output[:ice_flux]),
    "shortwave_net" => parent(output[:shortwave_net]),
    "longwave_net" => parent(output[:longwave_net])
)

matwrite("julia_synthetic_output.mat", output_dict)
println("\nOutput saved to julia_synthetic_output.mat for cross-validation")
