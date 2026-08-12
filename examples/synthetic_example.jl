# Simple example of running GEMB using synthetic climate forcing.
# Equivalent to MATLAB's GEMB_example_synthetic.m

using GEMB
using GEMB_ClimateForcing

## Set up the model and run it:

# Generate 3-hourly synthetic climate forcing data:
time_step_hours = 3
ds = simulate_climate_forcing("test_1", time_step_hours)
cf = GEMB.initialize_forcing(ds)

# Initialize model parameters:
mp = initialize_parameters(output_frequency=:daily)

# Create a climatological average time series:
cf_climatology = forcing_climatology(cf)

# Initialize a column:
profile = initialize_profile(mp, cf_climatology)

# Spin up a profile for up to 75 years of average forcing, exiting early once the
# depth-averaged density converges:
mp_spinup = initialize_parameters(output_frequency=:last)
profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup;
                             max_iterations=75, convergence_delta_density=0.01)

# The spun-up profile carries provenance: which climatology years were averaged
# and how the spinup converged (`metadata` is re-exported from DimensionalData).
prov = metadata(profile_spunup)
println("Spinup: ", prov[:spinup_cycles], " cycles, converged=", prov[:spinup_converged],
        ", climatology years=", prov[:climatology_n_years])

# Run GEMB with the spun-up profile. The provenance (including spinup_performed=true)
# is carried onto the output metadata; a profile run without spinup would instead
# record spinup_performed=false.
output = gemb(profile_spunup, cf, mp)

## Examine results:

# Get a 2D matrix of grid cell centers:
z_center = dz2z(parent(output[:dz]))

# Print summary statistics:
println("Simulation complete!")
println("  Time steps: ", size(output[:melt], 1))
println("  Profile layers: ", size(output[:temperature], 1))
println("  Mean surface albedo: ", round(mean(parent(output[:albedo_surface])), digits=3))
println("  Total melt: ", round(sum(parent(output[:melt])), digits=2), " kg/m²")
println("  Total runoff: ", round(sum(parent(output[:runoff])), digits=2), " kg/m²")

# To visualize results, use a plotting package:
# using CairoMakie
# fig = Figure()
# ax = Axis(fig[1,1], ylabel="Column height (m)")
# heatmap!(ax, 1:size(z_center,2), z_center[:,1], parent(output[:temperature]))
# display(fig)
