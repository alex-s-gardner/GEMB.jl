# Simple example of running GEMB using synthetic climate forcing.
# Equivalent to MATLAB's GEMB_example_synthetic.m

using GEMB

## Set up the model and run it:

# Generate 3-hourly synthetic climate forcing data:
time_step_hours = 3
cf = simulate_climate_forcing("test_1", time_step_hours)

# Initialize model parameters:
mp = ModelParameters(output_frequency="daily")

# Initialize a column:
profile = initialize_profile(mp, cf)

# Create a climatological average time series:
cf_climatology = forcing_climatology(cf)

# Spin up a profile for 75 years of average forcing:
mp_spinup = ModelParameters(output_frequency="last")
profile_spunup = gemb_spinup(profile, cf_climatology, mp_spinup, 75)

# Run GEMB with the spun-up profile:
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
