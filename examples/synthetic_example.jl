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

# Initialize a column. Its grid is sized to the depth this climate needs, which may
# be shallower than the configured column_depth_max; the spinup and gemb both take the
# column depth they hold fixed from the profile, so nothing has to be threaded.
profile = initialize_profile(mp, cf_climatology)

# Spin up a profile for up to 75 years of average forcing, exiting early once the
# depth-averaged density converges. gemb_spinup internally forces
# output_frequency=:last, so no separate :last params object is needed.
profile_spunup = gemb_spinup(profile, cf_climatology, mp;
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

# `output` is a DimStack: monolevel series over `Ti`, profile fields over `Z x Ti`.
# Displaying it in the REPL lists every layer with its dimensions and size, and the
# run's CF global attributes:
#
# julia> output
# ┌ 11688×222 DimStack ┐
# ├────────────────────┴──────────────────────────────────────────────────── dims ┐
#   ↓ Ti Sampled{DateTime} [1994-01-01T21:00:00, …, 2025-12-31T21:00:00] ForwardOrdered Irregular Points,
#   → Z  Sampled{Int64} 1:222 ForwardOrdered Regular Points
# ├───────────────────────────────────────────────────────────────────── layers ┤
#   :melt                     eltype: Float64 dims: Ti size: 11688
#   :runoff                   eltype: Float64 dims: Ti size: 11688
#   :refreeze                 eltype: Float64 dims: Ti size: 11688
#   ⋮
#   :temperature              eltype: Float64 dims: Z, Ti size: 222×11688
#   :density                  eltype: Float64 dims: Z, Ti size: 222×11688
#   :age                      eltype: Float64 dims: Z, Ti size: 222×11688
# ├─────────────────────────────────────────────────────────────────── metadata ┤
#   "source" => "GEMB.jl", "Conventions" => "CF-1.11", "spinup_performed" => true, …
# └───────────────────────────────────────────────────────────────────────────────┘

# Each layer carries its own CF metadata, so a value always travels with its units:
#
# julia> output[:melt]
# ┌ 11688-element DimArray{Float64, 1} melt ┐
# ├─────────────────────────────────────────┴─────────────────────────────── dims ┐
#   ↓ Ti Sampled{DateTime} [1994-01-01T21:00:00, …, 2025-12-31T21:00:00] ForwardOrdered Irregular Points
# ├─────────────────────────────────────────────────────────────────── metadata ┤
#   "units" => "kg m-2", "cell_methods" => "time: sum",
#   "long_name" => "meltwater produced", "standard_name" => "surface_snow_melt_amount"
# └───────────────────────────────────────────────────────────────────────────────┘
#  1994-01-01T21:00:00  11.4082
#  ⋮
#  2025-12-31T21:00:00   0.0

# Get a 2D matrix of grid cell centers (Z is a cell index, not a depth):
z_center = dz2z(parent(output[:dz]))

# Pull the surface row of any profile field:
#
# julia> surface_timeseries(output[:temperature])
#  1994-01-01T21:00:00  272.615
#  1994-01-02T21:00:00  272.934
#  ⋮
#  2025-12-31T21:00:00  261.531
T_surface = surface_timeseries(output[:temperature])

# Print summary statistics:
println("Simulation complete!")
println("  Time steps: ", size(output[:melt], 1))
println("  Profile layers: ", size(output[:temperature], 1))
println("  Mean surface albedo: ", round(mean(parent(output[:albedo_broadband])), digits=3))
println("  Mean surface temperature: ", round(mean(parent(T_surface)), digits=2), " K")
println("  Total melt: ", round(sum(parent(output[:melt])), digits=2), " kg/m²")
println("  Total runoff: ", round(sum(parent(output[:runoff])), digits=2), " kg/m²")
println("  Oldest firn in column: ",
        round(maximum(parent(output[:age])) / 365.25, digits=1), " yr")

## Plot the run:

# plot_output draws the whole run as a diagnostic dashboard: profile fields as
# heatmaps on the left, scalar time series on the right, sharing a linked time axis.
# It lives in a package extension, so it needs a Makie backend loaded first.
using CairoMakie
fig = plot_output(output; depthlims=(-10, 0))
save("gemb_diagnostics.png", fig; px_per_unit=2)
display(fig)
