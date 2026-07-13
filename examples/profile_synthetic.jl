# Comprehensive profiling script for GEMB.jl synthetic example
using GEMB
using Profile
using ProfileCanvas

println("="^70)
println("GEMB.jl Profiling Suite")
println("="^70)

# Step 1: Warmup run (compile everything)
println("\n[1/4] Warmup run (compiling)...")
time_step_hours = 3
cf_warmup = simulate_climate_forcing("test_1", time_step_hours)
mp_warmup = ModelParameters(output_frequency="daily")
profile_warmup = initialize_profile(mp_warmup, cf_warmup)
cf_clim_warmup = forcing_climatology(cf_warmup)
mp_spinup_warmup = ModelParameters(output_frequency="last")
profile_spunup_warmup = gemb_spinup(profile_warmup, cf_clim_warmup, mp_spinup_warmup, 75)
output_warmup = gemb(profile_spunup_warmup, cf_warmup, mp_warmup)
println("   Warmup complete!")

# Step 2: Profile full simulation
println("\n[2/4] Profiling full simulation...")
Profile.clear()
Profile.init(n=10^7, delay=0.001)  # Increase sample buffer

@profile begin
    cf = simulate_climate_forcing("test_1", time_step_hours)
    mp = ModelParameters(output_frequency="daily")
    profile = initialize_profile(mp, cf)
    cf_clim = forcing_climatology(cf)
    mp_spinup = ModelParameters(output_frequency="last")
    profile_spunup = gemb_spinup(profile, cf_clim, mp_spinup, 75)
    output = gemb(profile_spunup, cf, mp)
end

# Save flamegraph
html_file = "gemb_profile_full.html"
ProfileCanvas.html_file(html_file)
println("   Full profile saved to $html_file")

# Step 3: Profile spinup phase only (the bottleneck)
println("\n[3/4] Profiling spinup phase...")
Profile.clear()

cf_spinup = simulate_climate_forcing("test_1", time_step_hours)
mp_spinup_prof = ModelParameters(output_frequency="daily")
profile_spinup = initialize_profile(mp_spinup_prof, cf_spinup)
cf_clim_spinup = forcing_climatology(cf_spinup)
mp_sp = ModelParameters(output_frequency="last")

@profile begin
    profile_spunup_prof = gemb_spinup(profile_spinup, cf_clim_spinup, mp_sp, 75)
end

html_file_spinup = "gemb_profile_spinup.html"
ProfileCanvas.html_file(html_file_spinup)
println("   Spinup profile saved to $html_file_spinup")

# Step 4: Profile main simulation only
println("\n[4/4] Profiling main simulation...")
Profile.clear()

@profile begin
    output_main = gemb(profile_spunup_prof, cf_spinup, mp_spinup_prof)
end

html_file_main = "gemb_profile_main.html"
ProfileCanvas.html_file(html_file_main)
println("   Main simulation profile saved to $html_file_main")

println("\n" * "="^70)
println("PROFILING COMPLETE")
println("="^70)
println("Generated profiles:")
println("  1. $html_file - Full simulation flamegraph")
println("  2. $html_file_spinup - Spinup phase (64% of runtime)")
println("  3. $html_file_main - Main simulation (35% of runtime)")
println("\nFor allocation profiling, run:")
println("  julia --track-allocation=user --project=. examples/synthetic_example.jl")
println("  find src -name \"*.mem\" -exec grep -H \"[^0]\" {} \\;")
println("="^70)
