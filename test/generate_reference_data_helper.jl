"""
Helper script to generate reference data for MATLAB validation.

This script should be run AFTER generate_reference_data.m has been executed in MATLAB.
It verifies that all expected .mat files exist and can be loaded.
"""

using MAT
using Printf

const REF_DIR = joinpath(@__DIR__, "reference_data")
const EXPECTED_FILES = [
    "thermal_conductivity.mat",
    "turbulent_heat_flux.mat",
    "initialize_profile.mat",
    "calculate_shortwave_radiation.mat",
    "gemb_core.mat"
]

println("=" ^ 70)
println("MATLAB Reference Data Validation")
println("=" ^ 70)
println()

# Check if reference data directory exists
if !isdir(REF_DIR)
    println("❌ Reference data directory not found:")
    println("   $(REF_DIR)")
    println()
    println("Please run the following in MATLAB:")
    println("  cd /Users/gardnera/Documents/GitHub/GEMB.jl/test")
    println("  generate_reference_data")
    exit(1)
end

# Check each expected file
all_present = true
file_info = Dict()

for filename in EXPECTED_FILES
    filepath = joinpath(REF_DIR, filename)

    if isfile(filepath)
        try
            data = matread(filepath)
            nkeys = length(keys(data))
            filesize = filesize(filepath)

            println("✓ $(filename)")
            println("  - Keys: $(nkeys)")
            println("  - Size: $(round(filesize/1024, digits=2)) KB")
            println("  - Variables: $(join(keys(data), ", "))")
            println()

            file_info[filename] = (status=:ok, nkeys=nkeys)
        catch e
            println("⚠ $(filename) - Error reading file")
            println("  Error: $e")
            println()
            all_present = false
            file_info[filename] = (status=:error, error=e)
        end
    else
        println("❌ $(filename) - Not found")
        println()
        all_present = false
        file_info[filename] = (status=:missing,)
    end
end

println("=" ^ 70)

if all_present
    println("✅ All reference data files present and readable!")
    println()
    println("You can now run the test suite with MATLAB validation:")
    println("  julia --project=. -e 'using Pkg; Pkg.test()'")
    println()
    exit(0)
else
    println("❌ Some reference data files are missing or unreadable")
    println()
    println("To generate reference data, run in MATLAB:")
    println("  cd /Users/gardnera/Documents/GitHub/GEMB.jl/test")
    println("  generate_reference_data")
    println()
    println("Alternatively, if MATLAB licensing issues persist, the tests")
    println("will run with physics-based validation only (no MATLAB comparison)")
    println()
    exit(1)
end
