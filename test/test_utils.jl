"""
Test utilities for GEMB.jl testing framework.

Provides common functions for loading MATLAB reference data and comparing outputs.
"""

using Test

# Try to load MAT.jl for MATLAB reference data
const MAT_AVAILABLE = try
    using MAT
    true
catch
    false
end

"""
    load_matlab_reference(filename::String)

Load MATLAB reference data from test/reference_data/ directory.
Returns `nothing` if file doesn't exist or can't be loaded.
"""
function load_matlab_reference(filename::String)
    if !MAT_AVAILABLE
        return nothing
    end

    ref_file = joinpath(@__DIR__, "reference_data", filename)

    if !isfile(ref_file)
        return nothing
    end

    try
        return matread(ref_file)
    catch e
        @warn "Failed to load MATLAB reference data: $filename" exception=e
        return nothing
    end
end

"""
    compare_with_matlab(julia_val, matlab_val, name::String; rtol=1e-12, atol=1e-14)

Compare Julia output with MATLAB reference, printing informative messages.
Returns true if comparison passes, false otherwise.
"""
function compare_with_matlab(julia_val, matlab_val, name::String; rtol=1e-12, atol=1e-14)
    try
        if isapprox(julia_val, matlab_val; rtol=rtol, atol=atol)
            return true
        else
            max_diff = maximum(abs.(julia_val .- matlab_val))
            rel_diff = maximum(abs.((julia_val .- matlab_val) ./ (matlab_val .+ eps())))
            @warn "MATLAB validation mismatch for $name" max_abs_diff=max_diff max_rel_diff=rel_diff
            return false
        end
    catch e
        @warn "Error comparing with MATLAB reference for $name" exception=e
        return false
    end
end

"""
    matlab_validation_testset(testname::String, reffile::String, test_fn::Function)

Create a testset that validates against MATLAB if reference data is available.

# Arguments
- `testname`: Name of the test set
- `reffile`: Name of the .mat reference file
- `test_fn`: Function that receives the reference data dict and runs tests

# Example
```julia
matlab_validation_testset("thermal_conductivity", "thermal_conductivity.mat") do ref
    k_julia = calculate_something(ref["input"])
    @test k_julia ≈ ref["output"] rtol=1e-12
end
```
"""
function matlab_validation_testset(test_fn::Function, testname::String, reffile::String)
    ref_data = load_matlab_reference(reffile)

    if !isnothing(ref_data)
        @testset "$testname - MATLAB validation" begin
            test_fn(ref_data)
            @info "✓ MATLAB validation passed: $testname"
        end
    else
        @testset "$testname - MATLAB validation" begin
            @test_skip "MATLAB reference data not available: $reffile"
            @info "⊘ MATLAB reference not available for $testname"
            @info "  Run test/generate_reference_data.m in MATLAB to enable"
        end
    end
end

# Export functions
export load_matlab_reference, compare_with_matlab, matlab_validation_testset
