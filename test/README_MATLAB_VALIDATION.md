# MATLAB Cross-Validation System

This directory contains the infrastructure for validating GEMB.jl outputs against the reference MATLAB implementation.

## Overview

The test suite supports **dual-mode testing**:

1. **Standalone Mode** (Default): Tests validate using physics-based assertions
2. **MATLAB Validation Mode**: When reference data is available, tests additionally validate exact numerical agreement with MATLAB

## Quick Start

### Generate MATLAB Reference Data

```matlab
% In MATLAB, run:
cd /Users/gardnera/Documents/GitHub/GEMB.jl/test
generate_reference_data
```

This creates `.mat` files in `reference_data/`:
- `thermal_conductivity.mat`
- `turbulent_heat_flux.mat`
- `initialize_profile.mat`
- `calculate_shortwave_radiation.mat`
- `gemb_core.mat`

### Verify Reference Data

```bash
# In Julia, check that reference data was generated successfully
julia --project=. test/generate_reference_data_helper.jl
```

### Run Tests with MATLAB Validation

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Architecture

### Test Files Structure

Each physics module test file (e.g., `test_thermal_conductivity.jl`) contains:

1. **Physics-based tests**: Validate against known formulas and behavior
2. **MATLAB validation test**: Compares Julia output with MATLAB reference data

Example structure:
```julia
@testset "Physics validation" begin
    # Tests using formulas and expected behavior
end

# MATLAB validation (automatically skipped if data unavailable)
matlab_validation_testset("function_name", "reference_file.mat") do ref
    # Compare Julia vs MATLAB outputs
    @test julia_output ≈ ref["matlab_output"] rtol=1e-12
end
```

### Utility Functions

**File**: `test_utils.jl`

Provides:
- `load_matlab_reference(filename)` - Load .mat files
- `compare_with_matlab(julia_val, matlab_val, name)` - Compare with diagnostics
- `matlab_validation_testset(fn, name, file)` - Create MATLAB validation test block

### Reference Data Generation

**File**: `generate_reference_data.m`

MATLAB script that:
1. Loads MATLAB GEMB from `/Users/gardnera/Documents/GitHub/GEMB/src`
2. Runs each function with controlled inputs
3. Saves outputs to `.mat` files in `reference_data/`

## Tolerance Levels

MATLAB validation uses strict tolerances matching the CLAUDE.md specification:

- **Relative tolerance**: `1e-12` (0.0000000001%)
- **Absolute tolerance**: `1e-14`

This ensures near-perfect numerical agreement between implementations.

## Test Status

### Files with MATLAB Validation ✅

- ✅ `test_thermal_conductivity.jl`
- ✅ `test_turbulent_heat_flux.jl`
- ✅ `test_calculate_shortwave_radiation.jl`
- ✅ `test_initialize_profile.jl`
- ✅ `test_gemb_core.jl`

### Files Needing MATLAB Validation 🔄

- 🔄 `test_calculate_temperature.jl`
- 🔄 `test_calculate_albedo.jl`
- 🔄 `test_calculate_density.jl`
- 🔄 `test_calculate_accumulation.jl`
- 🔄 `test_calculate_melt.jl`
- 🔄 `test_calculate_grain_size.jl`
- 🔄 `test_manage_layers.jl`

Note: These functions can be validated once corresponding reference data is generated in MATLAB.

## Troubleshooting

### MAT.jl Not Available

The test system gracefully handles missing MAT.jl:
- MATLAB validation tests are automatically skipped
- Physics-based tests still run normally
- No errors are thrown

To enable MATLAB validation:
```bash
julia --project=test -e 'using Pkg; Pkg.add("MAT")'
```

### Reference Data Not Found

When `.mat` files aren't available:
- Tests display: `⊘ MATLAB reference not available`
- Tests are marked as `broken` (not failed)
- Instructions are provided to generate reference data

### MATLAB Licensing Issues

If MATLAB license fails:
1. Tests will run in standalone mode (physics-based validation only)
2. Reference data can be generated later when license is available
3. Tests can be re-run to enable MATLAB validation once data exists

### Validation Failures

If MATLAB validation fails (Julia ≠ MATLAB):
1. Check that MATLAB GEMB is up-to-date at `/Users/gardnera/Documents/GitHub/GEMB`
2. Verify reference data is current (regenerate if needed)
3. Inspect the specific discrepancy:
   ```julia
   # Tests print diagnostics on failure:
   # max_abs_diff = ...
   # max_rel_diff = ...
   ```

## Continuous Integration

### Local CI Workflow

```bash
# Step 1: Generate reference data (requires MATLAB)
matlab -batch "cd test; generate_reference_data"

# Step 2: Verify data
julia --project=. test/generate_reference_data_helper.jl

# Step 3: Run tests with validation
julia --project=. -e 'using Pkg; Pkg.test()'
```

### CI Without MATLAB

Tests are designed to run in CI environments without MATLAB:
- Physics-based tests provide comprehensive coverage
- MATLAB validation is optional enhancement
- No CI failures due to missing reference data

## Adding New MATLAB Validations

To add MATLAB validation for a new function:

### 1. Update `generate_reference_data.m`

```matlab
%% Test your_function
fprintf('Generating your_function reference data...\n')

% Set up inputs
input1 = ...;
input2 = ...;

% Call MATLAB function
output = your_function(input1, input2);

% Save reference data
save('reference_data/your_function.mat', 'input1', 'input2', 'output');
```

### 2. Add Validation to Test File

```julia
# At end of test_your_function.jl:
matlab_validation_testset("your_function", "your_function.mat") do ref
    # Extract inputs
    input1 = ref["input1"]
    input2 = ref["input2"]
    
    # Run Julia function
    output_julia = your_function(input1, input2)
    
    # Compare with MATLAB
    @test output_julia ≈ ref["output"] rtol=1e-12 atol=1e-14
end
```

### 3. Regenerate Reference Data

```matlab
cd /Users/gardnera/Documents/GitHub/GEMB.jl/test
generate_reference_data
```

### 4. Verify

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Performance

- Loading .mat files: ~10-50ms per file
- MATLAB validation adds: ~1-5% to test runtime
- Total test time with validation: <30 seconds

## References

- MATLAB GEMB: `/Users/gardnera/Documents/GitHub/GEMB`
- Julia GEMB: `/Users/gardnera/Documents/GitHub/GEMB.jl`
- Tolerance spec: `CLAUDE.md` (1e-12 relative error)
- MAT.jl docs: https://github.com/JuliaIO/MAT.jl

---

**Last Updated**: 2026-07-11  
**Status**: Infrastructure complete, partial validation active  
**Maintainer**: See `IMPLEMENTATION_STATUS.md`
