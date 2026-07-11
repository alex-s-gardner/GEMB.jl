# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GEMB.jl is a Julia implementation of the Glacier Energy and Mass Balance model. This is a physics-based snow/firn/ice model that simulates surface mass and energy balance processes including:
- Snow grain metamorphism
- Albedo evolution
- Shortwave radiation absorption with depth
- Temperature profile evolution
- Melt and refreezing
- Densification (compaction and wet compaction)
- Layer management

The model is a translation from MATLAB to Julia, maintaining fidelity to the original implementation while leveraging Julia's performance and type system.

## Development Commands

### Testing
```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a specific test file (from Julia REPL)
julia --project=.
using GEMB
using Test
include("test/test_thermal_conductivity.jl")
```

### Documentation
```bash
# Build documentation locally
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

### Package Management
```bash
# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Add a new dependency
julia --project=. -e 'using Pkg; Pkg.add("PackageName")'
```

## Architecture

### Core Structure

GEMB follows a modular physics-based architecture:

1. **Types** (`types.jl`): Defines three key structs:
   - `ModelParameters`: All model configuration (38 fields matching MATLAB implementation)
   - `ClimateForcing`: Time-series meteorological forcing with DimensionalData.jl `DimArray`s
   - `ClimateForcingStep`: Single timestep forcing values (plain struct for performance)

2. **Initialization** (`initialize_*.jl`):
   - `initialize_parameters()`: Creates `ModelParameters` with defaults or user overrides
   - `initialize_forcing()`: Loads/creates `ClimateForcing` time series
   - `initialize_profile()`: Sets up initial vertical profile of temperature, density, etc.

3. **Physics Modules** (`calculate_*.jl`): Individual physics processes called in sequence:
   - `calculate_grain_size()`: Snow grain metamorphism
   - `calculate_albedo()`: Surface albedo evolution
   - `calculate_shortwave_radiation()`: SW absorption profile
   - `calculate_temperature()`: Energy balance and temperature evolution
   - `calculate_accumulation()`: Snow/rain addition
   - `calculate_melt()`: Melt, runoff, and refreezing
   - `calculate_density()`: Non-melt densification
   - `manage_layers()`: Grid management (merging/splitting)

4. **Integration**:
   - `gemb_core()`: Single timestep integration calling all physics modules in sequence
   - `gemb()` in `gemb_driver.jl`: Main driver function that loops over time, accumulates output

5. **Utilities**:
   - `gemb_spinup()` in `spinup.jl`: Cycles forcing to reach equilibrium
   - `gemb_profile()`, `gemb_interp()` in `profile_extract.jl`: Extract profiles at specific times/depths
   - `forcing_climatology()`: Creates synthetic climatology from time series
   - `simulate_climate_forcing()` in `simulate/`: Generates synthetic forcing data

### Data Flow

```
Initialize Profile → Time Loop [ 
    Extract Forcing Step → 
    gemb_core() [
        1. Grain metamorphism
        2. Albedo
        3. Shortwave radiation
        4. Temperature evolution
        5. Accumulation
        6. Melt/refreeze
        7. Layer management
        8. Densification
    ] → 
    Accumulate Output
] → Return DimStack
```

### Key Design Principles

- **MATLAB Fidelity**: Function signatures and physics match original MATLAB implementation. Comments reference MATLAB line numbers.
- **DimensionalData.jl**: All input/output arrays use `DimArray` with explicit dimensions (`Ti` for time, `Z` for vertical)
- **Type Stability**: Physics calculations use plain Arrays, not DimArrays, for performance
- **Immutable Parameters**: `ModelParameters` is immutable; create new instance for modifications
- **Energy/Mass Conservation**: When `verbose=true`, `gemb_core()` validates conservation laws each timestep
- **Output Padding**: Profile outputs include padding (`output_padding` parameter) to accommodate column growth without reallocation

### Test Structure

Tests validate against MATLAB reference data stored in `test/` directory (MAT files). Each physics module has a corresponding test file that:
1. Loads MATLAB reference inputs
2. Runs Julia function
3. Compares outputs with reference (typically to ~1e-12 relative tolerance)

The test suite in `runtests.jl` runs all module tests in order of dependency.

### Constants

Physical constants are defined in `constants.jl`:
- `C_ICE`: Specific heat capacity of ice
- `LF`: Latent heat of fusion
- `CtoK`: Celsius to Kelvin offset
- etc.

These are used throughout the codebase and match MATLAB conventions.

## Testing Against MATLAB

When modifying physics functions, always run the corresponding test to ensure MATLAB equivalence is maintained. Tests use the MAT.jl package to load reference data. If you add new physics, generate reference data using `test/generate_reference_data.m` (MATLAB script).

## DimensionalData.jl Usage

- Always use `Ti` dimension for time indices
- Use `Z` dimension for vertical (depth) indices
- When extracting values from DimArray for computation, use `parent()` or `collect()` to get plain Arrays
- Results are packaged back into DimStack for return
