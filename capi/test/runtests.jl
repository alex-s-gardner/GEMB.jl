# Tests for the C API boundary.
#
# The compiled library passes numbers across as bare arrays in a declared order, so the contract
# is the order itself. These tests hold `capi/gemb.h`, `capi/gemb_capi.jl` and the Julia structs
# to the same one: a field added to `ClimateForcingStep` or to the flux set fails here rather than
# silently shifting a column of numbers in a host model.
#
# The end-to-end numerical check is `capi/test/test_capi.c`, driven by `capi/run_tests.jl`.
#
# Usage:
#   julia --project=capi capi/test/runtests.jl

using Test
using GEMB

include(joinpath(@__DIR__, "..", "gemb_capi.jl"))
using .GEMBCApi

const HEADER = read(joinpath(@__DIR__, "..", "gemb.h"), String)

"""
    header_define(name) -> Int

The value of a `#define` in `capi/gemb.h`.
"""
function header_define(name::AbstractString)
    m = match(Regex("^#define\\s+" * name * "\\s+\\(?(-?\\d+)\\)?\\s*(?:/\\*.*)?\$", "m"), HEADER)
    m === nothing && error("no #define for $name in capi/gemb.h")
    return parse(Int, m.captures[1])
end

@testset "C API boundary" begin
    @testset "array sizes agree with the header" begin
        @test header_define("GEMB_N_OPTIONS") == GEMBCApi.N_OPTIONS
        @test header_define("GEMB_N_NUMERIC") == GEMBCApi.N_NUMERIC
        @test header_define("GEMB_N_FLAGS") == GEMBCApi.N_FLAGS
        @test header_define("GEMB_N_FORCING") == GEMBCApi.N_FORCING
        @test header_define("GEMB_N_FLUX") == GEMBCApi.N_FLUX
    end

    @testset "status codes agree with the header" begin
        @test header_define("GEMB_OK") == GEMBCApi.GEMB_OK
        @test header_define("GEMB_ERR_BAD_ARGUMENT") == GEMBCApi.GEMB_ERR_BAD_ARGUMENT
        @test header_define("GEMB_ERR_BAD_HANDLE") == GEMBCApi.GEMB_ERR_BAD_HANDLE
        @test header_define("GEMB_ERR_BAD_OPTION") == GEMBCApi.GEMB_ERR_BAD_OPTION
        @test header_define("GEMB_ERR_PHYSICS") == GEMBCApi.GEMB_ERR_PHYSICS
        @test header_define("GEMB_ERR_UNEXPECTED") == GEMBCApi.GEMB_ERR_UNEXPECTED
    end

    @testset "forcing order matches ClimateForcingStep" begin
        # The `forcing` array is unpacked positionally into `ClimateForcingStep`, so the header's
        # GEMB_FRC_* indices must be that struct's field order.
        fields = fieldnames(GEMB.ClimateForcingStep)
        @test length(fields) == GEMBCApi.N_FORCING
        for (i, field) in enumerate(fields)
            @test header_define("GEMB_FRC_" * uppercase(string(field))) == i - 1
        end
    end

    @testset "flux order matches gemb_core" begin
        # `gemb_step` writes flux terms one by one, so the header's GEMB_FLUX_* indices must match
        # the order `gemb_core` returns them in.
        mp = GEMB.ModelParameters(; dt_divisors=GEMB.fast_divisors(108000000) ./ 10000)
        n = 8
        state = (temperature=fill(263.0, n), dz=fill(0.5, n), density=fill(400.0, n),
            water=zeros(n), grain_radius=fill(0.5, n), grain_dendricity=zeros(n),
            grain_sphericity=fill(0.5, n), age=zeros(n),
            evaporation_condensation=0.0, melt_surface=0.0)
        forcing = GEMB.ClimateForcingStep(10800.0, 263.0, 90000.0, 0.0, 5.0, 100.0, 250.0,
            300.0, 258.0, 5.0, 0.3, 2.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1)
        _, flux = GEMB.gemb_core(state, forcing, mp, false;
            thermal_workspace=GEMB.ThermalWorkspace())

        names = keys(flux)
        # `viscosity` is last and crosses through its own buffer, not `flux_out`.
        @test names[end] === :viscosity
        @test length(names) - 1 == GEMBCApi.N_FLUX
        for i in 1:GEMBCApi.N_FLUX
            @test header_define("GEMB_FLUX_" * uppercase(string(names[i]))) == i - 1
        end
    end

    @testset "spec arrays name real ModelParameters fields" begin
        M = GEMB.ModelParameters
        for (fields, T) in ((GEMBCApi.OPTION_FIELDS, Symbol),
                            (GEMBCApi.NUMERIC_FIELDS, Float64),
                            (GEMBCApi.FLAG_FIELDS, Bool))
            for field in fields
                @test field in fieldnames(M)
                @test fieldtype(M, field) === T
            end
        end
    end

    @testset "every option value is accepted by validate_parameters" begin
        # `surface_slope` is set because `runoff_method = :Darcy` requires a positive slope; it is
        # inert for every other value under test.
        for (i, field) in enumerate(GEMBCApi.OPTION_FIELDS)
            for value in GEMBCApi.OPTION_VALUES[i]
                mp = GEMB.ModelParameters(; field => value, surface_slope=0.01,
                    dt_divisors=[1.0])
                @test_nowarn GEMB.validate_parameters(mp)
            end
        end
    end

    @testset "every in-scope option field is declared" begin
        # Nothing a column step reads may be missing from the spec arrays. `output_frequency` and
        # `initialize_age` select output cadence and profile initialization, which a single step
        # does not consult; `thermal_solver` is a type parameter, not a code.
        out_of_scope = (:output_frequency, :initialize_age, :thermal_solver)
        declared = (GEMBCApi.OPTION_FIELDS..., GEMBCApi.FLAG_FIELDS...)
        for field in fieldnames(GEMB.ModelParameters)
            fieldtype(GEMB.ModelParameters, field) in (Symbol, Bool) || continue
            @test field in declared || field in out_of_scope
        end
    end

    @testset "option codes agree with the header" begin
        # Spot-check that the header's code constants index the same entries the shim resolves.
        @test header_define("GEMB_DENSIFICATION_ARTHERN") ==
              findfirst(==(:Arthern), GEMBCApi.OPTION_VALUES[1])
        @test header_define("GEMB_ALBEDO_GARDNERSHARP") ==
              findfirst(==(:GardnerSharp), GEMBCApi.OPTION_VALUES[14])
        @test header_define("GEMB_CONDUCTIVITY_CALONNE2019") ==
              findfirst(==(:Calonne2019), GEMBCApi.OPTION_VALUES[7])
        @test header_define("GEMB_RUNOFF_DARCY") ==
              findfirst(==(:Darcy), GEMBCApi.OPTION_VALUES[13])
    end

    @testset "option field indices agree with the header" begin
        for (i, field) in enumerate(GEMBCApi.OPTION_FIELDS)
            @test header_define("GEMB_OPT_" * uppercase(string(field))) == i - 1
        end
        for (i, field) in enumerate(GEMBCApi.NUMERIC_FIELDS)
            @test header_define("GEMB_NUM_" * uppercase(string(field))) == i - 1
        end
        for (i, field) in enumerate(GEMBCApi.FLAG_FIELDS)
            @test header_define("GEMB_FLAG_" * uppercase(string(field))) == i - 1
        end
    end
end
