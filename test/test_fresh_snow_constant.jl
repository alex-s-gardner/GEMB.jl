# `new_snow_method = :Constant`: a fresh-snow density of the caller's choosing (`new_snow_density`),
# with the fresh-grain properties of the fixed constants.

# The thin-snow forcing in a module of its own, so the regression pins' later `include` of the same file
# does not redefine its methods.
module ConstantSnowForcing
using GEMB
include(joinpath(@__DIR__, "synthetic_thin_snow_forcing.jl"))
end

@testset "fresh snow: :Constant density" begin
    for ρ in (50.0, 250.0, 300.0, 550.0)
        mp = GEMB.ModelParameters(new_snow_method=:Constant, new_snow_density=ρ)
        @test GEMB.fresh_snow_density(mp, 260.0, 200.0, 5.0, 250.0, 11.0) === ρ
        @test_nowarn GEMB.validate_parameters(mp)
    end
    # No other method reads `new_snow_density`.
    for method in (:Constant150, :Constant315, :Constant350, :Fausto, :FaustoFit, :Pahaut,
                   :Kaspers, :KuipersMunneke, :vanKampenhout2017)
        a = GEMB.ModelParameters(new_snow_method=method)
        b = GEMB.ModelParameters(new_snow_method=method, new_snow_density=250.0)
        @test GEMB.fresh_snow_density(a, 260.0, 200.0, 5.0, 250.0, 7.0) ===
              GEMB.fresh_snow_density(b, 260.0, 200.0, 5.0, 250.0, 7.0)
    end
end

@testset "fresh snow: :Constant density is validated" begin
    @test_throws "new_snow_density must be in [50.0, 550.0] kg m-3 (got 40.0)" GEMB.validate_parameters(
        GEMB.ModelParameters(new_snow_method=:Constant, new_snow_density=40.0))
    @test_throws "new_snow_density must be in [50.0, 550.0] kg m-3 (got 600.0)" GEMB.validate_parameters(
        GEMB.ModelParameters(new_snow_method=:Constant, new_snow_density=600.0))
end

@testset "fresh snow: :Constant cell has the default grain properties" begin
    n = 5
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, temperature_air=250.0, precipitation=50.0,
        wind_speed=5.0, temperature_air_mean=260.0, precipitation_mean=200.0, wind_speed_mean=5.0)
    cell(mp) = GEMB.calculate_accumulation(260.0 * ones(n), 0.1 * ones(n), 400.0 * ones(n),
        zeros(n), 0.5 * ones(n), 0.5 * ones(n), 0.5 * ones(n), zeros(n), cfs, mp, false)
    out = cell(GEMB.ModelParameters(new_snow_method=:Constant, new_snow_density=250.0, column_dzmin=0.05))
    @test out[3][1] == 250.0
    @test out[2][1] ≈ 50.0 / 250.0 rtol = 1e-14
    @test out[5][1] == GEMB.RE_NEW_SNOW
    @test out[6][1] == GEMB.GDN_NEW_SNOW
    @test out[7][1] == GEMB.GSP_NEW_SNOW
    # At 350 every returned array equals `:Constant350`'s.
    a = cell(GEMB.ModelParameters(new_snow_method=:Constant, new_snow_density=350.0, column_dzmin=0.05))
    b = cell(GEMB.ModelParameters(new_snow_method=:Constant350, column_dzmin=0.05))
    @test all(isequal.(a, b))
end

@testset "fresh snow: :Constant at 350 runs bit-identical to :Constant350" begin
    cf = ConstantSnowForcing.synthetic_thin_snow_forcing()
    run(mp) = gemb(initialize_profile(mp, cf; constant_density=true, constant_temperature=true), cf, mp)
    kw = (; output_frequency=:daily, albedo_ice=0.30)
    a = @test_logs (:warn, r"above the melt point") match_mode = :any run(ModelParameters(; kw...,
        new_snow_method=:Constant, new_snow_density=350.0))
    b = @test_logs (:warn, r"above the melt point") match_mode = :any run(ModelParameters(; kw...,
        new_snow_method=:Constant350))
    @test keys(a) == keys(b)
    @test all(k -> isequal(parent(a[k]), parent(b[k])), keys(a))
    # A different density moves the run, so the comparison above has teeth.
    c = @test_logs (:warn, r"above the melt point") match_mode = :any run(ModelParameters(; kw...,
        new_snow_method=:Constant, new_snow_density=250.0))
    @test !isequal(parent(a[:density]), parent(c[:density]))
end
