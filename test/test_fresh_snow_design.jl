using GEMB
using GEMB_ClimateForcing
using Test
using Dates
using GEMB: Ti, dims

# Whole-model checks of the three fresh-snow options, alone and together:
#   new_snow_method = :vanKampenhout2017       (deposition density)
#   destructive_metamorphism = true            (settling under :Crocus)
#   blowing_snow_method = :Keenan2021          (wind packing)
# `verbose=true` makes `gemb` check the mass and energy budgets at every step and over the run,
# and throw if either fails, so a completed run is the closure test.

isdefined(@__MODULE__, :CFM_SNOW) || include(joinpath(@__DIR__, "cfm_snow_configuration.jl"))
isdefined(@__MODULE__, :synthetic_maritime_forcing) ||
    include(joinpath(@__DIR__, "synthetic_maritime_forcing.jl"))

# The options are measured against the schemes they replace.
const _FS_OLD = CFM_SNOW

const _FS_CONFIGS = [
    "vK17 deposition" => (; new_snow_method=:vanKampenhout2017, densification_method=:Crocus),
    "destructive metamorphism" => (; densification_method=:Crocus, destructive_metamorphism=true),
    "Keenan 2021 packing" => (; blowing_snow_method=:Keenan2021),
    "all three" => (; new_snow_method=:vanKampenhout2017, densification_method=:Crocus,
        destructive_metamorphism=true, blowing_snow_method=:Keenan2021),
]

function _fs_forcings()
    ds = GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3)
    # Four years of the 3-hourly cold, windy test_1 record, by index.
    span = Ti(1:(8 * Dates.value(Date(1998, 1, 1) - Date(1994, 1, 1))))
    return ["test_1 (cold, windy, 4 yr)" => initialize_forcing(ds)[span],
            "maritime (melt and rain, 3 yr)" => synthetic_maritime_forcing(; years=3)]
end

@testset "fresh-snow options: multi-year mass and energy closure" begin
    for (site, cf) in _fs_forcings()
        mp0 = initialize_parameters(; output_frequency=:monthly, _FS_OLD...)
        profile = initialize_profile(mp0, cf)
        z_target = sum(profile[:dz])
        n_cells = length(profile[:dz])
        out0 = gemb(profile, cf, mp0)

        for (label, kw) in _FS_CONFIGS
            @testset "$site: $label" begin
                mp = initialize_parameters(; output_frequency=:monthly, _FS_OLD..., kw...)
                out = gemb(profile, cf, mp; verbose=true)
                # The packing and settling terms move no mass across the column boundary.
                @test all(iszero, out[:blowing_snow])
                @test size(out[:dz], 1) == n_cells
                for i in axes(out[:dz], 2)
                    @test sum(view(out[:dz], :, i)) ≈ z_target atol = 1e-9
                end
                ρ = parent(out[:density])
                @test all(x -> 1.0 <= x <= mp.density_ice, ρ)
                if startswith(site, "maritime") && label == "Keenan 2021 packing"
                    # The maritime gusts (at most ~8.8 m s-1 at 5 m) stay below the
                    # 9.1 m s-1 threshold of 350 kg m-3 dendritic snow, so nothing erodes
                    # and the run is bit-identical to the default.
                    @test parent(out[:density]) == parent(out0[:density])
                else
                    @test parent(out[:density]) != parent(out0[:density])
                end
            end
        end
    end
end

@testset "fresh-snow options: the new settings are inert when off" begin
    cf = synthetic_maritime_forcing(; years=1)
    mp0 = initialize_parameters(; output_frequency=:monthly, _FS_OLD...)
    profile = initialize_profile(mp0, cf)
    out0 = gemb(profile, cf, mp0)
    # Off-state values of the new settings that are only read when their option is on.
    out1 = gemb(profile, cf, initialize_parameters(; output_frequency=:monthly, _FS_OLD...,
        destructive_metamorphism=false, destructive_metamorphism_density=300.0,
        drift_fetch_length=3.0))
    for k in (:temperature, :density, :dz, :melt, :runoff, :refreeze)
        @test parent(out1[k]) == parent(out0[k])
    end
end

@testset "fresh-snow options: wind packing raises near-surface density" begin
    cf = initialize_forcing(GEMB_ClimateForcing.simulate_climate_forcing("test_1", 3))[
        Ti(1:(8 * 365))]
    base = (; output_frequency=:monthly, CFM_SNOW..., new_snow_method=:vanKampenhout2017,
        densification_method=:Crocus, destructive_metamorphism=true)
    mp_off = initialize_parameters(; base...)
    profile = initialize_profile(initialize_parameters(; output_frequency=:monthly, _FS_OLD...), cf)
    top_mean(out) = sum(view(parent(out[:density]), 1, :)) / size(out[:density], 2)
    ρ_off = top_mean(gemb(profile, cf, mp_off))
    ρ_on = top_mean(gemb(profile, cf, initialize_parameters(; base...,
        blowing_snow_method=:Keenan2021)))
    ρ_L3 = top_mean(gemb(profile, cf, initialize_parameters(; base...,
        blowing_snow_method=:Keenan2021, drift_fetch_length=3.0)))
    @test ρ_on > ρ_off
    # A shorter fetch reworks more mass.
    @test ρ_L3 > ρ_on
end
