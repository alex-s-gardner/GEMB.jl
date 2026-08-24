# Tests for calculate_shortwave_radiation

# Helper to create a ClimateForcingStep with specified shortwave fields
function _make_sw_cfs(; shortwave_downward=200.0, shortwave_downward_diffuse=50.0)
    return GEMB.ClimateForcingStep(
        86400.0,    # dt
        260.0,      # temperature_air
        80000.0,    # pressure_air
        0.0,        # precipitation
        5.0,        # wind_speed
        shortwave_downward,        # shortwave_downward
        250.0,      # longwave_downward
        300.0,      # vapor_pressure
        260.0,      # temperature_air_mean
        5.0,        # wind_speed_mean
        200.0,      # precipitation_mean
        2.0,        # temperature_observation_height
        10.0,       # wind_observation_height
        0.0,        # black_carbon_snow
        0.0,        # black_carbon_ice
        0.0,        # cloud_optical_thickness
        0.0,        # solar_zenith_angle
        shortwave_downward_diffuse,  # shortwave_downward_diffuse
        0.1,        # cloud_fraction
    )
end

@testset "Surface absorption basic (no penetration)" begin
    # No penetration (shortwave_subsurface_absorption = false)
    # Albedo method != "GardnerSharp" (standard calculation)
    n = 5
    dz = 0.1 * ones(n)
    density = 350.0 * ones(n)
    grain_radius = 0.5 * ones(n)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=false,
        albedo_method=:GreuellKonzelmann,
    )
    cfs = _make_sw_cfs()

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    expected_net = (1 - albedo_broadband) * 200.0

    @test shortwave_flux[1] ≈ expected_net atol = 1e-6
    @test sum(shortwave_flux[2:end]) ≈ 0.0 atol = 1e-6
end

@testset "Surface absorption GardnerSharp method" begin
    # No penetration, GardnerSharp albedo method
    # Separates diffuse and direct components
    n = 5
    dz = 0.1 * ones(n)
    density = 350.0 * ones(n)
    grain_radius = 0.5 * ones(n)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=false,
        albedo_method=:GardnerSharp,
    )
    cfs = _make_sw_cfs()

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    dsw_direct = 200.0 - 50.0
    expected = (1 - albedo_broadband) * dsw_direct + (1 - albedo_diffuse) * 50.0

    @test shortwave_flux[1] ≈ expected atol = 1e-6
end

@testset "Density override (ice at surface reverts to surface absorption)" begin
    # Penetration requested but surface is ICE -> revert to surface absorption
    n = 5
    dz = 0.1 * ones(n)
    density = 350.0 * ones(n)
    density[1] = 917.0  # Ice density at top
    grain_radius = 0.5 * ones(n)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=true,
        albedo_method=:GreuellKonzelmann,
    )
    cfs = _make_sw_cfs()

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    expected = (1 - albedo_broadband) * 200.0

    @test shortwave_flux[1] ≈ expected atol = 1e-6
    @test sum(shortwave_flux[2:end]) ≈ 0.0 atol = 1e-6
end

@testset "Penetration with BrunLefebre method" begin
    # Subsurface penetration with spectral BrunLefebre method
    n = 5
    dz = 0.1 * ones(n)
    density = 350.0 * ones(n)
    grain_radius = 0.5 * ones(n)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=true,
        albedo_method=:BrunLefebre,
    )
    cfs = _make_sw_cfs()

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    # Energy must penetrate below surface
    @test shortwave_flux[2] > 0.0

    # Absorption should decrease with depth
    @test shortwave_flux[1] > shortwave_flux[2]
    @test shortwave_flux[2] > shortwave_flux[3]

    # Total absorption should be less than incoming (albedo reflects some)
    @test sum(shortwave_flux) < 200.0

    # Should absorb a significant fraction
    @test sum(shortwave_flux) > 0.05 * 200.0
end

@testset "Penetration with standard method (energy conservation)" begin
    # Subsurface penetration with GreuellKonzelmann (density-dependent extinction)
    # Use a deep column to prevent flux from escaping the bottom
    n_deep = 100
    dz = 0.1 * ones(n_deep)
    density = 350.0 * ones(n_deep)
    grain_radius = 0.5 * ones(n_deep)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=true,
        albedo_method=:GreuellKonzelmann,
    )
    cfs = _make_sw_cfs()

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    # Conservation: sum of absorbed flux must equal surface net flux
    expected_total = (1 - albedo_broadband) * 200.0
    @test sum(shortwave_flux) ≈ expected_total atol = 1e-4

    # Top cell should absorb NIR band + penetrating UV/Vis part
    @test shortwave_flux[1] > 0.36 * expected_total

    # Energy should penetrate to second layer
    @test shortwave_flux[2] > 0.0
end

@testset "Zero flux (nighttime)" begin
    # Zero incoming SW should give zero absorption everywhere
    n = 5
    dz = 0.1 * ones(n)
    density = 350.0 * ones(n)
    grain_radius = 0.5 * ones(n)
    albedo_broadband = 0.7
    albedo_diffuse = 0.8

    mp = GEMB.ModelParameters(
        density_ice=917.0,
        shortwave_subsurface_absorption=true,
        albedo_method=:GreuellKonzelmann,
    )
    cfs = _make_sw_cfs(shortwave_downward=0.0, shortwave_downward_diffuse=0.0)

    shortwave_flux = GEMB.calculate_shortwave_radiation(dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

    @test sum(shortwave_flux) ≈ 0.0 atol = 1e-10
end


@testset "Cumulative extinction profile" begin
    # `_extinction_profile!` replaces `cumprod`, which cannot be statically compiled (see
    # `capi/README.md`). The replacement is only correct if it agrees to the last bit — a
    # tolerance here would let a real change in the absorption profile through.
    for n in (2, 3, 7, 20, 61)
        B = 200.0 .* rand(n)
        dz = 0.05 .+ rand(n)

        out = ones(n + 1)
        GEMB._extinction_profile!(out, B, dz)

        @test out == vcat([1.0], cumprod(exp.(-B .* dz)))
        @test out[1] == 1.0
        # Transmission decreases monotonically with depth and stays a fraction.
        @test all(0.0 .<= out .<= 1.0)
        @test issorted(out; rev=true)
    end

    # A uniform column has a closed form: transmission is exp(-B*dz) per cell.
    n = 10
    B = fill(3.0, n)
    dz = fill(0.25, n)
    out = ones(n + 1)
    GEMB._extinction_profile!(out, B, dz)
    for i in 1:n
        @test out[i+1] ≈ exp(-3.0 * 0.25 * i) rtol = 1e-14
    end

    # Mismatched lengths are a caller error, not a silently truncated profile.
    @test_throws DimensionMismatch GEMB._extinction_profile!(ones(5), rand(4), rand(3))
    @test_throws DimensionMismatch GEMB._extinction_profile!(ones(3), rand(4), rand(4))
end
