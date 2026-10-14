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

    # Every method applies the diffuse albedo to the diffuse stream.
    expected_net = (1 - albedo_broadband) * 150.0 + (1 - albedo_diffuse) * 50.0

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

    expected = (1 - albedo_broadband) * 150.0 + (1 - albedo_diffuse) * 50.0

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
    expected_total = (1 - albedo_broadband) * 150.0 + (1 - albedo_diffuse) * 50.0
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

@testset "Absorbed shortwave equals the surface net for every method" begin
    # Deep, graded column so nothing is transmitted through the base. For each albedo method,
    # with and without subsurface absorption, the column absorbs exactly
    # (1 − α)(SW↓ − SW↓_d) + (1 − α_d)·SW↓_d with the albedos `calculate_albedo` returned —
    # including over thin snow on ice, where those differ from the grain-size albedo.
    n = 400
    dz = fill(0.1, n)
    cfs = GEMB.ClimateForcingStep(; dt=3600.0, shortwave_downward=500.0,
        shortwave_downward_diffuse=200.0, solar_zenith_angle=40.0,
        cloud_optical_thickness=1.0, cloud_fraction=0.5)
    profiles = (collect(range(300.0, 917.0; length=n)),
        [fill(350.0, 1); fill(917.0, n - 1)])
    dz_profiles = (dz, [0.04; fill(0.1, n - 1)])
    for method in (:GardnerSharp, :BrunLefebre, :GreuellKonzelmann), sub in (false, true),
        (density, dzp) in zip(profiles, dz_profiles)
        mp = GEMB.ModelParameters(albedo_method=method, shortwave_subsurface_absorption=sub,
            albedo_ice=0.35)
        grain_radius = fill(0.5, n)
        α, α_d = GEMB.calculate_albedo(dzp, density, grain_radius, cfs, mp)
        sw = GEMB.calculate_shortwave_radiation(dzp, density, grain_radius, α, α_d, cfs, mp)
        @test sum(sw) ≈ (1 - α) * 300.0 + (1 - α_d) * 200.0 rtol = 1e-10
        @test GEMB.effective_albedo(α, α_d, cfs) ≈ 1 - sum(sw) / 500.0 rtol = 1e-10
    end
end

@testset "Brun-Lefebre spectral albedo" begin
    sF = (0.606, 0.301, 0.093)
    @test sum(sF) ≈ 1.0 rtol = 1e-15
    for r in (0.05, 0.3, 1.0, 2.5)
        gsz = 2r / 1000
        a = (min(0.98, 0.95 - 1.58 * sqrt(gsz)), max(0.0, 0.95 - 15.4 * sqrt(gsz)),
            max(0.127, 0.88 + 346.3 * gsz - 32.31 * sqrt(gsz)))
        @test GEMB._albedo_brun_lefebre(r) ≈ sum(sF .* a) rtol = 1e-12
    end
end

@testset "Density-branch extinction (Beer-Lambert, uniform B)" begin
    # ρ = 300 everywhere gives B = 10 m⁻¹. Cell i absorbs 0.64·SW_net·(e^{−Bz_{i−1}} − e^{−Bz_i}),
    # cell 1 additionally 0.36·SW_net, and the rest, 0.64·SW_net·e^{−Bz_m}, leaves the base.
    n = 12
    dz = [0.01 * 1.3^(i - 1) for i in 1:n]
    density = fill(300.0, n)
    mp = GEMB.ModelParameters(albedo_method=:GreuellKonzelmann,
        shortwave_subsurface_absorption=true)
    α = 0.6
    for sw_d in (0.0, 120.0)
        cfs = GEMB.ClimateForcingStep(; dt=3600.0, shortwave_downward=400.0,
            shortwave_downward_diffuse=sw_d)
        sw = GEMB.calculate_shortwave_radiation(dz, density, fill(0.5, n), α, α, cfs, mp)
        net = (1 - α) * 400.0
        z = cumsum(dz)
        z_top = [0.0; z[1:end-1]]
        expected = 0.64 * net .* (exp.(-10.0 .* z_top) .- exp.(-10.0 .* z))
        expected[1] += 0.36 * net
        @test sw ≈ expected rtol = 1e-12
        @test sum(sw) + 0.64 * net * exp(-10.0 * z[end]) ≈ net rtol = 1e-12
    end
end
