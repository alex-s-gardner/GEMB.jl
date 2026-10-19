# Eq. 9's grain-type gate on the Crocus viscosity (`crocus_f2_applies`). Crocus 3.0.1
# (`SNOWCROCOMPACTN`, snowcro.F90:2112-2121) raises the viscosity by eq. 9's f2 only in layers
# that are non-dendritic by its optical-diameter test, angular (sphericity below 51/99) and have
# been faceted without ever being wet. GEMB reads the present state: dendricity 0 and the same
# optical-diameter test, sphericity below 51/99, no liquid water.

function _f2_cfs(dt)
    return GEMB.ClimateForcingStep(; dt, temperature_air=263.15, pressure_air=80000.0,
        precipitation=0.0, wind_speed=0.0, vapor_pressure=100.0, temperature_air_mean=263.15,
        wind_speed_mean=0.0, precipitation_mean=300.0, temperature_observation_height=2.0,
        wind_observation_height=10.0)
end

@testset "crocus_f2_applies: the gate" begin
    s_max = GEMB.CROCUS_F2_SPHERICITY_MAX
    @test s_max == 51 / 99
    r = 0.3        # mm: gs = 0.6 mm, coarse enough for Crocus's optical-diameter test
    # Non-dendritic, angular, dry and coarse: the depth hoar eq. 9 is meant for.
    @test GEMB.crocus_f2_applies(0.0, 0.0, 0.0, r)
    @test GEMB.crocus_f2_applies(0.0, prevfloat(s_max), 0.0, r)
    # The sphericity bound is strict, as in Crocus (`PSNOWSPHERI < XVGRAN6`).
    @test !GEMB.crocus_f2_applies(0.0, s_max, 0.0, r)
    @test !GEMB.crocus_f2_applies(0.0, 1.0, 0.0, r)               # rounded grains
    # Dendritic snow, however little dendricity is left; the boundary is the one
    # `calculate_grain_size` classifies cells by.
    @test !GEMB.crocus_f2_applies(1.0, 0.0, 0.0, r)
    @test !GEMB.crocus_f2_applies(2 * GEMB.GDN_TOLERANCE, 0.0, 0.0, r)
    @test GEMB.crocus_f2_applies(GEMB.GDN_TOLERANCE, 0.0, 0.0, r)
    # Liquid water closes it, at the tolerance every other wet/dry test in GEMB uses.
    @test !GEMB.crocus_f2_applies(0.0, 0.0, 2 * GEMB.WATER_TOLERANCE, r)
    @test GEMB.crocus_f2_applies(0.0, 0.0, GEMB.WATER_TOLERANCE, r)
    # Crocus's non-dendritic test on the optical diameter, d = g·S + (1 − S)·(4e-4 + g)/2 >=
    # 1e-4·(4 − S) with g = 2r: it needs g >= 2e-4·(2 + S)/(1 + S), i.e. r >= 0.2 mm at S = 0 and
    # r >= 1/6 mm at S = 0.5.
    @test GEMB.crocus_f2_applies(0.0, 0.0, 0.0, 0.2)
    @test !GEMB.crocus_f2_applies(0.0, 0.0, 0.0, prevfloat(0.2))
    @test !GEMB.crocus_f2_applies(0.0, 0.0, 0.0, 0.15)            # fine faceted grains
    @test GEMB.crocus_f2_applies(0.0, 0.5, 0.0, 0.17)
    @test !GEMB.crocus_f2_applies(0.0, 0.5, 0.0, 0.16)
    # Inside the gate eq. 9 is therefore between e^1.32 = 3.74 and its cap of 4.
    gs_min = 2e-4 * (2 + s_max) / (1 + s_max)      # smallest admitted grain size 2r [m]
    @test exp((gs_min - 2e-4) / 1e-4) > 3.7
end

@testset "_crocus_viscosity: f2 only where the gate is open" begin
    ρ, T, D = 300.0, 263.15, 0.05
    base = 7.62237e6 * (ρ / 250.0) * exp(0.1 * (GEMB.CtoK - T) + 0.023 * ρ)
    # Closed: eq. 7 with f1, whatever the grain size.
    for r in (0.01, 0.05, 0.15, 0.3, 2.5)
        @test GEMB._crocus_viscosity(ρ, T, 0.0, D, r, false) ≈ base rtol = 1e-14
    end
    W = 1.0
    f1 = 1 / (1 + 60 * W / (1000 * D))
    @test GEMB._crocus_viscosity(ρ, T, W, D, 0.3, false) ≈ f1 * base rtol = 1e-14
    # Open: eq. 9, e at gs = 0.3 mm (radius 0.15 mm) and the cap of 4 by gs = 0.339 mm.
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.15, true) ≈ exp(1.0) * base rtol = 1e-14
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.3, true) ≈ 4.0 * base rtol = 1e-14
    # The default is open, so a bare five-argument call is eq. 9 as published (floored at 1).
    @test GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.3) == GEMB._crocus_viscosity(ρ, T, 0.0, D, 0.3, true)
end

@testset "calculate_density: one step against the closed form, per grain type" begin
    # Four 0.05 m cells of 0.3 mm radius (gs = 0.6 mm, so eq. 9 would give its cap of 4):
    #   1 dendritic (d 0.4, s 0.6)            gate closed, f2 = 1
    #   2 non-dendritic rounded (s 0.9)       gate closed, f2 = 1
    #   3 non-dendritic faceted (s 0.2), dry  gate open,   f2 = 4
    #   4 non-dendritic faceted, wet          gate closed, f2 = 1, f1 < 1
    dz = fill(0.05, 4)
    density = [200.0, 250.0, 300.0, 350.0]
    T = fill(263.15, 4)
    r = fill(0.3, 4)
    water = [0.0, 0.0, 0.0, 0.5]
    dend = [0.4, 0.0, 0.0, 0.0]
    sph = [0.6, 0.9, 0.2, 0.2]
    dt = 86400.0
    cfs = _f2_cfs(dt)
    mp = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)
    visc = fill(-1.0, 4)
    (dz_out, d_out) = GEMB.calculate_density(T, copy(dz), copy(density), r, water, cfs, mp;
        grain_dendricity=dend, grain_sphericity=sph, viscosity=visc)

    setprecision(BigFloat, 256) do
        f2 = big.([1, 1, 4, 1])
        load = big"0"
        for i in 1:4
            ρ, D, W = big(density[i]), big(dz[i]), big(water[i])
            η = f2[i] / (1 + 60 * W / (1000 * D)) * big"7.62237e6" * (ρ / 250) *
                exp(big"0.1" * (big"273.15" - big(T[i])) + big"0.023" * ρ)
            σ = (load + (ρ * D + W) / 2) * big"9.81"
            dz_e = D * exp(-σ / η * big(dt))
            @test dz_out[i] ≈ Float64(dz_e) rtol = 1e-12
            @test d_out[i] ≈ Float64(ρ * D / dz_e) rtol = 1e-12
            @test visc[i] ≈ Float64(η) rtol = 1e-12
            load += ρ * D + W
        end
    end
    # Only the faceted dry cell carries the factor of 4 in its viscosity.
    η_closed = [GEMB._crocus_viscosity(density[i], T[i], water[i], dz[i], r[i], false) for i in 1:4]
    @test visc ./ η_closed ≈ [1.0, 1.0, 4.0, 1.0] rtol = 1e-14
    @test dz_out .* d_out ≈ dz .* density rtol = 1e-14
end

@testset "calculate_density: the gate under :Crocus, and what the other schemes read" begin
    n = 3
    dz = fill(0.05, n)
    density = [200.0, 300.0, 400.0]
    T = fill(263.15, n)
    r = fill(0.3, n)
    cfs = _f2_cfs(86400.0)
    mp_c = GEMB.ModelParameters(density_ice=910.0, densification_method=:Crocus)
    mp_p = GEMB.ModelParameters(density_ice=910.0, densification_method=:CrocusPure)

    # Below the 450 handover `:Crocus` is the pure law, gate included.
    for (dend, sph) in ((zeros(n), zeros(n)), (zeros(n), ones(n)), (ones(n), fill(0.5, n)))
        a = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp_c;
            grain_dendricity=dend, grain_sphericity=sph)
        b = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp_p;
            grain_dendricity=dend, grain_sphericity=sph)
        @test a == b
    end
    # Rounded snow of the same size compacts faster than faceted snow: the gate is what
    # separates them now that the radius is the same.
    (_, d_fac) = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp_c;
        grain_dendricity=zeros(n), grain_sphericity=zeros(n))
    (_, d_rnd) = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp_c;
        grain_dendricity=zeros(n), grain_sphericity=ones(n))
    @test all(d_rnd[2:n] .> d_fac[2:n])
    @test all(d_rnd .>= d_fac)

    # The Crocus schemes fail fast without grain shape...
    for mp in (mp_c, mp_p)
        @test_throws "needs grain_dendricity and grain_sphericity" GEMB.calculate_density(
            T, copy(dz), copy(density), r, zeros(n), cfs, mp)
        @test_throws "needs grain_dendricity and grain_sphericity" GEMB.calculate_density(
            T, copy(dz), copy(density), r, zeros(n), cfs, mp;
            grain_dendricity=zeros(n - 1), grain_sphericity=zeros(n))
    end
    # ...and no other scheme reads it: with and without, bit-identical.
    for meth in (:Arthern, :ArthernB, :GSFC2020, :HerronLangway, :Barnola1991)
        mp = GEMB.ModelParameters(density_ice=910.0, densification_method=meth)
        a = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp)
        b = GEMB.calculate_density(T, copy(dz), copy(density), r, zeros(n), cfs, mp;
            grain_dendricity=zeros(n), grain_sphericity=zeros(n))
        @test a == b
    end
end
