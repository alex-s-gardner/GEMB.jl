# Tests for calculate_grain_size - translated from MATLAB test_calculate_grain_size.m

# Helper to create a ClimateForcingStep for grain size tests
function _make_grain_cfs(; dt=86400.0)
    return GEMB.ClimateForcingStep(
        dt,               # dt [s]
        265.0,            # temperature_air
        100000.0,         # pressure_air
        0.0,              # precipitation
        5.0,              # wind_speed
        200.0,            # shortwave_downward
        300.0,            # longwave_downward
        400.0,            # vapor_pressure
        260.0,            # temperature_air_mean
        5.0,              # wind_speed_mean
        200.0,            # precipitation_mean
        2.0,              # temperature_observation_height
        10.0,             # wind_observation_height
        0.0,              # black_carbon_snow
        0.0,              # black_carbon_ice
        0.0,              # cloud_optical_thickness
        0.0,              # solar_zenith_angle
        0.0,              # shortwave_downward_diffuse
        0.1,              # cloud_fraction
    )
end

@testset "Albedo method skip" begin
    n = 5
    temperature = 260.0 * ones(n)
    dz = 0.1 * ones(n)
    density = 300.0 * ones(n)
    water = zeros(n)
    grain_radius = 0.5 * ones(n)
    grain_dendricity = 0.5 * ones(n)
    grain_sphericity = 0.5 * ones(n)
    cfs = _make_grain_cfs()

    # Method "None" should skip grain evolution
    mp_none = GEMB.ModelParameters(albedo_method=:None)
    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp_none)

    @test gs_out == grain_radius
    @test gdn_out == grain_dendricity
    @test gsp_out == grain_sphericity

    # Method "GreuellKonzelmann" should also skip
    mp_gk = GEMB.ModelParameters(albedo_method=:GreuellKonzelmann)
    (gs_out2, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp_gk)

    @test gs_out2 == grain_radius
end

@testset "Dendritic dry low gradient" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [260.0, 260.0, 260.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)
    gr_before = copy(grain_radius)

    (gs_out, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Dendricity should decrease (decay)
    @test all(gdn_out .< gdn_before)
    # Sphericity should increase
    @test all(gsp_out .> gsp_before)
    # Grain radius should change
    @test gs_out != gr_before
end

@testset "Dendritic dry high gradient" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    # High temperature gradient (> 5 K/m)
    temperature = [260.0, 270.0, 280.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Under high gradient: dendricity decreases, sphericity decreases
    @test all(gdn_out .< gdn_before)
    @test all(gsp_out .< gsp_before)
end

@testset "Dendritic wet snow" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [250.0, 250.0, 250.0]
    water = [1.0, 1.0, 1.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.8, 0.8, 0.8]
    grain_sphericity = [0.2, 0.2, 0.2]

    gdn_before = copy(grain_dendricity)
    gsp_before = copy(grain_sphericity)

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Wet snow causes rapid rounding
    @test all(gdn_out .< gdn_before)
    @test all(gsp_out .> gsp_before)
end

@testset "Nondendritic dry (Marbouty)" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    # Moderate gradient (~20 K/m)
    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]  # must be < 400 for growth
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth in grain size
    @test all(gs_out .> grain_radius)
    # Dendricity stays at 0
    @test gdn_out == grain_dendricity
end

@testset "Marbouty density limit" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [250.0, 252.0, 254.0]
    dz = [0.1, 0.1, 0.1]
    density = [450.0, 450.0, 450.0]  # > 400 threshold
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # No growth expected above density threshold
    @test gs_out ≈ grain_radius atol = 1e-10
end

@testset "Nondendritic wet snow (Brun)" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]
    water = [1.5, 1.5, 1.5]
    grain_radius = [0.5, 0.5, 0.5]
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [0.5, 0.5, 0.5]

    (gs_out, gdn_out, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Expect growth via wet snow mechanism
    @test all(gs_out .> grain_radius)
    # Dendricity stays at 0
    @test gdn_out == grain_dendricity
end

@testset "Clamping limits" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    # Large dt to force dendricity toward 0
    cfs = _make_grain_cfs(dt=86400.0 * 100)

    temperature = [260.0, 260.0, 260.0]
    dz = [0.1, 0.1, 0.1]
    density = [200.0, 200.0, 200.0]
    water = [0.0, 0.0, 0.0]
    grain_radius = [0.2, 0.2, 0.2]
    grain_dendricity = [0.1, 0.1, 0.1]  # close to 0
    grain_sphericity = [0.2, 0.2, 0.2]

    (_, gdn_out, gsp_out) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    @test all(gdn_out .>= 0.0)
    @test all(gsp_out .<= 1.0)
end

@testset "Grain size cap" begin
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0 * 50)

    temperature = [GEMB.CtoK, GEMB.CtoK, GEMB.CtoK]
    dz = [0.1, 0.1, 0.1]
    density = [300.0, 300.0, 300.0]
    water = [2.0, 2.0, 2.0]
    grain_radius = [1.9, 1.9, 1.9]  # radius 1.9mm
    grain_dendricity = [0.0, 0.0, 0.0]
    grain_sphericity = [1.0, 1.0, 1.0]

    (gs_out, _, _) = GEMB.calculate_grain_size(
        temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, cfs, mp)

    # Radius should be capped at 1.0mm for spherical grains
    @test all(gs_out .<= 1.0 + 1e-10)
end

# ============================================================================
# PHYSICS-CORRECTNESS TESTS
#
# The tests above check qualitative behavior (increases / decreases). The tests
# below pin the EXACT published coefficients so any change that breaks fidelity
# to the source equations fails immediately.
#
# References (verified against Vionnet et al. 2012, GMD 5, 773-791 — the
# open-access CROCUS description that reproduces these laws — and the
# authoritative SURFEX/Crocus Fortran source SNOW3L_MARBOUTY):
#   - Brun et al. (1992), J. Glaciol., 38, 13-22   : dendricity/sphericity laws, size caps
#   - Brun (1989), Ann. Glaciol., 13, 22-26        : wet-snow grain volume growth
#   - Marbouty (1980), J. Glaciol., 26, 303-312    : dry-snow temperature-gradient growth
#
# Unit conventions (Vionnet 2012 Tables 1-2 & symbol table):
#   temperature T in KELVIN (exp(-6000/T)); gradient in K m^-1; time in DAYS;
#   liquid-water content theta as MASS PERCENT. The Brun (1989) volume-growth
#   rate is the sole exception: it is per SECOND (mm^3 s^-1).
# ============================================================================

@testset "PHYSICS: dendritic dry low-gradient coefficients (Brun 1992)" begin
    # Table 1, dendritic, G <= 5 K/m:
    #   d d/dt = -2e8 * exp(-6000/T)   [per day]
    #   d s/dt = +1e9 * exp(-6000/T)   [per day]
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    T = 260.0
    dt = 100.0                 # small dt so values stay away from the 0/1 clamps
    dt_days = dt / 86400.0
    cfs = _make_grain_cfs(dt=dt)

    # single cell => temperature gradient is exactly 0 (< 5 K/m branch)
    d0, s0 = 0.8, 0.2
    (_, gdn, gsp) = GEMB.calculate_grain_size(
        [T], [0.1], [200.0], [0.0], [0.25], [d0], [s0], cfs, mp)

    ex = exp(-6000.0 / T)
    @test gdn[1] ≈ d0 + (-2e8 * ex * dt_days) rtol = 1e-13
    @test gsp[1] ≈ s0 + ( 1e9 * ex * dt_days) rtol = 1e-13
end

@testset "PHYSICS: dendritic dry high-gradient coefficient (Brun 1992)" begin
    # Table 1, dendritic, G > 5 K/m:
    #   d d/dt = d s/dt = -2e8 * exp(-6000/T) * G^0.4   [per day]
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    dt = 10.0
    dt_days = dt / 86400.0
    cfs = _make_grain_cfs(dt=dt)

    # 3 cells, linear T; gradient at cell 1 = (T3-T1)/(dz1/2 + dz2 + dz3/2)
    T = [260.0, 270.0, 280.0]
    dz = [0.1, 0.1, 0.1]
    G1 = (T[3] - T[1]) / (dz[1]/2 + dz[2] + dz[3]/2)   # = 100 K/m
    @test G1 > 5

    d0, s0 = 0.8, 0.5
    (_, gdn, gsp) = GEMB.calculate_grain_size(
        T, dz, fill(200.0, 3), zeros(3), fill(0.25, 3),
        fill(d0, 3), fill(s0, 3), cfs, mp)

    C1 = (-2e8 * exp(-6000.0 / T[1]) * dt_days) * G1^0.4
    @test gdn[1] ≈ d0 + C1 rtol = 1e-12
    @test gsp[1] ≈ s0 + C1 rtol = 1e-12   # same coefficient applied to both
end

@testset "PHYSICS: dendritic wet-snow coefficient + 9% LWC mass cap (Brun 1992/1980)" begin
    # Table 2, dendritic wet:  d d/dt = -(1/16) theta^3,  d s/dt = +(1/16) theta^3
    # theta = liquid water content as MASS % = water / (density*dz) * 100,
    # capped at 9% by mass (Brun 1980).
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    dt = 100.0
    dt_days = dt / 86400.0
    cfs = _make_grain_cfs(dt=dt)

    # theta_raw = 2.0/(100*0.1)*100 = 20% -> capped to 9%
    d0, s0 = 0.8, 0.2
    (_, gdn, gsp) = GEMB.calculate_grain_size(
        [GEMB.CtoK], [0.1], [100.0], [2.0], [0.25], [d0], [s0], cfs, mp)

    D = (1.0 / 16.0) * 9.0^3 * dt_days     # uses capped theta = 9, NOT 20
    @test gdn[1] ≈ d0 - D rtol = 1e-12
    @test gsp[1] ≈ s0 + D rtol = 1e-12

    # theta = 9% exactly is the boundary and must NOT be treated as > cap
    # (water chosen so theta = 9): water/(density*dz)*100 = 9 -> water = 0.9
    (_, gdn2, _) = GEMB.calculate_grain_size(
        [GEMB.CtoK], [0.1], [100.0], [0.9], [0.25], [d0], [s0], cfs, mp)
    @test gdn2[1] ≈ d0 - (1.0/16.0) * 9.0^3 * dt_days rtol = 1e-12
end

@testset "PHYSICS: wet-snow grain volume growth (Brun 1989)" begin
    # Table 2, nondendritic spherical:  dV/dt = v0 + v1*theta^3  [mm^3 s^-1]
    #   v0 = 1.28e-8 mm^3/s,  v1 = 4.22e-10 mm^3/s,  theta in mass %.
    # New diameter from spherical-volume conservation:
    #   V_new = (4/3) pi (D0/2)^3 + dV/dt * dt_seconds
    #   D_new = 2 * (3/(4 pi) * V_new)^(1/3)
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    dt = 86400.0                     # 1 day
    cfs = _make_grain_cfs(dt=dt)

    r0 = 0.5                         # radius mm -> diameter 1.0 mm
    density, dz, water = 300.0, 0.1, 1.5
    theta = water / (density * dz) * 100.0     # = 5% mass, below the 9% cap

    (gs, _, _) = GEMB.calculate_grain_size(
        [GEMB.CtoK], [dz], [density], [water], [r0], [0.0], [1.0], cfs, mp)

    dV = (1.28e-8 + 4.22e-10 * theta^3) * dt   # mm^3 (dt in seconds)
    V0 = (4.0/3.0) * pi * (1.0/2)^3
    D_new = 2.0 * (3.0 / (4pi) * (V0 + dV))^(1/3)
    @test gs[1] ≈ D_new / 2 rtol = 1e-13
end

@testset "PHYSICS: effective grain size mapping (Vionnet 2012 Eq. 13)" begin
    # Dendritic optical size: gsz[mm] = 0.1*(d + (1-d)*(3s + 4(1-s)))  (= 0.1*(d+(1-d)*(4-s)))
    # GEMB normalizes d,s by 0.99. Use dt=0 so d,s are unchanged and only the
    # mapping is exercised. Returned value is radius = gsz/2.
    #
    # This mapping applies ONLY to dendritic cells (d > 0); for d == 0 the grain
    # size comes from the nondendritic growth path instead, so only d > 0 is
    # tested here.
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=0.0)

    for (d, s) in [(1.0, 0.0), (0.3, 1.0), (0.5, 0.5), (1.0, 1.0)]
        (gs, _, _) = GEMB.calculate_grain_size(
            [260.0], [0.1], [200.0], [0.0], [0.25], [d], [s], cfs, mp)
        dn, sn = d / 0.99, s / 0.99
        gsz = 0.1 * (dn + (1.0 - dn) * (sn * 3.0 + (1.0 - sn) * 4.0))
        @test gs[1] ≈ gsz / 2 rtol = 1e-13
    end
end

@testset "PHYSICS: Marbouty temperature coefficient F (source-code form)" begin
    # Marbouty (1980) Fig. 9 temperature factor, isolated by rho<150 (H=1) and
    # G=1 (dT >= 70 K/m). Q = F*H*G*E with E = 0.09 mm/day.
    # Endpoints: F(0C)=0.7, F(-6C)=1.0, F(-22C)=0.2, F(-40C)=0, F(<-40C)=0.
    #
    # NOTE: F(0C)=0.7 is the CORRECT value from the authoritative SURFEX source.
    # The published Vionnet (2012) Eq. B2 warm branch is a typo that would give
    # 1.3; GEMB (and this test) use the correct 0.7. Do not "fix" to 1.3.
    E = 0.09
    for (Tc, Fexp) in [(0.0, 0.7), (-6.0, 1.0), (-22.0, 0.2), (-40.0, 0.0), (-50.0, 0.0)]
        Q = GEMB._marbouty_Q(Tc + 273.15, 100.0, 100.0)   # H=1 (rho<150), G=1 (dT=100)
        @test Q / E ≈ Fexp atol = 1e-12
    end
end

@testset "PHYSICS: Marbouty density coefficient H (Brun 1992 / Marbouty 1980)" begin
    # Density factor: H=1 for rho<150, linear H=1-(rho-150)/250 to 0 at 400,
    # H=0 for rho>=400 (no depth-hoar growth in dense firn). Isolate with
    # Tc=0 (F=0.7) and G=1 (dT>=70).
    E = 0.09
    F0 = 0.7
    for (rho, Hexp) in [(100.0, 1.0), (150.0, 1.0), (275.0, 0.5), (400.0, 0.0), (450.0, 0.0)]
        Q = GEMB._marbouty_Q(273.15, rho, 100.0)
        @test Q / (E * F0) ≈ Hexp atol = 1e-12
    end
end

@testset "PHYSICS: Marbouty growth constant E = 0.09 mm/day (diameter)" begin
    # With F=H=G=1 the dry-snow grain-size (diameter) growth rate is exactly the
    # Marbouty growth constant E = 0.09 mm/day. Verify through the full routine:
    # nondendritic (d=0), dry, rho<150 (H=1), Tc=-6 (F=1), high gradient (G=1).
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    cfs = _make_grain_cfs(dt=86400.0)          # 1 day

    # 3 cells centered at Tc=-6C with a strong gradient so cell 2 sees G=1.
    T = [273.15 - 16.0, 273.15 - 6.0, 273.15 + 4.0]   # cell-2 gradient = 20/0.2 = 100 K/m
    r0 = 0.25
    (gs, gdn, _) = GEMB.calculate_grain_size(
        T, [0.1, 0.1, 0.1], fill(100.0, 3), zeros(3),
        fill(r0, 3), zeros(3), zeros(3), cfs, mp)

    gsz0 = r0 * 2
    r_expected = (gsz0 + 0.09 * 1.0) / 2       # diameter grows by E*1 day, then /2
    @test gs[2] ≈ r_expected rtol = 1e-12
    @test gdn == zeros(3)                       # dendricity stays 0 (nondendritic)
end

@testset "PHYSICS: grain-size diameter caps 2mm / 5mm (Brun 1992)" begin
    # Spherical grains (sphericity == 1): diameter <= 2 mm  => radius <= 1.0 mm
    # Non-spherical grains (sphericity < 1): diameter <= 5 mm => radius <= 2.5 mm
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)

    # Spherical: drive large wet growth, sphericity = 1
    cfs_long = _make_grain_cfs(dt=86400.0 * 50)
    (gs_sph, _, _) = GEMB.calculate_grain_size(
        [GEMB.CtoK], [0.1], [300.0], [2.0], [1.9], [0.0], [1.0], cfs_long, mp)
    @test gs_sph[1] ≈ 1.0 atol = 1e-10          # capped exactly at 2mm diameter

    # Non-spherical: drive very large dry growth, sphericity = 0
    cfs_vlong = _make_grain_cfs(dt=86400.0 * 1000)
    T = [273.15 - 21.0, 273.15 - 6.0, 273.15 + 9.0]   # strong gradient, Tc~-6 at cell 2
    (gs_ns, _, _) = GEMB.calculate_grain_size(
        T, [0.1, 0.1, 0.1], fill(100.0, 3), zeros(3),
        fill(3.0, 3), zeros(3), zeros(3), cfs_vlong, mp)
    @test all(gs_ns .<= 2.5 + 1e-10)            # capped at 5mm diameter
    @test gs_ns[2] ≈ 2.5 atol = 1e-10           # cell 2 actually reaches the cap
end

@testset "PHYSICS: temperature-gradient units (K/m) and central-difference stencil" begin
    # The gradient must be computed in K/m using grid-center spacing
    #   dz[i-1]/2 + dz[i] + dz[i+1]/2  (interior central difference).
    # A gradient of exactly 5 K/m is the low/high branch boundary. Build a case
    # straddling it and confirm the low-gradient (positive sphericity) law fires
    # rather than the high-gradient (negative) one — this pins the K/m unit
    # convention (a K/cm misread would be 100x off and flip the branch).
    mp = GEMB.ModelParameters(albedo_method=:GardnerSharp)
    dt = 100.0
    dt_days = dt / 86400.0
    cfs = _make_grain_cfs(dt=dt)

    # Interior cell gradient = (T3-T1)/(dz1/2+dz2+dz3/2). With dz=0.1 and
    # T=[259.8,260.0,260.2], gradient = 0.4/0.2 = 2 K/m  (< 5 -> low branch).
    T = [259.8, 260.0, 260.2]
    dz = [0.1, 0.1, 0.1]
    Gmid = (T[3] - T[1]) / (dz[1]/2 + dz[2] + dz[3]/2)
    @test Gmid ≈ 2.0 rtol = 1e-13
    @test Gmid < 5      # low-gradient branch

    d0, s0 = 0.8, 0.2
    (_, gdn, gsp) = GEMB.calculate_grain_size(
        T, dz, fill(200.0, 3), zeros(3), fill(0.25, 3),
        fill(d0, 3), fill(s0, 3), cfs, mp)

    # Low-gradient dendritic law at cell 2: sphericity INCREASES by +1e9*exp(-6000/T)*dt_days
    ex = exp(-6000.0 / T[2])
    @test gsp[2] ≈ s0 + 1e9 * ex * dt_days rtol = 1e-12
    @test gsp[2] > s0                    # low-gradient => sphericity up (would be down if >5)
    @test gdn[2] ≈ d0 - 2e8 * ex * dt_days rtol = 1e-12
end
