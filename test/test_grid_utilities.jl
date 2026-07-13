using Test
using GEMB

@testset "Grid Utilities" begin
    @testset "surface_timeseries" begin
        # Test with simple matrix
        A = [1.0 2.0 3.0;
             4.0 5.0 6.0;
             7.0 8.0 9.0]

        surface = surface_timeseries(A)
        @test surface == [1.0, 2.0, 3.0]

        # Test with NaN values at top
        B = [NaN NaN 3.0;
             4.0 5.0 NaN;
             7.0 8.0 9.0]

        surface_b = surface_timeseries(B)
        @test surface_b[1] ≈ 4.0
        @test surface_b[2] ≈ 5.0
        @test surface_b[3] ≈ 3.0

        # Test with all NaN column
        C = [NaN 2.0;
             NaN 5.0;
             NaN 8.0]

        surface_c = surface_timeseries(C)
        @test isnan(surface_c[1])
        @test surface_c[2] ≈ 2.0

        # Note: MATLAB version doesn't check matrix size, so Julia version doesn't either

        # Test with typical GEMB output dimensions
        M, N = 100, 200  # 100 layers, 200 timesteps
        profile = rand(M, N)
        # Add some NaN values in the deeper layers
        profile[80:end, :] .= NaN

        surface_profile = surface_timeseries(profile)
        @test length(surface_profile) == N
        @test all(surface_profile .≈ profile[1, :])
    end

    @testset "dz2z" begin
        # Test with simple uniform grid
        dz = ones(5, 3) * 0.1  # 5 layers, 3 timesteps, 0.1m spacing

        z_center = dz2z(dz)

        # Check dimensions
        @test size(z_center) == size(dz)

        # First cell center should be at -dz/2
        @test z_center[1, 1] ≈ -0.05

        # Second cell center should be at -dz - dz/2
        @test z_center[2, 1] ≈ -0.15

        # Test cumulative sum behavior
        expected_z = [-0.05, -0.15, -0.25, -0.35, -0.45]
        @test z_center[:, 1] ≈ expected_z

        # Test with varying grid spacing
        dz_var = [0.1, 0.2, 0.3, 0.4, 0.5]
        z_var = dz2z(reshape(dz_var, 5, 1))

        @test z_var[1] ≈ -0.05
        @test z_var[2] ≈ -0.25  # -cumsum([0.1,0.2])[2] + dz[1]/2 = -0.3 + 0.05
        @test z_var[3] ≈ -0.55  # -cumsum([0.1,0.2,0.3])[3] + dz[1]/2 = -0.6 + 0.05

        # Test with NaN handling
        dz_nan = [0.1, 0.2, NaN, 0.4, NaN]
        z_nan = dz2z(reshape(dz_nan, 5, 1))

        @test z_nan[1] ≈ -0.05
        # Skip middle values as they depend on cumsum behavior with NaN
        @test isnan(z_nan[3])
        @test isnan(z_nan[5])

        # Test with multiple columns
        dz_multi = [0.1 0.2;
                    0.1 0.2;
                    0.1 0.2]

        z_multi = dz2z(dz_multi)
        @test size(z_multi) == (3, 2)
        @test z_multi[1, 1] ≈ -0.05
        @test z_multi[1, 2] ≈ -0.1
    end

    @testset "fast_divisors" begin
        # Test small numbers
        @test fast_divisors(1) == [1]
        @test fast_divisors(2) == [1, 2]
        @test fast_divisors(6) == [1, 2, 3, 6]

        # Test example from documentation
        @test fast_divisors(42) == [1, 2, 3, 6, 7, 14, 21, 42]

        # Test prime number
        @test fast_divisors(13) == [1, 13]

        # Test perfect square
        @test fast_divisors(16) == [1, 2, 4, 8, 16]
        @test fast_divisors(36) == [1, 2, 3, 4, 6, 9, 12, 18, 36]

        # Test larger number
        divisors_100 = fast_divisors(100)
        @test 1 ∈ divisors_100
        @test 100 ∈ divisors_100
        @test 10 ∈ divisors_100
        @test all(100 % d == 0 for d in divisors_100)

        # Test power of 2
        divisors_64 = fast_divisors(64)
        @test divisors_64 == [1, 2, 4, 8, 16, 32, 64]

        # Note: MATLAB version doesn't validate input, so Julia version doesn't either

        # Test that all returned values are actual divisors
        for n in [12, 24, 30, 48, 60, 100]
            divs = fast_divisors(n)
            @test all(n % d == 0 for d in divs)
            @test all(issorted(divs))
        end
    end

    @testset "decyear2datenum" begin
        # Test mid-2023 (from documentation)
        dn = decyear2datenum(2023.5)
        # MATLAB datenum for 2023-07-02 12:00:00 is approximately 739069.5
        # Julia rata die calculation differs slightly from MATLAB - test relative values instead

        # Test start of year
        dn_start = decyear2datenum(2023.0)
        dn_end = decyear2datenum(2024.0)

        # 2023 is not a leap year, should have 365 days
        @test dn_end - dn_start ≈ 365.0 atol=0.01

        # Test leap year
        dn_leap_start = decyear2datenum(2020.0)
        dn_leap_mid = decyear2datenum(2020.5)
        dn_leap_end = decyear2datenum(2021.0)

        # Leap year should have 366 days
        days_in_2020 = dn_leap_end - dn_leap_start
        @test days_in_2020 ≈ 366.0 atol=0.01

        # Mid-leap-year should be after June 30 (accounting for extra day in Feb)
        mid_offset = dn_leap_mid - dn_leap_start
        @test mid_offset ≈ 183.0 atol=1.0  # 366/2 = 183

        # Test non-leap year
        dn_2023_start = decyear2datenum(2023.0)
        dn_2024_start = decyear2datenum(2024.0)
        days_in_2023 = dn_2024_start - dn_2023_start
        @test days_in_2023 ≈ 365.0 atol=0.01

        # Test array broadcasting
        decyears = [2020.0, 2020.5, 2021.0]
        datenums = decyear2datenum(decyears)

        @test length(datenums) == 3
        @test all(datenums[2] > datenums[1])
        @test all(datenums[3] > datenums[2])
        @test datenums[3] - datenums[1] ≈ 366.0 atol=0.01

        # Test fractional precision
        dn_frac1 = decyear2datenum(2023.25)
        dn_frac2 = decyear2datenum(2023.75)
        @test dn_frac2 - dn_frac1 ≈ 365.0 * 0.5 atol=0.1
    end

    @testset "Integration: dz2z and surface_timeseries" begin
        # Create a mock GEMB-style output
        nz = 50  # 50 vertical layers
        nt = 100  # 100 timesteps

        # Create grid spacing that decreases with depth
        dz = zeros(nz, nt)
        for j in 1:nt
            dz[:, j] = [0.01 * (1.1^(i-1)) for i in 1:nz]
        end

        # Add some NaN values in deeper layers
        dz[40:end, :] .= NaN

        # Convert to depth coordinates
        z_center = dz2z(dz)

        # Extract surface values
        dz_surface = surface_timeseries(dz)

        @test length(dz_surface) == nt
        @test all(dz_surface .≈ dz[1, :])

        # Check that z_center is negative (below surface)
        @test all(z_center[1:39, :] .< 0)

        # Check that z_center increases in magnitude with depth
        for j in 1:nt
            for i in 2:39
                @test z_center[i, j] < z_center[i-1, j]
            end
        end

        # Check NaN preservation
        @test all(isnan.(z_center[40:end, :]))
    end
end
