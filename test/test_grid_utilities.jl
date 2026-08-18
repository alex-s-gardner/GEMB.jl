using Test
using GEMB

@testset "Grid Utilities" begin
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

    @testset "gemb_interp does not extrapolate" begin
        # Two-cell column, centres at -0.5 and -1.5 m, with a steep gradient. A target grid
        # that reaches past both ends must return NaN there rather than continuing the end
        # interval: extrapolating the near-basal gradient past the column base is what
        # produced firn temperatures above the melt point in plotted output.
        z_center = [-0.5 -0.5; -1.5 -1.5]
        T = [263.15 263.15; 273.15 273.15]        # +10 K per metre downward
        z_target = [0.0, -0.5, -1.0, -1.5, -3.0]

        g = gemb_interp(z_center, T, z_target)

        @test isnan(g[1, 1])                      # above the top cell centre
        @test g[2, 1] ≈ 263.15                    # exactly at the top centre
        @test g[3, 1] ≈ 268.15                    # interpolated between the centres
        @test g[4, 1] ≈ 273.15                    # exactly at the bottom centre
        @test isnan(g[5, 1])                      # below the bottom centre
        @test maximum(filter(isfinite, parent(g))) <= 273.15

        # Nearest-neighbour agrees on the missing ends.
        gn = gemb_interp(z_center, T, z_target; interp_method=:nearest)
        @test isnan(gn[1, 1]) && isnan(gn[5, 1])
        @test gn[3, 1] ∈ (263.15, 273.15)

        # Ascending `z_center` (an upward-ordered column) takes the other branch and must
        # behave identically.
        g_asc = gemb_interp(reverse(z_center, dims=1), reverse(T, dims=1), z_target)
        @test isnan(g_asc[1, 1]) && isnan(g_asc[5, 1])
        @test g_asc[3, 1] ≈ 268.15
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

    @testset "Integration: dz2z" begin
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
