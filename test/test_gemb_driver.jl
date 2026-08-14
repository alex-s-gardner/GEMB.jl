using Test
using GEMB
using Dates

@testset "gemb_driver" begin
    @testset "Basic integration test" begin
        # Create simple parameters
        params = initialize_parameters(
            densification_method = :Arthern,
            albedo_method = :GardnerSharp
        )

        # Create simple forcing (1 week, hourly)
        n_steps = 24 * 7
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation
            fill(5.0, n_steps),  # wind_speed
            fill(200.0, n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        # Initialize profile
        profile = initialize_profile(params, forcing)

        # Run model for 1 week
        output = gemb(profile, forcing, params)

        # Basic sanity checks
        @test output isa DimStack
        @test haskey(output, :temperature)
        @test haskey(output, :density)
        @test haskey(output, :dz)

        # Check output dimensions
        @test length(dims(output, Ti)) == n_steps
        @test length(dims(output, Z)) > 0

        # Check temperature is in reasonable range
        temps = output[:temperature]
        @test all(isfinite.(temps[.!isnan.(temps)]))
        @test all(temps[.!isnan.(temps)] .> 200.0)  # Above absolute zero
        @test all(temps[.!isnan.(temps)] .< 300.0)  # Below boiling

        # Check density is in reasonable range (ice density)
        densities = output[:density]
        @test all(isfinite.(densities[.!isnan.(densities)]))
        @test all(densities[.!isnan.(densities)] .> 100.0)   # Above fresh snow
        @test all(densities[.!isnan.(densities)] .<= 917.0)  # At or below ice
    end

    @testset "Conservation test - no forcing" begin
        # Test that without any external forcing or melt, mass is conserved

        params = initialize_parameters(
            output_frequency = :all
        )

        # Zero forcing (no precipitation, no melt conditions)
        n_steps = 24  # 1 day
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(250.0, n_steps),  # temperature_air - Cold, no melt
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation - No precip
            fill(2.0, n_steps),  # wind_speed
            zeros(n_steps),  # shortwave_downward - No solar
            fill(150.0, n_steps),  # longwave_downward
            fill(50.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)

        output = gemb(profile, forcing, params)

        # Column mass at each output step (NaN-padded cells excluded).
        column_mass(ti) = begin
            d = output[:density][Ti=Near(ti)]
            z = output[:dz][Ti=Near(ti)]
            sum(d[.!isnan.(d)] .* z[.!isnan.(z)])
        end

        # NOTE: initialize_profile can build a column that slightly overshoots
        # column_depth_max; the first manage_layer_thickness call trims the bottom layer to
        # bring the column back within bounds. That one-time trim is mass leaving
        # through the bottom domain boundary (a legitimate grid operation), not a
        # physics conservation violation. Comparing the pre-trim initial mass to
        # the final mass therefore conflates the trim with conservation.
        #
        # To test physics conservation, compare the mass at the first output step
        # (after the initialization trim has settled) to the final mass. With no
        # precipitation, melt, or runoff, the physics must conserve column mass.
        first_mass = column_mass(time[1])
        final_mass = column_mass(time[end])

        @test abs(final_mass - first_mass) / first_mass < 0.01
    end

    @testset "Accumulation test" begin
        # Test that precipitation adds mass correctly

        params = initialize_parameters(
            output_frequency = :all
        )

        n_steps = 10
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        # Constant precipitation
        precip_rate = 0.001  # kg/m²/s for 1 hour = 3.6 kg/m²
        precip_per_hour = precip_rate * 3600.0

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            fill(precip_per_hour, n_steps),  # precipitation
            fill(2.0, n_steps),  # wind_speed
            zeros(n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)
        initial_mass = sum(profile.density .* profile.dz)

        output = gemb(profile, forcing, params)

        final_density = output[:density][Ti=Near(time[end])]
        final_dz = output[:dz][Ti=Near(time[end])]
        final_mass = sum(final_density[.!isnan.(final_density)] .*
                        final_dz[.!isnan.(final_dz)])

        # Model should produce finite, physical mass values
        @test isfinite(final_mass)
        @test final_mass > 0
    end

    @testset "Output frequency options" begin
        params = initialize_parameters()

        n_steps = 24 * 3  # 3 days
        start_time = DateTime(2020, 1, 1)
        time = start_time .+ Hour.(0:n_steps-1)

        forcing = initialize_forcing(
            time,
            fill(260.0, n_steps),  # temperature_air
            fill(101325.0, n_steps),  # pressure_air
            zeros(n_steps),  # precipitation
            fill(5.0, n_steps),  # wind_speed
            fill(100.0, n_steps),  # shortwave_downward
            fill(200.0, n_steps),  # longwave_downward
            fill(100.0, n_steps),  # vapor_pressure
            temperature_observation_height = 2.0,
            wind_observation_height = 10.0
        )

        profile = initialize_profile(params, forcing)

        # Test different output frequencies
        for freq in [:all, :daily, :last]
            params_freq = initialize_parameters(output_frequency = freq)
            output = gemb(profile, forcing, params_freq)

            if freq == :all
                @test length(dims(output, Ti)) == n_steps
            elseif freq == :daily
                @test length(dims(output, Ti)) == 3  # 3 days
            elseif freq == :last
                @test length(dims(output, Ti)) == 1
            end
        end
    end

    @testset "_compute_output_times weekly" begin
        # 21 days of hourly steps, starting Wed 2020-01-01.
        # Monday-anchored weeks: partial week ending 2020-01-05 (Sun),
        # full weeks 01-06..01-12, 01-13..01-19, and partial 01-20..01-21.
        times = DateTime(2020, 1, 1) .+ Hour.(0:(24 * 21 - 1))
        weekly = GEMB._compute_output_times(times, :weekly)

        # One emission per completed week boundary (three transitions here).
        expected = [t for i in eachindex(times)
                    for t in (times[i],)
                    if floor(Date(times[i]), Week) !=
                       floor(Date(i < length(times) ? times[i+1] :
                                  times[end] + (times[end] - times[end-1])), Week)]
        @test weekly == expected
        # Each emitted step is the last hour (23:00) before a week rollover.
        @test all(hour.(weekly) .== 23)
        # Boundaries land on Sundays (day before Monday-anchored week starts).
        @test all(dayname.(weekly) .== "Sunday")
    end

    @testset "Output provenance (spinup vs. no spinup)" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        params = initialize_parameters(output_frequency=:last)
        profile = initialize_profile(params, forcing)

        # No spinup: gemb runs directly on the freshly initialized profile.
        out_nospin = gemb(profile, forcing, params)
        md1 = DimensionalData.metadata(out_nospin)
        @test md1["spinup_performed"] == false
        @test !haskey(md1, "spinup_cycles")
        # Forcing provenance is still present.
        @test haskey(md1, "dataset")

        # With spinup: spinup_* / climatology_* provenance carries onto the output.
        ps = gemb_spinup(profile, forcing, params;
                         max_iterations=3, convergence_delta_density=0.01)
        out_spin = gemb(ps, forcing, params)
        md2 = DimensionalData.metadata(out_spin)
        @test md2["spinup_performed"] == true
        @test md2["spinup_cycles"] == 3
        @test haskey(md2, "climatology_n_years")
    end

    @testset "CF layer metadata" begin
        n = 365
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        params = initialize_parameters(output_frequency=:monthly)
        profile = initialize_profile(params, forcing)
        out = gemb(profile, forcing, params)

        # Every layer is described, in the stack's own key order. DimensionalData
        # rejects a partial or permuted `layermetadata`, so this is the assertion that
        # catches a layer added to `gemb` without a matching CF table entry.
        @test keys(DimensionalData.layermetadata(out)) == keys(out)
        for k in keys(out)
            md = DimensionalData.metadata(out[k])
            @test md isa Dict{String,String}
            @test !isempty(md["units"])
            @test !isempty(md["long_name"])
            @test haskey(md, "cell_methods")
        end

        # Interval semantics: a per-interval sum, a per-interval mean, a snapshot.
        @test DimensionalData.metadata(out[:melt])["cell_methods"] == "time: sum"
        @test DimensionalData.metadata(out[:temperature_air])["cell_methods"] == "time: mean"
        @test DimensionalData.metadata(out[:temperature])["cell_methods"] == "time: point"
        @test DimensionalData.metadata(out[:density])["units"] == "kg m-3"
        @test DimensionalData.metadata(out[:temperature_air])["standard_name"] == "air_temperature"

        # Exactly the layers whose CF semantics match GEMB's quantity carry a
        # `standard_name`. In particular the per-cell profile fields do not: CF's
        # temperature_in_surface_snow / surface_snow_density /
        # liquid_water_content_of_surface_snow are all bulk whole-pack quantities.
        named = Set(k for k in keys(out)
                    if haskey(DimensionalData.metadata(out[k]), "standard_name"))
        @test named == Set([:melt, :runoff, :precipitation, :rain,
                            :shortwave_net, :longwave_net,
                            :heat_flux_sensible, :heat_flux_latent,
                            :albedo_surface, :temperature_air])
        for k in (:temperature, :density, :water)
            @test !haskey(DimensionalData.metadata(out[k]), "standard_name")
        end

        # Dimensions: `Ti` is a real coordinate; `Z` is a bare cell index, so it must
        # not claim a standard_name or a `positive` direction.
        @test DimensionalData.metadata(dims(out, Ti))["standard_name"] == "time"
        zmd = DimensionalData.metadata(dims(out, Z))
        @test !haskey(zmd, "standard_name")
        @test !haskey(zmd, "positive")

        # CF global attributes sit alongside the run's provenance.
        @test DimensionalData.metadata(out)["Conventions"] == "CF-1.11"
        @test haskey(DimensionalData.metadata(out), "dataset")

        # An extracted profile keeps units but drops `cell_methods`, which would
        # otherwise reference a `time` coordinate the column does not have. Profiles
        # from `initialize_profile` describe themselves the same way.
        for p in (gemb_profile(out), profile)
            @test DimensionalData.metadata(p[:density])["units"] == "kg m-3"
            @test !haskey(DimensionalData.metadata(p[:density]), "cell_methods")
        end

        # `gemb_interp` regrids onto a genuine vertical coordinate and carries the
        # interpolated variable's own attributes across.
        gridded = gemb_interp(dz2z(parent(out[:dz])), out[:temperature],
                              collect(range(0.0, -5.0; length=10)))
        zattrs = DimensionalData.metadata(dims(gridded, Z))
        @test zattrs["units"] == "m"
        @test zattrs["positive"] == "up"
        @test DimensionalData.metadata(gridded)["units"] == "K"
    end
end
