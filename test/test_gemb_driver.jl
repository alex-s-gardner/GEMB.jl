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
                         simulation_years_maximum=3 * GEMB._years_per_cycle(forcing), convergence_delta_density=0.01)
        out_spin = gemb(ps, forcing, params)
        md2 = DimensionalData.metadata(out_spin)
        @test md2["spinup_performed"] == true
        @test md2["spinup_cycles"] == 3
        @test haskey(md2, "climatology_n_years")
    end

    @testset "Thermal-solve freezing diagnostics are interval sums" begin
        # Six 30-minute steps at the melting point, saturated air and 100 W m-2 of net longwave loss,
        # from a column at 273.15 K whose 6 mm surface cell holds more water than ice. The thermal
        # solve freezes that cell's water itself in the first step and never again.
        T0 = GEMB.CtoK
        n = 6
        time = DateTime(2020, 4, 1) .+ Minute.(30 .* (0:n-1))
        forcing = initialize_forcing(time, fill(T0, n), fill(100000.0, n), fill(0.0, n),
            fill(5.0, n), fill(0.0, n), fill(GEMB.SB * T0^4 - 100.0 / 0.98, n),
            fill(GEMB.saturation_vapor_pressure_water(T0), n);
            temperature_air_mean=260.0, wind_speed_mean=5.0, precipitation_mean=200.0,
            temperature_observation_height=2.0, wind_observation_height=10.0)
        profile = initialize_profile(initialize_parameters(), forcing)
        parent(profile[:temperature]) .= T0
        parent(profile[:dz])[1] = 0.0061
        parent(profile[:density])[1] = 93.4
        parent(profile[:water])[1] = 0.597
        each = gemb(profile, forcing, initialize_parameters(output_frequency=:all))
        whole = gemb(profile, forcing, initialize_parameters(output_frequency=:last))
        @test collect(each[:refreeze_thermal_count]) == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        @test each[:refreeze_thermal][1] > 0.0 && all(iszero, each[:refreeze_thermal][2:end])
        @test each[:refreeze_thermal][1] <= each[:refreeze][1]
        @test only(whole[:refreeze_thermal_count]) == sum(each[:refreeze_thermal_count])
        @test only(whole[:refreeze_thermal]) == sum(each[:refreeze_thermal])
        md = DimensionalData.metadata(each[:refreeze_thermal_count])
        @test md["units"] == "1" && md["cell_methods"] == "time: sum"
        @test DimensionalData.metadata(each[:refreeze_thermal])["units"] == "kg m-2"
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
                            :albedo_broadband, :temperature_air])
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
    @testset "restart layer schema" begin
        n = 90
        time = DateTime(2020, 1, 1) .+ Day.(0:n-1)
        forcing = initialize_forcing(
            time, fill(255.0, n), fill(85000.0, n), fill(0.5, n), fill(3.0, n),
            fill(50.0, n), fill(180.0, n), fill(80.0, n);
            temperature_air_mean=255.0, wind_speed_mean=3.0, precipitation_mean=182.6)
        mp = initialize_parameters(output_frequency=:monthly)
        profile = initialize_profile(mp, forcing)
        out = gemb(profile, forcing, mp)

        # `RESTART_LAYERS` is the column state, so it must be exactly what the two
        # constructors of a profile produce and exactly what `gemb` consumes. Pinning
        # both directions here is what makes the constant safe to read from outside the
        # package: a layer added to the state without joining the list fails on the
        # `gemb_profile` side, and a name left in the list after the state drops it
        # fails on the `initialize_profile` side.
        @test Set(keys(gemb_profile(out))) == Set(RESTART_LAYERS)
        @test Set(RESTART_LAYERS) ⊆ Set(keys(profile))

        # A profile short of a layer is rejected before the time loop, naming the gap.
        # Without the up-front check this surfaced as a `FieldError` on an internal
        # `state` field, from deep inside the run.
        incomplete = DimStack(NamedTuple(l => profile[l]
                                         for l in keys(profile) if l !== :age))
        err = try
            gemb(incomplete, forcing, mp)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("age", err.msg)
        @test occursin("RESTART_LAYERS", err.msg)

        # `dt_divisors` is derived inside `gemb`, so whatever the caller passes is
        # overwritten rather than honoured — that is why it is in `DERIVED_PARAMETERS`.
        @test :dt_divisors in DERIVED_PARAMETERS
        @test all(f -> f in fieldnames(ModelParameters), DERIVED_PARAMETERS)
    end
    @testset "heat_flux_basal output" begin
        # The conductive flux across the deepest interior face was already computed and
        # consumed by the verbose energy budget; Fix 6 reports it. It is an interval mean,
        # like the other three heat fluxes.
        mp = initialize_parameters(output_frequency=:daily)
        n_steps = 24 * 14
        time = DateTime(2020, 1, 1) .+ Hour.(0:n_steps-1)
        forcing = initialize_forcing(time, fill(255.0, n_steps), fill(101325.0, n_steps),
            fill(0.1, n_steps), fill(5.0, n_steps), fill(50.0, n_steps),
            fill(200.0, n_steps), fill(100.0, n_steps);
            temperature_air_mean=255.0, wind_speed_mean=5.0, precipitation_mean=300.0,
            temperature_observation_height=2.0, wind_observation_height=10.0)
        # The base of a deep initialized column can be isothermal to the last bit, so the
        # reservoir is set 1 K warmer than the cell above it to give the face a gradient.
        profile = initialize_profile(mp, forcing)
        parent(profile[:temperature])[end] = parent(profile[:temperature])[end-1] + 1.0
        out = gemb(profile, forcing, mp; verbose=true)

        bhf = parent(out[:heat_flux_basal])
        @test haskey(out, :heat_flux_basal)
        @test length(bhf) == length(parent(out[:heat_flux_sensible]))
        # Every interval must be filled — a NaN would mean an output step the accumulator
        # never reached.
        @test all(isfinite, bhf)
        # A real conductive flux, not an identically-zero placeholder, and bounded well
        # inside the range of any plausible geothermal flux confusion.
        @test any(!iszero, bhf)
        @test all(abs.(bhf) .< 100.0)

        md = GEMB.cf_attributes(:heat_flux_basal)
        @test md["units"] == "W m-2"
        @test md["cell_methods"] == "time: mean"
        # Diagnosed against a Dirichlet reservoir, so no CF geothermal standard name applies.
        @test !haskey(md, "standard_name")
    end

    @testset "viscosity output (output_viscosity)" begin
        n_steps = 24 * 14
        time = DateTime(2020, 1, 1) .+ Hour.(0:n_steps-1)
        forcing = initialize_forcing(time, fill(255.0, n_steps), fill(101325.0, n_steps),
            fill(0.1, n_steps), fill(5.0, n_steps), fill(50.0, n_steps),
            fill(200.0, n_steps), fill(100.0, n_steps);
            temperature_air_mean=255.0, wind_speed_mean=5.0, precipitation_mean=300.0,
            temperature_observation_height=2.0, wind_observation_height=10.0)

        _run(; kw...) = begin
            mp = initialize_parameters(; output_frequency=:daily, kw...)
            (mp, gemb(initialize_profile(mp, forcing), forcing, mp))
        end

        # Off by default, and absent rather than all-NaN — the layer costs a per-timestep
        # allocation, so the default hot path must not carry it.
        @test initialize_parameters().output_viscosity == false
        (_, out_off) = _run(densification_method=:CrocusPure)
        @test !haskey(out_off, :viscosity)

        # `:CrocusPure` applies the settling law at every density, so every cell gets a
        # viscosity.
        (mp_p, out_p) = _run(densification_method=:CrocusPure, output_viscosity=true)
        @test haskey(out_p, :viscosity)
        v = parent(out_p[:viscosity])
        @test size(v) == size(parent(out_p[:density]))
        @test all(isfinite, v)
        # Physically bounded: firn viscosity spans many orders of magnitude but is
        # positive everywhere and nowhere near zero or infinite.
        @test all(v .> 0)
        @test all(1e9 .< v .< 1e22)

        # Colder and denser firn is stiffer, which is the sign the Vionnet law carries in
        # both its exponentials.
        d = parent(out_p[:density])
        deepest = size(v, 1)
        @test v[deepest, end] > v[1, end]
        @test d[deepest, end] >= d[1, end]

        # Under the `:density` rule `:Crocus` hands every cell at or above its handoff density to
        # `:GSFC2020`, which forms no viscosity, so exactly those cells must be NaN. This is the
        # property that makes NaN meaningful rather than a missing-value artifact. (Under the
        # default `:rate_match` rule Crocus may still govern above the floor, and those cells keep
        # their viscosity.)
        (_, out_c) = _run(; densification_method=:Crocus, output_viscosity=true,
                          GEMB.CROCUS_PUBLISHED_HANDOFF...)
        vc = parent(out_c[:viscosity])
        dc = parent(out_c[:density])
        for i in eachindex(vc)
            if isfinite(vc[i])
                @test dc[i] < GEMB.CROCUS_HYBRID_DENSITY
            end
        end

        # Every non-Crocus scheme is a `c(rho_i - rho)` relaxation with no stress in it, so
        # the layer is present but entirely NaN rather than a reconstructed number.
        for method in (:Arthern, :HerronLangway, :GSFC2020, :Barnola1991)
            (_, out_n) = _run(densification_method=method, output_viscosity=true)
            @test all(isnan, parent(out_n[:viscosity]))
        end

        md = GEMB.cf_attributes(:viscosity)
        @test md["units"] == "Pa s"
        @test md["cell_methods"] == "time: point"
        @test !haskey(md, "standard_name")
    end
end

@testset "Radiation geometry: forcing layers, fallback, validation" begin
    n = 48
    time = DateTime(2020, 6, 1) .+ Hour.(0:n-1)
    sw = [max(0.0, 600.0 * sin(2π * (k - 6) / 24)) for k in 0:n-1]
    sza = [s > 0 ? 50.0 : 100.0 for s in sw]
    diffuse = 0.3 .* sw
    base = (time, fill(265.0, n), fill(80000.0, n), zeros(n), fill(3.0, n), sw,
        fill(250.0, n), fill(250.0, n))
    kw = (temperature_air_mean=265.0, wind_speed_mean=3.0, precipitation_mean=500.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)

    cf_layers = initialize_forcing(base...; kw..., solar_zenith_angle=sza,
        cloud_optical_thickness=fill(2.0, n), shortwave_downward_diffuse=diffuse)
    @test parent(cf_layers.solar_zenith_angle) == sza
    @test parent(cf_layers.shortwave_downward_diffuse) == diffuse

    mp = initialize_parameters(output_frequency=:all)
    # Ice column with 1 m of snow on top, so the grain-size albedo is in play.
    profile = initialize_profile(mp, cf_layers; constant_density=true,
        constant_temperature=true)
    k_snow = findfirst(>=(1.0), cumsum(parent(profile[:dz])))
    parent(profile[:density])[1:k_snow] .= 350.0
    parent(profile[:grain_radius])[1:k_snow] .= 0.3

    # The forcing layers reach the physics: the run differs from one with no geometry, and
    # an `mp` scalar does not override a supplied layer.
    out_layers = gemb(profile, cf_layers, mp; verbose=true)
    cf_none = initialize_forcing(base...; kw...)
    out_none = gemb(profile, cf_none, mp)
    @test parent(out_layers[:albedo_broadband]) != parent(out_none[:albedo_broadband])
    mp_scalar = initialize_parameters(output_frequency=:all, solar_zenith_angle=30.0,
        cloud_optical_thickness=1.0)
    @test parent(gemb(profile, cf_layers, mp_scalar)[:albedo_broadband]) ==
          parent(out_layers[:albedo_broadband])

    # Without layers, the `mp` scalars are used for every step: a constant 30° layer and the
    # `mp` fallback at 30° give the same run.
    cf_const = initialize_forcing(base...; kw..., solar_zenith_angle=fill(30.0, n),
        cloud_optical_thickness=fill(1.0, n))
    @test parent(gemb(profile, cf_none, mp_scalar)[:albedo_broadband]) ==
          parent(gemb(profile, cf_const, mp)[:albedo_broadband])

    # A DimStack forcing supplies them as layers.
    tdim = Ti(time)
    names = (:temperature_air, :pressure_air, :precipitation, :wind_speed,
        :shortwave_downward, :longwave_downward, :vapor_pressure)
    layers = NamedTuple{names}(Tuple(DimArray(v, (tdim,)) for v in base[2:end]))
    meta = Dict{String,Any}(String(k) => v for (k, v) in pairs(kw))
    stack = DimStack(merge(layers, (solar_zenith_angle=DimArray(sza, (tdim,)),
            shortwave_downward_diffuse=DimArray(diffuse, (tdim,)))); metadata=meta)
    cf_stack = initialize_forcing(stack)
    @test parent(cf_stack.solar_zenith_angle) == sza
    @test parent(cf_stack.shortwave_downward_diffuse) == diffuse
    @test all(iszero, cf_stack.cloud_optical_thickness)

    # Per-step range checks.
    bad_sza = copy(sza); bad_sza[5] = 181.0
    @test_throws "outside [0, 180]" initialize_forcing(base...; kw...,
        solar_zenith_angle=bad_sza)
    @test_throws "negative" initialize_forcing(base...; kw...,
        cloud_optical_thickness=fill(-0.1, n))
    bad_diffuse = copy(diffuse); bad_diffuse[13] = sw[13] + 1.0
    @test_throws "step 13" initialize_forcing(base...; kw...,
        shortwave_downward_diffuse=bad_diffuse)
    # A scalar diffuse flux larger than a night step's zero shortwave is a violation too.
    @test_throws ArgumentError initialize_forcing(base...; kw...,
        shortwave_downward_diffuse=10.0)
    @test_throws "timesteps" initialize_forcing(base...; kw..., solar_zenith_angle=sza[1:10])
end

@testset "Impurity and cloud-fraction layers must agree with ModelParameters" begin
    n = 48
    time = DateTime(2020, 6, 1) .+ Hour.(0:n-1)
    sw = [max(0.0, 600.0 * sin(2π * (k - 6) / 24)) for k in 0:n-1]
    base = (time, fill(265.0, n), fill(80000.0, n), zeros(n), fill(3.0, n), sw,
        fill(250.0, n), fill(250.0, n))
    kw = (temperature_air_mean=265.0, wind_speed_mean=3.0, precipitation_mean=500.0,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    mp = initialize_parameters(output_frequency=:all)
    cf = initialize_forcing(base...; kw...)
    # Ice column with 1 m of snow on top, so black carbon reaches the snow albedo.
    profile = initialize_profile(mp, cf; constant_density=true, constant_temperature=true)
    k_snow = findfirst(>=(1.0), cumsum(parent(profile[:dz])))
    parent(profile[:density])[1:k_snow] .= 350.0
    parent(profile[:grain_radius])[1:k_snow] .= 0.3
    out = gemb(profile, cf, mp)

    # A forcing value the run would discard is refused, naming the layer and the parameter.
    @test_throws "`black_carbon_snow` is 0.05 at step 1" gemb(profile,
        initialize_forcing(base...; kw..., black_carbon_snow=0.05), mp)
    @test_throws "`ModelParameters.black_carbon_ice` is 0.0" gemb(profile,
        initialize_forcing(base...; kw..., black_carbon_ice=0.1), mp)
    @test_throws "`cloud_fraction` is 0.3 at step 1" gemb(profile,
        initialize_forcing(base...; kw..., cloud_fraction=0.3), mp)
    @test_throws ArgumentError gemb(profile, cf,
        initialize_parameters(output_frequency=:all, cloud_fraction=0.3))

    # A layer that departs from the parameter at a single step is refused at that step.
    bc = fill(mp.black_carbon_snow, n); bc[17] = 0.2
    cf_bc = GEMB.ClimateForcing(cf.temperature_air, cf.pressure_air, cf.precipitation,
        cf.wind_speed, cf.shortwave_downward, cf.longwave_downward, cf.vapor_pressure,
        DimArray(bc, (Ti(time),)), cf.black_carbon_ice, cf.cloud_optical_thickness,
        cf.solar_zenith_angle, cf.shortwave_downward_diffuse, cf.cloud_fraction,
        cf.time_step, cf.temperature_air_mean, cf.wind_speed_mean, cf.precipitation_mean,
        cf.temperature_observation_height, cf.wind_observation_height)
    @test_throws "is 0.2 at step 17" gemb(profile, cf_bc, mp)

    # A value read from forcing-file metadata is checked the same way.
    tdim = Ti(time)
    names = (:temperature_air, :pressure_air, :precipitation, :wind_speed,
        :shortwave_downward, :longwave_downward, :vapor_pressure)
    layers = NamedTuple{names}(Tuple(DimArray(v, (tdim,)) for v in base[2:end]))
    meta = Dict{String,Any}(String(k) => v for (k, v) in pairs(kw))
    meta["black_carbon_snow"] = 0.02
    @test_throws "`black_carbon_snow` is 0.02" gemb(profile,
        initialize_forcing(DimStack(layers; metadata=meta)), mp)

    # Matching values run, and give the run the parameter alone gives.
    mp_bc = initialize_parameters(output_frequency=:all, black_carbon_snow=0.05,
        cloud_fraction=0.3)
    cf_bc = initialize_forcing(base...; kw..., black_carbon_snow=0.05, cloud_fraction=0.3)
    out_bc = gemb(profile, cf_bc, mp_bc)
    @test parent(out_bc[:albedo_broadband]) != parent(out[:albedo_broadband])
    @test parent(gemb(profile, cf_bc, mp_bc)[:melt]) == parent(out_bc[:melt])
end

@testset "Elevation identity: basal flux, SMB, water and firn air" begin
    # At fixed column depth, −Σ ice_flux = (SMB − ΔW)/ρᵢ + ΔFAC exactly, with SMB the surface
    # mass balance, ΔW the change in liquid water and ΔFAC the change in firn air content.
    n = 365
    time = DateTime(2003, 1, 1) .+ Day.(0:n-1)
    seasonal = cos.(2π .* (1:n) ./ 365.0 .- π)
    cf = initialize_forcing(time, 266.0 .+ 9.0 .* seasonal, fill(85000.0, n), fill(2.5, n),
        fill(4.0, n), 120.0 .+ 160.0 .* max.(0.0, seasonal), fill(290.0, n), fill(300.0, n);
        temperature_air_mean=266.0, wind_speed_mean=4.0, precipitation_mean=2.5 * 365.25,
        temperature_observation_height=2.0, wind_observation_height=10.0)
    mp = initialize_parameters(output_frequency=:all)
    profile = initialize_profile(mp, cf)
    out = gemb(profile, cf, mp; verbose=true)
    ρi = mp.density_ice
    total(k) = sum(parent(out[k]))
    smb = total(:precipitation) + total(:evaporation_condensation) - total(:runoff)
    ΔW = sum(parent(out[:water])[:, end]) - sum(parent(profile[:water]))
    fac(dz, ρ) = firn_air_content(dz, ρ, ρi)
    ΔFAC = fac(parent(out[:dz])[:, end], parent(out[:density])[:, end]) -
           fac(parent(profile[:dz]), parent(profile[:density]))
    @test total(:melt) > 0.0
    @test -total(:ice_flux) ≈ (smb - ΔW) / ρi + ΔFAC atol = 1e-9
end
