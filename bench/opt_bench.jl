# Optimization benchmark + numerical snapshot driver (skill: julia_optimize)
using GEMB
using GEMB_ClimateForcing  # synthetic forcing lives in the companion package
using BenchmarkTools
using Statistics
using Printf

const TIME_STEP_HOURS = 3

function build_inputs()
    ds = simulate_climate_forcing("test_1", TIME_STEP_HOURS)  # DimStack
    cf = GEMB.initialize_forcing(ds)                          # convert to ClimateForcing
    cf_clim = forcing_climatology(cf)
    # `initialize_profile` returns (profile, mp); the returned mp carries the
    # climate-derived column limits, so rebind it and use it for both the spinup
    # and the transient run. `gemb_spinup` forces
    # output_frequency=:last internally, so no separate spinup params object is needed.
    mp = ModelParameters(output_frequency=:daily)
    profile, mp = initialize_profile(mp, cf)
    return cf, mp, profile, cf_clim
end

# Run the representative hot path: 1994-2025 3-hourly simulation with 75-year spinup.
function run_hotpath(profile, cf_clim, cf, mp)
    profile_spunup = gemb_spinup(profile, cf_clim, mp; max_iterations=75)
    output = gemb(profile_spunup, cf, mp)
    return output
end

# Numerical fingerprint for equivalence checks across changes.
function snapshot(output)
    keys_sorted = sort(collect(keys(output)))
    d = Dict{Symbol,Tuple{Float64,Float64,Float64,Int}}()
    for k in keys_sorted
        a = parent(output[k])
        v = filter(x -> !isnan(x) && isfinite(x), vec(Float64.(a)))
        if isempty(v)
            d[k] = (NaN, NaN, NaN, 0)
        else
            d[k] = (sum(v), minimum(v), maximum(v), length(v))
        end
    end
    return d
end

function print_snapshot(d)
    for k in sort(collect(keys(d)))
        s, lo, hi, n = d[k]
        @printf("  %-32s sum=%.10e min=%.10e max=%.10e n=%d\n", k, s, lo, hi, n)
    end
end

function compare_snapshots(a, b)
    maxrel = 0.0
    worst = :none
    for k in keys(a)
        for i in 1:3
            av = a[k][i]; bv = b[k][i]
            (isnan(av) && isnan(bv)) && continue
            denom = max(abs(av), 1e-30)
            rel = abs(av - bv) / denom
            if rel > maxrel
                maxrel = rel; worst = k
            end
        end
    end
    return maxrel, worst
end

cf, mp, profile, cf_clim = build_inputs()

println("Warmup...")
out = run_hotpath(profile, cf_clim, cf, mp)
snap = snapshot(out)
println("Numerical snapshot:")
print_snapshot(snap)

# Save/compare snapshot against a reference file if present
const SNAP_FILE = joinpath(@__DIR__, "opt_bench_snapshot.txt")
if get(ENV, "SAVE_SNAPSHOT", "0") == "1"
    open(SNAP_FILE, "w") do io
        for k in sort(collect(keys(snap)))
            s, lo, hi, n = snap[k]
            @printf(io, "%s %.15e %.15e %.15e %d\n", k, s, lo, hi, n)
        end
    end
    println("Saved snapshot to $SNAP_FILE")
elseif isfile(SNAP_FILE)
    ref = Dict{Symbol,Tuple{Float64,Float64,Float64,Int}}()
    for line in eachline(SNAP_FILE)
        parts = split(line)
        ref[Symbol(parts[1])] = (parse(Float64, parts[2]), parse(Float64, parts[3]),
                                  parse(Float64, parts[4]), parse(Int, parts[5]))
    end
    maxrel, worst = compare_snapshots(ref, snap)
    @printf("Max rel diff vs saved snapshot: %.3e (field: %s)\n", maxrel, worst)
end

println("\nBenchmarking hot path (1994-2025 3-hourly simulation with 75-year spinup)...")
b = @benchmark run_hotpath($profile, $cf_clim, $cf, $mp) samples=7 seconds=120 evals=1
@printf("  min    = %.3f s\n", minimum(b.times)/1e9)
@printf("  median = %.3f s\n", median(b.times)/1e9)
@printf("  allocs = %d\n", b.allocs)
@printf("  memory = %.1f MiB\n", b.memory/2^20)
