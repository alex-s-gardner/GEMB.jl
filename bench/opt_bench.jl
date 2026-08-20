# Optimization benchmark + numerical snapshot driver (skill: julia_optimize)
#
#   julia --project=bench bench/opt_bench.jl                    # compare against the baseline
#   SAVE_SNAPSHOT=1 julia --project=bench bench/opt_bench.jl    # re-pin the baseline
#
# `opt_bench_snapshot.txt` is tracked, so an optimization is checked against the
# baseline everyone else measured. Re-pin only when a change is *meant* to move the
# numbers, and say what moved them in the commit message.
#
# The comparison is a few-ULP tolerance, not bit-equality, because the snapshot is
# compared across machines: a different Julia/LLVM version or a CPU that vectorizes
# reductions differently reassociates the `sum` over ~90k elements and lands a few
# representable steps away. That is not a change in the model. The script prints every
# differing statistic (not just the worst) and exits non-zero only when something exceeds
# `SNAPSHOT_ULP_TOLERANCE` or changes structurally — a field appearing, disappearing,
# changing element count, or turning NaN.
#
# So: a table of 1-3 ULP differences with "VERDICT: equivalent" is the expected output on a
# machine other than the one that pinned the file. Re-pinning to silence it would only move
# the same noise onto the next person, and would discard the baseline's cross-machine value.
using GEMB
using GEMB_ClimateForcing  # synthetic forcing lives in the companion package
using BenchmarkTools
using Statistics
using Printf

const TIME_STEP_HOURS = 3

function build_inputs()
    ds = simulate_climate_forcing("test_1", TIME_STEP_HOURS)  # DimStack
    cf = initialize_forcing(ds)                          # convert to ClimateForcing
    cf_clim = forcing_climatology(cf)
    # `initialize_profile` sizes the column to the forcing climate; the derived depth
    # rides along in the profile, so the same `mp` serves the spinup and the transient
    # run. `gemb_spinup` forces output_frequency=:last internally, so no separate
    # spinup params object is needed.
    mp = ModelParameters(output_frequency=:daily)
    profile = initialize_profile(mp, cf)
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

const STAT_NAMES = ("sum", "min", "max")

# Floating-point distance in representable steps. Only meaningful for two finite values of
# the same sign, which is the only case the caller uses it for.
function ulp_distance(a::Float64, b::Float64)
    (isfinite(a) && isfinite(b) && signbit(a) == signbit(b)) || return -1
    return abs(reinterpret(Int64, a) - reinterpret(Int64, b))
end

# A snapshot is compared across machines, so it cannot be held to bit-equality: a different
# Julia/LLVM version or a CPU that vectorizes reductions differently reassociates the `sum`
# over ~90k elements and lands a few steps away. Pairwise summation's error grows as
# O(log2(n)·eps), and log2(90_000) ≈ 16.5, so differently-associated sums of this length are
# expected to sit within tens of ULP. The measured spread of the tracked snapshot against
# Julia 1.12.5 on Apple silicon is 1-3 ULP across 19 of 35 fields.
#
# 64 ULP is therefore "same computation, different summation order" — about 4x the expected
# scale, and still ~1e-14 relative, far below any difference a physics or optimization change
# would produce. Anything above it is a real change in the numbers.
const SNAPSHOT_ULP_TOLERANCE = 64

# Compare every statistic of every field and return one entry per difference found, worst
# first. Reporting only the single largest relative difference (as this did originally) hides
# the shape of a regression: a change that moves 19 fields and one that moves a single field
# look identical when only the max is printed.
function compare_snapshots(ref, now)
    diffs = []          # (field, stat, refval, nowval, ulp, rel)
    structural = String[]

    for k in sort(collect(keys(ref)))
        if !haskey(now, k)
            push!(structural, "$k: in snapshot, absent from this run")
            continue
        end
        # Element count is not a float comparison — a change here means the run produced a
        # different amount of data, which no summation order can explain.
        if ref[k][4] != now[k][4]
            push!(structural, "$k: element count $(ref[k][4]) -> $(now[k][4])")
        end
        for i in 1:3
            av = ref[k][i]; bv = now[k][i]
            (isnan(av) && isnan(bv)) && continue
            if isnan(av) != isnan(bv)
                push!(structural, "$k $(STAT_NAMES[i]): NaN-ness changed ($av -> $bv)")
                continue
            end
            av === bv && continue
            rel = abs(av - bv) / max(abs(av), 1e-30)
            push!(diffs, (k, STAT_NAMES[i], av, bv, ulp_distance(av, bv), rel))
        end
    end
    for k in sort(collect(keys(now)))
        haskey(ref, k) || push!(structural, "$k: new field, absent from snapshot")
    end

    # Sort by ULP where it is defined, else by relative difference (which is what a sign
    # change or a non-finite value leaves us).
    sort!(diffs, by = d -> (d[5] >= 0 ? Float64(d[5]) : Inf, d[6]), rev=true)
    return diffs, structural
end

function report_comparison(diffs, structural)
    if isempty(diffs) && isempty(structural)
        println("Snapshot: bit-identical to the tracked baseline.")
        return true
    end

    for msg in structural
        println("  STRUCTURAL  ", msg)
    end

    over = filter(d -> d[5] < 0 || d[5] > SNAPSHOT_ULP_TOLERANCE, diffs)
    if !isempty(diffs)
        n_fields = length(unique(d[1] for d in diffs))
        @printf("Snapshot: %d statistics across %d fields differ (%d beyond %d ULP).\n",
                length(diffs), n_fields, length(over), SNAPSHOT_ULP_TOLERANCE)
        println("  field                            stat  ULP  rel        baseline                 this run")
        for (k, stat, av, bv, ulp, rel) in diffs
            @printf("  %-32s %-4s %4s %.2e  %.17e  %.17e\n",
                    k, stat, ulp >= 0 ? string(ulp) : "n/a", rel, av, bv)
        end
    end

    ok = isempty(structural) && isempty(over)
    if ok
        println("VERDICT: equivalent — every difference is within $(SNAPSHOT_ULP_TOLERANCE) ULP, " *
                "which is summation order, not a change in the numbers.")
    else
        println("VERDICT: NOT equivalent — the numbers changed. Re-pin only if that was the point.")
    end
    return ok
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
    diffs, structural = compare_snapshots(ref, snap)
    equivalent = report_comparison(diffs, structural)
    # Exit non-zero on a real divergence so this can gate CI, but only after the benchmark
    # below has printed — the timings are wanted either way.
    global SNAPSHOT_FAILED = !equivalent
end

println("\nBenchmarking hot path (1994-2025 3-hourly simulation with 75-year spinup)...")
b = @benchmark run_hotpath($profile, $cf_clim, $cf, $mp) samples=7 seconds=120 evals=1
@printf("  min    = %.3f s\n", minimum(b.times)/1e9)
@printf("  median = %.3f s\n", median(b.times)/1e9)
@printf("  allocs = %d\n", b.allocs)
@printf("  memory = %.1f MiB\n", b.memory/2^20)

# `allocs` and `memory` are exact integers, unlike the timings: quote them when reporting an
# optimization, since they are the part of this output that is not noise.

if isdefined(Main, :SNAPSHOT_FAILED) && SNAPSHOT_FAILED
    error("numerical snapshot diverged from $(basename(SNAP_FILE)) beyond " *
          "$(SNAPSHOT_ULP_TOLERANCE) ULP — see the table above")
end
