# Benchmark + numerical fingerprint for the forcing-cycle and convergence-criteria code added on
# the `calibrate-initial-guess` branch.
#
#   julia --project=bench bench/forcing_cycle_bench.jl
#   SAVE_SNAPSHOT=1 julia --project=bench bench/forcing_cycle_bench.jl    # re-pin
#
# `opt_bench.jl` remains the canonical whole-model driver, but it exercises none of this: it
# calls `forcing_climatology(cf)` (the averaged default) and `gemb_spinup` with no convergence
# criteria, so `method=:representative`, the `rank_by` scorers, the accumulation guard and the FAC
# criteria are all cold in it. This driver measures those paths specifically, and pins their
# numerical output so an optimization can be checked for equivalence.
#
# Three benchmarked stages, chosen because they are what a real run pays:
#
#   1. `forcing_climatology(:average)`   — the default cycle build.
#   2. `forcing_climatology(:representative)` — the same with year ranking and the guard, which
#      integrates one model-year per candidate year.
#   3. `gemb_spinup` under each convergence-criteria combination — the question being whether an
#      unrequested criterion costs anything.

using GEMB
using GEMB_ClimateForcing
using BenchmarkTools
using Statistics
using Printf

const DD = GEMB.DimensionalData
const TIME_STEP_HOURS = 3
const SPINUP_CYCLES = 40          # enough to amortize setup; short enough to benchmark repeatedly

function build_inputs()
    ds = simulate_climate_forcing("test_1", TIME_STEP_HOURS)
    cf = initialize_forcing(ds)
    mp = ModelParameters(output_frequency=:last)
    profile = initialize_profile(mp, cf)
    return cf, mp, profile
end

# Numerical fingerprint of a built cycle: the layer sums plus the selection provenance. Enough to
# catch a change in which years were chosen or in any forcing value carried through.
function cycle_snapshot(cycle)
    d = Dict{Symbol,Float64}()
    for k in sort(collect(keys(cycle)))
        v = filter(isfinite, vec(Float64.(parent(cycle[k]))))
        d[k] = isempty(v) ? NaN : sum(v)
    end
    md = DD.metadata(cycle)
    ys = get(md, :climatology_representative_year, nothing)
    d[:selected_years] = ys === nothing ? NaN : Float64(sum(Int[ys...]))
    d[:accumulation_error] = Float64(get(md, :climatology_accumulation_error, NaN))
    d[:melt_ratio] = Float64(get(md, :climatology_melt_ratio, NaN))
    return d
end

# Fingerprint of a spun-up column, plus the convergence measures the criteria produced.
function spinup_snapshot(g)
    md = DD.metadata(g)
    d = Dict{Symbol,Float64}(
        :cycles => Float64(md[:spinup_cycles]),
        :final_delta_density => md[:spinup_final_delta_density],
        :final_drift_density => md[:spinup_final_drift_density],
        :final_delta_fac => md[:spinup_final_delta_fac],
        :final_drift_fac => md[:spinup_final_drift_fac],
    )
    for k in (:density, :temperature, :dz)
        d[k] = sum(parent(g[k]))
    end
    return d
end

cf, mp, profile = build_inputs()

println("Warmup...")
avg_cycle = forcing_climatology(cf)
rep_cycle = forcing_climatology(cf; method=:representative, model_parameters=mp, verbose=false)
rep_smb = forcing_climatology(cf; method=:representative, model_parameters=mp,
                              rank_by=:smb, verbose=false)

# Spinup variants: the point is whether an unrequested criterion is free.
const SPINUP_VARIANTS = [
    "none"          => (;),
    "density only"  => (convergence_delta_density=1e-9, convergence_drift_density=1e-9),
    "fac only"      => (convergence_delta_fac=1e-12, convergence_drift_fac=1e-12),
    "both"          => (convergence_delta_density=1e-9, convergence_drift_density=1e-9,
                        convergence_delta_fac=1e-12, convergence_drift_fac=1e-12),
]

spin(kw) = gemb_spinup(profile, avg_cycle, mp; max_iterations=SPINUP_CYCLES,
                       thermal_workspace=GEMB.ThermalWorkspace(), kw...)

snap = Dict{Symbol,Float64}()
for (k, v) in cycle_snapshot(avg_cycle); snap[Symbol("avg_", k)] = v; end
for (k, v) in cycle_snapshot(rep_cycle); snap[Symbol("rep_", k)] = v; end
for (k, v) in cycle_snapshot(rep_smb); snap[Symbol("smb_", k)] = v; end
for (name, kw) in SPINUP_VARIANTS
    tag = replace(name, " " => "_")
    for (k, v) in spinup_snapshot(spin(kw)); snap[Symbol("spin_", tag, "_", k)] = v; end
end

println("\nNumerical snapshot:")
for k in sort(collect(keys(snap)))
    @printf("  %-44s %.12e\n", k, snap[k])
end

const SNAP_FILE = joinpath(@__DIR__, "forcing_cycle_snapshot.txt")
if get(ENV, "SAVE_SNAPSHOT", "0") == "1"
    open(SNAP_FILE, "w") do io
        for k in sort(collect(keys(snap)))
            @printf(io, "%s %.15e\n", k, snap[k])
        end
    end
    println("Saved snapshot to $SNAP_FILE")
elseif isfile(SNAP_FILE)
    ref = Dict{Symbol,Float64}()
    for line in eachline(SNAP_FILE)
        p = split(line)
        ref[Symbol(p[1])] = parse(Float64, p[2])
    end
    bad = String[]
    for k in sort(collect(keys(ref)))
        haskey(snap, k) || (push!(bad, "$k: absent from this run"); continue)
        a, b = ref[k], snap[k]
        (isnan(a) && isnan(b)) && continue
        # Relative tolerance rather than bit-equality: these sums run over ~10^4-10^5 elements
        # and are compared across machines, where reduction order differs.
        if isnan(a) != isnan(b) || !isapprox(a, b; rtol=1e-12)
            push!(bad, @sprintf("%s: %.15e -> %.15e", k, a, b))
        end
    end
    for k in sort(collect(keys(snap)))
        haskey(ref, k) || push!(bad, "$k: new field")
    end
    if isempty(bad)
        println("\nSnapshot: equivalent to the tracked baseline.")
    else
        println("\nSnapshot: NOT equivalent — the numbers changed:")
        foreach(m -> println("  ", m), bad)
        global SNAPSHOT_FAILED = true
    end
end

println("\nBenchmarking cycle construction...")
b_avg = @benchmark forcing_climatology($cf) samples=20 seconds=60
@printf("  :average          min %8.2f ms   allocs %9d   %8.1f MiB\n",
        minimum(b_avg.times) / 1e6, b_avg.allocs, b_avg.memory / 2^20)

b_rep = @benchmark forcing_climatology($cf; method=:representative, model_parameters=$mp,
                                       verbose=false) samples=5 seconds=120
@printf("  :representative   min %8.2f ms   allocs %9d   %8.1f MiB\n",
        minimum(b_rep.times) / 1e6, b_rep.allocs, b_rep.memory / 2^20)

# The per-criteria spinup sweep is opt-in: four 3-second spinups x 5 samples is ~15 minutes, and
# it has already answered its question — the four combinations measured 2.92-2.96 s, i.e. within
# noise of each other, so an unrequested criterion is free. Re-run it only when that could have
# changed (a new criterion, or a change to `_update!`/`_column_fac`).
if get(ENV, "BENCH_SPINUP_CRITERIA", "0") == "1"
    println("\nBenchmarking gemb_spinup ($SPINUP_CYCLES cycles) per criteria combination...")
    println("  (an unrequested criterion should cost nothing)")
    for (name, kw) in SPINUP_VARIANTS
        b = @benchmark spin($kw) samples=5 seconds=180
        @printf("  %-14s min %8.3f s   allocs %9d   %8.1f MiB\n",
                name, minimum(b.times) / 1e9, b.allocs, b.memory / 2^20)
    end
else
    println("\nSkipping the per-criteria spinup sweep (~15 min). " *
            "Set BENCH_SPINUP_CRITERIA=1 to run it.")
end

if isdefined(Main, :SNAPSHOT_FAILED) && SNAPSHOT_FAILED
    error("numerical snapshot diverged from $(basename(SNAP_FILE)) — see above")
end
