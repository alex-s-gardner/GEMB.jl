# Capture the pre-change `initialize_profile` baseline as the bit-identity gate
# for the universal-initialization work (plan verification steps 1-2).
#
# Run BEFORE any source edits to write the snapshot, then again AFTER to compare:
#
#   julia --project=test bench/baseline_init_snapshot.jl write
#   julia --project=test bench/baseline_init_snapshot.jl check
#
# `write` records the legacy `steady_state=false, depth_autoadjust=false` output;
# `check` asserts the new `constant_density=true, constant_temperature=true`
# escape hatch reproduces it exactly (`==`, not `isapprox`).

using GEMB
using GEMB_ClimateForcing
using Printf
using Serialization

const SNAPSHOT = joinpath(@__DIR__, "baseline_init_snapshot.jls")

const LAYERS = (:dz, :temperature, :density, :water, :grain_radius,
                :grain_dendricity, :grain_sphericity, :albedo, :albedo_diffuse)

"""Build the synthetic forcing the regression pin uses (test_1, 3-hourly)."""
function synthetic_forcing()
    ds = simulate_climate_forcing("test_1", 3)
    return GEMB.initialize_forcing(ds)
end

"""Capture the legacy pure-ice initialization."""
function capture_legacy()
    cf = synthetic_forcing()
    mp = GEMB.ModelParameters(output_frequency=:daily)
    profile, mp_out = GEMB.initialize_profile(mp, cf;
        steady_state=false, depth_autoadjust=false)
    return snapshot(profile, mp_out)
end

"""Capture the new escape-hatch initialization."""
function capture_hatch()
    cf = synthetic_forcing()
    mp = GEMB.ModelParameters(output_frequency=:daily)
    profile, mp_out = GEMB.initialize_profile(mp, cf;
        constant_density=true, constant_temperature=true)
    return snapshot(profile, mp_out)
end

snapshot(profile, mp) = (
    layers = NamedTuple(k => collect(profile[k]) for k in LAYERS),
    column_depth = mp.column_depth,
)

# The Herron–Langway steady-state cases exercised directly by
# test/test_gemb_spinup.jl:137-152. `herron_langway_steady_state` is being
# replaced by the universal marcher, so pin its output here: the new marcher with
# `densification_method = :HerronLangway` must reproduce these exactly.
const HL_CASES = (
    (label = "zc_default",  z = -collect(0.025:0.05:10.0),  T = 253.0, A = 300.0),
    (label = "deep_250m",   z = -collect(0.5:1.0:250.0),    T = 253.0, A = 300.0),
    (label = "zero_accum",  z = -collect(0.025:0.05:10.0),  T = 253.0, A = 0.0),
)

function capture_hl()
    mp = GEMB.ModelParameters()
    ρ0 = GEMB.fresh_snow_density(mp, 253.0, 300.0, 5.0)
    return NamedTuple(Symbol(c.label) =>
        GEMB.herron_langway_steady_state(c.z, c.T, c.A, ρ0, mp.density_ice)
        for c in HL_CASES)
end

function do_write()
    snap = capture_legacy()
    hl = capture_hl()
    serialize(SNAPSHOT, (profile = snap, hl = hl))
    n = length(snap.layers.dz)
    @printf("wrote %s\n  %d cells, column_depth = %.6f m\n", SNAPSHOT, n, snap.column_depth)
    for k in LAYERS
        v = snap.layers[k]
        @printf("  %-18s sum=%.17g  min=%.17g  max=%.17g\n", k, sum(v), minimum(v), maximum(v))
    end
    println("  Herron–Langway steady-state cases:")
    for c in HL_CASES
        v = hl[Symbol(c.label)]
        @printf("  %-18s n=%d sum=%.17g  min=%.17g  max=%.17g\n",
            c.label, length(v), sum(v), minimum(v), maximum(v))
    end
end

function do_check()
    isfile(SNAPSHOT) || error("no snapshot at $SNAPSHOT; run `write` on the pre-change code first")
    stored = deserialize(SNAPSHOT)
    ref = stored.profile
    new = capture_hatch()

    ok = true
    if new.column_depth != ref.column_depth
        @printf("MISMATCH column_depth: ref=%.17g new=%.17g\n", ref.column_depth, new.column_depth)
        ok = false
    end
    for k in LAYERS
        a, b = ref.layers[k], new.layers[k]
        if length(a) != length(b)
            @printf("MISMATCH %-18s length ref=%d new=%d\n", k, length(a), length(b))
            ok = false
        elseif a != b                      # exact equality, not isapprox
            d = maximum(abs.(a .- b))
            @printf("MISMATCH %-18s max|Δ|=%.17g\n", k, d)
            ok = false
        else
            @printf("ok       %-18s (%d cells, exact)\n", k, length(a))
        end
    end

    # The universal marcher with :HerronLangway must reproduce the old
    # `herron_langway_steady_state` exactly.
    if isdefined(GEMB, :steady_state_density)
        mp = GEMB.ModelParameters(densification_method=:HerronLangway)
        ρ0 = GEMB.fresh_snow_density(mp, 253.0, 300.0, 5.0)
        for c in HL_CASES
            a = stored.hl[Symbol(c.label)]
            b = GEMB.steady_state_density(c.z, c.T, c.A, ρ0, mp.density_ice, mp)
            if a != b
                @printf("MISMATCH hl:%-14s max|Δ|=%.17g\n", c.label,
                    length(a) == length(b) ? maximum(abs.(a .- b)) : NaN)
                ok = false
            else
                @printf("ok       hl:%-14s (%d cells, exact)\n", c.label, length(a))
            end
        end
    else
        println("note: GEMB.steady_state_density not defined yet — skipping H–L comparison")
    end

    if ok
        println("\nPASS: escape hatch is bit-identical to the legacy initialization.")
    else
        error("escape hatch is NOT bit-identical — this is a bug in the escape-hatch path, " *
              "not grounds to re-baseline the regression pin.")
    end
end

const mode = isempty(ARGS) ? "write" : ARGS[1]
mode == "write" ? do_write() :
mode == "check" ? do_check() :
error("unknown mode $mode (expected \"write\" or \"check\")")
