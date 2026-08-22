# Does selecting a single average-melt year work, on real Greenland forcing?
#
#   julia --project=bench --threads=auto bench/greenland_year_selection.jl
#
# `n_years=1` is the case the method was originally conceived as: pick the one real year whose
# melt is closest to the record mean. This isolates it against `:average` and against longer
# blocks.
#
# **Deliberately spinup-free.** An earlier attempt (`greenland_block_length.jl`) ran the same
# comparison with a full strict-criteria spinup per site per variant and was abandoned after an
# hour and 4000+ model-years without completing. But the two questions that decide whether a
# selection rule works — does the cycle reproduce the record's melt, and does it reproduce its
# accumulation — are properties of the *forcing*, answerable with one pass of `gemb` per cycle and
# no equilibration at all. That is what this measures, for a few seconds per site.
#
# Reuses the forcing cache written by `greenland_spinup_forcing.jl`.

include("greenland_common.jl")
const BLOCK_LENGTHS = (1, 3, 5)

function analyze(id, mp)
    cf = initialize_forcing(cached_forcing(id))
    mp_last = last_step_parameters(mp)
    profile = initialize_profile(mp, cf)
    cs = GEMB.initialize_climate_summary(cf, mp)
    melt_record = melt_rate(profile, cf, mp_last)

    cycles = Dict{String,Any}("avg" => forcing_climatology(cf))
    for n in BLOCK_LENGTHS
        # `melt` ranking (the shipped default) against `smb` ranking. SMB combines accumulation
        # and melt with the weighting the mass balance gives them, so the expectation is that it
        # trades a little melt fidelity for a lot of accumulation fidelity.
        cycles["m$n"] = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                            n_years=n, rank_by=:model, verbose=false)
        cycles["s$n"] = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                            n_years=n, rank_by=:smb, verbose=false)
    end

    out = Dict{String,Any}()
    for (k, c) in cycles
        cs_c = GEMB.initialize_climate_summary(c, mp)
        out[k] = (melt = melt_rate(profile, c, mp_last),
                  dacc = cs.accumulation > 0 ? cs_c.accumulation / cs.accumulation - 1 : NaN,
                  years_sel = get(DD.metadata(c), :climatology_representative_year, nothing))
    end
    return (id=id, accumulation=cs.accumulation, melt_record=melt_record,
            ratio=cs.accumulation > 0 ? melt_record / cs.accumulation : NaN, v=out)
end

const MP = GEMB.ModelParameters(output_frequency=:last)
const KEYS = ["avg", ("m$n" for n in BLOCK_LENGTHS)..., ("s$n" for n in BLOCK_LENGTHS)...]

@printf("Year-selection comparison (no spinup), %d sites, threads=%d\n",
        length(site_ids()), Threads.nthreads())

results = Any[]
lk = ReentrantLock()
Threads.@threads for id in site_ids()
    try
        r = analyze(id, MP)
        lock(lk) do; push!(results, r); end
    catch err
        lock(lk) do
            @printf("  %-12s FAILED: %s\n", id, first(split(sprint(showerror, err), '\n')))
        end
    end
end
sort!(results, by = r -> r.melt_record)

println("\n", "="^104)
println("MELT recovered as % of the real record")
println("="^104)
@printf("%-12s %9s | %s\n", "site", "melt_rec", join((lpad(k,7) for k in KEYS), " "))
for r in results
    pct(k) = r.melt_record > 1.0 ? 100 * r.v[k].melt / r.melt_record : NaN
    @printf("%-12s %9.1f | %s\n", r.id, r.melt_record,
            join((@sprintf("%6.0f%%", pct(k)) for k in KEYS), " "))
end

println("\n", "="^104)
println("ACCUMULATION offset of the cycle from the record")
println("="^104)
@printf("%-12s %9s | %s\n", "site", "accum", join((lpad(k,7) for k in KEYS), " "))
for r in results
    @printf("%-12s %9.1f | %s\n", r.id, r.accumulation,
            join((@sprintf("%+6.1f%%", 100 * r.v[k].dacc) for k in KEYS), " "))
end

println("\n", "="^104)
println("SUMMARY — melt over the 10 melting sites, accumulation over all 12")
println("="^104)
melters = [r for r in results if r.melt_record > 1.0]
@printf("%-5s %-38s | %s\n", "", "melt retention", "|accumulation offset|")
for k in KEYS
    rets = sort([100 * r.v[k].melt / r.melt_record for r in melters])
    daccs = [abs(100 * r.v[k].dacc) for r in results if isfinite(r.v[k].dacc)]
    @printf("%-5s median %6.0f%%  worst %5.0f%%  <10%%: %d/%-2d | median %5.1f%%  worst %5.1f%%\n",
            k, median(rets), minimum(rets), count(<(10), rets), length(rets),
            median(daccs), maximum(daccs))
end

# The trade the choice actually turns on: melt fidelity against accumulation fidelity. Reported
# as a joint score so a rule that wins on one and loses the other cannot look like a win.
println("\nJoint error (mean of |melt error| and |accumulation error|, melting sites only):")
for k in KEYS
    errs = [(abs(100 * r.v[k].melt / r.melt_record - 100) + abs(100 * r.v[k].dacc)) / 2
            for r in melters if isfinite(r.v[k].dacc)]
    @printf("  %-5s median %5.1f%%   mean %5.1f%%\n", k, median(errs), mean(errs))
end
