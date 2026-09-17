# Reference column run, stepped through `GEMB.gemb_core` in-process.
#
# `capi/test/test_capi.c` performs the same run through the compiled library. The two must agree
# to the last bit: the C path reaches the same physics with the same inputs, so any difference is
# a fault in the boundary, not a tolerance to be widened.
#
# Usage:
#   julia --project=capi capi/test/reference.jl > capi/test/reference.txt

using GEMB

include(joinpath(@__DIR__, "column.jl"))

"""
    run_reference() -> (state, flux_total)

Step the reference column `N_STEPS` times, accumulating each flux term over the run.
"""
function run_reference()
    state, forcing, mp = reference_setup()
    workspace = GEMB.ThermalWorkspace()
    flux_total = zeros(length(FLUX_NAMES))

    for _ in 1:N_STEPS
        state, flux = GEMB.gemb_core(state, forcing, mp, false; thermal_workspace=workspace)
        for (i, name) in enumerate(FLUX_NAMES)
            flux_total[i] += getfield(flux, name)
        end
    end
    return state, flux_total
end

state, flux_total = run_reference()

# `%.17g` round-trips a Float64 exactly, so the comparison is on the values rather than on how
# they were printed.
for name in (:temperature, :dz, :density, :water, :grain_radius, :grain_dendricity,
             :grain_sphericity, :age)
    column = getfield(state, name)
    for i in eachindex(column)
        println(rpad(string(name), 20), " ", i, " ", Printf.@sprintf("%.17g", column[i]))
    end
end
for (i, name) in enumerate(FLUX_NAMES)
    println(rpad("flux_" * string(name), 20), " 0 ", Printf.@sprintf("%.17g", flux_total[i]))
end
