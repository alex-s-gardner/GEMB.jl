#!/usr/bin/env julia
#
# Full C API gate: boundary tests, then the compiled library checked against Julia.
#
# Assumes the library is already built. Build first with:
#   julia capi/build.jl --trim=safe
#
# Usage:
#   julia capi/run_tests.jl [--lib <path>] [--rtol <x>]
#
# `--rtol` accepts a value that differs from the Julia reference by at most `x` relative to it, in
# place of bit identity. A library compiled for several CPUs (`JULIA_CPU_TARGET`) can run a
# different code path from the Julia session that writes the reference, which reorders
# floating-point operations; names, indices and NaNs must still agree exactly.

const CAPI_DIR = @__DIR__

lib = joinpath(CAPI_DIR, "build", "libgemb_safe")
rtol = 0.0
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--lib"
            i += 1
            i <= length(ARGS) || error("--lib requires an argument")
            global lib = ARGS[i]
        elseif ARGS[i] == "--rtol"
            i += 1
            i <= length(ARGS) || error("--rtol requires an argument")
            global rtol = parse(Float64, ARGS[i])
            rtol >= 0 || error("--rtol must be non-negative, got $rtol")
        else
            error("unrecognized argument `$(ARGS[i])`")
        end
        i += 1
    end
end

libdir = dirname(abspath(lib))
libname = basename(lib)
startswith(libname, "lib") || error("library name must start with `lib`: $libname")
linkname = libname[4:end]

isfile(lib * "." * Base.BinaryPlatforms.platform_dlext()) ||
    error("no library at $lib — run `julia capi/build.jl --trim=safe` first")

println("== boundary tests ==")
run(`$(Base.julia_cmd()) --startup-file=no --project=$CAPI_DIR $(joinpath(CAPI_DIR, "test", "runtests.jl"))`)

println("\n== reference run (Julia) ==")
reference = joinpath(CAPI_DIR, "test", "reference.txt")
open(reference, "w") do io
    run(pipeline(`$(Base.julia_cmd()) --startup-file=no --project=$CAPI_DIR $(joinpath(CAPI_DIR, "test", "reference.jl"))`;
        stdout=io))
end
println("wrote $reference")

println("\n== C driver (compiled library) ==")
cc = get(ENV, "CC", "cc")
driver = joinpath(tempdir(), "gemb_test_capi")
run(`$cc -O2 -Wall -o $driver $(joinpath(CAPI_DIR, "test", "test_capi.c"))
     -I$CAPI_DIR -L$libdir -l$linkname -Wl,-rpath,$libdir`)

actual = joinpath(tempdir(), "gemb_capi_out.txt")
open(actual, "w") do io
    run(pipeline(`$driver`; stdout=io))
end

# The C path reaches the same physics with the same inputs, so the two agree to the last bit or
# the boundary is at fault. `%.17g` round-trips a Float64, so comparing the text compares the
# values.
expected_lines = readlines(reference)
actual_lines = readlines(actual)

# A line is `<name> <index> <value>`. Within `rtol` it agrees when the name and index match and the
# values are equal, both NaN, or within `rtol` of the reference.
function within_rtol(e, a, rtol)
    fe, fa = split(e), split(a)
    length(fe) == length(fa) == 3 && fe[1:2] == fa[1:2] || return false
    ve, va = parse(Float64, fe[3]), parse(Float64, fa[3])
    (isnan(ve) && isnan(va)) || ve == va || abs(va - ve) <= rtol * abs(ve)
end

relative_difference(e, a) = (ve = parse(Float64, split(e)[3]); va = parse(Float64, split(a)[3]);
                             ve == va ? 0.0 : abs(va - ve) / abs(ve))

if expected_lines == actual_lines
    println("\nC API output is bit-identical to the Julia reference ",
            "($(length(expected_lines)) values)")
elseif rtol > 0 && length(expected_lines) == length(actual_lines) &&
       all(within_rtol(e, a, rtol) for (e, a) in zip(expected_lines, actual_lines))
    differing = [(e, a) for (e, a) in zip(expected_lines, actual_lines) if e != a]
    worst = maximum(relative_difference(e, a) for (e, a) in differing)
    println("\nC API output agrees with the Julia reference within rtol = $rtol: ",
            "$(length(differing)) of $(length(expected_lines)) values differ, by at most $worst relative")
else
    println("\nC API output differs from the Julia reference",
            rtol > 0 ? " by more than rtol = $rtol:" : ":")
    for (i, (e, a)) in enumerate(zip(expected_lines, actual_lines))
        (e == a || (rtol > 0 && within_rtol(e, a, rtol))) || println("  line $i\n    julia: $e\n    c:     $a")
    end
    length(expected_lines) == length(actual_lines) ||
        println("  line counts differ: $(length(expected_lines)) vs $(length(actual_lines))")
    exit(1)
end
