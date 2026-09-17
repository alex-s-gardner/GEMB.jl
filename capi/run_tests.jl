#!/usr/bin/env julia
#
# Full C API gate: boundary tests, then the compiled library checked against Julia.
#
# Assumes the library is already built. Build first with:
#   julia capi/build.jl --trim=safe
#
# Usage:
#   julia capi/run_tests.jl [--lib <path>]

const CAPI_DIR = @__DIR__

lib = joinpath(CAPI_DIR, "build", "libgemb_safe")
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--lib"
            i += 1
            i <= length(ARGS) || error("--lib requires an argument")
            global lib = ARGS[i]
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

if expected_lines == actual_lines
    println("\nC API output is bit-identical to the Julia reference ",
            "($(length(expected_lines)) values)")
else
    println("\nC API output differs from the Julia reference:")
    for (i, (e, a)) in enumerate(zip(expected_lines, actual_lines))
        e == a || println("  line $i\n    julia: $e\n    c:     $a")
    end
    length(expected_lines) == length(actual_lines) ||
        println("  line counts differ: $(length(expected_lines)) vs $(length(actual_lines))")
    exit(1)
end
