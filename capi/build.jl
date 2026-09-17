#!/usr/bin/env julia
#
# Build the GEMB C API shared library with JuliaC.
#
# Usage:
#   julia capi/build.jl [--trim=safe|unsafe-warn|no] [--out <path>]
#
# Writes the library next to the requested output path and the full compiler transcript to
# `<output>.log`. The log is the point of the exercise as much as the library is: under
# `--trim` the verifier reports every call it could not resolve, and that report is what says
# whether the compiled surface is still within scope.

const CAPI_DIR = @__DIR__
const ENTRY = joinpath(CAPI_DIR, "gemb_capi.jl")

# JuliaC drives the build but must not be a dependency of the environment being compiled — it
# would be pulled into the library's own dependency graph. `JULIAC_PROJECT` names an
# environment that has JuliaC installed; the default is the app environment
# `Pkg.Apps.add("JuliaC")` creates.
const JULIAC_PROJECT = get(ENV, "JULIAC_PROJECT",
    joinpath(homedir(), ".julia", "environments", "jcdrv"))

trim = "safe"
out = nothing
let i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if startswith(arg, "--trim=")
            global trim = split(arg, '=')[2]
        elseif arg == "--out"
            i += 1
            i <= length(ARGS) || error("--out requires an argument")
            global out = ARGS[i]
        else
            error("unrecognized argument `$arg`")
        end
        i += 1
    end
end

if out === nothing
    out = joinpath(CAPI_DIR, "build", "libgemb_$(replace(trim, "-" => "_"))")
end
mkpath(dirname(out))

# The `juliac` shim `Pkg.Apps.add("JuliaC")` installs into `~/.julia/bin` does the same thing;
# invoking the module directly keeps the build independent of whether that directory is on PATH.
cmd = `$(Base.julia_cmd()) --startup-file=no --project=$JULIAC_PROJECT -m JuliaC
       --project $CAPI_DIR
       --output-lib $out
       --export-abi $(out * ".abi.json")
       --compile-ccallable
       --trim=$trim
       --experimental
       --verbose
       $ENTRY`

logfile = out * ".log"
println("Building at --trim=$trim")
println("  output: $out")
println("  log:    $logfile")

open(logfile, "w") do io
    ok = success(pipeline(cmd; stdout=io, stderr=io))
    println(io)
    println(io, ok ? "BUILD SUCCEEDED" : "BUILD FAILED")
    ok || (println("build failed — see $logfile"); exit(1))
end

println("BUILD SUCCEEDED")
