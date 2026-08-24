# GEMB C API

A native shared library exposing GEMB's column physics to C, Fortran, and anything else that
can call a C ABI. Built by [JuliaC.jl](https://github.com/JuliaLang/JuliaC.jl) with `--trim`,
which strips everything not provably reachable from the declared entry points.

```bash
julia capi/build.jl --trim=safe     # -> capi/build/libgemb_safe.{dylib,so}
julia capi/run_tests.jl             # boundary tests + compiled library vs. Julia
```

The build needs Julia ≥ 1.12, a C compiler, and JuliaC in a separate environment:

```bash
julia --project=.juliac -e 'using Pkg; Pkg.add(name="JuliaC", version="0.3.9")'
JULIAC_PROJECT=.juliac julia capi/build.jl --trim=safe
```

`JULIAC_PROJECT` defaults to `~/.julia/environments/jcdrv`. JuliaC is deliberately kept out of
`capi/Project.toml`: it drives the build and must not become a dependency of the code being
compiled.

## Scope

**One column, one timestep.** The library exposes `GEMB.gemb_core` and the leaf physics it
calls. The host owns the time loop, the forcing, and all I/O.

Out of scope, and not compiled:

| Excluded | Why |
| --- | --- |
| `gemb`, `gemb_spinup` | Drive the time loop and build a `DimStack` of output. The host does the former; the latter is not expressible across a C ABI. |
| `initialize_profile`, `initialize_climate_summary`, `steady_state_profile` | Model setup. Pulls in DataInterpolations, untested under `--trim`. |
| `initialize_forcing`, `ClimateForcing` | DimensionalData at the API boundary. The C API takes a flat forcing array instead. |
| The output `DimStack` and CF metadata | `Dict{String,Any}`; no C equivalent. |
| `plot_output` | Makie, which is known not to trim. |

This split is what makes the library possible: `gemb_core` and every physics function beneath it
already operate on `Vector{Float64}` and concretely-typed `NamedTuple`s, with no DimensionalData
anywhere in the tree.

The scope is enforced rather than described. `capi/gemb_capi.jl` names every entry point
explicitly, `--trim=safe` compiles only what those reach, and the build fails if anything in that
tree becomes unresolvable.

## Using it

Take the defaults and override what you mean to change:

```c
#include "gemb.h"

int options[GEMB_N_OPTIONS];
double numeric[GEMB_N_NUMERIC];
int flags[GEMB_N_FLAGS];
gemb_defaults(options, numeric, flags);

options[GEMB_OPT_DENSIFICATION_METHOD] = GEMB_DENSIFICATION_LIGTENBERG;
flags[GEMB_FLAG_SHORTWAVE_SUBSURFACE_ABSORPTION] = 1;

int params, workspace;
gemb_params_create(10800.0, options, numeric, flags, &params);
gemb_workspace_create(&workspace);

for (int step = 0; step < n_steps; step++) {
    /* fill forcing[] for this step from your own data */
    if (gemb_step(n, temperature, dz, density, water, grain_radius, grain_dendricity,
                  grain_sphericity, age, scalars, forcing, params, workspace,
                  flux, NULL) != GEMB_OK) {
        char message[512];
        gemb_last_error(message, sizeof(message));
        fprintf(stderr, "GEMB: %s\n", message);
        break;
    }
}

gemb_workspace_destroy(workspace);
gemb_params_destroy(params);
```

`capi/test/test_capi.c` is a complete worked example. `capi/gemb.h` declares every index constant
and documents the units.

### Rules the host must honor

- **A workspace belongs to one column.** It is mutated in place. Threads stepping different
  columns concurrently need one workspace each. Parameters are immutable and may be shared.
- **Create parameters with the `dt` you will step with.** The thermal sub-step divisors are
  derived from it at creation. This is why parameters are built by a call rather than filled in
  as a struct.
- **The profile arrays are read and then overwritten** with the new state. The cell count and the
  total column depth are invariant, so `n` and `sum(dz)` are unchanged on return.
- **Errors come back as status codes**, never as an exception. A Julia exception must not unwind
  into C, so every entry point catches and returns a code; `gemb_last_error` retrieves what the
  code cannot carry.

### Handles, and the mutable state behind them

`ModelParameters` holds `Symbol` fields and a `Vector{Float64}`; `ThermalWorkspace` owns
growable buffers. Neither has a C layout, so both stay on the Julia side and are named across the
boundary by an opaque `int`.

The handle tables are mutable module-level state, which the model in `src/` deliberately has
none of — that is what lets one `ModelParameters` be read by any number of threads. Confining the
tables to `capi/`, which no normal `using GEMB` loads, is how both properties hold at once.

## Keeping it working

Two gates, both run by `julia capi/run_tests.jl`:

1. **`capi/test/runtests.jl`** holds `gemb.h`, `gemb_capi.jl` and the Julia structs to the same
   field orders. The library passes numbers as bare arrays, so the order *is* the contract: a
   field added to `ClimateForcingStep` or to the flux set fails this test rather than silently
   shifting a column of numbers in someone's model.
2. **`capi/test/test_capi.c`** steps a column 240 times through the compiled library and compares
   against the same run through `gemb_core` in Julia. Both reach the same physics with the same
   inputs, so they must agree **to the last bit** — this is not a tolerance to be widened. It also
   exercises the failure paths, since a host that mishandles a status code gets silent corruption.

`.github/workflows/juliac.yml` runs both on every pull request, over `--trim=safe` and
`--trim=no` on Linux and macOS.

### The CI job is non-blocking, for now

JuliaC is pre-1.0 and its trim behavior tracks upstream Julia, so a hard gate would fail pull
requests over churn that is not this package's. The verifier log is uploaded as an artifact
either way.

**Criterion for making it blocking:** `--trim=safe` green on both platforms for ten consecutive
main-branch runs with no intervening fix. Then delete `continue-on-error: true` from the workflow.

### What breaks trimming

Within `gemb_core`'s call tree, avoid:

- **`Base` functions that reach `dims`-style machinery.** `cumprod` on a vector was the one real
  obstacle found in the whole physics tree: it routes through `Base._accumulate!`, whose
  `CartesianIndices` construction the verifier cannot resolve. `GEMB._extinction_profile!` in
  `src/calculate_shortwave_radiation.jl` replaces it with an explicit running product.
- **Reflection over types.** `fieldnames(ModelParameters)` returns the `UnionAll`'s field list and
  makes the reconstructed value infer as abstract `ModelParameters`, which leaves every downstream
  call unresolvable. Name fields explicitly.
- **`findfirst` and friends over heterogeneous tuples.** Use an explicit loop.
- **`showerror` / `sprint` on a caught exception**, whose type is `Any`. `_record_error` in
  `gemb_capi.jl` reports the exception's type name instead.

What is fine, and was verified: `try`/`catch`, `error()` with string interpolation, `@warn` in the
hot path, `Dict` lookups keyed on a runtime `String` (`densification_lookup_M01`), and the
`Symbol` `if`/`elseif` option chains throughout the physics — `===` on a `Symbol` is a pointer
compare and each branch is concretely typed. New `Symbol` option values need no special care.

## Deployment

`--trim=safe` produces a library of about 2 MB; `--trim=no` produces about 310 MB. Both link
against `libjulia`. For a relocatable distribution, add `--bundle <dir>` to the `juliac` call in
`capi/build.jl`: it copies `libjulia`, the stdlibs and the artifacts alongside the output and sets
relative rpaths.
