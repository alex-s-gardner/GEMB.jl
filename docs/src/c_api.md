# C API

GEMB's column physics can be compiled to a native shared library and called from C, Fortran, or
any language with a C foreign-function interface. This is the route for embedding GEMB in a host
model — an ice-sheet model, a land-surface scheme, a coupled ESM component — without embedding
Julia's package ecosystem alongside it.

The library covers **one column, one timestep**. The host owns the time loop, the forcing, and
all I/O. Model setup (spinup, the climate-derived initial profile) and the `DimStack` output
remain Julia-side: use GEMB as a library for those, then hand the host a column to step.

## Building

Requires Julia 1.12 or later, a C compiler, and
[JuliaC.jl](https://github.com/JuliaLang/JuliaC.jl):

```bash
julia --project=.juliac -e 'using Pkg; Pkg.add(name="JuliaC", version="0.3.9")'
julia --project=capi -e 'using Pkg; Pkg.instantiate()'

JULIAC_PROJECT=.juliac julia capi/build.jl --trim=safe
```

This writes `capi/build/libgemb_safe.{so,dylib}` — about 2 MB — along with the trim verifier log
and a JSON dump of the exported ABI. `--trim=no` skips the reachability analysis and produces a
much larger library (about 310 MB) that does not depend on the analysis succeeding.

Verify the result:

```bash
julia capi/run_tests.jl
```

This checks the declared array orders against the Julia structs, then steps a column 240 times
through the compiled library and confirms the answer is bit-identical to the same run through
`GEMB.gemb_core`.

## Calling it

`capi/gemb.h` is the interface. Take the model defaults and override only what you mean to
change:

```c
#include "gemb.h"

int options[GEMB_N_OPTIONS];
double numeric[GEMB_N_NUMERIC];
int flags[GEMB_N_FLAGS];
gemb_defaults(options, numeric, flags);

options[GEMB_OPT_DENSIFICATION_METHOD] = GEMB_DENSIFICATION_LIGTENBERG;
numeric[GEMB_NUM_DENSITY_ICE] = 910.0;
flags[GEMB_FLAG_SHORTWAVE_SUBSURFACE_ABSORPTION] = 1;

int params, workspace;
gemb_params_create(10800.0, options, numeric, flags, &params);
gemb_workspace_create(&workspace);

for (int step = 0; step < n_steps; step++) {
    /* fill forcing[] for this step from your own data, in GEMB_FRC_* order */
    int status = gemb_step(n, temperature, dz, density, water, grain_radius,
                          grain_dendricity, grain_sphericity, age, scalars, forcing,
                          params, workspace, flux, NULL);
    if (status != GEMB_OK) {
        char message[512];
        gemb_last_error(message, sizeof(message));
        fprintf(stderr, "GEMB failed at step %d: %s\n", step, message);
        break;
    }
    /* the profile arrays now hold the new column state; flux[] holds this step's budget */
}

gemb_workspace_destroy(workspace);
gemb_params_destroy(params);
```

Link against the library:

```bash
cc -o model model.c -Icapi -Lcapi/build -lgemb_safe -Wl,-rpath,$PWD/capi/build
```

`capi/test/test_capi.c` is a complete, working example.

## What a caller has to know

- **The forcing timestep is fixed at parameter creation.** `gemb_params_create` takes `dt` because
  the thermal solver's sub-step divisors are derived from it. Stepping with a different `dt` than
  the parameters were built for is a mistake the API is shaped to prevent.
- **One workspace per column.** The workspace holds solver scratch and is mutated in place, so
  threads stepping different columns concurrently each need their own. Parameters are immutable
  and can be shared freely.
- **The profile arrays are in-out.** They are read, then overwritten with the new state. The cell
  count and the total column depth are invariant across a step, so `n` and the sum of `dz` are
  unchanged on return.
- **Failure is a return code, never an exception.** Julia exceptions cannot cross a C boundary, so
  every entry point catches and returns a status; `gemb_last_error` retrieves the message.
- **Depth ordering and units.** Index 1 is the surface. `capi/gemb.h` documents the unit of every
  forcing and flux entry.

## Scope, and why it is narrow

`--trim=safe` compiles only what is provably reachable from the declared entry points, which
rules out anything requiring runtime type resolution. GEMB's physics core satisfies that already —
`gemb_core` and every function beneath it operate on plain `Vector{Float64}` and concretely-typed
`NamedTuple`s. The layers above it do not: the `DimStack` output, the CF attribute dictionaries,
and `plot_output`'s Makie dependency are all outside what static compilation accepts.

Drawing the boundary at the column step rather than at the full model run is what makes the
library both possible and small. `capi/README.md` records the full exclusion list, the constructs
that break trimming, and how the boundary is defended in CI.
