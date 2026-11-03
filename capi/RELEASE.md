# libgemb

The column physics of the Glacier Energy and Mass Balance model (GEMB) as a C shared library: one
column, one timestep. The host owns the time loop, the forcing and all input and output. The bundle
runs on a machine without Julia.

| Path | Contents |
| --- | --- |
| `include/gemb.h` | The interface: every function, array index and unit |
| `lib/libgemb.{so,dylib}` | The library |
| `lib/libjulia*`, `lib/julia/`, `share/julia/` | The runtime the library loads; keep them beside it |
| `VERSION` | GEMB.jl version and commit, platform, and the Julia and JuliaC that built it |
| `build.log` | The trim verifier's report |

Requirements: glibc 2.17 or later on Linux, macOS 15 or later on Apple silicon. The library runs on
any CPU of its platform.

## Linking

Point the compiler at `include/`, the linker at `lib/`, and the run-time search path at `lib/`:

```bash
B=/path/to/libgemb-<version>-<platform>
cc -o model model.c -I"$B/include" -L"$B/lib" -lgemb -Wl,-rpath,"$B/lib"
```

The bundle can be moved as a whole; its libraries find each other by relative paths. If the host
is linked without an rpath, put `lib/` on `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH` (macOS)
instead.

## Calling

```c
#include "gemb.h"

int options[GEMB_N_OPTIONS];
double numeric[GEMB_N_NUMERIC];
int flags[GEMB_N_FLAGS];
gemb_defaults(options, numeric, flags);       /* the model defaults; override what you mean to */

int params, workspace;
gemb_params_create(10800.0, options, numeric, flags, &params);   /* dt of the forcing, s */
gemb_workspace_create(&workspace);                               /* one per column */

for (int step = 0; step < n_steps; step++) {
    /* fill forcing[] for this step, in GEMB_FRC_* order */
    if (gemb_step(n, temperature, dz, density, water, grain_radius, grain_dendricity,
                  grain_sphericity, age, scalars, forcing, params, workspace,
                  flux, NULL) != GEMB_OK) {
        char message[512];
        gemb_last_error(message, sizeof(message));
        fprintf(stderr, "GEMB: %s\n", message);
        break;
    }
    /* the profile arrays hold the new column; flux[] holds this step's budget */
}

gemb_workspace_destroy(workspace);
gemb_params_destroy(params);
```

- The forcing timestep is fixed when the parameters are created.
- Parameters can be shared by any number of threads; each concurrently stepped column needs its own
  workspace.
- The profile arrays are read and overwritten. Index 1 is the surface; the cell count and the total
  depth do not change.
- Every function returns a status code; `gemb_last_error` gives the message.

Model setup (spin-up and the initial profile) is not in the library; build the column with
GEMB.jl and hand it to the host. The C API page of the documentation describes the interface in
full: https://alex-s-gardner.github.io/GEMB.jl/stable/c_api
