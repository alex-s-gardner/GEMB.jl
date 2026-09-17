/*
 * Correctness driver for the GEMB C API.
 *
 * Steps the same column as capi/test/reference.jl, through the compiled library instead of
 * through Julia, and prints the result in the same format. capi/test/runtests.jl diffs the two:
 * both reach the same physics with the same inputs, so they must agree exactly.
 *
 * Also exercises the failure paths that a status-code API has to get right — bad handles, bad
 * option codes, NULL pointers — since a host that mishandles those gets silent corruption rather
 * than an error.
 *
 * Build:
 *   cc -o test_capi capi/test/test_capi.c -Icapi -Lcapi/build -lgemb_safe
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gemb.h"

/* Must match capi/test/column.jl. */
#define N_CELLS 20
#define N_STEPS 240
#define DT 10800.0

static int failures = 0;

static void check(int condition, const char *what) {
    if (!condition) {
        fprintf(stderr, "FAIL: %s\n", what);
        failures++;
    }
}

static void check_status(int status, int expected, const char *what) {
    if (status != expected) {
        fprintf(stderr, "FAIL: %s: got %d, expected %d\n", what, status, expected);
        failures++;
    }
}

/*
 * Print one value the way reference.jl does, so the two outputs diff cleanly. NaN is spelled as
 * Julia spells it — several diagnostics are NaN when the feature they describe is absent, and a
 * formatting difference there would read as a numerical one.
 */
static void emit(const char *name, int index, double value) {
    if (value != value) {
        printf("%-20s %d NaN\n", name, index);
    } else {
        printf("%-20s %d %.17g\n", name, index, value);
    }
}

int main(void) {
    int options[GEMB_N_OPTIONS];
    double numeric[GEMB_N_NUMERIC];
    int flags[GEMB_N_FLAGS];

    check_status(gemb_defaults(options, numeric, flags), GEMB_OK, "gemb_defaults");

    int params = 0, workspace = 0;
    check_status(gemb_params_create(DT, options, numeric, flags, &params), GEMB_OK,
                 "gemb_params_create");
    check_status(gemb_workspace_create(&workspace), GEMB_OK, "gemb_workspace_create");

    double temperature[N_CELLS], dz[N_CELLS], density[N_CELLS], water[N_CELLS];
    double grain_radius[N_CELLS], grain_dendricity[N_CELLS], grain_sphericity[N_CELLS];
    double age[N_CELLS];
    double scalars[2] = {0.0, 0.0};

    for (int i = 0; i < N_CELLS; i++) {
        temperature[i] = 263.0;
        dz[i] = 0.5;
        density[i] = 350.0 + 20.0 * (double)i;
        water[i] = 0.0;
        grain_radius[i] = 0.5;
        grain_dendricity[i] = 0.0;
        grain_sphericity[i] = 0.5;
        age[i] = 0.0;
    }

    double forcing[GEMB_N_FORCING] = {0.0};
    forcing[GEMB_FRC_DT] = DT;
    forcing[GEMB_FRC_TEMPERATURE_AIR] = 272.0;
    forcing[GEMB_FRC_PRESSURE_AIR] = 90000.0;
    forcing[GEMB_FRC_PRECIPITATION] = 0.5;
    forcing[GEMB_FRC_WIND_SPEED] = 5.0;
    forcing[GEMB_FRC_SHORTWAVE_DOWNWARD] = 350.0;
    forcing[GEMB_FRC_LONGWAVE_DOWNWARD] = 280.0;
    forcing[GEMB_FRC_VAPOR_PRESSURE] = 400.0;
    forcing[GEMB_FRC_TEMPERATURE_AIR_MEAN] = 258.0;
    forcing[GEMB_FRC_WIND_SPEED_MEAN] = 5.0;
    forcing[GEMB_FRC_PRECIPITATION_MEAN] = 300.0;
    forcing[GEMB_FRC_TEMPERATURE_OBSERVATION_HEIGHT] = 2.0;
    forcing[GEMB_FRC_WIND_OBSERVATION_HEIGHT] = 10.0;
    forcing[GEMB_FRC_CLOUD_FRACTION] = 0.1;
    /* The two refinements default to the scalars they refine, as the Julia tail constructor does. */
    forcing[GEMB_FRC_ACCUMULATION_MEAN] = 300.0;
    forcing[GEMB_FRC_TEMPERATURE_AIR_EFFECTIVE] = 258.0;

    double flux[GEMB_N_FLUX];
    double flux_total[GEMB_N_FLUX] = {0.0};

    for (int step = 0; step < N_STEPS; step++) {
        int status = gemb_step(N_CELLS, temperature, dz, density, water, grain_radius,
                               grain_dendricity, grain_sphericity, age, scalars, forcing,
                               params, workspace, flux, NULL);
        if (status != GEMB_OK) {
            char message[512];
            gemb_last_error(message, (int)sizeof(message));
            fprintf(stderr, "FAIL: gemb_step returned %d at step %d: %s\n",
                    status, step, message);
            return 1;
        }
        for (int i = 0; i < GEMB_N_FLUX; i++) {
            flux_total[i] += flux[i];
        }
    }

    static const char *const PROFILE_NAMES[8] = {
        "temperature", "dz", "density", "water",
        "grain_radius", "grain_dendricity", "grain_sphericity", "age"
    };
    const double *const profiles[8] = {
        temperature, dz, density, water,
        grain_radius, grain_dendricity, grain_sphericity, age
    };
    for (int p = 0; p < 8; p++) {
        for (int i = 0; i < N_CELLS; i++) {
            emit(PROFILE_NAMES[p], i + 1, profiles[p][i]);
        }
    }

    static const char *const FLUX_NAMES[GEMB_N_FLUX] = {
        "flux_albedo_broadband", "flux_shortwave_net", "flux_heat_flux_sensible",
        "flux_heat_flux_latent", "flux_heat_flux_basal", "flux_longwave_upward",
        "flux_rain", "flux_melt", "flux_runoff", "flux_refreeze", "flux_mass_added",
        "flux_mass_lateral", "flux_mass_blowing_snow", "flux_E_added",
        "flux_densification_from_compaction", "flux_densification_from_melt",
        "flux_percolation_depth", "flux_ice_slab_thickness", "flux_ice_slab_depth",
        "flux_aquifer_thickness", "flux_aquifer_depth"
    };
    for (int i = 0; i < GEMB_N_FLUX; i++) {
        emit(FLUX_NAMES[i], 0, flux_total[i]);
    }

    /* --- Failure paths. ------------------------------------------------------------------ */

    check_status(gemb_step(N_CELLS, temperature, dz, density, water, grain_radius,
                           grain_dendricity, grain_sphericity, age, scalars, forcing,
                           params + 9999, workspace, flux, NULL),
                 GEMB_ERR_BAD_HANDLE, "gemb_step rejects an unknown params handle");

    check_status(gemb_step(N_CELLS, temperature, dz, density, water, grain_radius,
                           grain_dendricity, grain_sphericity, age, scalars, forcing,
                           params, workspace + 9999, flux, NULL),
                 GEMB_ERR_BAD_HANDLE, "gemb_step rejects an unknown workspace handle");

    check_status(gemb_step(1, temperature, dz, density, water, grain_radius,
                           grain_dendricity, grain_sphericity, age, scalars, forcing,
                           params, workspace, flux, NULL),
                 GEMB_ERR_BAD_ARGUMENT, "gemb_step rejects n < 2");

    check_status(gemb_step(N_CELLS, NULL, dz, density, water, grain_radius,
                           grain_dendricity, grain_sphericity, age, scalars, forcing,
                           params, workspace, flux, NULL),
                 GEMB_ERR_BAD_ARGUMENT, "gemb_step rejects a NULL profile");

    /* An out-of-range option code must be caught at the boundary, not inside the physics. */
    int bad_options[GEMB_N_OPTIONS];
    memcpy(bad_options, options, sizeof(options));
    bad_options[GEMB_OPT_ALBEDO_METHOD] = 99;
    int bad_params = 0;
    check_status(gemb_params_create(DT, bad_options, numeric, flags, &bad_params),
                 GEMB_ERR_BAD_OPTION, "gemb_params_create rejects an out-of-range option code");

    char message[512];
    int length = gemb_last_error(message, (int)sizeof(message));
    check(length > 0, "gemb_last_error reports the rejected option");

    check_status(gemb_params_create(0.0, options, numeric, flags, &bad_params),
                 GEMB_ERR_BAD_ARGUMENT, "gemb_params_create rejects dt <= 0");

    check_status(gemb_params_destroy(params + 9999), GEMB_ERR_BAD_HANDLE,
                 "gemb_params_destroy rejects an unknown handle");

    /* Viscosity is written only when the parameter set asked for it. */
    int viscosity_options[GEMB_N_OPTIONS];
    double viscosity_numeric[GEMB_N_NUMERIC];
    int viscosity_flags[GEMB_N_FLAGS];
    gemb_defaults(viscosity_options, viscosity_numeric, viscosity_flags);
    viscosity_flags[GEMB_FLAG_OUTPUT_VISCOSITY] = 1;
    int viscosity_params = 0, viscosity_workspace = 0;
    check_status(gemb_params_create(DT, viscosity_options, viscosity_numeric, viscosity_flags,
                                    &viscosity_params), GEMB_OK,
                 "gemb_params_create with viscosity requested");
    check_status(gemb_workspace_create(&viscosity_workspace), GEMB_OK,
                 "gemb_workspace_create for the viscosity run");

    double viscosity[N_CELLS];
    for (int i = 0; i < N_CELLS; i++) {
        viscosity[i] = -1.0;
    }
    check_status(gemb_step(N_CELLS, temperature, dz, density, water, grain_radius,
                           grain_dendricity, grain_sphericity, age, scalars, forcing,
                           viscosity_params, viscosity_workspace, flux, viscosity),
                 GEMB_OK, "gemb_step with a viscosity buffer");
    int viscosity_written = 1;
    for (int i = 0; i < N_CELLS; i++) {
        if (viscosity[i] == -1.0) {
            viscosity_written = 0;
        }
    }
    check(viscosity_written, "gemb_step fills the viscosity buffer when requested");

    check_status(gemb_params_destroy(viscosity_params), GEMB_OK, "destroy viscosity params");
    check_status(gemb_workspace_destroy(viscosity_workspace), GEMB_OK,
                 "destroy viscosity workspace");
    check_status(gemb_params_destroy(params), GEMB_OK, "gemb_params_destroy");
    check_status(gemb_workspace_destroy(workspace), GEMB_OK, "gemb_workspace_destroy");

    if (failures > 0) {
        fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    fprintf(stderr, "all C API checks passed\n");
    return 0;
}
