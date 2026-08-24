/*
 * GEMB C API — Glacier Energy and Mass Balance model, one column, one timestep.
 *
 * The compiled surface is deliberately narrow: this library advances a single firn column by a
 * single forcing step. The host owns the time loop, the forcing, and all I/O. Model setup
 * (spinup, climate-derived profile initialization) and the DimensionalData output stack are not
 * part of it — see capi/README.md.
 *
 * The array orders declared below are the contract. They are checked against the Julia structs
 * by capi/test/runtests.jl, so a field added to ClimateForcingStep or to the flux set fails that
 * test rather than silently shifting a column of numbers.
 *
 * Build:  julia capi/build.jl --trim=safe
 */

#ifndef GEMB_H
#define GEMB_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Status codes. Every entry point returns one of these. ------------------------------- */

#define GEMB_OK                0
#define GEMB_ERR_BAD_ARGUMENT (-1)  /* NULL pointer, or n < 2, or dt <= 0 */
#define GEMB_ERR_BAD_HANDLE   (-2)  /* no live object for that handle */
#define GEMB_ERR_BAD_OPTION   (-3)  /* an option code is out of range for its field */
#define GEMB_ERR_PHYSICS      (-4)  /* the step failed; see gemb_last_error */
#define GEMB_ERR_UNEXPECTED   (-5)

/* --- Array sizes. ------------------------------------------------------------------------- */

#define GEMB_N_OPTIONS  15
#define GEMB_N_NUMERIC  25
#define GEMB_N_FLAGS     3
#define GEMB_N_FORCING  22
#define GEMB_N_FLUX     21

/* --- options[] indices, and the codes each accepts. Codes are 1-based. ------------------- */

#define GEMB_OPT_DENSIFICATION_METHOD         0
#define GEMB_OPT_DENSIFICATION_COEFFS_M01     1
#define GEMB_OPT_DENSIFICATION_ACCUMULATION   2
#define GEMB_OPT_MEAN_TEMPERATURE_METHOD      3
#define GEMB_OPT_NEW_SNOW_METHOD              4
#define GEMB_OPT_EMISSIVITY_METHOD            5
#define GEMB_OPT_THERMAL_CONDUCTIVITY_METHOD  6
#define GEMB_OPT_HEAT_CAPACITY_METHOD         7
#define GEMB_OPT_RAIN_HEAT_CAPACITY           8
#define GEMB_OPT_GRAIN_GROWTH_METHOD          9
#define GEMB_OPT_WATER_IRREDUCIBLE_METHOD    10
#define GEMB_OPT_MELT_GEOMETRY               11
#define GEMB_OPT_RUNOFF_METHOD               12
#define GEMB_OPT_ALBEDO_METHOD               13
#define GEMB_OPT_BLOWING_SNOW_METHOD         14

/* densification_method */
#define GEMB_DENSIFICATION_HERRON_LANGWAY  1
#define GEMB_DENSIFICATION_ARTHERN         2
#define GEMB_DENSIFICATION_ARTHERN_B       3
#define GEMB_DENSIFICATION_BARNOLA1991     4
#define GEMB_DENSIFICATION_CROCUS          5
#define GEMB_DENSIFICATION_CROCUS_PURE     6
#define GEMB_DENSIFICATION_GSFC2020        7
#define GEMB_DENSIFICATION_SIMONSEN2013    8
#define GEMB_DENSIFICATION_LIGTENBERG      9

/* densification_coeffs_M01 (consulted only by GEMB_DENSIFICATION_LIGTENBERG) */
#define GEMB_COEFFS_ANT_ERA5_GS_SW0     1
#define GEMB_COEFFS_ANT_ERA5V4_PAOLO23  2
#define GEMB_COEFFS_ANT_ERA5_BF_SW1     3
#define GEMB_COEFFS_ANT_RACMO_GS_SW0    4
#define GEMB_COEFFS_ANT_LIGTENBERG      5
#define GEMB_COEFFS_GRE_ERA5_GS_SW0     6
#define GEMB_COEFFS_GRE_RACMO_GS_SW0    7
#define GEMB_COEFFS_GRE_RACMO_GB_SW1    8
#define GEMB_COEFFS_GRE_KUIPERSMUNNEKE   9

/* densification_accumulation */
#define GEMB_ACCUMULATION_SNOWFALL       1
#define GEMB_ACCUMULATION_PRECIPITATION  2

/* mean_temperature_method */
#define GEMB_MEAN_TEMPERATURE_ARITHMETIC  1
#define GEMB_MEAN_TEMPERATURE_ARRHENIUS   2

/* new_snow_method */
#define GEMB_NEW_SNOW_CONSTANT150      1
#define GEMB_NEW_SNOW_CONSTANT315      2
#define GEMB_NEW_SNOW_CONSTANT350      3
#define GEMB_NEW_SNOW_FAUSTO           4
#define GEMB_NEW_SNOW_FAUSTO_FIT       5
#define GEMB_NEW_SNOW_PAHAUT           6
#define GEMB_NEW_SNOW_KASPERS          7
#define GEMB_NEW_SNOW_KUIPERSMUNNEKE   8

/* emissivity_method */
#define GEMB_EMISSIVITY_UNIFORM                    1
#define GEMB_EMISSIVITY_GRAIN_RADIUS_THRESHOLD     2
#define GEMB_EMISSIVITY_GRAIN_RADIUS_W_THRESHOLD   3

/* thermal_conductivity_method */
#define GEMB_CONDUCTIVITY_STURM          1
#define GEMB_CONDUCTIVITY_CALONNE        2
#define GEMB_CONDUCTIVITY_CALONNE2019    3
#define GEMB_CONDUCTIVITY_CALONNE2019AIR 4
#define GEMB_CONDUCTIVITY_MARCHENKO2019  5

/* heat_capacity_method */
#define GEMB_HEAT_CAPACITY_CONSTANT        1
#define GEMB_HEAT_CAPACITY_CUFFEYPATERSON  2

/* rain_heat_capacity */
#define GEMB_RAIN_HEAT_CAPACITY_WATER  1
#define GEMB_RAIN_HEAT_CAPACITY_ICE    2

/* grain_growth_method */
#define GEMB_GRAIN_GROWTH_MARBOUTY  1
#define GEMB_GRAIN_GROWTH_ARTHERN   2
#define GEMB_GRAIN_GROWTH_HYBRID    3

/* water_irreducible_method */
#define GEMB_IRREDUCIBLE_CONSTANT         1
#define GEMB_IRREDUCIBLE_COLEOULESAFFRE   2

/* melt_geometry */
#define GEMB_MELT_GEOMETRY_THICKNESS  1
#define GEMB_MELT_GEOMETRY_DENSITY    2

/* runoff_method. GEMB_RUNOFF_DARCY requires numeric[GEMB_NUM_SURFACE_SLOPE] > 0. */
#define GEMB_RUNOFF_INSTANTANEOUS  1
#define GEMB_RUNOFF_ZUOOERLEMANS   2
#define GEMB_RUNOFF_DARCY          3

/* albedo_method */
#define GEMB_ALBEDO_NONE               1
#define GEMB_ALBEDO_GARDNERSHARP       2
#define GEMB_ALBEDO_BRUNLEFEBRE        3
#define GEMB_ALBEDO_GREUELLKONZELMANN  4

/* blowing_snow_method */
#define GEMB_BLOWING_SNOW_NONE    1
#define GEMB_BLOWING_SNOW_CROCUS  2

/* --- numeric[] indices. --------------------------------------------------------------- */

#define GEMB_NUM_DENSITY_ICE                        0
#define GEMB_NUM_RAIN_TEMPERATURE_THRESHOLD         1
#define GEMB_NUM_EMISSIVITY                         2
#define GEMB_NUM_EMISSIVITY_GRAIN_RADIUS_LARGE      3
#define GEMB_NUM_EMISSIVITY_GRAIN_RADIUS_THRESHOLD  4
#define GEMB_NUM_SURFACE_ROUGHNESS_EFFECTIVE_RATIO  5
#define GEMB_NUM_HEAT_CAPACITY_ICE                  6
#define GEMB_NUM_WATER_IRREDUCIBLE_SATURATION       7
#define GEMB_NUM_IMPERMEABLE_DENSITY                8
#define GEMB_NUM_IMPERMEABLE_THICKNESS              9
#define GEMB_NUM_PORE_SATURATION_MAX               10
#define GEMB_NUM_ALBEDO_DENSITY_THRESHOLD          11
#define GEMB_NUM_ALBEDO_SNOW                       12
#define GEMB_NUM_ALBEDO_ICE                        13
#define GEMB_NUM_ALBEDO_FIXED                      14
#define GEMB_NUM_COLUMN_ZTOP                       15
#define GEMB_NUM_COLUMN_DZTOP                      16
#define GEMB_NUM_COLUMN_DZMIN                      17
#define GEMB_NUM_COLUMN_DZMAX                      18
#define GEMB_NUM_COLUMN_DEPTH_MAX                  19
#define GEMB_NUM_COLUMN_ZY                         20
#define GEMB_NUM_HORIZONTAL_STRAIN_RATE            21
#define GEMB_NUM_SURFACE_SLOPE                     22
#define GEMB_NUM_DRIFT_RATE                        23
#define GEMB_NUM_THERMAL_EXPLICIT_SAFETY_FACTOR    24

/* --- flags[] indices. Nonzero is true. ----------------------------------------------- */

#define GEMB_FLAG_SHORTWAVE_SUBSURFACE_ABSORPTION  0
#define GEMB_FLAG_OUTPUT_VISCOSITY                 1
#define GEMB_FLAG_BLOWING_SNOW_SUBLIMATION         2

/* --- forcing[] indices. --------------------------------------------------------------- */

#define GEMB_FRC_DT                              0  /* [s]        step length */
#define GEMB_FRC_TEMPERATURE_AIR                 1  /* [K]                    */
#define GEMB_FRC_PRESSURE_AIR                    2  /* [Pa]                   */
#define GEMB_FRC_PRECIPITATION                   3  /* [kg m-2]   over the step */
#define GEMB_FRC_WIND_SPEED                      4  /* [m s-1]                */
#define GEMB_FRC_SHORTWAVE_DOWNWARD              5  /* [W m-2]                */
#define GEMB_FRC_LONGWAVE_DOWNWARD               6  /* [W m-2]                */
#define GEMB_FRC_VAPOR_PRESSURE                  7  /* [Pa]                   */
#define GEMB_FRC_TEMPERATURE_AIR_MEAN            8  /* [K]        climatological */
#define GEMB_FRC_WIND_SPEED_MEAN                 9  /* [m s-1]    climatological */
#define GEMB_FRC_PRECIPITATION_MEAN             10  /* [kg m-2 yr-1]          */
#define GEMB_FRC_TEMPERATURE_OBSERVATION_HEIGHT 11  /* [m]                    */
#define GEMB_FRC_WIND_OBSERVATION_HEIGHT        12  /* [m]                    */
#define GEMB_FRC_BLACK_CARBON_SNOW              13  /* [ng g-1]               */
#define GEMB_FRC_BLACK_CARBON_ICE               14  /* [ng g-1]               */
#define GEMB_FRC_CLOUD_OPTICAL_THICKNESS        15  /* [-]                    */
#define GEMB_FRC_SOLAR_ZENITH_ANGLE             16  /* [rad]                  */
#define GEMB_FRC_SHORTWAVE_DOWNWARD_DIFFUSE     17  /* [W m-2]                */
#define GEMB_FRC_CLOUD_FRACTION                 18  /* [-]                    */
#define GEMB_FRC_ACCUMULATION_MEAN              19  /* [kg m-2 yr-1] snowfall only */
#define GEMB_FRC_TEMPERATURE_AIR_EFFECTIVE      20  /* [K]        Arrhenius-weighted */
#define GEMB_FRC_SNOW_DRIFT                     21  /* [kg m-2 yr-1] + is erosion */

/* --- flux[] indices, written by gemb_step. -------------------------------------------- */

#define GEMB_FLUX_ALBEDO_BROADBAND                0  /* [-]       */
#define GEMB_FLUX_SHORTWAVE_NET                   1  /* [W m-2]   */
#define GEMB_FLUX_HEAT_FLUX_SENSIBLE              2  /* [W m-2]   */
#define GEMB_FLUX_HEAT_FLUX_LATENT                3  /* [W m-2]   */
#define GEMB_FLUX_HEAT_FLUX_BASAL                 4  /* [W m-2]   */
#define GEMB_FLUX_LONGWAVE_UPWARD                 5  /* [W m-2]   */
#define GEMB_FLUX_RAIN                            6  /* [kg m-2]  */
#define GEMB_FLUX_MELT                            7  /* [kg m-2]  */
#define GEMB_FLUX_RUNOFF                          8  /* [kg m-2]  */
#define GEMB_FLUX_REFREEZE                        9  /* [kg m-2]  */
#define GEMB_FLUX_MASS_ADDED                     10  /* [kg m-2]  basal, signed */
#define GEMB_FLUX_MASS_LATERAL                   11  /* [kg m-2]  strain export */
#define GEMB_FLUX_MASS_BLOWING_SNOW              12  /* [kg m-2]  */
#define GEMB_FLUX_E_ADDED                        13  /* [J m-2]   */
#define GEMB_FLUX_DENSIFICATION_FROM_COMPACTION  14  /* [m]       */
#define GEMB_FLUX_DENSIFICATION_FROM_MELT        15  /* [m]       */
#define GEMB_FLUX_PERCOLATION_DEPTH              16  /* [m]       */
#define GEMB_FLUX_ICE_SLAB_THICKNESS             17  /* [m]       */
#define GEMB_FLUX_ICE_SLAB_DEPTH                 18  /* [m]       */
#define GEMB_FLUX_AQUIFER_THICKNESS              19  /* [m]       */
#define GEMB_FLUX_AQUIFER_DEPTH                  20  /* [m]       */

/* --- Entry points. --------------------------------------------------------------------- */

/*
 * Fill the three spec arrays with the model defaults. Take these and override only what you
 * mean to change, rather than assembling a spec from scratch.
 */
int gemb_defaults(int *options, double *numeric, int *flags);

/*
 * Build a parameter set for a forcing step of dt seconds; writes an opaque handle.
 *
 * dt belongs here because the thermal sub-step divisors are derived from it at creation. The
 * parameters are immutable once created, so one handle may be shared across threads. Release
 * with gemb_params_destroy.
 */
int gemb_params_create(double dt, const int *options, const double *numeric,
                       const int *flags, int *handle_out);
int gemb_params_destroy(int handle);

/*
 * Create thermal-solver scratch; writes an opaque handle.
 *
 * One workspace belongs to one column: it is mutated in place, so threads stepping different
 * columns concurrently must each hold their own. Its buffers size themselves to the column on
 * first use and are reused thereafter, so create one per column and keep it for the run.
 */
int gemb_workspace_create(int *handle_out);
int gemb_workspace_destroy(int handle);

/*
 * Advance one n-cell column by one forcing step.
 *
 * The eight profile arrays (each n doubles) and scalars are read, then overwritten with the new
 * state. The column holds its cell count and its total depth, so n and the sum of dz are
 * unchanged on return.
 *
 *   scalars       2 doubles: {evaporation_condensation, melt_surface}
 *   forcing       GEMB_N_FORCING doubles, in GEMB_FRC_* order
 *   flux_out      GEMB_N_FLUX doubles, in GEMB_FLUX_* order
 *   viscosity_out n doubles, or NULL. Written only when the parameter set set
 *                 GEMB_FLAG_OUTPUT_VISCOSITY.
 *
 * Returns GEMB_OK, or GEMB_ERR_PHYSICS when the step failed a conservation check — call
 * gemb_last_error for which.
 */
int gemb_step(int n,
              double *temperature,       /* [K]       */
              double *dz,                /* [m]       */
              double *density,           /* [kg m-3]  */
              double *water,             /* [kg m-2]  */
              double *grain_radius,      /* [mm]      */
              double *grain_dendricity,  /* [-]       */
              double *grain_sphericity,  /* [-]       */
              double *age,               /* [days]    */
              double *scalars,
              const double *forcing,
              int params_handle,
              int workspace_handle,
              double *flux_out,
              double *viscosity_out);

/*
 * Copy the most recent error message into buf, NUL-terminated. Returns its length in bytes
 * excluding the terminator; a longer message is truncated to len - 1.
 */
int gemb_last_error(char *buf, int len);

#ifdef __cplusplus
}
#endif

#endif /* GEMB_H */
