# Physical constants used throughout GEMB
# All values match the MATLAB implementation exactly

const CtoK = 273.15        # Celsius to Kelvin conversion [K]

# Specific heat capacity of snow/firn/ice [J kg-1 K-1]. This seeds
# `ModelParameters.heat_capacity_ice` and is *not* authoritative — the model reads
# `heat_capacity(mp, T)`, which honours `mp.heat_capacity_method`. The MATLAB value 2102
# is c_p at the melting point.
const HEAT_CAPACITY_ICE_DEFAULT = 2102.0

# Cuffey and Paterson (2010) eq. 9.1: c_p(T) = a + b·T, T in kelvin.
const HEAT_CAPACITY_CUFFEY_A = 152.5
const HEAT_CAPACITY_CUFFEY_B = 7.122

const HEAT_CAPACITY_AIR = 1005.0  # Specific heat capacity of air [J kg-1 K-1]

# Specific heat capacity of liquid water [J kg-1 K-1]. Carries the sensible heat of
# *above-freezing* liquid entering the column — i.e. rain — via the two-argument
# `specific_enthalpy_water(mp, T)`, under the `rain_heat_capacity = :water` default. Pore water
# in GEMB is carried strictly at the melting point, so this does not appear in the percolation
# or refreeze budgets; `heat_capacity(mp, T)` remains the ice-matrix value everywhere else.
# Value from the Community Firn Model (`solver.py`, `c_liq`).
const HEAT_CAPACITY_WATER = 4219.9

const LF = 0.3345e6        # Latent heat of fusion [J kg-1]
const LV = 2.495e6         # Latent heat of vaporization [J kg-1]
const LS = 2.8295e6        # Latent heat of sublimation [J kg-1]
const SB = 5.67e-8         # Stefan-Boltzmann constant [W m-2 K-4]
const GRAVITY = 9.81       # Gravitational acceleration [m s-2]
const R_GAS = 8.314        # Universal gas constant [J mol-1 K-1]
const DENSITY_WATER = 1000.0  # Density of water [kg m-3]
const VON_KARMAN = 0.4     # Von Karman constant [-]

const SECONDS_PER_YEAR = 365.25 * 86400.0  # Seconds in a (Julian) year [s]

# Year length assumed by the densification rate coefficients, in days. Deliberately 365, not
# `SECONDS_PER_YEAR / 86400` (365.25): every densification scheme's `c` is applied as
# `c / DAYS_PER_YEAR_DENSIFICATION * dt` with `dt` in days, so any scheme published with SI
# per-second coefficients must be converted against *this* value for the two to cancel.
# Changing it rescales every densification rate by 0.07%; it is also ~0.07% of the offset
# against the Community Firn Model, which uses a 365.25-day year.
const DAYS_PER_YEAR_DENSIFICATION = 365.0

# Snow/firn state values shared by the accumulation physics and the steady-state
# initial guess, so a freshly initialized column and freshly fallen snow agree.
const RE_NEW_SNOW = 0.05   # New snow effective grain radius [mm]
const GDN_NEW_SNOW = 1.0   # New snow dendricity [-]
const GSP_NEW_SNOW = 0.5   # New snow sphericity [-]
const GRAIN_DIAMETER_MAX = 5.0  # Non-spherical grain diameter cap [mm] (calculate_grain_size)
const GRAIN_RADIUS_ICE = GRAIN_DIAMETER_MAX / 2  # Grain radius of bare ice [mm]
const DENSITY_PORE_CLOSEOFF = 830.0  # Pore close-off density [kg m-3]

# Activation energy for grain growth [J mol-1], Arthern et al. (2010). Sets the
# `exp(Eg/(R·T_mean))` grain-growth factor in the accumulation-driven densification schemes
# (`_arthern_scaled_c`, `_gsfc2020_c` in `calculate_density`). Named rather than restated at
# each site so a densification scheme's assumed grain growth cannot drift from the others'.
const GRAIN_GROWTH_EG = 42400.0
# Grain-growth rate coefficient [m2 s-1], Arthern et al. (2010): dr²/dt = kgr·exp(-Eg/RT).
# Read by `calculate_grain_size` under `grain_growth_method` `:Arthern`/`:hybrid`. SI units, so
# a conversion to the module's mm² working units is applied at that call site.
const GRAIN_GROWTH_KGR = 1.3e-7
# Density above which Marbouty's (1980) temperature-gradient growth law returns zero [kg m-3],
# i.e. the ceiling of its calibration range. Read by the `H` coefficient in `_marbouty_Q`,
# which ramps to zero across [150, this].
const DENSITY_MARBOUTY_MAX = 400.0

# Density separating the two densification stages [kg m-3]. Every scheme in
# `calculate_density` switches rate coefficients here: it is the transition from
# grain-settling-dominated to creep-dominated compaction, and all the published two-stage
# coefficient pairs (Herron & Langway, Arthern, Ligtenberg, Simonsen, GSFC) are fitted
# against this same split.
const DENSITY_STAGE_TRANSITION = 550.0

# Reference temperature of the Calonne et al. (2019) conductivity regressions [K]. Their
# eq. 5 rescales each regime fit by the ice (and air) conductivity *relative to its value
# here*, so both ratios are 1 at this temperature and the fits return their published
# values. Deliberately not `CtoK`: this is a property of the regression, not the melting
# point, and the 3 K gap between them is not a rounding artifact.
const CONDUCTIVITY_T_REF = 270.15

# Reid (1966) thermal conductivity of air, `K = a·T^1.5/(b + T)` [W m-1 K-1]. Read only by
# `_thermal_conductivity_air`, which supplies the `K_air(T)/K_air(T_ref)` factor on the snow
# branch of `:Calonne2019Air`.
const K_AIR_REID_A = 2.334e-3
const K_AIR_REID_B = 164.54

# Fresh-snow density regression of Fausto et al. (2018), `ρ0 = a + b·(T - CtoK)`
# [kg m-3, T in K], calibrated for Greenland. Read by `fresh_snow_density` under
# `new_snow_method = :FaustoFit`. The older `:Fausto` option returns
# `DENSITY_NEW_SNOW_FAUSTO_CONSTANT`, which is this same fit evaluated at one fixed
# temperature (T ≈ 256.2 K) — see that function's docstring.
const DENSITY_NEW_SNOW_FAUSTO_A = 362.1
const DENSITY_NEW_SNOW_FAUSTO_B = 2.78
const DENSITY_NEW_SNOW_FAUSTO_CONSTANT = 315.0

# Surface roughness lengths (Bougamont, 2005), shared by the transient surface
# energy balance (`calculate_temperature`) and the initial-guess one
# (`_seb_annual_melt`), which blends between the dry-snow and ice values.
const Z0_SNOW_DRY = 0.00012  # 0.12 mm, dry snow [m]
const Z0_SNOW_WET = 0.0013   # 1.3 mm, wet snow [m]
const Z0_ICE = 0.0032        # 3.2 mm, ice [m]

# How close to `density_ice` counts as "firn has become ice" when sizing the initial
# column. Densification laws approach ice *asymptotically*, so no finite depth ever
# reaches it to within a float epsilon — testing against `D_TOLERANCE` would report
# "never reached ice" for every column. 0.5 kg m-3 is well below the uncertainty of
# any densification law and far below one output digit.
const DENSITY_FIRN_TOLERANCE = 0.5     # [kg m-3]

# Convergence tolerance on the albedo/melt fixed point in the initial guess. This is
# an initial state a spinup then relaxes, so resolving the albedo beyond ~1e-4 buys
# nothing and each iteration costs a full sweep of the forcing series.
const ALBEDO_FIXED_POINT_TOLERANCE = 1e-4

# Numerical boundary tolerances.
# These are float-safe offsets on branch comparisons (e.g. `x < threshold - D_TOLERANCE`)
# that reproduce the original MATLAB's exact branch decisions. They are load-bearing for
# MATLAB fidelity (the reference regression test diverges if they are removed), not
# arbitrary fudge factors. Values are centralized here so there is a single source of truth.
const D_TOLERANCE = 1e-11              # density / depth comparisons
const T_TOLERANCE = 1e-10              # temperature-gradient branch boundaries
const W_TOLERANCE = 1e-13              # surface (top-cell) water presence
const WATER_TOLERANCE = 1e-13          # pore water presence
const GDN_TOLERANCE = 1e-10            # grain dendricity / sphericity [0,1] clamps
const E_TOLERANCE = 1e-3               # energy-conservation check [J] (verbose only)
const M_TOLERANCE = 1e-3               # mass-conservation check [kg m-2] (verbose only)
const Z_TOLERANCE = 1e-9               # column-depth conservation check [m]
const T_BOTTOM_TOLERANCE = 1e-3        # Dirichlet bottom-cell drift check [K] (verbose only)
const T_MELT_SWITCH_TOLERANCE = 1e-4   # emissivity melt-switch temperature offset [K]

# Smallest stable sub-timestep the thermal solver will accept without warning [s]. Below
# this, the Von Neumann stability limit has collapsed, which in practice means a near-zero
# `dz` cell rather than a genuinely stiff column.
#
# `ExplicitThermal` only. `ImplicitThermal` is unconditionally stable, never consults
# `_max_safe_dt`, and so never warns however thin a cell becomes.
const DT_MIN_WARN = 1e-4

# --- ImplicitThermal solver ---
#
# Newton convergence on the surface energy balance. The iterate is the surface temperature, so
# the tolerance is in kelvin, and 1e-10 matches `T_TOLERANCE`: converging harder than the
# branch boundaries elsewhere in the model can resolve buys nothing. Convergence is quadratic
# once the residual is small, so this typically costs one iteration more than 1e-6 would.
#
# The iteration count is a safety net, not an expected path. Divergence is not possible in the
# usual sense — the linearized surface slope `Λ ≤ 0` strengthens an already diagonally dominant
# M-matrix — but the two discontinuities held outside the loop (the emissivity melt switch and
# the LV/LS latent-heat switch at `turbulent_heat_flux.jl:107`) can make the residual chatter
# at the melting point rather than settle. Stopping at 20 and accepting the last iterate is
# correct there: the flux actually applied is the true nonlinear flux at that iterate, so the
# column still conserves energy exactly, only the surface balance is left slightly unsatisfied.
const THERMAL_IMPLICIT_T_TOLERANCE = 1e-10
const THERMAL_IMPLICIT_MAX_ITERATIONS = 20

# One-sided finite-difference step [K] for the turbulent-flux derivatives `dQ_shf/dT` and
# `dQ_lhf/dT`, which enter the Newton diagonal. `_turbulent_heat_flux` runs through the
# Beljaars-Holtslag stability branches, so differentiating it by hand is not worth the
# maintenance; one extra call per iteration is cheaper than the accuracy would be worth.
#
# Measured, not assumed: the two flux evaluations together are 0.70 s of a 15.2 s run, so this
# call is not the cost of the implicit path. An analytic replacement that froze the stability
# coefficients was tried and reverted — it overshoots the true derivative up to 10x and more than
# doubles the iteration count. See `GEMB._surface_energy_balance_slope`.
#
# 1e-4 K rather than the usual `sqrt(eps)` ≈ 1.5e-8: the derivative only sets the convergence
# *rate*, never the answer (see `_thermal_solve!(::ImplicitThermal, ...)`), so a step large
# enough to stay well clear of cancellation in the flux difference beats a step chosen to
# minimize truncation error. It also matches `T_MELT_SWITCH_TOLERANCE`, so the difference does
# not straddle the melting-point latent-heat switch any more often than the melt switch does.
const THERMAL_BC_DERIVATIVE_STEP = 1e-4

# Target sub-step [s] for `ImplicitThermal`, and the ceiling on the resulting count.
#
# Backward Euler is unconditionally stable, so this is an *accuracy* control, not a stability
# one — the distinction that makes the implicit path cheap. The count is bounded by how fast the
# surface forcing changes and is independent of the cell count, of `dz`, and of the stiffest
# cell in the column, unlike the explicit path's ~29 sub-steps.
#
# Calibrated on a year of 3-hourly synthetic forcing (264-cell column, no spinup), comparing
# whole-model output against a finely sub-stepped implicit run at 168.75 s. RMS temperature
# differences [K] and the annual melt difference:
#
#   target [s]   n_sub(3h)   RMS surface   RMS interior   Δmelt      runtime
#      5400          2          1.164         0.143       -2.46%      —
#      2700          4          0.526         0.140       -1.66%      —
#      1800          6          0.419         0.126       -1.19%     3.6 s
#       900         12          0.341         0.141       -1.03%     5.8 s
#       450         24          0.272         0.119       -0.08%     9.6 s
#
# Runtimes re-measured after the static condensation of the interior rows (2.67x); the accuracy
# columns are unchanged by it, since it alters only how the same solve is factorized.
#
# 900 s is the knee: melt is within ~1% and the surface RMS has come down 3.4x from the
# two-sub-step case, while halving again costs 79% more runtime for 0.07 K.
#
# The interior column is *not* converging in this table — it floors at ~0.13 K regardless of
# target. That floor is not solver error: `manage_layer_thickness` merges and splits on discrete
# thickness thresholds, so a single flipped merge changes the deep discretization by more than
# the time integration does. Whole-model differencing therefore cannot resolve implicit accuracy
# below ~0.13 K, which is why the surface RMS and the melt total are what this is tuned on.
const THERMAL_IMPLICIT_DT_TARGET = 900.0
const THERMAL_IMPLICIT_SUBSTEPS_MAX = 64

# Diagnostic threshold, not a branch tolerance: how far above irreducible a cell must be for
# `aquifer_diagnostics` to call it saturated. Deliberately loose, and expressed as a fraction
# of the cell's pore space rather than as a mass.
#
# The reason is physical, not arithmetic. `calculate_melt` leaves a retaining cell at exactly
# its irreducible water, but `calculate_density` then compacts the cell in the same timestep,
# shrinking the pore space that irreducible water was computed against. The cell is left
# marginally over-saturated — a genuine excess of order 1e-4 kg m-2 per timestep, which
# accumulates between melt events. A mass threshold cannot separate that from a thin aquifer,
# because the two differ in saturation, not in mass: compaction residue sits ~1e-7 of pore
# space above irreducible, while ponded water reaches 1e-1 to 1. Nothing in the physics reads
# this; it only decides what `aquifer_diagnostics` reports.
const AQUIFER_TOLERANCE = 1e-3         # standing-water detection [fraction of pore space]

# The "not requested" sentinel for the opt-in per-cell viscosity diagnostic. A shared, empty
# `Vector{Float64}` rather than `nothing` so that `gemb_core`'s `flux` NamedTuple has one
# concrete type regardless of `mp.output_viscosity`; see the note in `gemb_core.jl`. Shared
# rather than allocated per call because it is never written to — every write site tests
# `isempty` first — so there is nothing to alias.
const NO_VISCOSITY = Float64[]

"""
    energy_tolerance(E_reference)

Absolute tolerance [J] for a verbose-only energy-conservation check whose budget totals
`E_reference`.

`E_TOLERANCE` alone is a *relative* 1e-14 against a deep column's ~1e11 J of enthalpy, which
only ever passed because both sides of the check were bit-identical arithmetic. Working in
enthalpy adds real round-trip noise (~3e-6 J per deep-cell merge), so every check is
relative-with-floor: the floor keeps small columns strict, the relative term scales with the
budget being checked.
"""
@inline energy_tolerance(E_reference::Real) = max(E_TOLERANCE, 1e-12 * abs(E_reference))
