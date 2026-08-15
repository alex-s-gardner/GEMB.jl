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

# Density separating the two densification stages [kg m-3]. Every scheme in
# `calculate_density` switches rate coefficients here: it is the transition from
# grain-settling-dominated to creep-dominated compaction, and all the published two-stage
# coefficient pairs (Herron & Langway, Arthern, Ligtenberg, Simonsen, GSFC) are fitted
# against this same split.
const DENSITY_STAGE_TRANSITION = 550.0

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
const T_MELT_SWITCH_TOLERANCE = 1e-4   # emissivity melt-switch temperature offset [K]

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
