# Physical constants used throughout GEMB
# All values match the MATLAB implementation exactly

const CtoK = 273.15        # Celsius to Kelvin conversion [K]
const C_ICE = 2102.0       # Specific heat capacity of snow/ice [J kg-1 K-1]
const C_AIR = 1005.0       # Specific heat capacity of air [J kg-1 K-1]
const LF = 0.3345e6        # Latent heat of fusion [J kg-1]
const LV = 2.495e6         # Latent heat of vaporization [J kg-1]
const LS = 2.8295e6        # Latent heat of sublimation [J kg-1]
const SB = 5.67e-8         # Stefan-Boltzmann constant [W m-2 K-4]
const GRAVITY = 9.81       # Gravitational acceleration [m s-2]
const R_GAS = 8.314        # Universal gas constant [J mol-1 K-1]
const DENSITY_WATER = 1000.0  # Density of water [kg m-3]
const VON_KARMAN = 0.4     # Von Karman constant [-]

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
