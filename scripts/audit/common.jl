# Shared by the audit scripts in this directory: version detection and wrappers for the few
# internal functions whose signatures differ between `main` (0cda823) and `fix/physics-audit`,
# so that every script runs unchanged against either checkout.
using GEMB, Dates, Printf, Statistics

# `apply_vapor_exchange!` exists only on the branch.
const ON_BRANCH = isdefined(GEMB, :apply_vapor_exchange!)
println("GEMB at ", pkgdir(GEMB), " — ", ON_BRANCH ? "fix/physics-audit" : "main (0cda823)")

# calculate_albedo: main takes (dz, density, water, grain_radius, melt_surface, cfs, mp).
albedo(dz, density, grain_radius, cfs, mp) = ON_BRANCH ?
    GEMB.calculate_albedo(dz, density, grain_radius, cfs, mp) :
    GEMB.calculate_albedo(dz, density, zeros(length(dz)), grain_radius, 0.0, cfs, mp)

# calculate_melt: main takes the step's rain before `mp`. Returns
# (temperature, dz, density, water, ..., melt, melt_surface, runoff, refreeze, percolation_depth).
melt(T, dz, ρ, w, gr, gd, gs, age, rain, mp) = ON_BRANCH ?
    GEMB.calculate_melt(T, dz, ρ, w, gr, gd, gs, age, mp, true) :
    GEMB.calculate_melt(T, dz, ρ, w, gr, gd, gs, age, rain, mp, true)

# calculate_temperature with the column's pore water; main has no `water` keyword and treats
# only `water_surface` (cell 1) for roughness.
thermal(T, dz, ρ, water, gr, sw, cfs, mp) = ON_BRANCH ?
    GEMB.calculate_temperature(T, dz, ρ, water[1], gr, sw, cfs, mp, true; water) :
    GEMB.calculate_temperature(T, dz, ρ, water[1], gr, sw, cfs, mp, true)

# A forcing step built by keyword (both versions have this constructor).
step(; kw...) = GEMB.ClimateForcingStep(; kw...)
