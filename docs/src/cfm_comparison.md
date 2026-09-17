# GEMB and the Community Firn Model

The [Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) (CFM;
Stevens et al., 2020) is an open-source firn model built for ice-core interpretation and
firn densification research. GEMB is built for surface mass and energy balance over ice
sheets. The two overlap in the firn column and differ in what surrounds it.

This page summarizes how the two treat the physics they share, and records what GEMB
adopted while reading CFM. It is a description of two designs, not a ranking of them.

## Scope

| | GEMB | CFM |
|---|---|---|
| Surface energy balance | Solved: albedo, shortwave penetration with depth, turbulent fluxes, energy-coupled sublimation | Albedo and turbulent fluxes are prescribed inputs |
| Densification schemes | 9 selectable | ~24 selectable |
| Grain size | Prognostic: temperature-gradient (Marbouty) and wet (Brun) metamorphism, Arthern (2010) normal growth | Prognostic: Arthern (2010) normal growth |
| Percolation | Tipping bucket, with Darcy or Zuo–Oerlemans lateral drainage | Tipping bucket, Richards equation, dual-domain preferential flow |
| Vertical grid | Fixed cell count and fixed total depth, maintained by merge/split | Multi-resolution subgrids |
| Thermal solve | Explicit finite volume (default) or backward Euler | Implicit, shared with gas and isotope transport |
| Gas transport, isotope diffusion | Not represented | Firn air and isotope modules |

Both models use a tipping-bucket percolation scheme as their default, share the van
Genuchten / Yamaguchi / Calonne parameter set for Darcy drainage, and arrived independently
at the same impermeable-layer runoff thresholds (830 kg m⁻³ over 0.1 m).

## Physics GEMB shares with CFM

Where the two implement the same law, GEMB's version is cross-checked against CFM's. The
Arthern (2010) densification scheme (`:ArthernB`) agrees to 1 part in 10⁸; `:GSFC2020`,
`:Simonsen2013`, `:Barnola1991` and `:Crocus` agree to a 0.1% offset arising from a
different seconds-per-year convention. Full details are in the
[Physics notes](@ref "Physics notes").

Both models solve the subsurface heat equation as a finite volume problem with
harmonic-mean face conductivity and apply the surface energy balance as a flux. Fourteau et
al. (2024) place both in the same class in their taxonomy of surface-energy-balance
couplings.

## What GEMB adopted

Six changes came out of reading CFM against GEMB. Each is selected by a `ModelParameters`
option, so earlier behavior remains reproducible, and each is recorded in the
[Physics notes](@ref "Physics notes").

**Densification driven by snowfall rather than total precipitation.** Rain does not bury the
column, so accumulation-driven schemes now compact against a snow-only accumulation mean.
CFM's `GSFC2020` excludes instantaneous accumulation forcing for a related reason.
`densification_accumulation`.

**Arrhenius terms evaluated at an effective temperature.** `exp(E/RT)` is convex in `1/T`,
so evaluating it at the arithmetic mean temperature differs from the mean of the evaluated
term — by tens of percent at a site with a ±15 K seasonal cycle. CFM computes an effective
temperature for this; GEMB now offers one. `mean_temperature_method`.

**Arthern (2010) grain growth in dry firn.** GEMB's seasonal-snow metamorphism laws are
inactive above 400 kg m⁻³ and where the temperature gradient vanishes, so grain radius did
not evolve in deep firn. CFM integrates `dr²/dt = k_gr·exp(−E_g/RT)`, which is
density-independent. GEMB gained the same law, either alone or handed over from Marbouty at
400 kg m⁻³. `grain_growth_method`.

**The initialized age profile.** `steady_state_profile` integrates age along its march;
that age is now carried into the initial column instead of being reset to zero, which makes
a pore-close-off age diagnostic available without a multi-century spinup. `initialize_age`.

**Rain sensible heat at the water heat capacity.** Rain above the melting point now carries
its sensible heat at 4219.9 rather than 2102 J kg⁻¹ K⁻¹. About 1% of the latent heat term
for rain at 275 K, growing with air temperature. `rain_heat_capacity`.

**Basal heat flux and effective viscosity as outputs.** Both were computed internally and
discarded. `heat_flux_basal` is the conductive flux across the deepest interior face — a
diagnosed flux against a Dirichlet reservoir, not a prescribed geothermal flux. `viscosity`
is optional (`output_viscosity`) and populated only under the Crocus settling law, the one
GEMB scheme that forms a viscosity; other schemes write `NaN`.

## Differences retained

**Liquid water in the cell heat capacity.** CFM's enthalpy solver includes pore water in the
volumetric heat capacity, which its solver can do because phase change happens inside the
solve. GEMB is split-operator: phase change happens in `calculate_melt`, and pore water is
carried at the melting point. Giving that water a sensible heat capacity in GEMB would let a
wet cell store energy the melt step cannot see. Adopting CFM's treatment means adopting its
enthalpy-with-phase-change solver.

**Liquid water in thermal conductivity.** CFM weights conductivity by liquid fraction; GEMB's
four methods are functions of density and temperature only. This part is separable from the
phase-change question and is a real gap. It needs the full water vector inside
`calculate_temperature`, a hot-loop signature change deferred to its own benchmark pass.

**Richards-equation and preferential-flow percolation.** GEMB is bucket-only, citing RetMIP
(Vandecrux et al., 2020, Sect. 5.2), which found that the added complexity did not
systematically improve agreement with observations at their sites. See
[Meltwater percolation scheme](architecture.md#Meltwater-percolation-scheme).

**Percolation ordering within a cell.** A bucket scheme may compute the downward flux either
after refreezing (Crocus's order, and GEMB's) or before it. Neither ordering is observationally
constrained. On GEMB's synthetic forcing the two are indistinguishable: the impermeability
re-test after refreeze diverts water 14 times per run, and in every case the cell below was
already a barrier, so run totals are identical. No option was added. Worth revisiting with
real forcing at an ice-slab site, where lenses form within permeable firn.

**Densification schemes not ported.** `Morris2014`, `Goujon2003`, `Breant2017`,
`KuipersMunneke2015`, `Brils2022`, `Max2018` and the `Yamazaki1993` fresh-snow stage. Out of
scope for this pass; the cross-validation path for adding one exists and works.

**Grid and spinup structure.** CFM's multi-resolution subgrids and its analytic
Herron–Langway spin-up both have GEMB counterparts that work differently: a fixed-count,
fixed-depth column with two merge/split controllers, and a steady-state march that applies to
every densification scheme rather than one.

**Firn air transport and isotope diffusion.** Ice-core capabilities outside a surface
mass-balance model's scope. Both ride CFM's shared implicit solver.

**Strain softening.** GEMB's `apply_horizontal_strain!` thins layers at constant density;
CFM's `horizontal_divergence` rescales mass, and additionally softens the densification rate.
The mass-versus-thickness difference is unresolved.

**Implicit thermal solver.** GEMB has one (`ImplicitThermal`), added while this comparison was
open. The interior rows follow the standard construction; GEMB's nonlinear surface energy
balance has no analogue in CFM's tridiagonal and is handled by Newton iteration. It is 2.4×
slower than the explicit default on a well-conditioned column, and its value is that its cost
does not depend on the thinnest cell. See [Thermal solvers](@ref "Thermal solvers").

## Convention differences

- Gravity: CFM 9.8, GEMB 9.81 m s⁻².
- Seconds per year: CFM 31 557 600; GEMB uses 365.0 days for densification, a ~0.1% standing
  offset recorded in the [Physics notes](@ref "Physics notes").
- Fresh-snow grain size: CFM uses Linow (2012); GEMB uses a fixed value.

## References

- Stevens, C. M., Verjans, V., Lundin, J. M. D., Kahle, E. C., Horlings, A. N., Horlings,
  B. I., and Waddington, E. D. (2020). The Community Firn Model (CFM), v1.0.
  *Geoscientific Model Development*, 13, 4355–4377.
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010).
  In situ measurements of Antarctic snow compaction compared with predictions of models.
  *Journal of Geophysical Research*, 115, F03011.
- Marbouty, D. (1980). An experimental study of temperature-gradient metamorphism.
  *Journal of Glaciology*, 26(94), 303–312.
- Vandecrux, B., et al. (2020). The firn meltwater Retention Model Intercomparison Project
  (RetMIP). *The Cryosphere*, 14, 3785–3810.
- Oraschewski, F. M., and Grinsted, A. (2022). Modeling enhanced firn densification due to
  strain softening. *The Cryosphere*, 16, 2683–2700.
- Fourteau, K., Brondex, J., Brun, F., and Dumont, M. (2024). A novel numerical
  implementation for the surface energy budget of melting snowpacks and glaciers.
  *Geoscientific Model Development*, 17, 1903–1929.
