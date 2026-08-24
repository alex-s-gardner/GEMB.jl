# GEMB and IMAU-FDM

[IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM) (Ligtenberg et al., 2011;
Brils et al., 2022) is the Utrecht Institute for Marine and Atmospheric Research firn
densification model, configured per ice sheet domain. It takes surface temperature and
drifting snow as prescribed boundary conditions from a regional climate model, where GEMB
solves a surface energy balance and carries grain-size state.

This page summarizes how the two treat the physics they share, and records what GEMB adopted
while reading IMAU-FDM. Most of the capability surface follows from the difference in scope
and is not tabulated; the interesting content is where both models implement the same law.

## Scope

| | GEMB | IMAU-FDM |
|---|---|---|
| Surface boundary | Nonlinear surface energy balance solved for skin temperature | Prescribed surface temperature (Dirichlet) |
| Albedo, shortwave penetration, turbulent fluxes | Solved | Not represented |
| Grain size | Prognostic | Not carried |
| Drifting snow | Prescribed forcing layer, or the Crocus scheme | Prescribed from RACMO |
| Vertical grid | Fixed cell count and depth | Grows and shrinks in 100–200 layer blocks |
| Thermal solve | Explicit finite volume (default) or backward Euler | θ-weighted implicit, Thomas algorithm |
| Ligtenberg densification calibration | 9 sets, Antarctica and Greenland | 2 sets, one per domain |

The two agree on the heat capacity `c_p = 152.5 + 7.122·T`, on gravity, on the activation
energies for creep and grain growth, on ice density, and on the melting point. IMAU-FDM uses
365.25 days per year against GEMB's 365.0 for densification, a ~0.07% standing offset
recorded in the [Physics notes](@ref "Physics notes").

## What GEMB adopted

Three findings, each selected by an option or an optional forcing layer.

**The air-conductivity ratio in Calonne et al. (2019) eq. 5.** The equation blends a snow and
a firn conductivity fit and rescales each by its constituent conductivities relative to their
values at 270.15 K. The two branches carry different scalings: the firn branch describes a
connected ice skeleton and scales with ice alone, while the snow branch describes grains in
air and scales with both. IMAU-FDM carries both factors, using Reid (1966) for air. GEMB's
`:Calonne2019` applies only the ice ratio, following CFM. Because air conducts less well as
it cools, omitting the air factor raises conductivity in cold low-density snow — by ~20% at
220 K and 150–300 kg m⁻³, falling to zero at the reference temperature and above ~550 kg m⁻³
where the blend hands over to the air-free branch.

`thermal_conductivity_method = :Calonne2019Air` selects the fuller form. Since the air term
enters only as a ratio, Reid's absolute calibration cancels and only its temperature shape
matters.

**Fausto et al. (2018) as a temperature-dependent regression.** GEMB's `:Fausto` fresh-snow
density returns a constant 315 kg m⁻³, which is the published regression
`362.1 + 2.78·(T − 273.15)` evaluated at 256.2 K. IMAU-FDM implements the regression.
`new_snow_method = :FaustoFit` does the same in GEMB, reading instantaneous air temperature
where it is available and the climatological mean in the steady-state march. The regression
is unbounded below, so the transient call site clamps the result.

**Drifting snow as a prescribed surface flux.** IMAU-FDM reads a drift field alongside
precipitation and applies it with a sign-dependent density: erosion removes mass at the
surface layer's own density, and deposition adds it at the fresh-snow density. That asymmetry
matters — taking the cell density in both directions would let an erode-and-redeposit cycle
change surface density with no net mass exchange.

GEMB now has the same path: an optional `snow_drift` forcing layer (kg m⁻² yr⁻¹, positive
for erosion), with `drift_rate` as a constant-rate shorthand. Deposited mass enters at age
zero. Three differences from IMAU-FDM: the layer is a rate rather than a per-timestep mass,
so it is independent of the forcing interval; it is read independently of
`blowing_snow_method`, so a climate-model drift field can drive GEMB with no internal scheme;
and it is zero by default.

For an internally computed blowing-snow scheme the reference is Crocus rather than IMAU-FDM,
whose drift physics is upstream in RACMO. That scheme is available as
`blowing_snow_method = :Crocus`.

## Where the two implementations differ

Neither model was changed on GEMB's side for these; they are recorded so the comparison does
not have to be redone.

**Irreducible water content.** IMAU-FDM offers two forms of Coleou and Lesaffre (1998), one
of which divides by `(1 − W_m)`. GEMB's `irreducible_saturation` carries that factor.

**Layer merging.** IMAU-FDM merges two cells by mass-weighting temperature. With a
temperature-dependent `c_p` that is not exactly energy-conserving; the merged enthalpy differs
from the sum of the parts at second order in ΔT. GEMB merges through specific enthalpy, which
is exact for any `c_p(T)`. The same distinction appears in IMAU-FDM's refreeze routines, which
evaluate `c_p` at the layer temperature but apply the resulting increment with `c_p` at the
melting point. GEMB carries enthalpy as the budget currency throughout.

**Densification reference temperature.** IMAU-FDM evaluates the grain-growth Arrhenius term at
the bottom layer's temperature. Arthern et al. (2010) specifies mean annual surface
temperature, which is what CFM and GEMB use. The bottom-layer choice is an equilibrium
identity — a deep column's base approaches the mean annual surface temperature — but it
couples densification to the basal boundary condition during transients.

**Spinup convergence.** IMAU-FDM iterates on the change in surface elevation and firn air
content between cycles. GEMB's `gemb_spinup` tests a per-cycle delta and, optionally, a
least-squares drift slope over a trailing window, since a column creeping steadily at just
under a delta tolerance passes a delta-only test.

**Domain-specific densification recalibration.** IMAU-FDM recalibrates its densification
multiplier per domain. GEMB's Ligtenberg `M0`/`M1` lookup does the same thing and carries nine
calibration sets. The coefficients are not interchangeable between the two: different forcing
products, and a different Arrhenius reference temperature per the point above.

**Column growth.** IMAU-FDM adds and deletes layers at the base in blocks. GEMB's column holds
a fixed cell count and total depth, maintained by two merge/split controllers.

**Ice-shelf buoyancy.** IMAU-FDM diagnoses floating-ice elevation. GEMB represents the
ice-dynamic coupling through `apply_horizontal_strain!` only.

**Window-averaged temperature for fresh-snow density.** IMAU-FDM can drive fresh-snow density
from a temperature averaged over a recent-snowfall window. GEMB does not; the
`fresh_snow_density` signature admits either input.

## The thermal solve

IMAU-FDM's implicit solver is compact — θ-weighted, Thomas algorithm, following Versteeg and
Malalasekera — and its interior rows transfer directly to GEMB. Its surface row builds a
prescribed Dirichlet surface temperature, which has no analogue in GEMB, where the surface row
is a nonlinear energy balance and shortwave radiation is a distributed source. The bottom
boundary does match: both models hold a fixed bottom cell.

GEMB's `ImplicitThermal` linearizes the surface flux into the matrix diagonal and iterates
with Newton. A lagged surface flux diverges for a centimetre-scale surface cell at any
timestep of interest. Because the residual is built from the true nonlinear flux at each
iterate, the linearized terms cancel at convergence, so the slope affects the convergence rate
and not the answer. The nonlinearity being confined to row 1 also permits static condensation:
the interior is eliminated once per sub-step and Newton iterates on a single scalar equation.

`ExplicitThermal` remains the default. Fourteau et al. (2024) derive the same
Newton-with-Schur-complement construction independently, and quantify what a prescribed
Dirichlet surface temperature costs in their framework: a spurious flux of −14.5 W m⁻² for a
melting glacier surface, changing ablation by 40%. See [Thermal solvers](@ref "Thermal solvers")
and [Surface energy balance numerics](architecture.md#Surface-energy-balance-numerics).

## References

- Brils, M., Kuipers Munneke, P., van de Berg, W. J., and van den Broeke, M. (2022). Improved
  representation of the contemporary Greenland ice sheet firn layer by IMAU-FDM v1.2G.
  *Geoscientific Model Development*, 15, 7121–7138.
- Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R. (2011). An improved
  semi-empirical model for the densification of Antarctic firn. *The Cryosphere*, 5, 809–819.
- Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F., and Geindreau,
  C. (2019). Thermal conductivity of snow, firn, and porous ice from 3-D image-based
  computations. *Geophysical Research Letters*, 46, 13079–13089.
- Reid, R. C., Prausnitz, J. M., and Sherwood, T. K. (1966). *The Properties of Gases and
  Liquids*. McGraw-Hill.
- Yen, Y.-C. (1981). *Review of thermal properties of snow, ice and sea ice*. CRREL Report
  81-10.
- Fausto, R. S., Box, J. E., Vandecrux, B., van As, D., Steffen, K., MacFerrin, M. J.,
  Machguth, H., and Colgan, W. (2018). A snow density dataset for improving surface boundary
  conditions in Greenland ice sheet firn modeling. *Frontiers in Earth Science*, 6, 51.
- Coleou, C., and Lesaffre, B. (1998). Irreducible water saturation in snow: experimental
  results in a cold laboratory. *Annals of Glaciology*, 26, 64–68.
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R. (2010).
  In situ measurements of Antarctic snow compaction compared with predictions of models.
  *Journal of Geophysical Research*, 115, F03011.
- Vionnet, V., Brun, E., Morin, S., Boone, A., Faroux, S., Le Moigne, P., Martin, E., and
  Willemet, J.-M. (2012). The detailed snowpack scheme Crocus and its implementation in
  SURFEX v7.2. *Geoscientific Model Development*, 5, 773–791.
- Lenaerts, J. T. M., van den Broeke, M. R., Déry, S. J., van Meijgaard, E., van de Berg,
  W. J., Palm, S. P., and Sanz Rodrigo, J. (2012). Modeling drifting snow in Antarctica with a
  regional climate model: 1. Methods and model evaluation. *Journal of Geophysical Research*,
  117, D05108.
- Versteeg, H. K., and Malalasekera, W. (2007). *An Introduction to Computational Fluid
  Dynamics: The Finite Volume Method*, 2nd ed. Pearson.
- Fourteau, K., Brondex, J., Brun, F., and Dumont, M. (2024). A novel numerical implementation
  for the surface energy budget of melting snowpacks and glaciers. *Geoscientific Model
  Development*, 17, 1903–1929.
