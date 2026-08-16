# GEMB and IMAU-FDM

[IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM) (Ligtenberg et al.,
2011; Brils et al., 2022) is the Utrecht Institute for Marine and Atmospheric
Research firn densification model — ~3,900 lines of Fortran in `source/`,
configured per domain through TOML (`FGRN055` for Greenland, `ANT27` for
Antarctica). It is the third open firn model with a physics surface comparable to
GEMB's, after the [Community Firn Model comparison](cfm_comparison.md).

This page records a full read of IMAU-FDM's physics-bearing sources against
GEMB's `src/`, in two directions:

1. **What GEMB could take from it** — two findings, both adopted.
2. **Whether GEMB's existing physics agrees with an independent implementation** —
   which produced more than the first direction did, including one place where
   IMAU-FDM's own source comments flag its version as a bug GEMB never had.

IMAU-FDM is the **narrower** model of the two. It takes surface temperature as a
prescribed Dirichlet boundary rather than solving an energy balance, and has no
albedo, no shortwave penetration, no turbulent fluxes, and no grain-size state.
Most of the capability surface is therefore GEMB-has/FDM-lacks and is not
tabulated below; the interesting content is in the places where both models
implement the *same* physics and disagree about it.

Citations point at IMAU-FDM as `firn_physics.f90:247` and at GEMB as
`src/thermal_conductivity.jl`. Fortran line numbers are from the state of `main`
read in August 2026.

## Summary

| Physics | Verdict | Notes |
|---|---|---|
| Calonne (2019) air-conductivity ratio | **gap, adopted** | Finding 1 — new `:Calonne2019Air` method |
| Fausto (2018) fresh-snow density | **gap, adopted** | Finding 2 — new `:FaustoFit` method |
| Coleou irreducible water content | **GEMB correct** | GEMB never had the bug IMAU-FDM's own comment documents |
| Layer merging | **GEMB correct** | GEMB merges on enthalpy; IMAU-FDM mass-weights temperature |
| Refreeze energy accounting | **GEMB correct** | IMAU-FDM mixes `cp(T)` and `cp(T_melt)` in one budget |
| Densification reference temperature | **GEMB correct** | IMAU-FDM uses the *bottom* layer; the paper, CFM and GEMB use mean annual surface |
| Spinup convergence criterion | **GEMB richer** | GEMB tests drift slope as well as per-cycle delta |
| Heat capacity, physical constants | **agreement** | same `152.5 + 7.122·T`, `g`, `Ec`, `Eg`, `ρ_i`, `T_melt` |
| Ligtenberg `M0`/`M1` calibration | **GEMB ahead** | 9 calibration sets vs IMAU-FDM's 2 |
| Implicit tridiagonal thermal solver | **gap, deferred** | interior rows transfer; the surface row does not — see below |

## What GEMB adopted

Both findings land as **new `Symbol` option values**, leaving the defaults of the day
bit-identical. Neither appears in `README.md`'s **Physics notes** section, because neither
changed a default when it was added. (`:Calonne2019` later *became* the default in v2.0.0,
for the separate reasons recorded there.)

### Finding 1 — `:Calonne2019` omits the air-conductivity ratio of eq. 5

Calonne et al. (2019) eq. 5 blends a snow-regime and a firn-regime conductivity
fit with a logistic weight in density, and rescales each by the constituent
conductivities at the working temperature relative to their values at the
temperature the fits were built at (270.15 K). The two regimes carry **different**
scalings because they describe different microstructures: the firn fit describes a
connected ice skeleton and scales with ice alone, while the snow fit describes
grains suspended in air and scales with *both* constituents.

IMAU-FDM implements both factors (`firn_physics.f90:247-268`, `Thermal_Cond`):

```fortran
kair = (2.334E-3*Temp**(3./2.))/(164.54 + Temp)          ! Reid (1966)
kair_ref = (2.334E-3*270.15**(3./2.))/(164.54 + 270.15)
ki = (1.-theta)*kice/kice_ref*kair/kair_ref*kcal + theta*kice/kice_ref*kf
```

GEMB's `_thermal_conductivity_calonne2019` applied only the ice ratio, to both
branches. Its docstring was explicit about this, because the implementation had
been cross-checked against the CFM — which hardcodes `K_air = kref_a` in
`diffusion.py` with the standing comment *"at some point find equation for
T-dependence of air"*. GEMB inherited the CFM's simplification; **IMAU-FDM is the
only one of the three models that carries the term.**

Because air conducts less well as it cools, omitting the factor makes GEMB
**over-conduct cold snow**. Measured ratio of `:Calonne2019` to the fuller form:

| T [K] | `kair/kair_ref` | ρ=150 | ρ=300 | ρ=450 | ρ=550 | ρ=700 |
|---|---|---|---|---|---|---|
| 220.0 | 0.8307 | +20.4% | +20.4% | +9.9% | +0.28% | 0% |
| 233.15 | 0.8764 | +14.1% | +14.1% | +7.0% | +0.20% | 0% |
| 253.15 | 0.9440 | +5.9% | +5.9% | +3.1% | +0.09% | 0% |
| 270.15 | 1.0000 | 0% | 0% | 0% | 0% | 0% |

The bias vanishes at the reference temperature and above ρ ≈ 550, where the
logistic weight hands over to the air-free firn branch. It is largest exactly in a
cold, low-density polar winter surface layer — that is, exactly where the seasonal
cold wave propagates — so it damps the winter surface signal.

`mp.thermal_conductivity_method = :Calonne2019Air` selects the fuller form.
`:Calonne2019` remains the default and is unchanged. Since the air term enters only
as a *ratio*, Reid's absolute calibration cancels and only its temperature shape
matters; continuity into ice is preserved because the factor is weighted by
`(1-θ)`, which is ~1e-8 at ρ = 917.

One convention difference is worth recording: IMAU-FDM normalizes by the
*computed* Yen (1981) value at 270.15 K, where GEMB and the CFM both use the
published literal `2.107`. That is a uniform 1.15e-4 relative offset, it predates
this addition, and it applies equally to `:Calonne2019`. GEMB keeps the literal.

### Finding 2 — `:Fausto` is the published fit frozen at one temperature

GEMB's `fresh_snow_density` returned a bare `315.0` for `:Fausto`, with no
temperature dependence — unlike every other temperature-dependent option in the
same function (`:Kaspers`, `:KuipersMunneke`), which all carry theirs. IMAU-FDM
implements Fausto et al. (2018) as the actual regression
(`initialise_model.f90:96`, `112`, `116`):

```fortran
Rho0FM(step) = 362.1 + 2.78*(TempFM(step) - const%Tmelt)   ! Fausto et al. 2018
```

`362.1 + 2.78·(T − 273.15) = 315.0` at **T ≈ 256.2 K**. GEMB's constant is
demonstrably this same fit evaluated at one plausible Greenland annual-mean
temperature — not a different parameterization.

`mp.new_snow_method = :FaustoFit` selects the regression;  `:Fausto` still returns
315.0. Two implementation notes:

- `fresh_snow_density` gained an **optional fifth argument**, the *instantaneous*
  air temperature, which only `:FaustoFit` reads. The function is shared by
  `calculate_accumulation` (which has `cfs.temperature_air`) and
  `steady_state_profile` (which has only a climatological mean), so the argument
  defaults to the mean: the transient site passes the instantaneous temperature
  while the steady-state march still gets a well-defined climatological ρ₀, and the
  four older methods are bit-identical either way.
- The regression is **unbounded below** — it reaches zero at T ≈ 143 K — so the
  transient call site now clamps to `[1, ρ_i]`, mirroring the clamp
  `steady_state_profile` already applied. IMAU-FDM caps its *non*-Greenland fit at
  470 but leaves the Fausto branch uncapped, so there is no upstream cap to copy.
- `:FaustoFit` also takes the Crocus (Vionnet et al., 2012) wind-dependent
  fresh-snow grain properties that `:Fausto` takes. That coupling is orthogonal to
  the density fit, and omitting it would have silently changed fresh-snow
  dendricity.

## Consistency findings — where GEMB is correct or better

These need no code changes. They are independent corroboration of GEMB's physics,
recorded here so they are not re-litigated.

### Coleou irreducible water — GEMB never had IMAU-FDM's bug

`water_physics.f90:190-198` offers two variants of the Coleou & Lesaffre (1998)
irreducible water content. The one named `Coleou1998_corr` divides by `(1 − Wm)`,
with the source comment *"fixes known bug in FDM v1p2"*; the uncorrected
`Coleou1998_1p2` "underestimates LWC". GEMB's `irreducible_saturation`
(`src/calculate_melt.jl`) has always had the `wmi/(1-wmi)` correction. Two
implementations arrived independently at the same wet-mass-versus-dry-mass
subtlety, and GEMB was on the right side of it from the start.

### Layer merging — GEMB conserves energy, IMAU-FDM does not

`grid_routines.f90:54` merges two cells by mass-weighting *temperature*:

```fortran
T(k-1) = (T(k-1)*M(k-1) + T(k)*M(k))/(M(k-1)+M(k))
```

That is not energy-conserving given IMAU-FDM's own temperature-dependent
`c_p = 152.5 + 7.122·T`: the merged cell's enthalpy differs from the sum of its
parts at second order in ΔT. GEMB's `manage_layer_thickness` merges through
`specific_enthalpy`, which is exact for any `c_p(T)`.

The same defect appears in IMAU-FDM's refreeze routines
(`water_physics.f90:34-45`, `151-158`), which compute the available energy with
`cp` evaluated at the layer temperature but apply the resulting temperature
increment with `cp0` evaluated at the melting point. GEMB carries enthalpy — never
`T·c_p` — as the budget currency throughout, and its `verbose=true` per-timestep
conservation checks would not survive the mixed-`c_p` form.

### Densification reference temperature — IMAU-FDM is the outlier

`Densific` evaluates the grain-growth Arrhenius term at `T(1)`, the **bottom**
layer of the column (`firn_physics.f90:323`, `329`). Arthern et al. (2010)
specifies mean annual *surface* temperature; the CFM uses `self.T_mean`
(surface, `physics.py:334`); GEMB uses `tam` (`src/calculate_density.jl`). Two of
three implementations agree with the paper.

IMAU-FDM's choice is defensible as a steady-state identity — at deep-column
equilibrium the bottom temperature approaches the mean annual surface temperature —
but it is not the published form, and it couples densification to the basal
boundary condition during transients.

### Spinup convergence — GEMB's criterion is richer

IMAU-FDM iterates on the squared change in surface elevation *and* in firn air
content between cycles (`time_loop.f90:57`, `129-130`), with bounds 0.0001
(FGRN055) / 0.004 (ANT27) and 3–70 or 3–200 cycles. GEMB's `gemb_spinup`
(`src/spinup.jl`) tests **both** a per-cycle delta *and* a least-squares drift
slope over a window — precisely because, as its own docstring says, "a column
creeping steadily at just under the tolerance passes" a delta-only test. That is
the exact weakness of IMAU-FDM's criterion.

This also resolves an open question left by the CFM pass: the benchmark's
`converged = false` after 75 cycles is because `bench/opt_bench.jl` passes **no**
convergence criteria at all, so `checking = false` and the flag is never set. It
was never a convergence failure.

### Heat capacity and constants — agreement

GEMB's `:CuffeyPaterson` heat capacity is the same `152.5 + 7.122·T` IMAU-FDM
hardcodes throughout. `g = 9.81`, `Ec = 60000`, `Eg = 42400`, `ρ_i = 917`,
`T_melt = 273.15` all match (`settings/*/constants.toml`). IMAU-FDM's
`days_per_year = 365.25` against GEMB's deliberate
`DAYS_PER_YEAR_DENSIFICATION = 365.0` is the ~0.07% standing offset already
recorded in `README.md` deviation 4 — no new action.

## Deliberately not adopted

### Implicit tridiagonal thermal solver — the interior transfers, the surface does not

The CFM pass flagged an implicit thermal solver as GEMB's highest-value numerics
opportunity, and IMAU-FDM's `Solve_Temp_Imp` (`firn_physics.f90:169-241`) is
attractively compact: ~70 lines, θ-weighted, Thomas algorithm, following Versteeg
& Malalasekera. GEMB currently uses an explicit, Von Neumann sub-stepped enthalpy
solve.

**It is not a drop-in.** Line 219 builds the surface boundary as
`Su = 2·kip·Ts/DZ` — a *prescribed* Dirichlet surface temperature. GEMB solves a
nonlinear surface energy balance (longwave, turbulent fluxes) in
`calculate_temperature.jl`, plus shortwave penetrating to depth, so the surface row
would need linearization and outer iteration, and the shortwave source term has no
analogue in IMAU-FDM's tridiagonal at all. The bottom boundary does match: both
models hold a fixed bottom cell.

The useful half of the finding is that the **interior** rows transfer directly.
That lowers the cost of the eventual implementation below what the CFM assessment
assumed, but it does not make it small.

### `numSnow` window-averaged temperature

IMAU-FDM optionally drives fresh-snow density from a temperature averaged over a
recent-snowfall window (`initialise_model.f90:104-134`). That is a smoothing
choice, not new physics, and the `fresh_snow_density` signature change in Finding 2
makes either input expressible if it is ever wanted.

### Domain-specific `MO` recalibration

`Densific` (`firn_physics.f90:288-315`) recalibrates its densification
multiplier per domain, e.g. `MO_low = 0.7522 − 0.0178·log(acav)` for Greenland.
This is the direct analogue of GEMB's Ligtenberg `M0`/`M1`
(`src/densification_lookup.jl`), which already carries **nine** calibration sets
across Antarctica and Greenland against IMAU-FDM's two. GEMB is ahead here, and
the fits are not interchangeable — different forcing products, and a different
Arrhenius reference temperature per the point above — so cross-importing
coefficients would be unsound.

### Bottom-of-column layer add/delete

`Add_Layers` / `Delete_Layers` grow and shrink IMAU-FDM's column in 100- or
200-layer blocks. GEMB's fixed-count, fixed-depth column with its two `grid_ops.jl`
controllers is a deliberate and documented alternative, not a gap.

### Ice-shelf buoyancy

`vbouy` (`firn_physics.f90:106-113`) is a floating-ice elevation diagnostic,
outside a column mass-balance model's remit. GEMB's `apply_horizontal_strain!`
covers the ice-dynamic coupling it does model.

## References

- Brils, M., Kuipers Munneke, P., van de Berg, W. J., and van den Broeke, M.
  (2022). Improved representation of the contemporary Greenland ice sheet firn
  layer by IMAU-FDM v1.2G. *Geoscientific Model Development*, 15, 7121–7138.
- Ligtenberg, S. R. M., Helsen, M. M., and van den Broeke, M. R. (2011). An
  improved semi-empirical model for the densification of Antarctic firn.
  *The Cryosphere*, 5, 809–819.
- Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L., Flin, F.,
  and Geindreau, C. (2019). Thermal conductivity of snow, firn, and porous ice
  from 3-D image-based computations. *Geophysical Research Letters*, 46,
  13079–13089.
- Reid, R. C., Prausnitz, J. M., and Sherwood, T. K. (1966). *The Properties of
  Gases and Liquids*. McGraw-Hill.
- Yen, Y.-C. (1981). *Review of thermal properties of snow, ice and sea ice*.
  CRREL Report 81-10.
- Fausto, R. S., Box, J. E., Vandecrux, B., van As, D., Steffen, K., MacFerrin,
  M. J., Machguth, H., and Colgan, W. (2018). A snow density dataset for improving
  surface boundary conditions in Greenland ice sheet firn modeling. *Frontiers in
  Earth Science*, 6, 51.
- Coleou, C., and Lesaffre, B. (1998). Irreducible water saturation in snow:
  experimental results in a cold laboratory. *Annals of Glaciology*, 26, 64–68.
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.
  (2010). In situ measurements of Antarctic snow compaction compared with
  predictions of models. *Journal of Geophysical Research*, 115, F03011.
- Vionnet, V., Brun, E., Morin, S., Boone, A., Faroux, S., Le Moigne, P., Martin,
  E., and Willemet, J.-M. (2012). The detailed snowpack scheme Crocus and its
  implementation in SURFEX v7.2. *Geoscientific Model Development*, 5, 773–791.
- Versteeg, H. K., and Malalasekera, W. (2007). *An Introduction to Computational
  Fluid Dynamics: The Finite Volume Method*, 2nd ed. Pearson.
