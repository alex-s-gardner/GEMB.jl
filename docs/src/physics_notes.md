# Physics notes

```@meta
CurrentModule = GEMB
```

Where GEMB's implementation departs from, or interprets, the published law it is cited to,
and the reasoning behind choices a user of the model would otherwise have to read the source
to find. Organized by subject. Every scheme's full citation is in the docstring of the
function that implements it.

## Densification

**`:Arthern`** (Arthern et al., 2010) is the default: `dρ/dt = c(ρᵢ − ρ)` in two stages
either side of 550 kg m⁻³, driven by mean accumulation and an Arrhenius factor in mean
annual temperature.

**`:ArthernB`** is the stress-and-grain-size form of the same paper (eq. B1),
`dρ/dt ∝ σ/r²`. Overburden is the integral `Σ(ρⱼdzⱼ + waterⱼ)g` over the cells above, not
the nearest cell's density applied to the whole overlying depth. Cross-validated against the
Community Firn Model's `Arthern2010T` to 1 part in 10⁸. Because the rate goes as `1/r²`,
this is the scheme most sensitive to grain radius, and therefore to `grain_growth_method`.
The steady-state initial guess falls back to `:Arthern`, which spinup then relaxes.

**`:GSFC2020`** is Medley et al. (2022) GSFC-FDM v1.2.1 eq. 18, the recalibrated successor
to `:Arthern`: the same relaxation form with mean accumulation raised to a fitted exponent
(α₀ = 0.91, α₁ = 0.644) and a per-stage activation energy (59 500 and 56 870 J mol⁻¹ against
Arthern's single 60 000). Because α ≠ 1 the accumulation units are load-bearing rather than
absorbed into a prefactor: accumulation is kg m⁻² yr⁻¹, the convention the exponents were
fitted against. Cross-checked against the Community Firn Model to a 0.102% gravity offset.

**`:Simonsen2013`** is Simonsen et al. (2013), `:Arthern`'s form and activation energies
retuned for Greenland: a constant factor 0.8 below 550 kg m⁻³, and 1.25·γ above with
`γ = 61.7/√A · exp(−3800/RT̄)`, γ scaling the second stage only. The form follows Lundin et
al. (2017) FirnMICE eqs. A36–A37, which publishes no numeric tuning scalars, so 0.8 and 1.25
are the Community Firn Model's values, against which the implementation is cross-checked to
the same gravity offset. As for `:GSFC2020` the accumulation units are load-bearing, since
`γ ∝ A^(−1/2)`.

**`:Crocus` and `:CrocusPure`** are the viscous settling of Vionnet et al. (2012) eqs. 5–9,
`dD/D = −σ/η·dt` with `η = f1·f2·η₀(ρ/c_η)·exp(a_η(T_fus − T) + b_η·ρ)`. This is the only
scheme in which liquid water affects densification: `f1 = (1 + 60·W_liq/(ρ_w·D))⁻¹` reduces
viscosity by up to ~61× in wet cells, and
`f2 = min(4, exp(min(0.4 mm, gs − 0.2 mm)/0.1 mm))` raises it for coarse angular grains.
Overburden is the `Σ(ρdz + water)·g` integral taken at the cell midpoint — the paper's
half-own-weight rule for the surface layer, applied uniformly. Integration uses
`D·exp(−σ/η·dt)`, eq. 5's exact solution at constant σ and η.

Two departures from the published law:

- `f2` is floored at 1. Eq. 9 is bounded above but not below, so it *softens* fine-grained
  snow — 0.368, a 2.7× speedup, at GEMB's 0.05 mm fresh-snow radius. Crocus never reaches
  that regime because eq. 9 applies only to non-dendritic snow, whose grain size the paper
  puts at 0.3–0.4 mm where eq. 9 gives 2.7–4. GEMB carries one grain radius for both
  regimes, so the floor keeps `f2` a stiffening correction over the domain eq. 9 covers and
  inert below it.
- `:Crocus` hands cells at or above 450 kg m⁻³ to `:GSFC2020`. Eq. 7 is fitted to a 1–2 m
  alpine snowpack and `exp(b_η·ρ)` saturates in firn: unblended, the law gives 0.7–0.8×
  `:Arthern`'s compaction rate in the top few metres but only 0.02–0.09× below 20 m.
  450 kg m⁻³ is the threshold the Community Firn Model uses for the same purpose. The
  handover is a branch rather than a weighted blend, since neither paper prescribes a
  blending function. `:CrocusPure` applies eq. 5 at every density.

The paper's `c_η = 250` and its `f2` are used, where the Community Firn Model hardcodes
`f2 = 4` and adopts van Kampenhout et al. (2017)'s retuned `c_η = 358`.

**`:Barnola1991`** is Herron and Langway (1980) stage 1 below 550 kg m⁻³ and Barnola et al.
(1991) pressure sintering above, `dρ/dt = ρ·A₀·exp(−Q/RT)·f(ρ)·σ³` with
`A₀ = 2.54e4 MPa⁻³ s⁻¹` and `Q = 60 kJ mol⁻¹`. It is the only scheme here that treats the
firn–ice transition mechanistically rather than extrapolating a snow/firn fit: `f(ρ)` is a
polynomial fitted to Pimienta and Duval (1987) below 800 kg m⁻³ and the analytic
isolated-spherical-pore form `(3/16)(1 − ρ/ρᵢ)/(1 − (1 − ρ/ρᵢ)^⅓)³` above, which vanishes
with porosity, so the law self-limits as ρ → ρᵢ rather than being driven to zero by a
`(ρᵢ − ρ)` factor. Overburden and integration follow `:Crocus` and forward Euler
respectively. Stage 1 is bit-for-bit `:HerronLangway`, sharing one kernel. Cross-checked
against the Community Firn Model, which agrees on both branches and all four polynomial
coefficients.

Four constraints and departures:

- The self-limiting behavior holds only if the closed-pore branch is reached, and the
  800 kg m⁻³ handover is absolute while the polynomial carries no `ρᵢ`. `:Barnola1991`
  therefore requires `density_ice >= 900`, enforced at construction, against `[800, 950]`
  elsewhere.
- The paper's driving stress is overburden *minus bubble pressure*; only the overburden is
  applied, so effective stress is overstated above 830 kg m⁻³ where pores close. The paper
  gives no expression for the bubble term. Goujon et al. (2003) eqs. A11–A12 do, so it
  belongs in a separate scheme.
- Only the closed-pore branch scales with `density_ice`, and the branches meet in value at
  ρᵢ = 919.96 and in slope at 920.06 against the paper's stated C¹ matching, so the fit
  assumes ρᵢ ≈ 920. At GEMB's default 910 the rate steps down 14% crossing 800 kg m⁻³ — a
  discontinuity in `dρ/dt`, not in ρ.
- The 550 kg m⁻³ handover is also a rate discontinuity, by construction in the paper: it
  joins a depth-independent accumulation-driven rate to a σ³ one. On a 250 K,
  300 kg m⁻² yr⁻¹ column, sintering is ~1.3e4× slower than Herron–Langway just above 550 at
  1 m depth and overtakes it only near 13 m.

`n = 3` throughout, as in the paper, whose `A₀` was fitted against it over
0.55–0.8 g cm⁻³. The paper's `n = 1` remark scopes itself to bulk ice below close-off
(Pimienta and Duval torsion-tested 0.85 g cm⁻³ ice and report no densification rate), so it
is not an unimplemented firn branch. Calibrated over −14 to −57 °C and
2.2–65 g cm⁻² yr⁻¹, with no liquid-water term.

**Accumulation forcing.** Accumulation-driven schemes (`:Arthern`, `:Ligtenberg`,
`:Simonsen2013`, `:GSFC2020`, `:HerronLangway`, `:Crocus`, and `:Barnola1991` below its
stage transition) compact against `accumulation_mean`, which partitions precipitation on
`rain_temperature_threshold`, since rain does not bury the column. This is the same split
`initialize_climate_summary` applies, so the transient run and the initializer use the same
accumulation flux. `densification_accumulation = :precipitation` drives compaction with
total precipitation instead; on a site with 3.2% rain the difference in 32-year final mean
column density is ~1.1 kg m⁻³.

**Mean temperature.** `exp(E/RT)` is convex in `1/T`, so `exp(E/R⟨T⟩) ≠ ⟨exp(E/RT)⟩` — at a
site with a ±15 K seasonal amplitude about 250 K the two differ by tens of percent, in the
direction of under-densification. `mean_temperature_method = :arrhenius` uses the effective
temperature

```
T_eff = E_g / (R · log(mean(exp(E_g / (R·T_air)))))
```

with `E_g = 42 400 J mol⁻¹`. One `T_eff` is exact for all three consumers
(`:ArthernB`, `:GSFC2020`, and Simonsen's γ) only because they share that activation energy;
a scheme with a different `E_g` needs its own. The default `:arithmetic` uses the plain mean.

**Domain calibration.** `:Ligtenberg`'s `M0`/`M1` multipliers carry nine calibration sets
across Antarctica and Greenland (`densification_lookup.jl`). They are not interchangeable
with other models' fits, which were regressed against different forcing products and a
different Arrhenius reference temperature.

## Grain growth

Grain metamorphism runs on every timestep, independently of `albedo_method`. Four schemes
read `grain_radius` — `:ArthernB`, `:Crocus` and `:CrocusPure` densification, and the
`:grain_radius_threshold`/`:grain_radius_w_threshold` emissivity methods — and each selects
independently of the albedo scheme.

`grain_growth_method` governs the dry non-dendritic branch:

- **`:Marbouty`** — Marbouty (1980) temperature-gradient metamorphism. Its density factor is
  identically zero above 400 kg m⁻³ and its gradient factor is meaningless where `∇T → 0`,
  so it does not evolve grain radius in deep firn.
- **`:Arthern`** — Arthern et al. (2010) normal growth, `dr²/dt = k_gr·exp(−E_g/RT)`,
  density-independent. This is the default and what the Community Firn Model uses.
- **`:hybrid`** — Marbouty below 400 kg m⁻³ and Arthern at or above it, keeping the
  temperature-gradient physics in seasonal snow while grains continue to coarsen in firn.

The dendritic and wet (Brun, 1989) branches and the `GRAIN_DIAMETER_MAX` cap are common to
all three. Because `:ArthernB` densification goes as `1/r²`, faster grain growth densifies
*more* slowly: on synthetic forcing with `:ArthernB`, final mean column density is
902.6 kg m⁻³ under `:Marbouty` (median grain radius 0.41 mm), 891.4 under `:hybrid`
(0.589 mm) and 889.3 under `:Arthern` (0.592 mm).

## Fresh snow

`new_snow_method` selects the fresh-snow density:

| Option | Form |
|---|---|
| `:Constant350` (default), `:Constant315`, `:Constant150` | Bare constants; 350 is what the Community Firn Model ships |
| `:Fausto` | 315 kg m⁻³, the Fausto et al. (2018) regression at 256.2 K |
| `:FaustoFit` | `362.1 + 2.78·(T_air − 273.15)`, the regression itself, clamped to `[1, ρᵢ]` since it reaches zero near 143 K |
| `:Pahaut` | `max(50, 109 + 6·(T_air − 273.15) + 26·√U)` |
| `:Kaspers`, `:KuipersMunneke` | Temperature-dependent polar fits |

`:Pahaut` is Pahaut (1975) as implemented in Crocus, via Lafaysse et al. (2026) eq. 35. It is
the only option fitted to alpine seasonal snow rather than to a polar ice sheet, and it is
much lighter than the polar fits over their common range — 107 kg m⁻³ at −10 °C and 5 m s⁻¹
against `:FaustoFit`'s 334 — because alpine snowfall is warmer, wetter, and far less
wind-packed than the katabatic-scoured surfaces the Greenland and Antarctic fits were
regressed on. Prefer it for temperate and mid-latitude glaciers, where the polar fits
overestimate fresh-snow density, suppressing new-snow albedo and speeding its burial. The
published 50 kg m⁻³ floor is retained; it binds below about −10 °C in calm air.

`:Fausto`, `:FaustoFit` and `:Pahaut` also select the Crocus (Vionnet et al., 2012)
wind-dependent fresh-grain dendricity and sphericity, those being what Crocus itself pairs
these densities with. `fresh_snow_density` takes instantaneous air temperature and wind
speed where they are available and the climatological means in the steady-state march, so the
initial guess is well defined at every setting.

## Thermal physics

**Heat capacity.** `heat_capacity_method = :constant` holds `heat_capacity_ice`
(2102 J kg⁻¹ K⁻¹, the value at the melting point), matching the Community Firn Model's
melt-enabled path. `:CuffeyPaterson` is `c_p(T) = 152.5 + 7.122·T` (Cuffey and Paterson,
2010, eq. 9.1); the constant form overstates the cold content of firn at 240 K by 12.9% and
at 210 K by 27.5%.

**Internal energy is the enthalpy integral `∫c_p dT`**, not `M·T·c_p`. The latter is valid
only for constant `c_p`; under `:CuffeyPaterson` it overstates enthalpy by `(b/2)T²`, which
is 0.79·`L_f` at the melting point. Cell mixing, the thermal solver, and every energy budget
work in enthalpy. Absolute column enthalpy is not comparable between the two methods
(574 061 against 307 385 J kg⁻¹ at 273.15 K), so reported thermal energy shifts when
switching.

**Thermal conductivity.** `thermal_conductivity_method = :Calonne2019` is the default,
Calonne et al. (2019) eq. 5 — a sigmoid blend at ρ = 450 kg m⁻³ between the Calonne (2011)
snow quadratic and a firn branch, scaled by the temperature-dependent conductivity of ice.
It is continuous into ice by construction and so is deliberately not short-circuited at
`density_ice`; the `917` in its firn branch is the paper's fitted pure-ice density, not
`mp.density_ice`. RetMIP Sect. 5.1 recommends it by name, and it is what both the Community
Firn Model and IMAU-FDM ship.

`:Calonne2019Air` additionally carries eq. 5's air-conductivity ratio on the snow branch,
using Reid (1966) for air. See [the IMAU-FDM comparison](imau_fdm_comparison.md) for the
magnitude. `:Marchenko2019` is Marchenko et al. (2019) eq. 30, `k = 0.301e-2·ρ − 0.724`,
fitted over ρ = 350–900 and floored at the Calonne (2011) value below their crossing near
ρ = 321, where the bare fit heads negative. Both `:Calonne2019` and `:Marchenko2019` are
higher than `:Sturm`/`:Calonne` in firn, the direction RetMIP's cold bias at Summit and
Dye-2 implies.

Face conductivity is the harmonic mean `1/(dz[i+1]/2K[i+1] + dz[i]/2K[i])` — the series
resistance of the two half-cells, which is the exact flux for a piecewise-constant
conductivity where an arithmetic mean is not.

**The explicit sub-step limit.** The explicit solve is stable when each cell's
own-temperature coefficient stays non-negative, `dt ≤ ρᵢcᵢdzᵢ/(Gᵢ + Gᵢ₋₁)`, where `Gᵢ` is
the harmonic-mean face conductance the flux loop already uses. The Dirichlet bottom cell is
excluded, its enthalpy never being updated. This is the limit for the graded grid GEMB
actually solves on; the textbook uniform-grid form `0.5·ρᵢcᵢdzᵢ²/Kᵢ` substitutes `2Kᵢ/dzᵢ`
for `Gᵢ + Gᵢ₋₁`, which agrees exactly only when `dz` and `K` are uniform. Because
`Gᵢ + Gᵢ₋₁` ranges over `(0, 4Kᵢ/dzᵢ]`, the uniform form errs in both directions — measured
per cell on GEMB's own column, the ratio spans 0.66 to 1.82. `thermal_explicit_safety_factor`
(default 0.8) is not slack against an imprecise limit: `_max_safe_dt` bounds diffusion only,
while the surface cell also carries the surface energy balance's feedback `Λ = dQ/dT₁ ≤ 0`
in the same coefficient, and nothing in the model bounds `Λ`.

**Implicit surface-row damping.** `ImplicitThermal`'s Newton step weight is halved whenever a
step fails to shrink, and the least-residual iterate is retained at the iteration cap.
Without damping, 39% of sub-step solves reach the cap in limit cycles of ~0.02 K median
amplitude straddling the latent-heat switch at 273.15 K, rather than in divergence; damping
cuts that to 0.9% at no measurable runtime cost. See
[Thermal solvers](@ref "Thermal solvers").

## Turbulent fluxes

**Integrated stability functions.** Paulson's (1970) closed forms on the unstable branch take
the *inverse* profile function as their argument, so the exponents are positive:
`x_m = (1 − 19ζ)^¼` inverts Högström's (1988) `φ_m = (1 − 19ζ)^(−¼)`, and likewise
`x_h = (1 − 11.6ζ)^½`. Since `Ψ(ζ) = ∫₀^ζ (1 − φ)/z dz` must vanish at neutral stability,
Högström's `κ_H/κ_M = 0.95` ratio belongs to the neutral transfer coefficient and not inside
the integral — carried inside, it puts `Ψ_h(0) = −0.0999` and leaves the turbulent fluxes
discontinuous exactly where `T_surface` crosses `T_air`. That crossing is where the implicit
solver's Newton iteration converges, and where a jump can leave the surface energy balance
with no solution in `T_surface` (Fourteau et al., 2024, Appendix D, who move their own branch
point to `Ri_b = 0` for the same reason).

Both branches — this one and the stable Beljaars and Holtslag (1991) eqs. 28 and 32 — are
checked against the definition directly, by numerically integrating their own published `φ`
and comparing to the closed form, agreeing to 6 decimal places at `ζ = ±0.01, ±0.1, ±1, ±5`.

**`ζ` is bounded below by −100 on the unstable branch.** This is a numerical guard on the
bulk formulation, not a change to the physics. `ζ` is diagnosed from the bulk Richardson
number, which carries `wind_speed⁻²`, so at the `min_wind_speed = 0.01 m s⁻¹` floor over a
melting surface under very cold air the synthetic forcing reaches `ζ ≈ −4.9e5`. That is the
wind floor showing through rather than a stability regime, and Monin–Obukhov theory has no
observational support anywhere near it (Högström's fits span `|ζ| ≲ 2`). Unbounded, `Ψ_h`
grows past `log(z_T/z_Q)`, so `coefHT` crosses zero and the flux diverges and changes sign.
−100 is two orders of magnitude beyond the calibration range, so it never binds in physically
meaningful conditions, and it keeps the transfer coefficients positive for every roughness
GEMB uses. Bounding `ζ` rather than clamping the coefficients keeps the fluxes a continuous,
monotone function of `T_surface`, which the implicit solver needs.

## Meltwater

**Irreducible saturation.** `water_irreducible_method = :ColeouLesaffre` is the default:
`S_wi = wmi/(1 − wmi)·ρᵢρ/(ρ_w(ρᵢ − ρ))` with `wmi = 0.057(ρᵢ − ρ)/ρ + 0.017`, Coléou and
Lesaffre (1998) eq. 3 via Langen et al. (2017) eq. 4. The `1/(1 − wmi)` factor converts the
paper's wet-mass fraction to a dry-mass basis. Retention is zero at and above pore
close-off. Both the Community Firn Model and IMAU-FDM use this form, and RetMIP Sect. 5.4
attributes part of the multi-model under-retention in the percolation zone to a flat 0.07.
`:constant` holds `water_irreducible_saturation` at every density, at every retention site.

**Melt geometry.** `melt_geometry = :thickness` (default, as in Crocus) holds density fixed
and shrinks `dz`; `:density` holds `dz` fixed and lowers density, as SNOWPACK does and as
Fourteau et al. (2026) Sect. 2.3 argues for, on the grounds that the phase change occurs
within the microstructure and that the high density of wet snow is better explained by its
low viscosity under overburden. Refreezing is at constant thickness under both, so only
`:density` returns a cell's geometry along with its mass over a melt–refreeze cycle. The
choice affects melting cells' `dz` and `density`, and through the irreducible-retention
capacity their water and runoff.

**Impermeability.** Water is routed to runoff at a contiguous run of cells at or above
`impermeable_density` (default 830 kg m⁻³) thicker than `impermeable_thickness` (0.1 m).
RetMIP recommends bucket schemes adopt such a criterion and finds model spread at ice-slab
sites dominated by this pair; the participating models span 810 to 917 kg m⁻³. The
`:ColeouLesaffre` retention gate stays on the constant `DENSITY_PORE_CLOSEOFF` — where
capillary retention ceases for want of connected pore space is a different question from
where a lens stops conducting flow — so lowering the flow criterion does not also change
retention.

**Runoff timescale.** `runoff_method = :instantaneous` is the default: all blocked water
leaves within the timestep, which is what all three RetMIP bucket lineages and both
comparison models do. `:ZuoOerlemans` drains `min(1, Δt/τ)` of the excess with
`τ = c₁ + c₂·exp(−c₃·S)`, `c₁ = 0.33 d`, `c₂ = 25 d`, `c₃ = 140` (Zuo and Oerlemans, 1996,
eqs. 21–22, coefficients via Lefebre et al., 2003 / Langen et al., 2017). `:Darcy` drains
`min(excess, ρ_w·Δt·K_sat·K_rel·S)` with `K_sat` from Calonne et al. (2012) eq. 6 and `K_rel`
from van Genuchten (1980) with the Yamaguchi et al. (2012) parameterization. Both read
`surface_slope` as the hydraulic gradient. These are the schemes used by the two RetMIP
models with the lowest firn-temperature error at the ice-slab site KAN\_U (−1.6 °C and
+0.6 °C, against a spread reaching +4.7 °C).

Two published lineages disagree on `c₁` — the Community Firn Model uses 1.5 d while citing
the same source. GEMB follows Langen et al., the lineage RetMIP evaluated; the difference is
immaterial at ice-sheet slopes (21 d either way at `S = 0.01`).

**Saturated pore space.** Under `:instantaneous`, pore water is clamped to irreducible at
every retention site, so a saturated cell cannot exist and firn aquifers are structurally
unrepresentable — RetMIP Sect. 5.4 dropped their aquifer site from the retention evaluation
for this reason. When `runoff_method !== :instantaneous`, water blocked at a barrier or the
column base backs up into the pore space above it, filling each cell to `pore_saturation_max`
of its capacity from the barrier upward; only what reaches past the surface cell runs off.
Cells at or above `impermeable_density` are skipped as having no connected pore space. The
upward pass does not refreeze — those cells' cold content was consumed during percolation —
and a cold cell left holding liquid is resolved by the refreeze-pore-water step on the next
timestep.

**Rain sensible heat** above the melting point is carried at
`HEAT_CAPACITY_WATER = 4219.9 J kg⁻¹ K⁻¹` through the two-argument
`specific_enthalpy_water(mp, T)`. Pore water is pinned at the melting point, where the
one-argument form applies exactly. `rain_heat_capacity = :ice` carries it at
`heat_capacity_ice` instead, which understates it by about a factor of two — ~3.9 kJ kg⁻¹
against 334.5 kJ kg⁻¹ of latent heat for rain at 275 K, growing with `T_air − 273.15`. On the
synthetic site, 32-year melt differs by 48.7 kg m⁻² (0.43%).

**Percolation order.** The downward flux is computed after refreezing, as in Crocus, and
impermeability is re-tested after refreeze so a cell that has just become a barrier blocks
water in the same timestep. See
[the CFM comparison](cfm_comparison.md#Differences-retained) for the alternative ordering and
why no option is offered.

## Blowing and drifting snow

**`blowing_snow_method = :Crocus`** [default `:none`] is the SURFEX/Crocus `SNOWDRIFT`
scheme: a per-layer mobility index (Lafaysse et al., 2026, eqs. 59–61, the `B92` branch,
which is the one written in dendricity, sphericity and grain size as GEMB's state is) and
driftability index (eq. 62), an effective driftability decayed exponentially by the overburden
above each layer (eq. 63), wind-slab densification toward 350 kg m⁻³ on a 2 d timescale
(eq. 64), and grain fragmentation toward smaller, rounder, non-dendritic grains (eqs. 65–66).
The layer loop stops at the first non-drifting layer, as the Crocus source does and eq. 63
does not show. It affects density, `dz`, grain radius, dendricity and sphericity, and through
them albedo, conductivity, densification and percolation. It is **exactly mass-conserving** —
`dz` shrinks in proportion to the density rise — so it appears in no mass budget.

Two departures from Crocus: GEMB has no wet-snow history variable, so the eq. 61 mobility cap
uses instantaneous liquid water rather than "has ever been wet", which leaves
refrozen-but-now-dry layers more driftable; and melt-style thickness bookkeeping is not
applied, the compaction being density-only.

**`blowing_snow_sublimation`** [default `false`] sublimates the suspended fraction after
Gordon et al. (2006) as Crocus implements it (eqs. 67–68), read only under
`blowing_snow_method = :Crocus`. Off is Crocus's own default. This is the only part of the
computed scheme that moves mass across the surface: it removes mass from the surface cell,
capped at half of it per timestep, and reports it in `blowing_snow` and in the mass budget.
**No latent heat is charged to the column** — the particles are suspended in the air and draw
their heat of sublimation from it, unlike surface sublimation, which is charged. The Crocus
paper (Lafaysse et al., 2026, Table K3) and its source disagree on which exponent belongs to
the temperature ratio and which to the wind ratio; the source is followed, being the
evaluated code, and the discrepancy is recorded in `src/blowing_snow.jl`.

**`drift_rate`** [kg m⁻² yr⁻¹, default 0, bounded ±2000] and the optional `snow_drift`
forcing layer are the prescribed path, positive for erosion. Drift is a surface mass
source/sink with a sign-dependent density: eroded mass leaves at the surface cell's own
density, deposited mass arrives at `fresh_snow_density` and at age zero, diluting the surface
cell's age as condensation does. Taking the cell density in both directions would let an
erode-and-redeposit cycle change surface density with no net mass exchange. Both are read
independently of `blowing_snow_method`, so a regional-climate-model drift field can drive
GEMB with no internal scheme; the forcing layer takes precedence where present.

## The vertical grid

The column holds a **constant cell count** and a **constant total depth**, both fixed at
initialization and enforced every timestep by two independent controllers: count is restored
by exactly conservative merge and split operations at depth, and depth by a continuous signed
adjustment to the bottom cell, the model's only basal mass and energy flux. Mass and energy
are conserved throughout — checked every timestep under `verbose=true`, with the whole-run
budget closing to ~1e-11 kg m⁻².

Consequences for consumers of the output:

- **Profile arrays are sized exactly to the column and top-justified**, the surface always at
  row 1. There is no padding and no need to scan for `NaN`.
- **`column_depth_max` is a ceiling on the constructed column depth**, not the depth itself.
  The depth is derived by `initialize_profile` from the climate and then held.
- **Albedo is not carried state.** All four methods (`:None`, `:GardnerSharp`,
  `:BrunLefebre`, `:GreuellKonzelmann`) are diagnostic functions of the current column, so
  `calculate_albedo` returns two scalars rather than per-cell vectors. The surface series
  `albedo_broadband` is recorded from the value used in that timestep's shortwave balance.
- **`ice_flux` is a signed per-interval basal flux**, so `cumsum(ice_flux)` is the cumulative
  flux exactly. Its negation is surface elevation change against a datum fixed in the ice;
  see [`trim_bottom!`](@ref) for the identity.

**`horizontal_strain_rate`** [yr⁻¹, default 0] is the trace of the horizontal strain-rate
tensor `ε̇_xx + ε̇_yy`. By incompressibility it thins (positive, divergence) or thickens
(negative, convergence) every cell at constant density by `exp(−D·dt)`. The mass it exports
leaves laterally and is reported as `strain_thinning`, separately from the basal flux.

## Initialization and spinup

**Deep temperature.** The initialized deep temperature is the mean *surface* temperature,
accumulated from the skin temperatures `_seb_annual_melt` already solves for, not the mean air
temperature. The column is coupled to the atmosphere only through the surface energy balance,
so its deep mean tends to the surface mean, and the two differ by the whole radiative and
turbulent budget. This is not only an initialization detail: the deepest cell is a Dirichlet
reservoir, so the initialized value is a boundary condition no length of spinup relaxes.
Across a 21-site synthetic fleet the air-temperature value is biased +11.6 K warm on average
and never cold (worst +24.5 K) against the self-consistent value, and the surface-temperature
form cuts the mean temperature jump across that frozen cell from 0.79 K to 0.37 K. It also
*increases* cycles to convergence (309 to 455 on that fleet): a slower-densifying colder
column starts further from its density attractor. Refreezing's latent warming decays over the
annual accumulation layer rather than being applied at every depth.

**No initialized cell starts above the melt point.** `initialize_profile` clamps the
temperature it fills to 273.15 K, with a warning. This binds only at sites whose mean annual
air temperature is above freezing, where the column would otherwise begin as ice above its
melting point holding no water — enthalpy the column has no state to carry.

**Age.** `initialize_age = :steady_state` (default) carries the age `steady_state_profile`
integrates along its march, the same march temperature, density and all three grain variables
are initialized from. `:zero` starts every cell at zero, and the `constant_density` escape
hatch forces `:zero` regardless, since it discards the march's density and the marched age
describes a firn column that flag replaces with solid ice.

The `age` profile field is the mass-weighted mean age, in decimal days, of all mass in a cell
(matrix plus pore water), measured from column initialization. Snowfall, rain and vapour
deposition enter at age 0; proportional mass removal (melt, sublimation, runoff, basal trim
under accumulation, horizontal strain) is age-neutral; merges mass-weight it and splits
duplicate it. Meltwater carries the age of the firn it melted from and mass-weights that age
into the cell where it refreezes, so a heavy melt year does not read as artificially young
firn at depth. Age accumulates across `gemb_spinup` cycles, so for a spun-up run the epoch is
the start of spinup. No physics reads it. Under sustained ablation the deepest cell's age is a
lower bound, basal accretion inheriting that cell's own age — the same treatment `density` and
`temperature` get at that site.

**Convergence.** Spinup convergence is judged on the mass-weighted mean density of the whole
column. Because the column depth is fixed for the run, the full column is the same domain at
every cycle and the mean is exact without interpolation. `convergence_delta_density` bounds
the change between consecutive cycles; `convergence_drift_density` bounds the least-squares
slope of column-mean density against cycle over the trailing `drift_window` cycles (default
10), in kg m⁻³ per cycle. The trend test is the stricter claim, since a column creeping
steadily at just under the delta tolerance passes the step test while still densifying. When
both are given, both must hold.

**Spinup forcing.** `forcing_climatology(method=:average)` (default) averages complete years
into a one-year cycle. It preserves each field's mean but shrinks its variance by
~`1/n_years`, and melt is *rectified* — zero until the surface energy balance reaches the melt
point — so averaging cancels the warm excursions that carry it.
`forcing_climatology(method=:representative)` returns a block of `n_years` consecutive real
years scored against the record instead.

Measured on ERA5-Land forcing at 12 glacierized Greenland sites (2000–2019), melt retained on
the averaged cycle has a median of 12% of the record's, with 5 of 10 melting sites below 10%,
against a median 102% for a 3-year representative block. The failure tracks elevation and
absolute melt magnitude, not the melt/accumulation ratio: 0% retention at Saddle (2456 m) and
DYE-2 (2094 m), 86–115% at the three lowest sites (271–577 m). High sites melt from the tail
of the distribution, which averaging removes; low sites melt from the mean, which it
preserves.

Averaging remains the default because it is right where melt is negligible: it preserves
accumulation to 0.1–0.5% at the dry sites, against 15–22% for a melt-ranked block, and
integrates less than half the model-years. `rank_by` selects the scoring variable (`:model`
melt, `:smb` surface mass balance, `:estimate` the cheap surface-energy-balance melt), and
`fallback_accumulation_tolerance` (default 0.05) re-selects on surface mass balance when a
melt-ranked block's accumulation is too far off the record — measured to cut the worst
accumulation error from 30.2% to 11.3% and the median from 6.9% to 2.5% across those sites.

Two limits are recorded rather than hidden: every figure above is *forcing* fidelity, since
the equilibrated-column comparison could not be completed (6 of 12 sites do not converge in
800 cycles), and `:representative` costs about 2.2× the model-years of `:average`.

## Diagnostics

These take no part in the mass or energy budget.

- **`firn_air_content`** is metres of air, `Σ dz(1 − ρ/ρ_ice)`. `firn_air_content_10m` and
  `firn_air_content_20m` are the same integral to a fixed depth, which is how published
  values are reported (Vandecrux et al., 2019).
- **`percolation_depth`** is an interval maximum, matching what upward-looking radar measures
  (Heilig et al., 2018).
- **`ice_slab_thickness`, `ice_slab_depth`, `aquifer_thickness`, `aquifer_depth`** are
  instantaneous, so a `NaN` meaning "no slab" or "no standing water" cannot poison an interval
  mean. The slab terms are scanned after the grid controllers run, so recomputing slab depth
  from the output `dz` and `density` reproduces them. Aquifer detection is thresholded at
  `AQUIFER_TOLERANCE`, 1e-3 of a cell's pore space, rather than at the physics tolerance:
  `calculate_melt` leaves a retaining cell at exactly irreducible, but `calculate_density`
  then compacts it within the same timestep, shrinking the pore space retention was computed
  against and leaving a genuine ~1e-7-of-pore-space excess that accumulates between melt
  events.
- **`close_off_age`** is the age at the shallowest cell reaching pore close-off
  (830 kg m⁻³), `NaN` for an open column.
- **`heat_flux_basal`** is the conductive flux across the deepest interior face. It is a flux
  diagnosed against a Dirichlet reservoir, not a prescribed geothermal flux, and it tends to
  zero under a converged spinup — hence the name.
- **`viscosity`** [Pa s] is an optional profile layer under `output_viscosity = true`. It is
  populated by the `:Crocus`/`:CrocusPure` settling law and is `NaN` under every other scheme,
  which form no effective viscosity: they are `ρ̇ = c(ρᵢ − ρ)` relaxations with no stress in
  them, so a reconstructed `η = σ/(ρ̇/ρ)` would report a number derived from a stress the
  scheme never used.
