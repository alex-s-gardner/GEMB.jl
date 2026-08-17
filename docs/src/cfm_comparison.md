# GEMB and the Community Firn Model

The [Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel)
(CFM; Stevens et al., 2020, *Geoscientific Model Development* 13, 4355–4377) is
the other open-source firn model with a physics surface comparable to GEMB's.
This page records a capability-by-capability comparison of the two, made from a
full read of CFM's `CFM_main/` against GEMB's `src/`, together with what GEMB
adopted from it and — more often — what it deliberately did not.

The two models are largely **complementary rather than redundant**. They were
built for different questions: CFM for ice-core interpretation and firn
densification research, GEMB for surface mass and energy balance over ice sheets.
The differences follow from that.

Citations below point at CFM files as `physics.py:1508` and at GEMB as
`src/calculate_density.jl`. CFM line numbers are from the state of `main` read in
August 2026.

## Summary

| Capability | GEMB status | Notes |
|---|---|---|
| **Surface energy balance** | | |
| Albedo model | **have**; CFM lacks | GEMB has three schemes (`src/calculate_albedo.jl`); CFM takes albedo as prescribed input |
| Shortwave penetration with depth | **have**; CFM lacks | `src/calculate_shortwave_radiation.jl`; CFM has only a commented-out `rad_pen` sketch in `diffusion.py` |
| Turbulent fluxes (Monin–Obukhov / bulk Richardson) | **have**; CFM lacks | GEMB computes them; CFM takes QH/QL as *inputs* |
| Energy-coupled sublimation | **have**; CFM lacks | |
| **Densification** | | |
| Scheme count | 9 enabled (2 further gated) vs CFM's ~24 | see [Densification schemes](#densification-schemes) |
| Stress/microstructure family | partly (`:ArthernB`, `:Barnola1991`, `:Crocus`) | `Morris2014`, `Goujon2003`, `Breant2017` etc. not adopted this pass |
| Accumulation forcing for densification | **gap, adopted** | Fix 1 below |
| Arrhenius-weighted mean temperature | **gap, adopted** | Fix 2 below |
| Strain softening | **gap, deferred** | mass-vs-thickness discrepancy noted below |
| **Grain growth** | | |
| Arthern (2010) `dr²/dt` normal-growth law | **gap, adopted** | Fix 3 below |
| Temperature-gradient (Marbouty) metamorphism | **have**; CFM lacks | `src/calculate_grain_size.jl` |
| Linow (2012) surface grain size | **not adopted** | GEMB uses a fixed `RE_NEW_SNOW` |
| **Hydrology** | | |
| Bucket scheme | **have**, both | |
| Richards-equation percolation | **gap, deferred by decision** | see below |
| Dual-domain preferential flow | **gap, deferred by decision** | see below |
| Darcy lateral drainage | **have**, both, same parameter set | `src/hydrology.jl`, cross-checked against CFM's `hydrconducsat_Calonne` |
| Ponding / impermeable-layer runoff | **have**, both, independently identical defaults | see below |
| Percolation ordering within a cell | **open question, not adopted** | see below |
| Liquid water in thermal conductivity | **gap, deferred** | needs a hot-loop signature change |
| Liquid water in cell heat capacity | **not applicable** | architecturally incompatible; see below |
| **Numerics** | | |
| Implicit tridiagonal thermal solver | **implemented, opt-in** (`ImplicitThermal`) | see below |
| Multi-resolution subgrids | **not applicable** | GEMB solves this differently |
| Analytic Herron–Langway spin-up | **not a gap** | GEMB's march is more general |
| **Ice-core capabilities** | | |
| Firn air / gas transport | **not applicable** | outside an SMB model's remit |
| Isotope diffusion | **not applicable** | as above |

## What GEMB adopted

Reading CFM against GEMB surfaced a small set of findings that are *not* about
CFM's extra capabilities but about places where GEMB's own physics was internally
inconsistent, or where a cheap, well-cited parameterization was missing. Those
were fixed. Each is gated on a `ModelParameters` `Symbol` so the prior behaviour
remains reproducible bit-for-bit, and each appears in the **Physics notes**
section of `README.md`.

### Fix 1 — densification driven by mean snowfall, not mean precipitation

`calculate_density` compacted against `cfs.precipitation_mean`, which **includes
rain**. Rain does not bury the column, so every accumulation-driven scheme
(`:Arthern`, `:Ligtenberg`, `:Simonsen2013`, `:GSFC2020`) overstated the burial
rate at any raining site. GEMB already knew better in one place:
`initialize_climate_summary.jl` splits snow from rain on
`mp.rain_temperature_threshold` and hands `steady_state_profile` an
`accumulation_effective`. Initialization and the transient run therefore used
different, inconsistent accumulation forcings — and the initializer's was the
correct one. CFM's `GSFC2020` (`physics.py:1171`) rejects `bdot_type='instant'`
for the same class of reason.

`ClimateForcing` metadata gained `accumulation_mean`, derived by the same
snow/rain split, as a *fraction* of `precipitation_mean` so a caller-supplied
climatological mean keeps its own scale. Gated on
`mp.densification_accumulation ∈ (:accumulation, :precipitation)`, default
`:accumulation`.

### Fix 2 — Arrhenius terms evaluated at the arithmetic mean temperature

`exp(E/RT)` is strongly convex in `1/T`, so `exp(E/R⟨T⟩) ≠ ⟨exp(E/RT)⟩`. GEMB fed
a plain arithmetic mean air temperature to three Arrhenius grain-growth factors
(`_arthern_scaled_c`, `_gsfc2020_c`, Simonsen's `γ`). At a site with a ±15 K
seasonal amplitude about 250 K the two differ by tens of percent, biasing toward
under-densification. CFM handles this explicitly with `effectiveT()` in
`RCMpkl_to_spin.py`.

`ClimateForcing` metadata gained `temperature_air_effective`:

```
T_eff = Eg / (R · log(mean(exp.(Eg ./ (R .* temperature_air)))))
```

with `Eg = GRAIN_GROWTH_EG = 42400`. A single `T_eff` is exact for all three
consumers **only** because they share that activation energy; a scheme with a
different `Eg` would need its own. Gated on
`mp.mean_temperature_method ∈ (:arithmetic, :arrhenius)`, default `:arithmetic`;
`:arrhenius` is the physically preferable setting.

### Fix 3 — grain growth stopped dead at 400 kg m⁻³

GEMB evolved `grain_radius` only through seasonal-snow parameterizations:
dendritic metamorphism, Marbouty (1980) dry, Brun (1989) wet. Marbouty's density
factor is zero above 400 kg m⁻³ and its temperature-gradient factor is
meaningless at depth where `∇T → 0`, so grain radius was effectively **frozen
throughout the firn column** below a few metres. Meanwhile `:ArthernB`
densification goes as `σ/r²`, making deep densification maximally sensitive to
exactly the quantity that had stopped evolving.

CFM evolves `r²` with Arthern (2010), `dr²/dt = kgr·exp(-Eg/RT)`,
`kgr = 1.3e-7` m² s⁻¹ (`physics.py:1508-1577`, `GrGrowPhysics`) — density-
independent, and the standard firn law. GEMB gained
`mp.grain_growth_method ∈ (:Marbouty, :Arthern, :hybrid)`, governing the
non-dendritic *dry* branch only; dendritic and wet (Brun) branches and the
`GRAIN_DIAMETER_MAX` cap are unchanged. `:hybrid` uses Marbouty where it is
defined (ρ < 400) and Arthern at or above it, so seasonal snow keeps the
temperature-gradient physics and firn grains keep growing. Default is
`:Marbouty` (bit-identical to prior behaviour); **`:hybrid` is the recommended
setting**, particularly with `:ArthernB`.

Measured effect, synthetic forcing, 75-cycle climatological spinup then a 32-year
3-hourly transient with `densification_method=:ArthernB` — final mean column
density, against the 877.7 kg m⁻³ reference in `README.md` deviation 1:

| `grain_growth_method` | mean column density | median grain radius |
|---|---|---|
| `:Marbouty` | 902.6 | 0.41 mm |
| `:hybrid` | 891.4 | 0.589 mm |
| `:Arthern` | 889.3 | 0.592 mm |

`:hybrid` roughly halves the discrepancy. Note the direction: because `:ArthernB`
goes as `1/r²`, *faster* grain growth densifies *more slowly* — the frozen-radius
column sits above the reference, so growing the grains is corrective.

### Fix 4 — the steady-state march computed an age profile, then discarded it

`initialize_profile` set `age .= 0.0` on every path, even though
`steady_state_profile` integrates age along the march and temperature, density,
grain radius, dendricity and sphericity are all initialized from that same march.
Consequence: any age or residence-time diagnostic was meaningless until a
multi-century spinup had flushed the column, making a close-off-age diagnostic
unobtainable in practice.

The marched age is now used, behind `mp.initialize_age ∈ (:zero, :steady_state)`,
default `:steady_state`. Nothing in the physics reads `age`, so budgets are
unchanged. A `close_off_age` diagnostic was added alongside — age at the
shallowest cell reaching `DENSITY_PORE_CLOSEOFF`, `NaN` when the column never
closes off, following `ice_slab_depth`'s convention.

### Fix 5 — rain sensible heat carried at the ice heat capacity

Rain entered at air temperature via `mix_temperature_liquid`, whose `:constant`
path reduces to carrying the rain's sensible heat above 273.15 K at
`c_ice = 2102` rather than ~4200 J kg⁻¹ K⁻¹ — understating it by about a factor
of two. There was no liquid-water heat capacity constant anywhere in `src/`. The
verbose budget used the same convention on both sides, so the conservation check
could not catch it.

`HEAT_CAPACITY_WATER = 4219.9` (CFM's `c_liq`, `solver.py`) was added, together
with a two-argument `specific_enthalpy_water(mp, T)`. The one-argument form
(melting-point water) is unchanged: it is used correctly throughout the
melt/percolation budgets, where pore water *is* isothermal at `CtoK`. Gated on
`mp.rain_heat_capacity ∈ (:water, :ice)`, default `:water`. For rain at 275 K
this is ≈3.9 kJ kg⁻¹ against 334.5 kJ kg⁻¹ of latent heat — about 1%, small but
systematic and growing with `T_air − 273.15`.

### Fix 6 — exposing basal heat flux and effective viscosity

Both quantities were already computed and discarded.

`heat_flux_basal` is the diagnosed conductive flux across the deepest interior
face. It was returned to `gemb_core` but used only inside the `verbose` energy
check. It is now a first-class output. It is deliberately *not* named `ghf`: it
is a flux diagnosed against a Dirichlet reservoir, not a prescribed geothermal
flux, and it tends to zero under a converged spinup.

Effective viscosity — `_crocus_viscosity` is the only one GEMB forms anywhere —
is now available as an optional profile output via `mp.output_viscosity`
(default `false`, so the hot path is untouched). CFM makes viscosity a
first-class output and consumes it in `strain.py`. GEMB populates it **only**
under `:Crocus`/`:CrocusPure`, and writes `NaN` elsewhere. This is a deliberate
narrowing of the plan's original intent: every other GEMB scheme is a
`ρ̇ = c(ρᵢ − ρ)` relaxation with no stress in it, so reconstructing
`η = σ/(ρ̇/ρ)` would report a number derived from a stress the scheme never
used, in units that invite comparison with a quantity it is not. `NaN` says "not
computed here", which is honest.

## Deliberately not adopted

### Liquid water in the cell heat capacity — architecturally incompatible

CFM's enthalpy solver uses `c_vol = (g_ice·c_ice + g_liq·c_liq)·ρ_tot`
(`solver.py`, `transient_solve_EN`). GEMB's cell mass counts matrix only
(`src/calculate_temperature.jl`), documented as intentional in `src/types.jl`.
**This should stay.** GEMB is split-operator: phase change happens in
`calculate_melt`, not in the thermal solve, and pore water is carried strictly at
`CtoK`. Giving that water a sensible heat capacity would let a wet cell bank
energy as water warmed above the melting point — a fictitious reservoir that
`melt_mass_from_excess` (which scales by matrix `ρ·dz`) cannot see, breaking the
melt handoff. CFM can do this only because latent heat is *inside* its solve.
Adopting it requires adopting CFM's enthalpy-with-phase-change solver wholesale;
it is not a drop-in.

The Von Neumann stability limit uses the same water-free mass, which is
*conservative* — more thermal mass would relax it.

### Liquid water in thermal conductivity — a genuine gap, deferred

CFM uses `K_eff = g_liq·K_liq + g_ice·K_firn` with
`K_liq = 0.55575·(ρ_liq/1000)^1.885`. This half **is** separable from the
phase-change question and is a real gap: all four GEMB methods are functions of
`(ρ, T)` only. It needs the full `water` vector, and `calculate_temperature`
currently receives only `water_surface::Float64`, used solely to pick a roughness
length. Deferred as a follow-up — not because it is wrong, but because the
signature change touches the hot loop and deserves its own benchmark pass.

### Richards-equation and preferential-flow percolation

CFM's `re_snowpack.py`, `prefflow_snowpack.py`, `darcy_funcs.py`. GEMB's own
docstring already justifies bucket-only, citing RetMIP
(`src/calculate_melt.jl`, Vandecrux et al. 2020 §5.2). CFM's own implementations
print "still in development" and only engage after a hardcoded model year 1980.
GEMB's `:Darcy` lateral-drainage rate law already shares CFM's van
Genuchten/Yamaguchi/Calonne parameter set (`src/hydrology.jl`).

### Percolation ordering within a cell — an open question, not adopted

Lafaysse et al. (2026) §2.4.18 notes that a bucket scheme may compute the downward
flux `Fᵢ` either *after* refreezing `rᵢ` (Crocus's order, and GEMB's) or *before*
it, treating ice layers as impermeable first (Fourteau et al. 2024, "meant to avoid
liquid water percolation through an entire glacier"). Neither is validated: Fourteau
never tests the ordering — their §5 concedes their own bucket scheme runs all melt
off without percolation — and Lafaysse hedges, flagging only "potential impacts on
glaciers simulations".

Note that Crocus itself has no sub-ice-density barrier at all: in its Algorithm 1 the
only special case is `ρ = ρ_I`, and even there the excess grows the ice layer's
thickness rather than being blocked. GEMB pairs Crocus's ordering with an
`impermeable_density` criterion inherited from the MATLAB model, a combination
neither upstream model has. That is what made the stale-density barrier test a defect
rather than a choice, fixed by re-testing impermeability after refreeze.

With that fixed, the two orderings are indistinguishable on GEMB's synthetic forcing:
the re-test diverts water 14 times per run, and in every case the cell below was
already a barrier, so the water ran off one cell deeper and the run totals are
bit-identical. An option is therefore not added — the settings would differ only in
code, not in output. Worth revisiting with real forcing at an ice-slab site (KAN_U,
Dye-2), where lenses form *within* permeable firn rather than adjacent to existing
ones, which is the regime in which the orderings can diverge.

### Firn air / gas transport and isotope diffusion

`firn_air.py`, `isotopeDiffusion.py`. Ice-core capabilities, outside a surface
mass-balance model's remit. Both ride CFM's shared `transient_solve_TR`, so they
become comparatively cheap **if** an implicit tridiagonal solver ever lands in
GEMB — the two are coupled decisions.

### Additional densification schemes

`Morris2014`/HL-dynamic with its advected `Hx`, `Goujon2003`, `Breant2017`,
`KuipersMunneke2015`, `Brils2022`, `Max2018`/`MaxSP`, `Yamazaki1993` fresh-snow
stage. Out of scope this pass. Worth recording that GEMB already has
`:Barnola1991`, an `Arthern2010T`-equivalent (`:ArthernB`), `:GSFC2020`,
`:Simonsen2013`, `:Ligtenberg`, `:Crocus` and `:CrocusPure`, and that GEMB's
`:ArthernB` was validated against CFM's `Arthern2010T` to 1e-8 (`README.md`
deviation 3) — the cross-validation path exists and works.

### Implicit tridiagonal thermal solver — implemented, opt-in

Closed. `ImplicitThermal` (`src/calculate_temperature.jl`) is a backward-Euler
scheme on the same tridiagonal system, selected by dispatch on
`mp.thermal_solver`. `ExplicitThermal` remains the default and is bit-identical
to before, so existing runs are unaffected.

The blocker CFM's structure did not resolve was the nonlinear surface row, which
IMAU-FDM's tridiagonal also lacks an analogue for; see
[the IMAU-FDM page](imau_fdm_comparison.md). GEMB linearizes it with Newton into
the matrix diagonal (`Λ ≤ 0`, which strengthens diagonal dominance). Because the
residual is built from the *true* nonlinear flux at each iterate, the linearized
terms cancel at convergence — so the slope affects only the convergence rate, never
the answer.

Two properties were kept that a naive tridiagonal solve would lose:

- **Bit-level energy conservation.** The solved field is used as a predictor of the
  implicit face temperatures, then applied as one flux per face (`+F`/`−F` to the
  adjacent cells), exactly as the explicit path does. Pairwise cancellation is
  exact, so the column conserves independently of the solve residual, the Newton
  tolerance, and `c_p(T)`. No tolerance was widened.
- **The Dirichlet bottom cell is untouched by construction**, being excluded from
  the unknowns rather than solved and restored.

The expected speedup did **not** materialize, which is the honest result. Measured
over a year of 3-hourly forcing on a 264-cell column: explicit 2.44 s, implicit
5.75 s at the calibrated default. Removing the stability constraint cuts sub-steps
only ~2.4× (≈26 explicit vs 12 implicit for accuracy), and the 4.19 Newton
iterations per sub-step more than consume that. Profiling also shows why the
explicit path is hard to beat here: its inner loop runs at 1.03× a pure
memory-streaming floor, and the binding stability cell is not an outlier — the
second-tightest limit is only 1.08× the tightest, so there is no pathological cell
to rescue.

Both properties are the ones Fourteau et al. (2024) identify as load-bearing in the
finite-volume framing of this coupling — conservative flux application rather than a
Dirichlet surface temperature, and harmonic-mean face conductivity. Their taxonomy
also places GEMB and CFM in the same class, and leaves an explicit surface degree of
freedom as a live third option; see
[Surface energy balance numerics](index.md#Surface-energy-balance-numerics).

Where the implicit scheme does win is the case the explicit limit cannot bound: a
single 1e-4 m refrozen lens drops the explicit limit from 903 s to 3.98 s (227×)
while the implicit count does not move. That robustness, not throughput, is the
reason to have it.

The cheaper lever noted here originally still stands independently: the limit is
quadratic in `dz`, so raising `column_dzmin` from 0.025 → 0.05 m is ~4×.

### Strain softening

Oraschewski & Grinsted (2022), `strain.py`. GEMB's `apply_horizontal_strain!`
(`src/grid_ops.jl`) corresponds to CFM's `horizontal_divergence`, but **CFM
rescales mass where GEMB rescales thickness**, and CFM additionally softens the
densification rate itself. The mass-versus-thickness discrepancy is worth a
deliberate decision even though no code changed here.

### Multi-resolution subgrids

`regrid.py`, `regrid22`. GEMB's two-controller fixed-count grid
(`src/grid_ops.jl`) achieves the same goal differently, and its docstrings
document at length why the obvious alternatives are unstable.

### Analytic Herron–Langway spin-up

`hl_analytic.py`. GEMB's `steady_state_profile` march is more general — it works
for every densification scheme, not just Herron–Langway — and Fix 4 makes it
strictly more useful. Not a gap.

## Incidental differences, no action

- **`GRAVITY`**: CFM 9.8, GEMB 9.81.
- **Seconds per year**: CFM `31557600`; GEMB uses a deliberate
  `DAYS_PER_YEAR_DENSIFICATION = 365.0` (`src/constants.jl`). The resulting ~0.1%
  standing offset is already recorded in `README.md` deviation 4.
- **Surface grain size**: CFM uses Linow (2012) (`physics.py:1508`); GEMB uses a
  fixed `RE_NEW_SNOW`.
- **Ponding thresholds**: CFM's `RhoImp = 830` and `ThickImp = 0.1` match GEMB's
  `impermeable_density`/`impermeable_thickness` defaults exactly. Independent
  agreement, worth recording.

## References

- Stevens, C. M., Verjans, V., Lundin, J. M. D., Kahle, E. C., Horlings, A. N.,
  Horlings, B. I., and Waddington, E. D. (2020). The Community Firn Model (CFM),
  v1.0. *Geoscientific Model Development*, 13, 4355–4377.
- Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and Thomas, E. R.
  (2010). In situ measurements of Antarctic snow compaction compared with
  predictions of models. *Journal of Geophysical Research*, 115, F03011.
- Marbouty, D. (1980). An experimental study of temperature-gradient
  metamorphism. *Journal of Glaciology*, 26(94), 303–312.
- Vandecrux, B., et al. (2020). The firn meltwater Retention Model
  Intercomparison Project (RetMIP). *The Cryosphere*, 14, 3785–3810.
- Oraschewski, F. M., and Grinsted, A. (2022). Modeling enhanced firn
  densification due to strain softening. *The Cryosphere*, 16, 2683–2700.
- Fourteau, K., Brondex, J., Brun, F., and Dumont, M. (2024). A novel numerical
  implementation for the surface energy budget of melting snowpacks and glaciers.
  *Geoscientific Model Development*, 17, 1903–1929.
