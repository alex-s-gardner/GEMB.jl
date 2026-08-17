---
---

# Internals: Support {#Internals:-Support}

Non-exported constants, metadata tables, interpolation, and remaining helpers. See [Internals](/internals#Internals) for the caveat on depending on these.
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.AQUIFER_TOLERANCE' href='#GEMB.AQUIFER_TOLERANCE'><span class="jlbinding">GEMB.AQUIFER_TOLERANCE</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
AQUIFER_TOLERANCE
```


Diagnostic threshold, not a branch tolerance: how far above irreducible a cell must be for `aquifer_diagnostics` to call it saturated. Deliberately loose, and expressed as a fraction of the cell's pore space rather than as a mass.

The reason is physical, not arithmetic. `calculate_melt` leaves a retaining cell at exactly its irreducible water, but `calculate_density` then compacts the cell in the same timestep, shrinking the pore space that irreducible water was computed against. The cell is left marginally over-saturated — a genuine excess of order 1e-4 kg m-2 per timestep, which accumulates between melt events. A mass threshold cannot separate that from a thin aquifer, because the two differ in saturation, not in mass: compaction residue sits ~1e-7 of pore space above irreducible, while ponded water reaches 1e-1 to 1. Nothing in the physics reads this; it only decides what `aquifer_diagnostics` reports.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.DT_MIN_WARN' href='#GEMB.DT_MIN_WARN'><span class="jlbinding">GEMB.DT_MIN_WARN</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
DT_MIN_WARN
```


Smallest stable sub-timestep the thermal solver will accept without warning [s]. Below this, the Von Neumann stability limit has collapsed, which in practice means a near-zero `dz` cell rather than a genuinely stiff column.

`ExplicitThermal` only. `ImplicitThermal` is unconditionally stable, never consults `_max_safe_dt`, and so never warns however thin a cell becomes.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.HEAT_CAPACITY_WATER' href='#GEMB.HEAT_CAPACITY_WATER'><span class="jlbinding">GEMB.HEAT_CAPACITY_WATER</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
HEAT_CAPACITY_WATER
```


Specific heat capacity of liquid water [J kg-1 K-1]. Carries the sensible heat of _above-freezing_ liquid entering the column — i.e. rain — via the two-argument `specific_enthalpy_water(mp, T)`, under the `rain_heat_capacity = :water` default. Pore water in GEMB is carried strictly at the melting point, so this does not appear in the percolation or refreeze budgets; `heat_capacity(mp, T)` remains the ice-matrix value everywhere else. Value from the Community Firn Model (`solver.py`, `c_liq`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.THERMAL_BC_DERIVATIVE_STEP' href='#GEMB.THERMAL_BC_DERIVATIVE_STEP'><span class="jlbinding">GEMB.THERMAL_BC_DERIVATIVE_STEP</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
THERMAL_BC_DERIVATIVE_STEP
```


One-sided finite-difference step [K] for the turbulent-flux derivatives `dQ_shf/dT` and `dQ_lhf/dT`, which enter the Newton diagonal. `_turbulent_heat_flux` runs through the Beljaars-Holtslag stability branches, so differentiating it by hand is not worth the maintenance; one extra call per iteration is cheaper than the accuracy would be worth.

Measured, not assumed: the two flux evaluations together are 0.70 s of a 15.2 s run, so this call is not the cost of the implicit path. An analytic replacement that froze the stability coefficients was tried and reverted — it overshoots the true derivative up to 10x and more than doubles the iteration count. See `GEMB._surface_energy_balance_slope`.

1e-4 K rather than the usual `sqrt(eps)` ≈ 1.5e-8: the derivative only sets the convergence _rate_, never the answer (see `_thermal_solve!(::ImplicitThermal, ...)`), so a step large enough to stay well clear of cancellation in the flux difference beats a step chosen to minimize truncation error. It also matches `T_MELT_SWITCH_TOLERANCE`, so the difference does not straddle the melting-point latent-heat switch any more often than the melt switch does.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR' href='#GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR'><span class="jlbinding">GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
THERMAL_EXPLICIT_SAFETY_FACTOR
```


Fraction of `_max_safe_dt` the explicit sub-step may use. Default of `mp.thermal_explicit_safety_factor`; `ExplicitThermal` only.

This is a genuine margin, not a hedge against the diffusive limit being imprecise. That limit is exact for the stencil it is derived from (see `_max_safe_dt`), but it accounts for _diffusion only_. The surface cell also carries the explicitly-evaluated surface energy balance, whose linearization `Λ = dQ_sfc/dT₁ ≤ 0` (see `_surface_energy_balance_slope`) enters cell 1's own-temperature coefficient exactly as a face conductance does:

```julia
coef₁ = 1 − dt·(G₁ + |Λ|) / (ρ₁c₁dz₁)
```


`_max_safe_dt` omits `|Λ|`, so the true surface constraint is always stricter than the limit returned, by a factor that the model does not bound: `|Λ|` grows with wind speed (turbulent transfer) and with `T₁³` (longwave), while `ρ₁c₁dz₁` falls as the surface cell thins toward `column_dzmin`. Unlike the graded-grid error the face-based form fixed, this one is one-directional and unbounded, which is why a fixed fraction is the right instrument.

0.8 was inherited from the MATLAB model without a stated derivation. Measured over a year of 3-hourly synthetic forcing at `test_1` (2920 steps, spun-up column, `mp` defaults), the surface amplification factor `dt(G₁+|Λ|)/(ρ₁c₁dz₁)` actually realized was:

```julia
factor   median   p99     max     frac > 1   frac > 2   mean n_sub
0.8      0.196    0.609   0.955   0.0000     0.0        24.20
0.95     0.235    0.731   1.146   0.0065     0.0        20.20
1.0      0.235    0.731   1.146   0.0065     0.0        20.17
```


At 0.8 the surface coefficient stays non-negative on every step of that year; at 0.95 it goes negative on 0.65% of them. Negative-but-above-`−1` is damped ringing rather than divergence — outright blow-up needs a factor above 2, which this column never reaches — so 0.8 is not the difference between a run that works and one that does not _here_. It is the margin that keeps the scheme monotone, and the cost of holding it is ~17% more sub-steps on this column (the divisor grid is coarse, so the two often land on the same `dt` regardless).

Raising it toward 1.0 is therefore a defensible speed/monotonicity trade on a column resembling the one measured, not a free one, and 1.0 exactly is not safe in general: nothing in the model bounds `|Λ|/G₁`, so a thin, warm, windy surface cell can put the factor past 2. Runs that want the limit removed rather than tightened should use `ImplicitThermal`, which is unconditionally stable and does not consult this at all.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.THERMAL_IMPLICIT_DT_TARGET' href='#GEMB.THERMAL_IMPLICIT_DT_TARGET'><span class="jlbinding">GEMB.THERMAL_IMPLICIT_DT_TARGET</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
THERMAL_IMPLICIT_DT_TARGET
```


Target sub-step [s] for `ImplicitThermal`, and the ceiling on the resulting count.

Backward Euler is unconditionally stable, so this is an _accuracy_ control, not a stability one — the distinction that makes the implicit path cheap. The count is bounded by how fast the surface forcing changes and is independent of the cell count, of `dz`, and of the stiffest cell in the column, unlike the explicit path's ~29 sub-steps.

Calibrated on a year of 3-hourly synthetic forcing (264-cell column, no spinup), comparing whole-model output against a finely sub-stepped implicit run at 168.75 s. RMS temperature differences [K] and the annual melt difference:

target [s]   n_sub(3h)   RMS surface   RMS interior   Δmelt      runtime      5400          2          1.164         0.143       -2.46%      —      2700          4          0.526         0.140       -1.66%      —      1800          6          0.419         0.126       -1.19%     3.6 s       900         12          0.341         0.141       -1.03%     5.8 s       450         24          0.272         0.119       -0.08%     9.6 s

Runtimes re-measured after the static condensation of the interior rows (2.67x); the accuracy columns are unchanged by it, since it alters only how the same solve is factorized.

900 s is the knee: melt is within ~1% and the surface RMS has come down 3.4x from the two-sub-step case, while halving again costs 79% more runtime for 0.07 K.

The interior column is _not_ converging in this table — it floors at ~0.13 K regardless of target. That floor is not solver error: `manage_layer_thickness` merges and splits on discrete thickness thresholds, so a single flipped merge changes the deep discretization by more than the time integration does. Whole-model differencing therefore cannot resolve implicit accuracy below ~0.13 K, which is why the surface RMS and the melt total are what this is tuned on.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.energy_tolerance-Tuple{Real}' href='#GEMB.energy_tolerance-Tuple{Real}'><span class="jlbinding">GEMB.energy_tolerance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
energy_tolerance(E_reference)
```


Absolute tolerance [J] for a verbose-only energy-conservation check whose budget totals `E_reference`.

`E_TOLERANCE` alone is a _relative_ 1e-14 against a deep column's ~1e11 J of enthalpy, which only ever passed because both sides of the check were bit-identical arithmetic. Working in enthalpy adds real round-trip noise (~3e-6 J per deep-cell merge), so every check is relative-with-floor: the floor keeps small columns strict, the relative term scales with the budget being checked.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_height_attributes-Tuple{}' href='#GEMB.cf_height_attributes-Tuple{}'><span class="jlbinding">GEMB.cf_height_attributes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cf_height_attributes()
```


CF attributes for a genuine vertical coordinate, as produced by [`gemb_interp`](/api#GEMB.gemb_interp) when it regrids onto fixed heights. [`dz2z`](/api#GEMB.dz2z) returns heights measured negative-downward from the surface, hence `positive = "up"`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_layer_index_attributes-Tuple{}' href='#GEMB.cf_layer_index_attributes-Tuple{}'><span class="jlbinding">GEMB.cf_layer_index_attributes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cf_layer_index_attributes()
```


CF attributes for the `Z` dimension of GEMB output, which is a bare 1-based cell index rather than a coordinate. It therefore carries no `standard_name` and no `positive`: claiming either would present an index as a physical height. Depths come from the `dz` layer via [`dz2z`](/api#GEMB.dz2z).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_layermetadata-Tuple{NamedTuple}' href='#GEMB.cf_layermetadata-Tuple{NamedTuple}'><span class="jlbinding">GEMB.cf_layermetadata</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cf_layermetadata(layers::NamedTuple; time_axis=true)
```


Build the `layermetadata` NamedTuple for a `DimStack` whose layers are `layers`.

Mapping over `keys(layers)` rather than a fixed order is deliberate: DimensionalData requires the `layermetadata` keys to match the layer keys both in membership and in order, and the extracted-profile stack orders its layers differently from the `gemb` output stack.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cf_time_attributes-Tuple{}' href='#GEMB.cf_time_attributes-Tuple{}'><span class="jlbinding">GEMB.cf_time_attributes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cf_time_attributes()
```


CF attributes for the `Ti` dimension of GEMB output.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._interp1-Tuple{AbstractVector, AbstractVector, AbstractVector}' href='#GEMB._interp1-Tuple{AbstractVector, AbstractVector, AbstractVector}'><span class="jlbinding">GEMB._interp1</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_interp1(x, y, xi)
```


1D linear interpolation, **without** extrapolation: query points outside `[min(x), max(x)]` return `NaN` rather than a continuation of the end interval. `x` must be sorted (ascending or descending). See [`gemb_interp`](/api#GEMB.gemb_interp) for why the ends are missing rather than extrapolated.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._interp1_nearest-Tuple{AbstractVector, AbstractVector, AbstractVector}' href='#GEMB._interp1_nearest-Tuple{AbstractVector, AbstractVector, AbstractVector}'><span class="jlbinding">GEMB._interp1_nearest</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_interp1_nearest(x, y, xi)
```


1D nearest-neighbor interpolation, **without** extrapolation: query points outside `[min(x), max(x)]` return `NaN`, matching [`_interp1`](/internals_support#GEMB._interp1-Tuple{AbstractVector,%20AbstractVector,%20AbstractVector}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.air_density-Tuple{Real, Real}' href='#GEMB.air_density-Tuple{Real, Real}'><span class="jlbinding">GEMB.air_density</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
air_density(pressure, temperature) -> ρ_air [kg m-3]
```


Density of air from the ideal gas law, `ρ = 0.029·P/(R·T)`, with the 0.029 kg mol-1 mean molar mass of dry air the MATLAB implementation uses.

Shared by the transient surface energy balance ([`calculate_temperature`](/internals_physics#GEMB.calculate_temperature-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool})) and the initial-guess one (`_seb_annual_melt`), so both use one expression.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

