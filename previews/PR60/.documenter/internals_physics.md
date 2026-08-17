---
---

# Internals: Physics {#Internals:-Physics}

Non-exported physics kernels, called once per timestep from `gemb_core`. See [Internals](/internals#Internals) for the caveat on depending on these.
<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_accumulation-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters, Bool}' href='#GEMB.calculate_accumulation-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters, Bool}'><span class="jlbinding">GEMB.calculate_accumulation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_accumulation(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)
```


Add precipitation and deposition to the model column.

Precipitation is classified as snow or rain based on `mp.rain_temperature_threshold`. Snow is added as a new layer (if depth &gt; dzmin) or merged into the top cell. Rain is added by increasing the mass and temperature of the top grid cell, with temperature adjusted to account for latent heat of fusion.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, rain)`. Arrays grow by one cell (via [`open_slot!`](/internals_grid#GEMB.open_slot!-Tuple{NamedTuple,%20Int64})) when snow depth &gt; dzmin; the extra cell is reclaimed later in the timestep by [`manage_layer_thickness`](/internals_grid#GEMB.manage_layer_thickness-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters,%20Bool})'s count controller.

This is the model's principal age source. All three paths add mass at age 0: a fresh cell is set to 0 outright, and the two merge-into-cell-1 paths dilute `age[1]` by the arriving mass via [`dilute_age`](/internals_grid#GEMB.dilute_age-Tuple{Float64,%20Float64,%20Float64}).

Albedo is absent: it is diagnosed from the column at the top of the _next_ timestep (see [`calculate_albedo`](/internals_physics#GEMB.calculate_albedo-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ClimateForcingStep,%20ModelParameters})), which reads the fresh-snow grain size and density this function sets, so new snow brightens the surface without an albedo being stored here.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._albedo_gardner-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vararg{Float64, 4}}' href='#GEMB._albedo_gardner-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vararg{Float64, 4}}'><span class="jlbinding">GEMB._albedo_gardner</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_albedo_gardner(grain_radius, dz, density, c1, c2, SZA, t)
```


Broadband albedo parameterization from Gardner and Sharp (2010). Accounts for grain size, soot loading, solar zenith angle, and cloud optical thickness. Two-layer parameterization is applied when an ice layer exists below snow.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_albedo-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, ClimateForcingStep, ModelParameters}' href='#GEMB.calculate_albedo-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, ClimateForcingStep, ModelParameters}'><span class="jlbinding">GEMB.calculate_albedo</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_albedo(dz, density, water, grain_radius, melt_surface, cfs::ClimateForcingStep, mp::ModelParameters)
```


Calculate snow, firn, and ice albedo using one of several methods:
2. "GardnerSharp": function of effective grain radius (Gardner & Sharp, 2010)
  
3. "BrunLefebre": function of effective grain radius (Lefebre et al., 2003)
  
4. "GreuellKonzelmann": function of density and cloud amount (Greuell & Konzelmann, 1994)
  

Returns `(albedo, albedo_diffuse)` as scalars — the broadband surface albedo for the direct-beam and for the diffuse stream.

Every method here is **diagnostic**: the albedo is a function of the current column state alone, so nothing is carried across a timestep boundary and nothing is stored per cell. `albedo_diffuse` differs from `albedo` only under `:GardnerSharp` (recomputed at a fixed 50° effective zenith angle), which is also the only method whose consumer, [`calculate_shortwave_radiation`](/internals_physics#GEMB.calculate_shortwave_radiation-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20ClimateForcingStep,%20ModelParameters}), splits the direct and diffuse streams; the other methods return it equal to `albedo`.

**References**
- Gardner, A. S. and Sharp, M. J. (2010). J. Geophys. Res., 115, F01009.
  
- Lefebre, F., et al. (2003). J. Geophys. Res., 108, 4231.
  
- Greuell, W. and Konzelmann, T. (1994). Global Planet. Change, 9, 91-114.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.DensificationCoeffs' href='#GEMB.DensificationCoeffs'><span class="jlbinding">GEMB.DensificationCoeffs</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
DensificationCoeffs(mp::ModelParameters, pm, tam)
```


Loop-invariant densification terms for `mp.densification_method`, hoisted out of the per-cell / per-age-step loop. Shared by [`calculate_density`](/internals_physics#GEMB.calculate_density-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}) and [`steady_state_density`](/api#GEMB.steady_state_density) so both derive `c` from identical inputs.

`k0` and `k1` are the hoisted per-stage scalars (below and above `DENSITY_STAGE_TRANSITION`). Every scheme that reaches them needs exactly two, so their _meaning_ is per-method while their arity is not:

|                                                                   method |                            `k0`, `k1` |
| ------------------------------------------------------------------------:| -------------------------------------:|
| `:Arthern`, `:ArthernB`, `:Barnola1991`, `:CrocusPure`, `:HerronLangway` |                          unused (1.0) |
|                                                            `:Ligtenberg` | the accumulation-dependent `M0`, `M1` |
|                                                          `:Simonsen2013` |                 `F0 = 0.8` and `F1·γ` |
|                                                   `:GSFC2020`, `:Crocus` | `pm^α · g` with α0 = 0.91, α1 = 0.644 |


`pm` is mean annual accumulation [kg m-2 yr-1] and `tam` mean annual air temperature [K].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.BARNOLA_CLOSEOFF' href='#GEMB.BARNOLA_CLOSEOFF'><span class="jlbinding">GEMB.BARNOLA_CLOSEOFF</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
BARNOLA_CLOSEOFF
```


Density [kg m-3] at which `_barnola_f` switches from the fitted polynomial to the analytic closed-pore expression.

800, not GEMB's `DENSITY_PORE_CLOSEOFF = 830`: this is the density at which Barnola's _polynomial fit_ is superseded by the closed-pore geometry, which is a property of that fit rather than of pore close-off as GEMB's other code means it. The polynomial is fitted to Pimienta and Duval's data over roughly 0.55-0.8 g cm-3 and diverges from the analytic form beyond it — at ρ = 900 it gives 0.0432 against the analytic 0.0087, a factor of 5 — which is why the handover exists rather than running the polynomial to ρᵢ.

The handover is meant to be smooth. The paper states that the polynomial "was calculated in order to make the two functions, `fe(ρ)` and `fs(ρ)`, and their first derivatives equal for ρ = 0.8 g cm-3" — so both value and slope should match there, and any step is a mismatch rather than an intended feature of the law.

The size of that step depends on `mp.density_ice`, because only the closed-pore branch carries it. The two branches cross at ρᵢ = 919.96, and their first derivatives cross at 920.06 — the paper's C¹ statement is what fixes both, and it identifies the fit's own ice density as ~920 kg m-3. The branches agree to 0.06% at 920 against 4.3% at 917 and 14.1% at GEMB's default 910, so the gap at any other `density_ice` measures the mismatch between the configured value and the fit's. At the default 910 the rate steps _down_ by 14% crossing 800, which is a discontinuity in `dρ/dt` but not in ρ, and is small against the factor-5 error the alternative would introduce deeper in the column.

The sign of that step is not fixed: it is downward below the fit's ~920 and upward above it (+7.6% at 925, +49.5% at 950). See [`BARNOLA_MIN_DENSITY_ICE`](/internals_physics#GEMB.BARNOLA_MIN_DENSITY_ICE) for why the low side is gated and the high side is not.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.BARNOLA_MIN_DENSITY_ICE' href='#GEMB.BARNOLA_MIN_DENSITY_ICE'><span class="jlbinding">GEMB.BARNOLA_MIN_DENSITY_ICE</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
BARNOLA_MIN_DENSITY_ICE
```


Lowest `mp.density_ice` [kg m-3] that `:Barnola1991` accepts, enforced in `validate_parameters`. The rest of the model permits `[800, 950]`.

The gate exists because `BARNOLA_CLOSEOFF` is an absolute density while the polynomial branch below it carries no `density_ice`. As the configured ice density falls toward 800 the handover moves into firn whose true porosity is nearly zero, and the polynomial — which cannot be rescaled, its coefficients being fitted against the fit's own ρᵢ ≈ 920 — overstates `f` there: by 5.9x at ρ = 700 and 15.9x at ρ = 799.9 for `density_ice = 820`. At exactly 800 the polynomial covers the entire column, `f` never reaches 0, and the scheme loses the self-limiting behaviour that motivates it, leaving only the ice-density clamp.

900 admits GEMB's default 910 (a -14% step at the handover) and the common 917 (-4%) while refusing the range where the law degrades: the step is -27% at 900 and -94% at 820.

Deliberately one-sided. Above ~920 the step inverts rather than growing without bound (+7.6% at 925, +49.5% at 950) and the polynomial stays inside its fitted porosity range, so those configurations are merely inconsistent with the fit rather than unphysical. They are documented on [`BARNOLA_CLOSEOFF`](/internals_physics#GEMB.BARNOLA_CLOSEOFF) instead of rejected.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.BARNOLA_N' href='#GEMB.BARNOLA_N'><span class="jlbinding">GEMB.BARNOLA_N</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
BARNOLA_N
```


Stress exponent in `dρ/dt = ρ·A0·exp(-Q/RT)·f(ρ)·σⁿ`.

Fixed at 3, matching Barnola et al. (1991) for the whole firn range this scheme covers: "_n_ taken to be equal to 3, as the effective stress in the firn is rapidly higher than 0.1 MPa". Their `A0` was fitted against `n = 3` over 0.55-0.8 g cm-3, so the exponent and the prefactor are one calibration.

The paper's adjacent remark that "the exponent _n_ is 1 when the effective stress is lower than 0.1 MPa" is _not_ an unimplemented firn branch, though it reads like one. It scopes itself to bulk ice: "**below the close-off**, a good fit to the experimental data is obtained by taking successively _n_ = 3 and _n_ = 1 (Pimienta, 1987)", citing Doake and Wolff (1985) and Pimienta and Duval (1987). Those are shear-creep studies of dense polar ice — Pimienta and Duval torsion-tested 0.85 g cm-3 samples and analysed Dye 3 inclinometry at 1000-1784 m, and report no densification rate at all — so they supply no volumetric prefactor for porous firn, and none is published here to pair with `n = 1`.

Three reasons not to add the branch, beyond the missing prefactor:
- Pimienta and Duval's own exponent is not 1. Their abstract and conclusion both give "smaller than 2"; regressing their Table 1 (South Pole ice, -15 °C) gives n ≈ 1.55. Barnola's "_n_ is 1" rounds a sub-2 shear result.
  
- Matching an `n = 1` branch continuously at 0.1 MPa forces `A1 = A0·(1e5)²`, so the rate would rise by `(1e5/σ)²` below it — ×2 at 70 kPa, ×41 at 15 kPa, ×100 at 10 kPa.
  
- Barnola reproduced observed 0.55-0.83 g cm-3 profiles at Vostok and other Antarctic and Greenland sites with `n = 3` throughout that range. A 40× shallow acceleration would break the agreement the coefficients exist to produce.
  

The Community Firn Model reaches the same conclusion, carrying the switch as dead code (`# nBa[sigmaEff<1.0e5]=1.0`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.CROCUS_HYBRID_DENSITY' href='#GEMB.CROCUS_HYBRID_DENSITY'><span class="jlbinding">GEMB.CROCUS_HYBRID_DENSITY</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



```julia
CROCUS_HYBRID_DENSITY
```


Density [kg m-3] at which `:Crocus` hands a cell over to `:GSFC2020`.

Eq. 7 is fitted to a 1-2 m alpine snowpack, and `exp(bη·ρ)` reaches ~1e9 at ρ = 900: on a Greenland firn profile the pure law gives 0.7-0.8 of `:Arthern`'s compaction rate in the top few metres but only 0.02-0.09 below 20 m. That is the published law behaving as fitted, not a unit error, so `:Crocus` is a composite — Crocus viscosity where liquid water and grain type govern settling, `:GSFC2020` where creep does.

450 kg m-3 is the threshold the Community Firn Model's `Crocus` scheme uses for the same purpose (its `RHO_CC`), and it sits between the two laws' calibration ranges rather than inside either. The blend is a discontinuous handover, not a weighted average: nothing in either paper prescribes a blending function. `:CrocusPure` applies eq. 5 at every density.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._barnola_f-Tuple{Float64, Float64}' href='#GEMB._barnola_f-Tuple{Float64, Float64}'><span class="jlbinding">GEMB._barnola_f</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_barnola_f(ρ, density_ice) -> f [-]
```


Densification factor of Barnola et al. (1991), the geometric term in their eq. 2,

```julia
dρ/dt = ρ · A0 · exp(-Q/RT) · f(ρ) · ΔPⁿ
```


Two regimes, meeting at [`BARNOLA_CLOSEOFF`](/internals_physics#GEMB.BARNOLA_CLOSEOFF):
- Open pores, `ρ <= 800`: the paper's eq. 3, `fe(ρ)`, an empirical fit `log10 f = α(ρ/1000)³ + β(ρ/1000)² + δ(ρ/1000) + γ` "empirically deduced for the 0.55-0.8 g cm-3 density range", parameterizing the Pimienta and Duval (1987) pressure-sintering data. Note the polynomial is in g cm-3, so the coefficients are only meaningful against `ρ/1000`.
  
- Closed pores, `ρ > 800`: the paper's eq. 4, `fs(ρ)`, the Wilkinson and Ashby (1975) spherical pore model, `f = (3/16)·(1 − ρ/ρᵢ) / (1 − (1 − ρ/ρᵢ)^(1/3))³`, which the paper describes as "already valid below ρ = 0.8 g cm-3".
  

`ρᵢ` is `mp.density_ice`, so the closed-pore branch follows the configured ice density. The polynomial branch cannot: its coefficients were fitted with a fixed ice density (~920 kg m-3; see [`BARNOLA_CLOSEOFF`](/internals_physics#GEMB.BARNOLA_CLOSEOFF)), and rescaling them is not something the paper licenses.

Below `DENSITY_STAGE_TRANSITION` the scheme does not use this factor at all — it uses Herron & Langway stage 1, matching the paper ("from the surface to ρ = 0.55 g cm-3 the Herron and Langway (1980) model was used"), and that is exactly GEMB's `_hl_c0`. The identity was verified against the Community Firn Model's `Barnola1991` zone 1 to a ratio of 1.0 (12 significant figures) across three (T, ρ, accumulation) triples, so `:Barnola1991` and `:HerronLangway` share the one kernel rather than restating it.

Cross-checked against the Community Firn Model's `Barnola1991`, which agrees on both branches and on all four polynomial coefficients. CFM applies the full overburden `σ` where this uses the cell-midpoint value, matching the convention `:Crocus` already uses here; and CFM's zone 1 carries an `A^aHL` accumulation term with `aHL = 1`, which is the linear accumulation dependence already inside `_hl_c0`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._crocus_viscosity-NTuple{5, Float64}' href='#GEMB._crocus_viscosity-NTuple{5, Float64}'><span class="jlbinding">GEMB._crocus_viscosity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_crocus_viscosity(ρ, T, W_liq, D, grain_radius) -> η [kg s-1 m-1]
```


Snow viscosity of Vionnet et al. (2012) eq. 7, with the two microstructure correction factors of eqs. 8 and 9:
- `f1 = 1/(1 + 60·W_liq/(ρ_w·D))` (eq. 8) reduces viscosity in the presence of liquid water, bottoming out around 1/61 for a saturated cell. This is the term GEMB has no other expression of: every accumulation-driven scheme here is blind to `water`, and `:ArthernB` reads it only as overburden.
  
- `f2 = min(4, exp(min(g1, gs − g2)/g3))` (eq. 9) _raises_ viscosity for coarse angular grains, suppressing compaction in depth hoar.
  

`W_liq` is cell water [kg m-2], `D` cell thickness [m], `grain_radius` [mm]; `gs` in eq. 9 is a grain size, taken here as the grain diameter `2·grain_radius` in metres.

`f2` is floored at 1, which the published eq. 9 does not do. Two facts force the choice:
2. Eq. 9 is bounded above (`min(4.0, ...)`) but not below. It crosses 1 at `gs = g2 = 0.2 mm` and decays as `exp((gs − g2)/g3)` beneath that, so a _fine_-grained cell gets `f2 < 1` — a softening, from a factor the paper introduces to "account for [...] the increase of viscosity with angular grains". At GEMB's fresh-snow radius (`RE_NEW_SNOW = 0.05 mm`, so gs = 0.1 mm) it is `exp(-1) = 0.368`, i.e. new snow compacts 2.7× faster than the unmodified eq. 7 rate.
  
3. In Crocus that regime is unreachable, so the paper never has to bound it. Eq. 9 is a _non-dendritic_ relation, and Crocus tracks dendricity explicitly: fresh snow is dendritic (Sect. 3.3) and is described by `d` and `s`, not by `gs`. A layer only acquires a `gs` once `d` reaches 0, and the paper puts non-dendritic `gs` in the 0.3-0.4 mm range — where eq. 9 gives 2.7 to 4 (its cap binds at gs = 0.339 mm). Every `gs` the paper's own model can present to eq. 9 therefore yields `f2 >= 2.7`.
  

GEMB has no such gate: it carries one `grain_radius` for dendritic and non-dendritic snow alike, and initializes it at 0.05 mm — squarely inside the interval eq. 9 leaves undefined. Passing that radius through unfloored would apply a 1.6-2.7× softening to the fresh snow at the top of the column, inverting the factor's intent for the cells it was never meant to score. Flooring at 1 makes `f2` a pure stiffening correction, which is what eq. 9 does across its whole valid domain, and leaves it inert (`f2 = 1`) below it rather than extrapolating a relation off its domain.

`grain_dendricity` is available in the column state, so a closer port could gate eq. 9 on `d == 0` as Crocus does. That would need a defensible `f2` for dendritic snow — the paper supplies none, since the question does not arise there — so the floor is the conservative reading, not the only one.

Deviates from the Community Firn Model's `Crocus`, which was used to cross-check eq. 7 and the `dρ/dt = ρσ/η` form: CFM hardcodes `f2 = 4.0` (eq. 9's ceiling, so it stiffens every cell including fresh fine-grained snow) and adopts van Kampenhout et al. (2017) eq. 8's retuned `cη = 358`. The paper's `f2` and `cη = 250` are used here.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._densification_rate-Tuple{GEMB.DensificationCoeffs, Float64, Float64}' href='#GEMB._densification_rate-Tuple{GEMB.DensificationCoeffs, Float64, Float64}'><span class="jlbinding">GEMB._densification_rate</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_densification_rate(p::DensificationCoeffs, ρ, T) -> c [yr-1]
```


Densification rate coefficient for the method in `p`, dispatching on `p.method`. Used by [`steady_state_density`](/api#GEMB.steady_state_density), where the branch is taken once per age step rather than once per cell. `calculate_density`'s hot loop calls the per-method functions directly, so it takes the branch once per call.

The stress-driven schemes cannot be evaluated here: the march carries neither overburden stress, grain radius, nor water, and it only produces an initial guess that the spinup then relaxes. Each therefore falls back to the nearest accumulation-driven law — `:Crocus` to `:GSFC2020`, which is already its own above-threshold branch and so governs most of the column the march builds; `:Barnola1991` to `:HerronLangway`, whose stage 1 _is_ its own below-threshold branch, so the two agree exactly below `DENSITY_STAGE_TRANSITION` and the fallback is only an approximation above it; `:ArthernB` and `:CrocusPure` to `:Arthern`. The transient run uses the real law in every case.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._densify_cell!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, Float64, Float64, Float64}' href='#GEMB._densify_cell!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, Float64, Float64, Float64}'><span class="jlbinding">GEMB._densify_cell!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_densify_cell!(density, dz_in, dz_out, i, c, dt, density_ice)
```


Apply the densification increment for cell `i` given rate coefficient `c`: update `density[i]` in place, clamp it to the density of ice, and write the mass-conserving new grid-cell length to `dz_out[i]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_density-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters}' href='#GEMB.calculate_density-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters}'><span class="jlbinding">GEMB.calculate_density</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_density(temperature, dz, density, grain_radius, water, cfs::ClimateForcingStep, mp::ModelParameters)
```


Compute the densification of snow/firn using one of several models:
- `:HerronLangway`: Herron and Langway (1980)
  
- `:Arthern`: semi-empirical model of Arthern et al. (2010) [default]
  
- `:ArthernB`: physical model from Appendix B of Arthern et al. (2010)
  
- `:Crocus`: viscous settling of Vionnet et al. (2012) below `CROCUS_HYBRID_DENSITY`, `:GSFC2020` above; the only scheme here in which liquid water weakens the matrix
  
- `:CrocusPure`: Vionnet et al. (2012) at every density, without the firn handover
  
- `:GSFC2020`: recalibrated Arthern of Medley et al. (2022), GSFC-FDM v1.2.1
  
- `:Simonsen2013`: Arthern retuned for Greenland by Simonsen et al. (2013)
  
- `:Barnola1991`: Herron-Langway stage 1 below `DENSITY_STAGE_TRANSITION`, pressure sintering above; the only scheme here that models the firn-ice transition mechanistically
  
- `:LiZwally`: empirical model of Li and Zwally (2004)
  
- `:Helsen`: modified empirical model by Helsen et al. (2008)
  
- `:Ligtenberg`: semi-empirical model of Ligtenberg et al. (2011)
  

Returns `(dz, density)`. `density` is updated in place; `dz` is returned as a new array (recomputed from the conserved cell mass).

Pass a length-`m` `viscosity` vector to have the per-cell effective viscosity [Pa s] written into it. It is filled with `NaN` first and written only by the `:Crocus`/`:CrocusPure` settling branch, the one place a viscosity is genuinely formed; under every other scheme, and for the cells `:Crocus` hands to `:GSFC2020`, it stays `NaN`. Omitting the keyword leaves the hot path untouched: the default is the empty `NO_VISCOSITY`, and an empty vector means "not requested" (a concrete type on both settings — see the note in `gemb_core.jl`).

`water` [kg m-2] is read only by the stress-driven schemes: `:ArthernB` and `:Barnola1991` (as part of the overburden) and `:Crocus`/`:CrocusPure` (overburden, plus eq. 8's viscosity reduction); the accumulation-driven schemes proxy overburden with mean accumulation and ignore it.

The accumulation-driven branches are scalar-loop implementations that are numerically identical, element by element, to the reference vectorized MATLAB translation, but avoid the mask / gather / broadcast temporaries (`mass_init`, `idx`, `H`, `c0`, `c1`, ...) the vectorized form allocated per call. `:ArthernB` deviates from the reference — see the comment on that branch.

**References**
- Arthern, R. J., et al. (2010). J. Geophys. Res., 115, F03011.
  
- Herron, M. and Langway, C. (1980). J. Glaciol., 25, 373-385.
  
- Li, J. and Zwally, H. (2004). Ann. Glaciol., 38, 309-313.
  
- Helsen, M. M., et al. (2008). Science, 320, 1626-1629.
  
- Ligtenberg, S. R. M., et al. (2011). The Cryosphere, 5, 809-819.
  
- Medley, B., et al. (2022). The Cryosphere, 16, 3971-4011.
  
- Simonsen, S. B., et al. (2013). J. Glaciol., 59, 545-558.
  
- Vionnet, V., et al. (2012). Geosci. Model Dev., 5, 773-791.
  
- Lundin, J. M. D., et al. (2017). J. Glaciol., 63, 401-422. (FirnMICE; eqs. A36-A37)
  
- Pimienta, P. and Duval, P. (1987). J. Phys. Colloques, 48, C1-243-C1-248.
  
- Barnola, J.-M., Pimienta, P., Raynaud, D., and Korotkevich, Y. S. (1991). Tellus, 43B, 83-90.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._grain_gradient-Tuple{Vector{Float64}, Vector{Float64}, Int64, Int64}' href='#GEMB._grain_gradient-Tuple{Vector{Float64}, Vector{Float64}, Int64, Int64}'><span class="jlbinding">GEMB._grain_gradient</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_grain_gradient(temperature, dz, i, m)
```


Absolute temperature gradient [degC m-1] at cell `i` of a column of `m` cells, using the same finite-difference stencil as the reference implementation. The grid-point-center separations reduce to local `dz` combinations, so no cumulative depth array is needed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._grain_lwc-Tuple{Float64, Float64, Float64}' href='#GEMB._grain_lwc-Tuple{Float64, Float64, Float64}'><span class="jlbinding">GEMB._grain_lwc</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_grain_lwc(water, density, dz)
```


Liquid-water content as a mass fraction [%] for a single cell, capped at 9% (Brun, 1980).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._marbouty_Q-Tuple{Float64, Float64, Float64}' href='#GEMB._marbouty_Q-Tuple{Float64, Float64, Float64}'><span class="jlbinding">GEMB._marbouty_Q</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_marbouty_Q(temperature, density, dT)
```


Grain-growth rate coefficient Q [mm d-1] for a single cell per Fig. 9 of Marbouty (1980). No grain growth above `DENSITY_MARBOUTY_MAX` (H = 0). Scalar equivalent of the reference `_Marbouty`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_grain_size-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters}' href='#GEMB.calculate_grain_size-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters}'><span class="jlbinding">GEMB.calculate_grain_size</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_grain_size(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, cfs::ClimateForcingStep, mp::ModelParameters)
```


Model the evolution of effective snow grain size, dendricity, and sphericity.

Accounts for different physical processes depending on snow state:
- Dendritic Snow (fresh): Evolves based on temperature gradients and liquid water content.
  
- Nondendritic Dry Snow: Temperature gradient metamorphism using Marbouty (1980), or normal grain growth after Arthern et al. (2010), selected by `mp.grain_growth_method`.
  
- Wet Snow: Rapid grain growth due to liquid water using Brun (1989).
  

`mp.grain_growth_method` governs the non-dendritic _dry_ branch only; the dendritic and wet branches, and the sphericity caps, are the same on every setting.
- `:Marbouty` (default, MATLAB's behaviour) — Marbouty (1980) everywhere. Its density factor `H` vanishes at `DENSITY_MARBOUTY_MAX`, so grain radius is frozen throughout the firn column below a few metres.
  
- `:Arthern` — `dr²/dt = kgr·exp(-Eg/RT)` (`GRAIN_GROWTH_KGR`, `GRAIN_GROWTH_EG`) at every density, dropping the temperature-gradient dependence entirely.
  
- `:hybrid` — Marbouty below `DENSITY_MARBOUTY_MAX`, Arthern at or above it. Seasonal snow keeps the temperature-gradient physics; firn grains keep growing, which matters because `:ArthernB` densification goes as `1/r²`.
  

Runs unconditionally. Metamorphism is not contingent on model configuration — grains coarsen in a real snowpack whichever albedo scheme is selected — and `grain_radius` is read by four schemes across three modules ([`calculate_albedo`](/internals_physics#GEMB.calculate_albedo-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ClimateForcingStep,%20ModelParameters}), [`calculate_shortwave_radiation`](/internals_physics#GEMB.calculate_shortwave_radiation-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20ClimateForcingStep,%20ModelParameters}), `:ArthernB`/`:Crocus`/`:CrocusPure` densification in [`calculate_density`](/internals_physics#GEMB.calculate_density-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}), and the grain-radius emissivity methods in [`calculate_temperature`](/internals_physics#GEMB.calculate_temperature-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool})). MATLAB skips this function unless `albedo_method` is `:GardnerSharp` or `:BrunLefebre`, which covers only the first two, so the other configurations silently ran with grain size frozen at the initial profile. Enumerating the consumers to skip the work costs 11% on the runs that can skip it and re-breaks silently whenever a new consumer is added, so the work is simply always done.

Returns `(grain_radius, grain_dendricity, grain_sphericity)`. **All three are updated in place**, and the returned `grain_radius` is the same array that was passed in — pass a copy if the incoming values are still needed. `grain_radius` used to be returned as a fresh array; it now shares the mutation convention of the other two, which removed the largest allocation in the function (see the `gsz` comment below).

This is a scalar-loop implementation that is numerically identical, element by element, to the reference vectorized MATLAB translation, but avoids the ~30 mask / gather / broadcast temporaries the vectorized form allocated per call.

**References**
- Brun, E., et al. (1992). J. Glaciol., 38, 13-22.
  
- Marbouty, D. (1980). J. Glaciol., 26, 303-312.
  
- Brun, E. (1989). Ann. Glaciol., 13, 22-26.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.aquifer_diagnostics-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}' href='#GEMB.aquifer_diagnostics-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}'><span class="jlbinding">GEMB.aquifer_diagnostics</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
aquifer_diagnostics(dz, density, water, mp) -> (thickness, depth)
```


Summarize the standing (super-irreducible) water in a column.

`thickness` [m] is the total thickness of cells holding more water than capillary forces can retain — the saturated thickness a borehole or a firn core would find wet. `depth` [m] is the depth to the top of the shallowest such cell, `NaN` when there are none, so "dry column" stays distinguishable from "water table at the surface".

Both are zero and `NaN` under `mp.runoff_method === :instantaneous`, where no cell can hold more than its irreducible water. They are what makes the ponding in [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) observable in output, and the quantity RetMIP (Vandecrux et al., 2020, Sect. 5.4) evaluates aquifer representation against.

A cell counts as saturated only once it exceeds irreducible by [`AQUIFER_TOLERANCE`](/internals_support#GEMB.AQUIFER_TOLERANCE) of its pore space, rather than by the arithmetic `WATER_TOLERANCE` the physics uses. `calculate_melt` leaves a retaining cell at exactly irreducible, but `calculate_density` then compacts it within the same timestep, shrinking the pore space that retention was computed against and leaving the cell genuinely — if negligibly — over-saturated. That residue accumulates between melt events and would otherwise read as a metres-thick water table in a run that cannot have one. It is separable from real ponding by saturation but not by mass: the residue sits ~1e-7 of pore space above irreducible against 1e-1 or more for ponded water.

Diagnostic only: nothing in the physics reads either value. Called from [`gemb_core`](/internals_grid#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool}) after the grid controllers, alongside [`ice_slab_diagnostics`](/internals_physics#GEMB.ice_slab_diagnostics-Tuple{Vector{Float64},%20Vector{Float64},%20ModelParameters}), so the column scanned is the one the profile output records at that timestamp.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_melt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, ModelParameters, Bool}' href='#GEMB.calculate_melt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, ModelParameters, Bool}'><span class="jlbinding">GEMB.calculate_melt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_melt(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, rain, mp::ModelParameters, verbose::Bool)
```


Compute meltwater generation, percolation, refreezing, and runoff using a tipping bucket approach.

Processes:
2. Initial Refreeze: Existing pore water in cold layers is refrozen.
  
3. Melt Generation: Excess energy above 0 C is converted to liquid meltwater.
  
4. Percolation: Liquid water percolates downward, refreezing in cold layers, being retained as pore water, or running off at impermeable ice lenses.
  

Water is routed to runoff at a contiguous run of cells at or above `mp.impermeable_density` thicker than `mp.impermeable_thickness` (defaults 830 kg m-3 and 0.1 m, as in MATLAB; across the RetMIP models the density criterion spans 810–917 kg m-3 — see `ModelParameters`).

**Ponding above a barrier**

What happens to water that reaches a barrier depends on `mp.runoff_method`:
- `:instantaneous` (the default, and MATLAB's behaviour) — it leaves the column in the same timestep. Every cell is capped at its irreducible saturation, so the column can never hold standing water.
  
- `:ZuoOerlemans` / `:Darcy` — it _backs up_. After the downward percolation loop, an upward pass fills the cells above each barrier to their pore capacity `mp.pore_saturation_max·(ρᵢ − ρ)·dz` [kg m-2], deepest first; only what is left after cell 1 is full genuinely leaves as runoff, the column then being saturated to the surface. Water in excess of irreducible instead drains laterally, over a finite timescale, in [`apply_lateral_drainage!`](/internals_physics#GEMB.apply_lateral_drainage!-Tuple{NamedTuple,%20Float64,%20ModelParameters}). Barrier cells are excluded from the fill: at or above `mp.impermeable_density` there is no connected pore space to pond in, which is the premise of their being impermeable in the first place.
  The capacity uses the same `(ρᵢ − ρ)·S·dz` form as [`irreducible_saturation`](/internals_physics#GEMB.irreducible_saturation-Tuple{ModelParameters,%20Float64}) with `S = mp.pore_saturation_max`, so `water/capacity` is exactly the saturation `S_wi` is defined against and a cell at `S_wi = 1` reads as full rather than over-full.
  The upward pass does not refreeze what it deposits, even into a cell that is still cold (a cell can pass water down while retaining cold content if its refreeze was capped by `d_max` rather than by cold content). Any such cell is resolved at the top of the _next_ timestep by the REFREEZE PORE WATER block above, so the only cost is a one-step delay; mass and energy are conserved either way, because pore water and runoff carry the same melt-point enthalpy.
  Firn aquifers then form bottom-up with no further machinery: water still drains to irreducible everywhere it passes, piles up at the deepest barrier, and backs up from there. RetMIP (Vandecrux et al., 2020, Sect. 5.4) excluded the firn-aquifer site from its retention evaluation for exactly this reason — models that cannot hold water above irreducible saturation cannot represent an aquifer at all.
  This deviates from MATLAB, which has only the `:instantaneous` behaviour.
  

**No preferential-flow domain**

This is a single-domain bucket scheme: water descends cell by cell through the matrix, and there is no second, fast, heterogeneous ("piping") flow path, nor a Richards-equation matrix solver. That is deliberate. RetMIP (Vandecrux et al., 2020) found that the three of nine models carrying an explicit deep- or preferential-percolation scheme (CFM-Cr, CFM-KM, UppsalaUniDeepPerc) did worse than the bucket schemes at three of the four Greenland sites: they infiltrated too deeply and ran warm — firn-temperature mean error +3.6 to +6.2 °C at Dye-2 and +1.8 to +4.7 °C at KAN_U, plus a warm bias at the near-melt-free Summit site — and at Dye-2 in 2016 they percolated to 10 m against the 2.5 m the upward-looking radar observed, growing near-surface ice slabs several metres thick where none exist. Their one advantage was the firn-aquifer site, which only they recharged. RetMIP's conclusion (their Sect. 5.2) is that until preferential flow in firn is better constrained observationally, the more complex schemes do not necessarily give better results than simple bucket schemes; the productive levers for a bucket scheme are instead the impermeability criterion and the runoff timescale above — the two best-performing models at the ice-slab site KAN_U were bucket schemes that delay runoff rather than models that percolate deeper.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, melt_total, melt_surface, runoff_total, freeze_total, percolation_depth)`.

**Age transport**

`age` [d] is the mass-weighted mean age of all mass in a cell, matrix and pore water together, so a phase change that stays inside one cell (the initial pore-water refreeze) does not move it at all. Melt removal is likewise age-neutral: it takes a _fraction_ of the cell, which leaves the mean age of the remainder unchanged.

What does move age is percolation. Meltwater carries the age of the firn it came from, so the percolation loop transports the age _moment_ `mass × age` alongside `flux_dn`, in a companion `flux_age` vector holding the mass-weighted mean age of the water in transit. At each cell the water generated there and any water the cell expels (both at `age[i]`) mix with the through-flow from above (at `flux_age[i]`) into one pool; whatever refreezes enters the cell at the pool's age, and whatever continues down or runs off leaves at it. Refreezing therefore makes a cell _younger_ only to the extent that genuinely younger water reached it, rather than resetting it.

`percolation_depth` [m] is a diagnostic that takes no part in the mass or energy budget: the base of the deepest cell water reached this timestep, 0 when no water moved. Note that the percolation loop can walk past cells water never entered (its exit also waits on the deepest cell holding pre-existing excess pore water), so this tracks water arrival per cell rather than the loop index. The companion slab diagnostics are taken in [`gemb_core`](/internals_grid#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool}) via [`ice_slab_diagnostics`](/internals_physics#GEMB.ice_slab_diagnostics-Tuple{Vector{Float64},%20Vector{Float64},%20ModelParameters}), after the grid controllers run.

Arrays may shrink (cells deleted when mass=0).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.ice_slab_diagnostics-Tuple{Vector{Float64}, Vector{Float64}, ModelParameters}' href='#GEMB.ice_slab_diagnostics-Tuple{Vector{Float64}, Vector{Float64}, ModelParameters}'><span class="jlbinding">GEMB.ice_slab_diagnostics</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ice_slab_diagnostics(dz, density, mp::ModelParameters) -> (thickness, depth)
```


Summarize the ice slabs in a column under the same criterion the percolation scheme uses.

`thickness` [m] is the total thickness of all cells at or above `mp.impermeable_density`, whether or not they form a flow-blocking run — the quantity a firn core measures.

`depth` [m] is the depth to the top of the shallowest run that would actually block percolation, matching both clauses of the impermeable branch of [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}): a contiguous run thicker than `mp.impermeable_thickness`, or any single cell at `mp.density_ice`, which blocks unconditionally however thin it is. It is `NaN` when nothing qualifies, so "no slab" stays distinguishable from "slab at the surface".

Both thresholds are capped at `mp.density_ice`, since `mp.impermeable_density` may be set above it (see [`ice_slab_diagnostics`](/internals_physics#GEMB.ice_slab_diagnostics-Tuple{Vector{Float64},%20Vector{Float64},%20ModelParameters}) source) and solid ice must always count.

Diagnostic only: nothing in the physics reads either value. Called from [`gemb_core`](/internals_grid#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool}) after the grid controllers, not from [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}), so the column scanned is the one the profile output records at that timestamp.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.irreducible_saturation-Tuple{ModelParameters, Float64}' href='#GEMB.irreducible_saturation-Tuple{ModelParameters, Float64}'><span class="jlbinding">GEMB.irreducible_saturation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
irreducible_saturation(mp::ModelParameters, density) -> S_wi [-]
```


Irreducible (capillary-held) water saturation of the pore space, per `mp.water_irreducible_method`. The retention of a cell is `(ρᵢ − ρ) · S_wi · dz` [kg m-2] at every site that uses it, so only `S_wi` varies between methods.
- `:constant` — `mp.water_irreducible_saturation` at every density (Colbeck, 1974).
  
- `:ColeouLesaffre` — the default. Coléou and Lesaffre (1998) eq. 3 via Langen et al. (2017) eq. 4, where irreducible water mass per unit total mass is `wmi = 0.057(ρᵢ − ρ)/ρ + 0.017` and
  
  ```
  S_wi = wmi/(1 − wmi) · ρᵢ·ρ / (ρ_w(ρᵢ − ρ))
  ```
  
  `mp.water_irreducible_saturation` is ignored. Retention rises with density — ~0.069 at ρ = 300 against ~0.163 at ρ = 800 (ρᵢ = 917) — which is where the constant value under-retains most, in the percolation zone.
  

`:ColeouLesaffre` is gated at `min(DENSITY_PORE_CLOSEOFF, mp.density_ice)`, where there is no connected pore space left to hold water, returning zero there as the Community Firn Model and [`_irreducible_water`](/internals_grid#GEMB._irreducible_water-Tuple{Real,%20Real,%20Real,%20ModelParameters}) do. `:constant` is deliberately _not_ gated, so it stays a flat saturation at every density. The gate also removes the `ρ → ρᵢ` singularity in `S_wi`, which is why it takes the `min`: with `mp.density_ice` configured below pore closeoff the singularity would otherwise sit below the constant threshold.

The gate is the _constant_ `DENSITY_PORE_CLOSEOFF`, deliberately not the tunable `mp.impermeable_density`. The two thresholds answer different questions: this one is where capillary retention ceases for want of connected pore space, while `impermeable_density` is where a lens stops conducting flow at the scale the model represents. Tying them together would make lowering the flow criterion silently change retention too.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.pond_blocked_water!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, ModelParameters}' href='#GEMB.pond_blocked_water!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, ModelParameters}'><span class="jlbinding">GEMB.pond_blocked_water!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pond_blocked_water!(M, density, water, water_delta, age, runoff, flux_dn, flux_age, Xi, mp)
    -> runoff_blocked [kg m-2]
```


Back water that percolation could not pass — the runoff each cell generated, plus the basal outflow `flux_dn[Xi]` — up into the pore space above, and return only the mass that still leaves the column.

Called from [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) after the percolation loop and before `water_delta` is applied, so it works on `water .+ water_delta` as the current pore water. A no-op under `mp.runoff_method === :instantaneous`, where it returns `flux_dn[Xi]` with `runoff` untouched, leaving the total identical to MATLAB's `sum(runoff) + flux_dn[Xi]`.

Otherwise the blocked water is pooled (mass-weighted by [`mix_age`](/internals_grid#GEMB.mix_age-NTuple{4,%20Float64}), since each contribution leaves its cell at that cell's outflow age `flux_age[i+1]`), the contributing `runoff` entries are zeroed, and the pool is spent filling cells to [`pore_capacity`](/internals_physics#GEMB.pore_capacity-Tuple{ModelParameters,%20Float64,%20Float64}) from the deepest contributing cell upward. Cells at or above `mp.impermeable_density` are skipped: at pore close-off the pore space is disconnected, which is the premise of their blocking flow. Whatever the column cannot hold is returned as runoff — at that point every cell up to the surface is at capacity, so it genuinely leaves.

Mass is conserved by construction (every kilogram is either placed in a cell or returned), and so is energy: pore water and runoff carry the same melt-point enthalpy [`specific_enthalpy_water`](/internals_physics#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real}), so moving mass between them changes neither budget term's per-kilogram value. Nothing refreezes here; see the ponding section of [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) on why the next timestep handles that.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.pore_capacity-Tuple{ModelParameters, Float64, Float64}' href='#GEMB.pore_capacity-Tuple{ModelParameters, Float64, Float64}'><span class="jlbinding">GEMB.pore_capacity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pore_capacity(mp, density, dz) -> capacity [kg m-2]
```


Maximum water a cell can hold, `mp.pore_saturation_max·(ρᵢ − ρ)·dz`.

Deliberately the same `(ρᵢ − ρ)·S·dz` form as the retention in [`irreducible_saturation`](/internals_physics#GEMB.irreducible_saturation-Tuple{ModelParameters,%20Float64}), with `S = mp.pore_saturation_max` in place of `S_wi`, so that `water/capacity` is exactly the saturation `S_wi` is defined against: a cell holding its irreducible water sits at `S_wi/pore_saturation_max` of capacity, and one at `S_wi = 1` reads as full rather than over-full. The physically exact pore mass `ρ_w·dz·(1 − ρ/ρᵢ)` would disagree with the retention expression unless `ρ_w = ρᵢ`, and mixing the two conventions would let a cell exceed capacity while still under saturation 1.

Returns 0 for a cell at or above `mp.density_ice`, which has no pore space at all.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.remove_melt!-Tuple{ModelParameters, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, Float64}' href='#GEMB.remove_melt!-Tuple{ModelParameters, Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, Float64}'><span class="jlbinding">GEMB.remove_melt!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
remove_melt!(mp, M, density, dz, i, melt_i) -> dz_cell [m]
```


Take `melt_i` [kg m-2] of ice out of cell `i` and return the cell's thickness afterwards, the value every downstream pore-space expression in the percolation loop of [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) is written against.

Which of the cell's two geometric variables absorbs the loss is `mp.melt_geometry`:
- `:thickness` (the default) — density is held fixed and the thickness shrinks to `M/ρ`. This is Crocus's treatment (D'Amboise et al., 2017) and GEMB's historical behaviour.
  
- `:density` — the thickness is held fixed at `dz[i]` and `density[i]` is lowered to `M/dz`. This is SNOWPACK's treatment (Bartelt and Lehning, 2002) and the one Fourteau et al. (2026) Sect. 2.3 argues is physically correct: melting occurs _within_ the snow microstructure, removing ice from the pore walls rather than collapsing the layer, and the observed high density of wet snow is better attributed to the low viscosity of wet snow under overburden — a densification process — than to melt-driven thinning.
  

`:density` also removes an asymmetry. Refreezing in this loop is already at constant thickness (it fills pore space, `density[i] = M[i]/dz_0`), so under `:thickness` a melt-refreeze cycle that restores a cell's mass does not restore its geometry, while under `:density` it does.

Both settings conserve mass — only the (`dz`, `ρ`) split of the same `M` differs — so neither affects the budget checks. `dz` itself is rebuilt as `M ./ density` at the end of [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}), which is why this returns the thickness rather than writing it.

**References**
- Fourteau, K., Brondex, J., Cancès, C., and Dumont, M. (2026). Numerical strategies for representing Richards' equation and its couplings in snowpack models. _Geosci. Model Dev._, 19, 3193–3212. Sect. 2.3.
  
- Bartelt, P. and Lehning, M. (2002). A physical SNOWPACK model for the Swiss avalanche warning: Part I. _Cold Reg. Sci. Technol._, 35, 123–145.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_shortwave_radiation-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, ClimateForcingStep, ModelParameters}' href='#GEMB.calculate_shortwave_radiation-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, ClimateForcingStep, ModelParameters}'><span class="jlbinding">GEMB.calculate_shortwave_radiation</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_shortwave_radiation(dz, density, grain_radius, albedo_broadband, albedo_diffuse, cfs::ClimateForcingStep, mp::ModelParameters)
```


Distribute absorbed shortwave radiation vertically within snow/ice.

Depending on model configuration:
2. Surface Absorption: All net shortwave energy is absorbed by the top grid cell (`shortwave_subsurface_absorption = false`).
  
3. Subsurface Penetration: Shortwave energy penetrates and is absorbed by deeper layers (`shortwave_subsurface_absorption = true`), using either:
  - Density-dependent extinction (Bassford, 2002)
    
  - Spectral-dependent extinction for "BrunLefebre" (Lefebre et al., 2003)
    
  

Returns `shortwave_flux` vector [W m-2] of absorbed shortwave radiation per grid cell.

**References**
- Lefebre, F., et al. (2003). J. Geophys. Res., 108, 4231.
  
- Greuell, W. and Konzelmann, T. (1994). Global Planet. Change, 9, 91-114.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._ThermalSurface' href='#GEMB._ThermalSurface'><span class="jlbinding">GEMB._ThermalSurface</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
_ThermalSurface
```


Sub-timestep-invariant surface quantities, computed once per `calculate_temperature` call and handed to the thermal solver. `isbits`, so it is stack-allocated and the field reads fold away.

Carries only the surface energy balance's own inputs — roughness lengths, air density, the hoisted turbulent-flux factors, and the emissivity state. Diagnostics travel separately.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._emissivity_initialize-Tuple{Float64, ModelParameters}' href='#GEMB._emissivity_initialize-Tuple{Float64, ModelParameters}'><span class="jlbinding">GEMB._emissivity_initialize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_emissivity_initialize(grain_radius_surface, mp::ModelParameters)
```


Initialize emissivity based on surface grain radius and model parameters. Returns `(emissivity, emissivity_melt_switch)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._find_dt_divisor-Tuple{Float64, Vector{Float64}}' href='#GEMB._find_dt_divisor-Tuple{Float64, Vector{Float64}}'><span class="jlbinding">GEMB._find_dt_divisor</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_find_dt_divisor(dt_target, dt_divisors)
```


Find the largest dt_divisor that is &lt;= dt_target. Allocation-free.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._max_safe_dt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}' href='#GEMB._max_safe_dt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}'><span class="jlbinding">GEMB._max_safe_dt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_max_safe_dt(temperature, dz, density, K, mp)
```


Largest explicit sub-step [s] that keeps the thermal solve stable. Allocation-free.

The scheme in [`calculate_temperature`](/internals_physics#GEMB.calculate_temperature-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters,%20Bool}) updates cell `i` as `H[i] += F_i - F_{i-1}` with `F_i = G_i·dt·(T_{i+1} - T_i)`, so the coefficient of `T_i` in the equivalent temperature update is `1 - dt·(G_i + G_{i-1})/(ρᵢcᵢdzᵢ)`. Requiring it to stay non-negative — the positivity condition that is Von Neumann stability for this stencil — gives

```julia
dt ≤ ρᵢ cᵢ dzᵢ / (Gᵢ + Gᵢ₋₁),    Gᵢ = 1 / (dz[i+1]/(2K[i+1]) + dz[i]/(2K[i]))
```


evaluated here by reusing the very harmonic-mean face conductances the solve uses, so the limit and the scheme cannot drift apart.

This replaced the textbook uniform-grid form `0.5·ρᵢcᵢdzᵢ²/Kᵢ`, which substitutes `2Kᵢ/dzᵢ` for `Gᵢ + Gᵢ₋₁`. The two agree exactly when `dz` and `K` are uniform, but GEMB's grid is graded, and the substitution is not conservative there: since `Gᵢ ≤ 2Kᵢ/dzᵢ` with equality only in the limit of a vanishing neighbour resistance, `Gᵢ + Gᵢ₋₁` can reach `4Kᵢ/dzᵢ`, so the old form could overestimate the true limit by up to a factor of two. Measured per cell on GEMB's own column the ratio spanned 0.66 to 1.82 — an error in both directions, the low end of which the 0.8 safety factor does not cover. It is also _looser_ at the cell that binds, so the correct limit is both sound and cheaper.

Cell `m` is excluded: it is the Dirichlet reservoir, its enthalpy is never updated (`Q_sw[m] = 0`, no flux is drawn from it), so it imposes no stability constraint. The old form included it.

::: warning Diffusion only — the surface cell is under-constrained

Cell 1 also carries the explicitly-evaluated surface energy balance, whose slope `Λ = dQ_sfc/dT₁ ≤ 0` enters its own-temperature coefficient exactly as a face conductance does: `coef₁ = 1 − dt·(G₁ + |Λ|)/(ρ₁c₁dz₁)`. The `|Λ|` is not included here, so the value returned is always an _over_estimate of the true surface constraint, by a factor the model does not bound. That is what `mp.thermal_explicit_safety_factor` covers, and why it is not simply slack — see [`THERMAL_EXPLICIT_SAFETY_FACTOR`](/internals_support#GEMB.THERMAL_EXPLICIT_SAFETY_FACTOR) for the measured margin. Folding `|Λ|` in here directly is the open option; it would need the surface state, which this function deliberately does not take.

:::

**References**
- Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_, Ch. 4 (positivity of the explicit finite-volume coefficients).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._surface_energy_balance-Tuple{Float64, Float64, GEMB._ThermalSurface, ClimateForcingStep}' href='#GEMB._surface_energy_balance-Tuple{Float64, Float64, GEMB._ThermalSurface, ClimateForcingStep}'><span class="jlbinding">GEMB._surface_energy_balance</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_surface_energy_balance(T1, emissivity, sfc, cfs)
```


Surface energy balance at surface-cell temperature `T1` [K]. Returns `(longwave_upward, heat_flux_sensible, heat_flux_latent, latent_heat, T_surface)` with the three fluxes in W m-2 (signed as energy _into_ the column, so `longwave_upward` is negative) and `T_surface = min(CtoK, T1)`.

The clamp is what makes the balance non-smooth at the melting point: above it the surface stops responding to `T1` entirely, so every flux — and every derivative — goes flat. Shared by the implicit solver's Newton iteration and its final applied-flux evaluation so the two cannot drift.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._surface_energy_balance_slope-Tuple{Float64, Float64, GEMB._ThermalSurface, ClimateForcingStep, Float64, Float64, Float64}' href='#GEMB._surface_energy_balance_slope-Tuple{Float64, Float64, GEMB._ThermalSurface, ClimateForcingStep, Float64, Float64, Float64}'><span class="jlbinding">GEMB._surface_energy_balance_slope</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_surface_energy_balance_slope(T1, emissivity, sfc, cfs, longwave_upward, heat_flux_sensible, heat_flux_latent)
```


Slope `dQ_surface/dT1` [W m-2 K-1] of [`_surface_energy_balance`](/internals_physics#GEMB._surface_energy_balance-Tuple{Float64,%20Float64,%20GEMB._ThermalSurface,%20ClimateForcingStep}) at `T1`, given the fluxes already evaluated there. Returns a value `≤ 0`.

The longwave term is analytic, `-4σεT³`. The two turbulent terms come from a single extra `_turbulent_heat_flux` call, differenced with [`THERMAL_BC_DERIVATIVE_STEP`](/internals_support#GEMB.THERMAL_BC_DERIVATIVE_STEP): the Beljaars-Holtslag stability branches make hand-differentiation a maintenance liability for no benefit, because this slope sets only the Newton convergence _rate_ and never the converged answer (see [`_thermal_solve!`](/internals_physics#GEMB._thermal_solve!-Tuple{ExplicitThermal,%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Vector{Float64},%20GEMB._ThermalSurface,%20ClimateForcingStep,%20ModelParameters,%20Bool})).

Three deliberate crudenesses, each licensed by that:
- **Above the melting point the slope is exactly zero.** The `min(CtoK, T1)` clamp makes `Q` independent of `T1` there, so this is not an approximation.
  
- **The difference is taken away from the clamp.** A forward step that would cross `CtoK` is replaced by a backward one; `Q` is smooth below the melting point, so either is equally good.
  
- **Each component is clamped at `≤ 0`.** A positive slope would weaken the Newton diagonal instead of strengthening it, and it can arise spuriously — the `LV`/`LS` latent-heat switch at `turbulent_heat_flux.jl` is a genuine jump in `heat_flux_latent` right at the melting point, where wet columns sit. Clamping costs convergence rate, never correctness.
  

::: tip An analytic turbulent slope was tried and reverted

The extra `_turbulent_heat_flux` call is the dominant cost of the implicit path, so an analytic replacement was measured. `dQ_shf/dT` and `dQ_lhf/dT` are both algebraic in quantities the flux evaluation already computes — the transfer-coefficient products and the saturation vapour pressure — _if_ the coefficients themselves are held fixed. They cannot be: `coefM`, `coefHT`, `coefHQ` depend on `T_surface` through the bulk Richardson number, and the term dropped by freezing them substantially cancels the rest. Measured over 675 sampled surface states, the frozen-coefficient slope overshot the true derivative by up to 10x, with 9% of states outside a factor of two. Overshoot under-relaxes Newton: the iteration count more than doubled, hit `THERMAL_IMPLICIT_MAX_ITERATIONS`, and whole-model runtime _rose_ 43% (15.2 s to 21.8 s) while the unconverged output drifted 15 K. Correct in principle — `Λ` really does affect only the rate — but the rate is exactly what it cost. Differentiating the stability branches too would be sound, and is the open option here.

:::


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_solve!-Tuple{ExplicitThermal, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, GEMB._ThermalSurface, ClimateForcingStep, ModelParameters, Bool}' href='#GEMB._thermal_solve!-Tuple{ExplicitThermal, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, GEMB._ThermalSurface, ClimateForcingStep, ModelParameters, Bool}'><span class="jlbinding">GEMB._thermal_solve!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_solve!(::ExplicitThermal, temperature, dz, density, K, shortwave_flux, water_surface, grain_radius, sfc, cfs, mp, verbose)
```


Advance the column temperature over one forcing timestep with the explicit finite-volume scheme, sub-stepped to the Von Neumann stability limit (see [`_max_safe_dt`](/internals_physics#GEMB._max_safe_dt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters})).

`temperature` is updated in place. Returns the forcing-step averages `(longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal, evaporation_condensation)`; the first four are W m-2, the last is kg m-2 accumulated over the step.

`water_surface` and `grain_radius` are carried for the verbose diagnostic only.

This is one implementation of the [`AbstractThermalSolver`](/api#GEMB.AbstractThermalSolver) interface; see there for the invariants every scheme must honour.

**References**
- Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_, Ch. 3-4.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_solve!-Tuple{ImplicitThermal, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, GEMB._ThermalSurface, ClimateForcingStep, ModelParameters, Bool}' href='#GEMB._thermal_solve!-Tuple{ImplicitThermal, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, GEMB._ThermalSurface, ClimateForcingStep, ModelParameters, Bool}'><span class="jlbinding">GEMB._thermal_solve!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_solve!(::ImplicitThermal, temperature, dz, density, K, shortwave_flux, water_surface, grain_radius, sfc, cfs, mp, verbose)
```


Advance the column temperature over one forcing timestep with a backward-Euler scheme on a tridiagonal system, solved by the Thomas algorithm. Unconditionally stable, so it takes no stability sub-steps: [`_max_safe_dt`](/internals_physics#GEMB._max_safe_dt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters}), [`_find_dt_divisor`](/internals_physics#GEMB._find_dt_divisor-Tuple{Float64,%20Vector{Float64}}), `mp.dt_divisors` and [`DT_MIN_WARN`](/internals_support#GEMB.DT_MIN_WARN) are all unused on this path, and a thin refrozen lens costs nothing.

`temperature` is updated in place; the return is the [`AbstractThermalSolver`](/api#GEMB.AbstractThermalSolver) 5-tuple. `water_surface` and `grain_radius` are carried for the verbose diagnostic only.

**The system**

Unknowns are cells `1 … m-1`. Cell `m` is the Dirichlet reservoir and is never an unknown, so its temperature is returned bit-unchanged by construction; it enters row `m-1` as the known term `A_face[m-1]·temperature[m]` on the right-hand side. Its shortwave share is absorbed by cell `m-1`, exactly as on the explicit path.

Row `i` is the backward-Euler enthalpy balance `M_i·(h(T_i) − h(T_i^old)) = Q_sw,i + F_i − F_{i-1}` with the diffusive fluxes evaluated at the _new_ temperatures, giving diagonal `M_i·c_i + A_face[i-1] + A_face[i]` and off-diagonals `−A_face[i]`. The face conductances `A_face` are the same harmonic means the explicit path uses, so the two schemes discretize identical conductances and differ only in time integration.

`c_i` is the [`chord_heat_capacity`](/internals_physics#GEMB.chord_heat_capacity-Tuple{ModelParameters,%20Real,%20Real}) between the old and new temperatures — an exact identity for `h = aT + (b/2)T²`, not a linearization. This is what lets enthalpy stay the conserved quantity while the system solved is linear in temperature. It also makes the system mildly nonlinear under `:CuffeyPaterson`, resolved by the same iteration as the surface row; under the default `:constant` it is a true constant and the iteration only sees the surface.

`M_i·c_i > 0` always, so the matrix is an irreducibly diagonally dominant M-matrix: non-singular, **no pivoting required**, and `A⁻¹ ≥ 0`, which is the discrete maximum principle — the solve cannot manufacture a new extremum.

**The surface row: Newton, not lagged Picard**

The surface flux `Q(T_1)` is strongly nonlinear and must be taken implicitly. Lagging it — a Picard iteration on the previous iterate's flux — _diverges_ here: the fixed-point gain `|dQ/dT_1|·Δt/(M_1 c_1)` exceeds 1 for a centimetre-scale surface cell at any Δt of interest. (That is precisely why the explicit path's lagged sub-steps are stable: their tiny `dt` drives the same gain far below 1.)

So the flux is linearized into the diagonal, `Q(T_1) ≈ Q_k + Λ_k(T_1 − T_{1,k})` with `Λ_k ≤ 0` from [`_surface_energy_balance_slope`](/internals_physics#GEMB._surface_energy_balance_slope-Tuple{Float64,%20Float64,%20GEMB._ThermalSurface,%20ClimateForcingStep,%20Float64,%20Float64,%20Float64}), which _strengthens_ the already dominant diagonal. Because `Q_k` is the true nonlinear flux at the iterate, the two `Λ` terms cancel at convergence and the residual satisfied is the true nonlinear surface balance. **`Λ` therefore affects only the convergence rate, never the answer** — the licence for its finite-difference turbulent terms and its clamp.

The two discontinuities are kept as far _outside_ the iteration as they can be: the emissivity melt switch is evaluated once per sub-step, as on the explicit path, and the `LV`/`LS` latent-heat switch is absorbed by the `Λ` clamp.

That is not sufficient on its own. Measured over a year of 3-hourly synthetic forcing (3.77M sub-step solves), 39% of solves reached `THERMAL_IMPLICIT_MAX_ITERATIONS` without converging, and the iterate traces show **limit cycles** rather than divergence: 2- and 3-cycles of ~0.02 K median amplitude, many straddling 273.15 K, where the latent-heat switch still enters through `Q_k` itself even though `Λ` is clamped. The iteration is therefore **damped**: the step weight is halved whenever a step fails to shrink, which contracts a cycle onto its own mean, and the least-residual iterate is retained so a cycle still running at the cap is not left on its worse phase. This cut non-convergence 44x, to 0.9% of solves, at no measurable runtime cost, and leaves converging solves bit-identical — on the quadratic path the residual falls every iteration, so the weight never leaves 1 and the update is exactly the Newton step. See `THERMAL_IMPLICIT_DAMPING_FLOOR` for the measurement and for why sub-step halving — the remedy Fourteau et al. (2026) Sect. 3.1.3 pairs with their own Newton solve — was tried and rejected here.

**Conservation**

A Thomas sweep leaves a round-off residual, so the solved field is used as a _predictor of the implicit face temperatures_ rather than written to `temperature` directly. The enthalpy update is then applied as one flux per face, `+F` to cell `i` and `−F` to cell `i+1` — structurally the same pass the explicit path uses. The pairwise cancellation is exact to the last bit, so the column total conserves independently of the solve residual, the Newton tolerance, and `c_p(T)`, and the verbose budget check is the explicit path's with no tolerance moved.

The surface flux applied in that pass is the _true_ nonlinear flux at the converged iterate, not its linearization, and the returned averages are those same applied values — so `gemb_core`'s energy budget and the `evaporation_condensation` mass budget close by construction.

**Sub-stepping**

Sub-steps here buy **accuracy, not stability**: backward Euler is first-order and strongly damping, and the emissivity melt switch only resolves between sub-steps. The count follows [`THERMAL_IMPLICIT_DT_TARGET`](/internals_support#GEMB.THERMAL_IMPLICIT_DT_TARGET) and is independent of the cell count, of `dz`, and of the stiffest cell in the column — the property the explicit path lacks.

**Static condensation: Newton on a scalar, not on the column**

The nonlinearity is confined to row 1 — the surface energy balance. Rows `2 … n` are linear in the iterate under the default `:constant` heat capacity, so re-eliminating them once per Newton iteration is wasted work. Instead the interior is condensed **bottom-up, once per sub-step**, each row reduced to an affine function of the cell above it (`T_i = rhs[i] + sup[i]·T_{i-1}`). Newton then iterates on a single scalar equation in `T_1`, at O(1) per iteration, and one forward substitution propagates the answer back down the column.

Under `:CuffeyPaterson` the chord capacities depend on the iterate, so the condensation is repeated in an outer loop until they settle; under `:constant` the outer loop provably runs once.

This is what the exponential-integrator idea in the design notes reduces to in cheap form: the operator is constant across all Newton iterations of a sub-step, so factorize it once and reuse it. Measured: **15.20 s to 5.70 s (2.67x)** on the benchmark below, with whole-model output agreeing to 1.2e-6 K and annual melt to 3.2e-8 kg m-2 — the same converged answer to round-off, as the `Λ`-independence property requires.

**Performance: still slower than explicit, and why**

Measured on a year of 3-hourly synthetic forcing over a 264-cell column, whole-model runtime:

|            Solver |                   `DT_TARGET` | Runtime |
| -----------------:| -----------------------------:| -------:|
| `ExplicitThermal` |                             — |  2.44 s |
| `ImplicitThermal` |          1800 s (6 sub-steps) |  3.59 s |
| `ImplicitThermal` | 900 s (12 sub-steps, default) |  5.75 s |
| `ImplicitThermal` |          450 s (24 sub-steps) |  9.62 s |


So 1.5x at the coarsest usable accuracy and 2.4x at the calibrated default. Unconditional stability is the reason this scheme exists; throughput is not, and on a column with no stiff cell it does not beat the explicit path.

What remains is the Newton iteration itself: **4.19 iterations per sub-step solve**, measured over 1.12M solves. Capping the iteration at 1 takes the run to 4.98 s, so the iteration is most of the remaining cost, and it is irreducible — surface-only and whole-profile convergence criteria give identical counts (4.190 vs 4.190), i.e. cell 1 is always the last cell to converge.

Two things are measured _not_ to be the lever, recorded so they are not retried:
- **The turbulent-flux calls.** Direct timing puts `_surface_energy_balance` at 65 ns and `_surface_energy_balance_slope` at 69 ns, so all flux evaluations together are 0.70 s of the pre-condensation 15.20 s. The earlier attribution of the cost to doubled flux calls was wrong; the cost was the O(n) sweep per iteration, which is what condensation removed.
  
- **The convergence tolerance.** Relaxing `THERMAL_IMPLICIT_T_TOLERANCE` from 1e-10 to 1e-3 — seven orders of magnitude, far looser than the model's other branch tolerances — cut the pre-condensation runtime only from 15.2 s to 10.8 s.
  

Where this scheme wins is the case the explicit path cannot bound: its cost is independent of the stiffest cell. A single 1e-4 m refrozen lens drops the explicit stability limit from 903 s to 3.98 s on this column (measured, `test_calculate_temperature.jl`) — a 227x sub-step increase for one thin cell — while the implicit count does not move. That is the reason to keep it available.

**References**
- Patankar, S. V. (1980). _Numerical Heat Transfer and Fluid Flow_, Ch. 4.
  
- Versteeg, H. K. & Malalasekera, W. (2007). _An Introduction to Computational Fluid Dynamics_, Ch. 8.
  
- Thomas, L. H. (1949). _Elliptic Problems in Linear Difference Equations over a Network_.
  
- Beljaars, A. C. M. & Holtslag, A. A. M. (1991). Flux parameterization over land surfaces. _J. Appl. Meteorol._ 30, 327-341.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.calculate_temperature-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters, Bool}' href='#GEMB.calculate_temperature-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}, ClimateForcingStep, ModelParameters, Bool}'><span class="jlbinding">GEMB.calculate_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
calculate_temperature(temperature, dz, density, water_surface, grain_radius, shortwave_flux, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool)
```


Compute new temperature profile accounting for energy absorption and thermal diffusion.

Solves the 1D heat transfer equation using a finite-volume explicit scheme (Patankar, 1980). Accounts for:
- Surface energy balance (turbulent fluxes, radiative fluxes)
  
- Subsurface thermal diffusion
  
- Shortwave penetration as a source term
  
- Thermal conductivity updates (Sturm, 1997)
  

Sub-time steps are determined by Von Neumann stability analysis.

Returns `(temperature, longwave_upward, heat_flux_sensible, heat_flux_latent, heat_flux_basal, evaporation_condensation)`.

**References**
- Bougamont, M., et al. (2005). (Surface roughness).
  
- Foken, T. (2008). Micrometeorology. (Roughness lengths).
  
- Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow.
  
- Sturm, M., et al. (1997). (Thermal conductivity).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.DENSIFICATION_COEFFS_M01' href='#GEMB.DENSIFICATION_COEFFS_M01'><span class="jlbinding">GEMB.DENSIFICATION_COEFFS_M01</span></a> <Badge type="info" class="jlObjectType jlConstant" text="Constant" /></summary>



Calibrated densification coefficients for the Ligtenberg model.

Each entry maps a region/calibration identifier to a matrix of coefficients:   [M0_550_offset M0_550_slope M0_830_offset M0_830_slope] or (for regions with bare-ice calibration):   [M0_550_offset M0_550_slope M0_830_offset M0_830_slope;    M1_550_offset M1_550_slope M1_830_offset M1_830_slope]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.densification_lookup_M01-Tuple{Symbol}' href='#GEMB.densification_lookup_M01-Tuple{Symbol}'><span class="jlbinding">GEMB.densification_lookup_M01</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
densification_lookup_M01(densification_coeffs_M01::Symbol)
```


Return calibrated densification coefficients for the Ligtenberg model. Matches MATLAB's `densification_lookup_M01.m`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.add_energy_temperature-Tuple{ModelParameters, Real, Real, Real}' href='#GEMB.add_energy_temperature-Tuple{ModelParameters, Real, Real, Real}'><span class="jlbinding">GEMB.add_energy_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
add_energy_temperature(mp::ModelParameters, T, M, Q) -> T [K]
```


Temperature of mass `M` [kg] at temperature `T` [K] after absorbing `Q` joules: `h⁻¹(h(T) + Q/M)`.

For constant `c_p` this is evaluated as `Q/M/c + T` rather than via an `h`/`h⁻¹` round-trip, so the default path is bit-identical to the reference arithmetic. Returns `T` unchanged when `M <= 0`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.chord_heat_capacity-Tuple{ModelParameters, Real, Real}' href='#GEMB.chord_heat_capacity-Tuple{ModelParameters, Real, Real}'><span class="jlbinding">GEMB.chord_heat_capacity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
chord_heat_capacity(mp::ModelParameters, T1, T2) -> c_chord [J kg-1 K-1]
```


Heat capacity that makes the enthalpy difference between two temperatures an exact product:

```julia
specific_enthalpy(mp, T2) - specific_enthalpy(mp, T1) == c_chord * (T2 - T1)
```


For `c_p = a + bT` the enthalpy is `h = aT + (b/2)T²`, so the secant slope of `h` between `T1` and `T2` is `a + b·(T1 + T2)/2` — the heat capacity evaluated at the _midpoint_. This is an algebraic identity, not a linearization: `h` is quadratic, so its chord slope is exact for any `T1`, `T2`, however far apart. Under `:constant` it is `mp.heat_capacity_ice`, and it reduces to [`heat_capacity`](/internals_physics#GEMB.heat_capacity-Tuple{ModelParameters,%20Real})`(mp, T)` when `T1 == T2`.

This is what lets an implicit thermal solver carry enthalpy as the conserved quantity while solving a _linear_ system in temperature. Writing `M·c_p(T_old)·(T_new − T_old)` instead — the obvious lagged form — injects a spurious `(b/2)(T_new − T_old)²` of energy per step, which is the error the [`specific_enthalpy`](/internals_physics#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real}) warning describes in its per-cell form. Only the midpoint evaluation cancels it.

`c_p > 0` over the physical range, so `c_chord > 0` and the diagonal it scales stays positive.

See [`specific_enthalpy`](/internals_physics#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real}) for the enthalpy this is the chord slope of.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cold_content_mass-Tuple{ModelParameters, Real, Real, Real}' href='#GEMB.cold_content_mass-Tuple{ModelParameters, Real, Real, Real}'><span class="jlbinding">GEMB.cold_content_mass</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cold_content_mass(mp::ModelParameters, T, density, dz) -> [kg]
```


As above, for a cell given by its density [kg m-3] and thickness [m] rather than its mass.

Not merely `cold_content_mass(mp, T, density*dz)`: the reference multiplies `ρ` and `dz` into the product separately at this site and pre-multiplied at the other, and floating-point multiplication is not associative. Keeping both groupings is what makes the `:constant` path bit-identical at _both_ call sites.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.cold_content_mass-Tuple{ModelParameters, Real, Real}' href='#GEMB.cold_content_mass-Tuple{ModelParameters, Real, Real}'><span class="jlbinding">GEMB.cold_content_mass</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
cold_content_mass(mp::ModelParameters, T, M) -> [kg]
```


Mass of melting-point water that mass `M` [kg] of ice matrix at temperature `T` [K] can refreeze before reaching the melting point: `M·(h(CtoK) − h(T))/LF`. Zero for `T >= CtoK`.

This is the `freeze_max` of the melt equations. Making it the enthalpy difference rather than `M·c·(CtoK−T)` is what keeps it consistent with [`refreeze_temperature`](/internals_physics#GEMB.refreeze_temperature-Tuple{ModelParameters,%20Real,%20Real,%20Real}) — if the two used different energy maps, refreezing the full cold content would not land the cell exactly at the melting point.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.column_enthalpy-Tuple{ModelParameters, AbstractVector{Float64}, AbstractVector{Float64}}' href='#GEMB.column_enthalpy-Tuple{ModelParameters, AbstractVector{Float64}, AbstractVector{Float64}}'><span class="jlbinding">GEMB.column_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
column_enthalpy(mp::ModelParameters, M, temperature, water) -> E [J]
```


Total enthalpy [J] of a column given per-cell matrix mass `M` [kg], `temperature` [K], and pore water `water` [kg m-2]: `Σ Mᵢh(Tᵢ) + (Σwaterᵢ)·h_water`.

Written as an explicit loop rather than a generator. The verbose energy budgets sit inside large functions (`calculate_melt` especially), and a generator's closure there was enough to push inference over its limit and infer _every_ local in the enclosing function as `Any` — boxing the whole hot path for an 18x allocation regression, even with `verbose=false`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.enthalpy_temperature_scale-Tuple{ModelParameters}' href='#GEMB.enthalpy_temperature_scale-Tuple{ModelParameters}'><span class="jlbinding">GEMB.enthalpy_temperature_scale</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
enthalpy_temperature_scale(mp::ModelParameters) -> s
```


Divisor that [`temperature_from_scaled_enthalpy`](/internals_physics#GEMB.temperature_from_scaled_enthalpy-Tuple{ModelParameters,%20Real}) expects to have been applied to a specific enthalpy already: `mp.heat_capacity_ice` under `:constant`, `1.0` otherwise.

Exists so a caller with a per-cell mass reciprocal can fold `1/c` into it once and turn the inverse map into a multiply. See [`temperature_from_scaled_enthalpy`](/internals_physics#GEMB.temperature_from_scaled_enthalpy-Tuple{ModelParameters,%20Real}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.excess_specific_enthalpy-Tuple{ModelParameters, Real}' href='#GEMB.excess_specific_enthalpy-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.excess_specific_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
excess_specific_enthalpy(mp::ModelParameters, T_excess) -> [J kg-1]
```


Specific enthalpy held above the melting point by ice that was `T_excess` [K] warmer than `CtoK`: `h(CtoK + T_excess) − h(CtoK)`. This is the energy available to melt the cell, so comparing it directly against `LF` replaces the reference's `LF/c_p` temperature threshold — a quantity that only exists because `h` was linear.

Evaluated in the cancellation-free form `a·ΔT + (b/2)·ΔT·(ΔT + 2·CtoK)` rather than as a difference of two large enthalpies.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.excess_temperature_from_specific_enthalpy-Tuple{ModelParameters, Real}' href='#GEMB.excess_temperature_from_specific_enthalpy-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.excess_temperature_from_specific_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
excess_temperature_from_specific_enthalpy(mp::ModelParameters, e) -> ΔT [K]
```


Inverse of [`excess_specific_enthalpy`](/internals_physics#GEMB.excess_specific_enthalpy-Tuple{ModelParameters,%20Real}): the temperature excess above `CtoK` that holds `e` [J kg-1] of excess specific enthalpy. Uses the cancellation-free root of `(b/2)ΔT² + (a + b·CtoK)ΔT − e = 0`.

Needed because the melt equations branch on kelvin thresholds (`T_TOLERANCE` is in kelvin) while the arithmetic is in joules, so the two representations must be interconvertible.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.full_melt_excess_temperature-Tuple{ModelParameters}' href='#GEMB.full_melt_excess_temperature-Tuple{ModelParameters}'><span class="jlbinding">GEMB.full_melt_excess_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
full_melt_excess_temperature(mp::ModelParameters) -> ΔT [K]
```


Temperature excess above the melting point at which a cell holds exactly enough enthalpy to melt its entire mass — the `LF/c_p` of the reference, generalized. A cell hotter than this has surplus energy to pass downward.

Kept in kelvin because the melt equations' branch tests use a kelvin tolerance (`T_tolerance = 1e-10`); reformulating those tests in joules would silently rescale them by ~2100.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.heat_capacity-Tuple{ModelParameters, Real}' href='#GEMB.heat_capacity-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.heat_capacity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
heat_capacity(mp::ModelParameters, T) -> c_p [J kg-1 K-1]
```


Specific heat capacity of the ice matrix at temperature `T` [K], per `mp.heat_capacity_method`:
- `:constant` — `mp.heat_capacity_ice` (default 2102, the melting-point value). MATLAB behaviour.
  
- `:CuffeyPaterson` — Cuffey and Paterson (2010) eq. 9.1, `c_p = 152.5 + 7.122·T`. Gives 2097.9 at the melting point but only 1862 at 240 K and 1648 at 210 K, so the constant value over-permits refreezing in cold firn by +12.6% and +27.5% respectively.
  

Applies to snow, firn, and ice alike: the load-bearing matrix is ice, and pore air and pore water contribute nothing to this term.

For anything that accumulates or balances joules, use [`specific_enthalpy`](/internals_physics#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real}) — not `T * heat_capacity(mp, T)`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.melt_mass_from_excess-Tuple{ModelParameters, Real, Real, Real}' href='#GEMB.melt_mass_from_excess-Tuple{ModelParameters, Real, Real, Real}'><span class="jlbinding">GEMB.melt_mass_from_excess</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
melt_mass_from_excess(mp::ModelParameters, T_excess, density, dz) -> [kg]
```


Mass of ice [kg] that the enthalpy held above the melting point can melt, for a cell of `density` [kg m-3] and thickness `dz` [m] whose temperature exceeded `CtoK` by `T_excess` [K]: `h_excess·ρ·dz/LF`.

Grouped to match the reference's `T_excess·ρ·dz·c/LF` on the `:constant` path.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.mix_temperature-Tuple{ModelParameters, Vararg{Real, 4}}' href='#GEMB.mix_temperature-Tuple{ModelParameters, Vararg{Real, 4}}'><span class="jlbinding">GEMB.mix_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mix_temperature(mp::ModelParameters, T1, M1, T2, M2) -> T [K]
mix_temperature(mp::ModelParameters, h1, M1, h2, M2, ::Val{:enthalpy}) -> T [K]
```


Temperature of the mixture of two masses, conserving enthalpy: solve `M₁h(T₁) + M₂h(T₂) = (M₁+M₂)·h(T)` for `T`.

This is the only correct way to combine cells of different temperature under a temperature-dependent `c_p`. The mass-weighted mean temperature `(M₁T₁ + M₂T₂)/(M₁+M₂)` is _not_ equivalent: `h` is convex, so the mean loses `(b/2)·M·Var_M(T)` joules — about 2 225 J for an equal-mass pair 5 K apart, against an energy-conservation tolerance of `E_TOLERANCE = 1e-3` J.

The second form takes specific enthalpies directly, for callers that already hold them (the mixing partner's enthalpy is often `specific_enthalpy_water` or `h(T_air) + LF` rather than a temperature).

Under `:constant` the first form mixes in temperature space. That is not an approximation: for constant `c_p` the `c` cancels out of the enthalpy balance, so the mass-weighted mean _is_ the exact enthalpy-conserving mixture. Doing it directly rather than via a `h`/`h⁻¹` round-trip keeps the default path bit-identical to the MATLAB arithmetic, which would otherwise drift at ~1e-13 per merge.

Returns `T1` when the total mass is zero.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.mix_temperature_liquid-Tuple{ModelParameters, Vararg{Real, 4}}' href='#GEMB.mix_temperature_liquid-Tuple{ModelParameters, Vararg{Real, 4}}'><span class="jlbinding">GEMB.mix_temperature_liquid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mix_temperature_liquid(mp::ModelParameters, T_liquid, M_liquid, T_solid, M_solid) -> T [K]
```


Temperature after mixing `M_liquid` [kg] of liquid water at `T_liquid` [K] into `M_solid` [kg] of ice matrix at `T_solid` [K], conserving enthalpy. The liquid carries the latent heat of fusion, so its specific enthalpy is [`specific_enthalpy_water`](/internals_physics#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real})`(mp, T_liquid)`.

Two callers, and they sit on opposite sides of the melting point:
- **Rain on snow** (`calculate_accumulation`), with `T_liquid = T_air > CtoK`. The rain's sensible heat above the melting point is carried at `mp.rain_heat_capacity` — see [`specific_enthalpy_water`](/internals_physics#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real}).
  
- **Refreezing** ([`refreeze_temperature`](/internals_physics#GEMB.refreeze_temperature-Tuple{ModelParameters,%20Real,%20Real,%20Real})), with `T_liquid = CtoK` exactly. The sensible term vanishes and this reduces to the melting-point form.
  

Under `:constant` the `c` cancels out of the enthalpy balance, so both are evaluated as a weighted mean of temperatures with the liquid at an effective temperature `h_liquid/c` — `T_liquid + LF/c` when the liquid's heat capacity is the matrix's, which is the grouping the MATLAB reference uses. Keeping that grouping is what makes the refreeze path (and `rain_heat_capacity = :ice`) bit-identical to the reference rather than drifting through an `h`/`h⁻¹` round-trip.

Returns `T_solid` when the total mass is zero.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.refreeze_temperature-Tuple{ModelParameters, Real, Real, Real}' href='#GEMB.refreeze_temperature-Tuple{ModelParameters, Real, Real, Real}'><span class="jlbinding">GEMB.refreeze_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
refreeze_temperature(mp::ModelParameters, T, M_new, freeze) -> T [K]
```


Temperature after `freeze` [kg] of melting-point pore water refreezes into ice matrix that was at temperature `T` [K]. `M_new` is the matrix mass _after_ the refrozen mass is added, matching the melt equations' bookkeeping order.

Releases both the latent heat and the sensible heat of cooling the new ice from `CtoK` to the mixture temperature — i.e. it is [`mix_temperature_liquid`](/internals_physics#GEMB.mix_temperature_liquid-Tuple{ModelParameters,%20Vararg{Real,%204}}) with the liquid at the melting point, written in the reference's incremental form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.specific_enthalpy-Tuple{ModelParameters, Real}' href='#GEMB.specific_enthalpy-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.specific_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
specific_enthalpy(mp::ModelParameters, T) -> h [J kg-1]
```


Specific enthalpy of the ice matrix at temperature `T` [K], measured from a 0 K reference: `h(T) = ∫₀ᵀ c_p dT′`.
- `:constant` — `c·T` exactly, so the default path is bit-identical to the MATLAB `M·T·C_ICE` accounting.
  
- `:CuffeyPaterson` — `a·T + (b/2)T²`.
  

::: warning Warning

Never use `T * heat_capacity(mp, T)` in an energy budget. For `c_p = a + bT` that product is `a·T + b·T²`, which overstates the enthalpy by `(b/2)T² = 2.66e5 J kg-1` at the melting point — 0.79 × the latent heat of fusion. It is a dominant error, not a correction. `test/test_heat_capacity.jl` guards against reintroducing it.

:::

The 0 K reference makes absolute enthalpies incomparable _between_ methods (574 061 vs 307 385 J kg-1 at the melting point), so a reported `E_thermal` jumps when `heat_capacity_method` changes. Only differences are physically meaningful.

See [`temperature_from_specific_enthalpy`](/internals_physics#GEMB.temperature_from_specific_enthalpy-Tuple{ModelParameters,%20Real}) for the inverse.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.specific_enthalpy_water-Tuple{ModelParameters, Real}' href='#GEMB.specific_enthalpy_water-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.specific_enthalpy_water</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
specific_enthalpy_water(mp::ModelParameters, T) -> h [J kg-1]
```


Specific enthalpy of liquid water at temperature `T` [K]: melting-point water plus its sensible heat above the melting point, `h_water + c_water·(T − CtoK)`.

Rain is the only liquid GEMB carries above `CtoK`, and it is the only caller. Everything else — pore water, refreezing, runoff — is isothermal at the melting point and uses the one-argument form, which this reduces to exactly at `T == CtoK`.

`mp.rain_heat_capacity` selects the heat capacity of that sensible term:
- `:water` (default) — [`HEAT_CAPACITY_WATER`](/internals_support#GEMB.HEAT_CAPACITY_WATER) ≈ 4220 J kg-1 K-1, the physical value.
  
- `:ice` — `heat_capacity(mp, CtoK)` ≈ 2102, which is what MATLAB and GEMB.jl before this option used. It understates the rain's sensible heat by about a factor two: ~3.9 kJ kg-1 for rain at 275 K, against 334.5 kJ kg-1 of `LF`. Small (~1% of the rain's energy) but systematic, and growing linearly with `T_air − CtoK`.
  

Below `CtoK` this returns the melting-point value unchanged rather than extrapolating: liquid colder than the melting point is not a state the column represents, and the callers reach this only on the rain branch, which is gated on `T_air > mp.rain_temperature_threshold`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.specific_enthalpy_water-Tuple{ModelParameters}' href='#GEMB.specific_enthalpy_water-Tuple{ModelParameters}'><span class="jlbinding">GEMB.specific_enthalpy_water</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
specific_enthalpy_water(mp::ModelParameters) -> h [J kg-1]
```


Specific enthalpy of liquid water at the melting point: the enthalpy of ice at `CtoK` plus the latent heat of fusion. This is the `LF + CtoK * C_ICE` of the MATLAB budgets, generalized to a temperature-dependent `c_p`.

GEMB carries pore water only at the melting point, so no liquid heat capacity is needed.

See the two-argument form for liquid _above_ the melting point, which only rain is.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.surplus_energy-Tuple{ModelParameters, Real, Real, Real}' href='#GEMB.surplus_energy-Tuple{ModelParameters, Real, Real, Real}'><span class="jlbinding">GEMB.surplus_energy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
surplus_energy(mp::ModelParameters, T_excess, T_surplus, M) -> Q [J]
```


Energy held by mass `M` [kg] _beyond_ what melting its entire mass would consume, given a temperature excess `T_excess` [K] above the melting point and the corresponding surplus `T_surplus = max(0, T_excess − full_melt_excess_temperature(mp))`.

Zero when `T_surplus` is zero, so this agrees exactly with the kelvin-space mask the melt equations branch on. Under `:constant` it evaluates as `T_surplus·c·M` (the reference's grouping); otherwise as `(h_excess(T_excess) − LF)·M`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.temperature_from_scaled_enthalpy-Tuple{ModelParameters, Real}' href='#GEMB.temperature_from_scaled_enthalpy-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.temperature_from_scaled_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
temperature_from_scaled_enthalpy(mp::ModelParameters, hs) -> T [K]
```


[`temperature_from_specific_enthalpy`](/internals_physics#GEMB.temperature_from_specific_enthalpy-Tuple{ModelParameters,%20Real}) for a caller that has already divided the specific enthalpy by [`enthalpy_temperature_scale`](/internals_physics#GEMB.enthalpy_temperature_scale-Tuple{ModelParameters}). Under `:constant` that division _is_ the whole inverse, so this is the identity; under `:CuffeyPaterson` the scale is 1 and this is the plain inverse.

The point is to get the division off an inner loop. The thermal solver converts enthalpy to temperature twice per cell per sub-step; folding `1/c` into the mass reciprocal it already precomputes replaces both divisions with multiplies, which is worth ~40% of the solver's runtime — division is the one arithmetic op on that loop that does not pipeline.

Folding the two constants into one reciprocal is a reassociation, so results differ from `h/M/c` in the last bit or two. That is not gated by the bench fingerprint: the 75-year non-converged spinup it drives amplifies any 1-ULP perturbation to O(1e-2) in the melt fields, as a control reassociation in the _unmodified_ solver confirms.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.temperature_from_specific_enthalpy-Tuple{ModelParameters, Real}' href='#GEMB.temperature_from_specific_enthalpy-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB.temperature_from_specific_enthalpy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
temperature_from_specific_enthalpy(mp::ModelParameters, h) -> T [K]
```


Inverse of [`specific_enthalpy`](/internals_physics#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real}). `h` is in J kg-1 from a 0 K reference.

For `:CuffeyPaterson` this inverts the monotone quadratic `h = aT + (b/2)T²` using the cancellation-free root

```julia
T = 2h / (a + √(a² + 2bh))
```


which is well-conditioned for small `h`, unlike the textbook quadratic formula. `c_p > 0` over the physical range, so `h` is strictly increasing and the root is unique.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.apply_lateral_drainage!-Tuple{NamedTuple, Float64, ModelParameters}' href='#GEMB.apply_lateral_drainage!-Tuple{NamedTuple, Float64, ModelParameters}'><span class="jlbinding">GEMB.apply_lateral_drainage!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
apply_lateral_drainage!(cols, dt_seconds, mp) -> runoff [kg m-2]
```


Drain water held above irreducible saturation out of the column laterally, and return the mass removed.

Runs after [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) in [`gemb_core`](/internals_grid#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool}), on the water that [`pond_blocked_water!`](/internals_physics#GEMB.pond_blocked_water!-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters}) left standing. Only `cols.water` is touched — not `dz`, `density`, or the cell count — so both grid invariants are untouched and the caller needs only to add the returned mass to its runoff total. No energy term is needed either: the runoff enthalpy in [`gemb_core`](/internals_grid#GEMB.gemb_core-Tuple{Any,%20ClimateForcingStep,%20ModelParameters,%20Bool}) is derived from the runoff mass at [`specific_enthalpy_water`](/internals_physics#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real}), the same value the water carried as pore water.

Per cell, with `excess = water - irreducible` and by `mp.runoff_method`:
- `:instantaneous` — returns 0 and touches nothing. Blocked water has already left inside [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}), and no cell holds any excess to drain.
  
- `:ZuoOerlemans` — `excess·min(1, dt/τ)` with `τ` from [`runoff_timescale`](/internals_physics#GEMB.runoff_timescale-Tuple{ModelParameters}). Slope enters only through `τ`. The `min(1, ·)` matters at long steps: at the 25 d flat-terrain timescale a daily step drains 4% of the excess, but a monthly forcing step would otherwise ask for 120% of it.
  
- `:Darcy` — `ρ_w·dt·K_sat·K_rel·S`, capped at the excess. Darcy's law for a lateral flux under a hydraulic gradient equal to the surface slope, with conductivity from [`hydraulic_conductivity_saturated`](/internals_physics#GEMB.hydraulic_conductivity_saturated-Tuple{Float64,%20Float64}) scaled by [`relative_permeability`](/internals_physics#GEMB.relative_permeability-Tuple{Float64,%20Float64,%20Float64}). Slope, grain size, density, and saturation all enter, so unlike the timescale law this drains a coarse saturated layer orders of magnitude faster than a fine barely-wet one.
  

Both are strictly capped at the excess, so a cell is never drained below irreducible: the scheme removes standing water, never capillary-held water.

`age` is untouched. GEMB carries one mean age per cell, over matrix and pore water together, so the water leaving carries `age[i]` by construction and the mean of what remains is unchanged — the same argument that makes melt removal age-neutral in [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.hydraulic_conductivity_saturated-Tuple{Float64, Float64}' href='#GEMB.hydraulic_conductivity_saturated-Tuple{Float64, Float64}'><span class="jlbinding">GEMB.hydraulic_conductivity_saturated</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
hydraulic_conductivity_saturated(grain_radius, density) -> K_sat [m s-1]
```


Saturated hydraulic conductivity of snow/firn, Calonne et al. (2012) eq. 6:

```julia
K_sat = 3r²·ρ_w·g/μ·exp(-0.013ρ)
```


with `r` the grain radius **in metres** and `μ` the dynamic viscosity of water. GEMB carries `grain_radius` in millimetres (see `RE_NEW_SNOW`), so this converts; the conductivity goes as `r²`, making a missed conversion a factor of 10⁶.

Cross-checked against the Community Firn Model's `hydrconducsat_Calonne` (`CFM_main/darcy_funcs.py`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.relative_permeability-Tuple{Float64, Float64, Float64}' href='#GEMB.relative_permeability-Tuple{Float64, Float64, Float64}'><span class="jlbinding">GEMB.relative_permeability</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_permeability(grain_radius, density, saturation_effective) -> K_rel [-]
```


Relative hydraulic conductivity from the van Genuchten (1980) model with the Yamaguchi et al. (2012) parameterization, as Hirashima et al. (2010) eq. 10 writes it:

```julia
n     = 1 + 2.7e-3·(ρ/2r)^0.61        (Yamaguchi 2012 eq. 7, r in metres)
m     = 1 - 1/n
K_rel = θ_e^½·(1 - (1 - θ_e^(1/m))^m)²
```


`saturation_effective` is `θ_e`, the water content above irreducible as a fraction of the drainable pore space, clamped away from both endpoints by the caller. `K_rel` rises from 0 at `θ_e = 0` to 1 at `θ_e = 1`, so a barely-wet cell drains far more slowly than a saturated one at the same conductivity — the nonlinearity that distinguishes this from a fixed timescale.

`grain_radius` is in millimetres, as everywhere in GEMB; the Yamaguchi fit needs `ρ/2r` in kg m-4, so the conversion happens here.

Cross-checked against the Community Firn Model's `vG_Yama_params`/`krel_vG` (`CFM_main/darcy_funcs.py`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.runoff_timescale-Tuple{ModelParameters}' href='#GEMB.runoff_timescale-Tuple{ModelParameters}'><span class="jlbinding">GEMB.runoff_timescale</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
runoff_timescale(mp) -> tau [s]
```


Characteristic lateral runoff timescale of Zuo and Oerlemans (1996) eq. 22, `τ = c1 + c2·exp(-c3·S)`, with `S = mp.surface_slope` [m m-1].

Steeper terrain drains faster, asymptotically to `c1`; the timescale saturates at `c1 + c2` (≈ 25 d) on flat terrain rather than diverging, so a zero slope still drains. Independent of grain size, density, and saturation — the whole slope dependence of the scheme is here.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_conductivity_air-Tuple{Real}' href='#GEMB._thermal_conductivity_air-Tuple{Real}'><span class="jlbinding">GEMB._thermal_conductivity_air</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_conductivity_air(T) -> K [W m-1 K-1]
```


Thermal conductivity of air, `a*T^1.5/(b + T)` (Reid, 1966), with coefficients `K_AIR_REID_A`/`K_AIR_REID_B`. Read only as the ratio `K_air(T)/K_air(T_ref)` inside [`_thermal_conductivity_calonne2019`](/internals_physics#GEMB._thermal_conductivity_calonne2019-Tuple{Real,%20Real})'s `:Calonne2019Air` form, so its absolute calibration never enters the result — see the note there.

`T^1.5` is written `T*sqrt(T)`, which the hardware does in one instruction where the generic `pow` takes ~15x longer — measured at 14.9 ns vs 1.0 ns per element, which was the whole of a 3.5x slowdown on this function when the air form was added. The two agree exactly for three quarters of the doubles in 150-320 K and differ by 1 ulp (~2e-16 relative) for the rest, twelve orders of magnitude below the fit's own accuracy.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_conductivity_calonne2019-Tuple{Real, Real}' href='#GEMB._thermal_conductivity_calonne2019-Tuple{Real, Real}'><span class="jlbinding">GEMB._thermal_conductivity_calonne2019</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_conductivity_calonne2019(T, ρ; air_factor::Bool) -> K [W m-1 K-1]
```


Calonne et al. (2019) eq. 5. A logistic weight `θ` in density blends two regime-specific fits, each rescaled by the constituent conductivities at `T` relative to their values at the temperature the fits were calibrated at (`CONDUCTIVITY_T_REF` = 270.15 K):

```julia
θ      = 1/(1 + exp(-2a(ρ - ρ_t))),      a = 0.02, ρ_t = 450 kg m-3
k_snow = 0.024 - 1.23e-4ρ + 2.5e-6ρ^2   (the :Calonne 2011 fit, ρ -> 0 gives air)
k_firn = 2.107 + 0.003618(ρ - 917)       (ρ -> 917 gives ice)
K      = (1-θ)·(K_ice/k_i)·(K_air/k_a)·k_snow + θ·(K_ice/k_i)·k_firn
```


The two regimes carry different scalings because they are physically different: the firn fit describes a connected ice skeleton and so scales with ice alone, while the snow fit describes grains in air and scales with both constituents. The reference values are internally consistent with that reading — `k_i = 2.107` is `_thermal_conductivity_ice` at 270.15 K to within 1.2e-4 relative, and `k_a = 0.024` is `_thermal_conductivity_air` there to within 6.6e-3 — so the ratios are 1 at the reference temperature by construction and the fits return their published values.

`air_factor` selects which of the two published readings of the snow term is used:
- `false` (`:Calonne2019`) pins `K_air` at its reference value, so the snow prefactor collapses to `K_ice/k_i` and the whole blend takes a single factor. This matches the Community Firn Model, whose `diffusion.py` hardcodes `K_air = kref_a` with a standing TODO to find the temperature dependence, and it is the form GEMB has always used.
  
- `true` (`:Calonne2019Air`) carries `K_air(T)/K_air(T_ref)` as well. This matches IMAU-FDM (`firn_physics.f90`, `Thermal_Cond`). Since it enters only as a ratio, Reid's absolute calibration cancels and only its temperature _shape_ matters — which is why the 0.7% gap between `K_air(T_ref)` and `k_a` is irrelevant here.
  

Air conducts more poorly as it cools, so the ratio is below 1 below the reference temperature (0.83 at 220 K, 0.88 at 233 K, 0.94 at 253 K) and `:Calonne2019` therefore returns a _higher_ conductivity than `:Calonne2019Air` for cold snow — by up to ~20% at 220 K for ρ ≲ 300, ~10% at ρ = 450, and under 0.3% by ρ = 550 where `θ → 1` and the air-free firn branch takes over. The default is `:Calonne2019`; `:Calonne2019Air` is the fuller form of the published equation.

Both forms are continuous into ice: at ρ = 917 the weight is 1 to within 1e-8 and `k_firn = k_i`, so `K = K_ice(T)` to the same relative precision and no ice branch is needed. The air factor cannot disturb this because it is weighted by `(1-θ)`.

Unlike most density constants in GEMB this uses the literal 917, not `mp.density_ice`: it is a fitted coefficient of the published regression, not a configurable property of the column (the same reasoning as the `:Barnola1991` note in `initialize_parameters.jl`).

Cross-checked against the Community Firn Model's `Calonne2019` conductivity in `CFM_main/diffusion.py` and IMAU-FDM's `Thermal_Cond` in `source/firn_physics.f90`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_conductivity_ice-Tuple{Real}' href='#GEMB._thermal_conductivity_ice-Tuple{Real}'><span class="jlbinding">GEMB._thermal_conductivity_ice</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_conductivity_ice(T) -> K [W m-1 K-1]
```


Thermal conductivity of pure ice, `9.828*exp(-5.7e-3*T)` (Yen, 1981). Used both for the ρ &gt;= `density_ice` branch and as the temperature scaling inside the [`_thermal_conductivity_calonne2019`](/internals_physics#GEMB._thermal_conductivity_calonne2019-Tuple{Real,%20Real}) blend.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._thermal_conductivity_marchenko2019-Tuple{Real}' href='#GEMB._thermal_conductivity_marchenko2019-Tuple{Real}'><span class="jlbinding">GEMB._thermal_conductivity_marchenko2019</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_thermal_conductivity_marchenko2019(ρ) -> K [W m-1 K-1]
```


Marchenko et al. (2019) eq. 30, `K = 0.301e-2ρ - 0.724`, calibrated by fitting modelled to observed subsurface temperatures at Lomonosovfonna, Svalbard, over ρ ≈ 350–900 kg m-3. It gives higher firn conductivity than `:Sturm` or `:Calonne`, which is the direction the multi-model cold bias in RetMIP implies.

The fit is linear and crosses zero at ρ ≈ 241, so it cannot be used at low density. It is floored by the `:Calonne` (2011) snow fit, which the two curves cross near ρ ≈ 321 — just below the calibration range — so the floor takes over continuously and exactly where the regression stops being supported by data. Temperature-independent, so the ice branch (`_thermal_conductivity_ice`) still applies at ρ &gt;= density_ice.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.thermal_conductivity-Tuple{AbstractVector, AbstractVector, ModelParameters}' href='#GEMB.thermal_conductivity-Tuple{AbstractVector, AbstractVector, ModelParameters}'><span class="jlbinding">GEMB.thermal_conductivity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
thermal_conductivity(temperature, density, mp::ModelParameters)
```


Compute thermal conductivity for snow/firn/ice based on density and temperature. Matches MATLAB's `thermal_conductivity.m` for the two methods MATLAB carries.

For snow/firn (density &lt; density_ice), by `mp.thermal_conductivity_method`:
- `:Sturm` — Sturm et al. (1997): `K = 0.138 - 1.01e-3*ρ + 3.233e-6*ρ^2`
  
- `:Calonne` — Calonne et al. (2011): `K = 0.024 - 1.23e-4*ρ + 2.5e-6*ρ^2`
  
- `:Calonne2019` — Calonne et al. (2019) eq. 5, a temperature-dependent sigmoid blend of a snow and a firn regime, with `K_air` pinned at its reference value (matching the Community Firn Model); see [`_thermal_conductivity_calonne2019`](/internals_physics#GEMB._thermal_conductivity_calonne2019-Tuple{Real,%20Real})
  
- `:Calonne2019Air` — the same equation carrying the `K_air(T)/K_air(T_ref)` factor on its snow branch as well (matching IMAU-FDM). Gives lower conductivity than `:Calonne2019` for cold low-density snow — up to ~20% at 220 K — and is identical to it at and above ρ ≈ 550 and at the 270.15 K reference temperature
  
- `:Marchenko2019` — Marchenko et al. (2019) eq. 30, a linear firn fit floored by `:Calonne`; see [`_thermal_conductivity_marchenko2019`](/internals_physics#GEMB._thermal_conductivity_marchenko2019-Tuple{Real})
  

For ice (density &gt;= density_ice):
- K = 9.828 * exp(-5.7e-3 * T)
  

The two 2019 forms are the exception: each is continuous into ice by construction (both return `K_ice(T)` at ρ = 917) and so is evaluated at all densities rather than short-circuiting to the ice branch.

`:Sturm` and `:Calonne` are the MATLAB-equivalent options and either reproduces the reference to 1e-12. The `:Calonne2019` pair and `:Marchenko2019` are additions with no MATLAB counterpart, all recommended over the older fits by Vandecrux et al. (2020, RetMIP Sect. 5.1), who attribute part of a multi-model cold bias at Summit and Dye-2 to conductivity parameterizations.

Returns vector of thermal conductivities [W m-1 K-1].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.thermal_conductivity-Tuple{Real, Real, ModelParameters}' href='#GEMB.thermal_conductivity-Tuple{Real, Real, ModelParameters}'><span class="jlbinding">GEMB.thermal_conductivity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
thermal_conductivity(temperature::Real, density::Real, mp) -> K [W m-1 K-1]
```


Scalar method, for callers with a single (temperature, density) pair — the steady-state initial guess marches one parcel and would otherwise allocate a one-element array per step.

The vector method above loops over this one, so the branch and the calibrated coefficients have a single definition.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._turbulent_heat_flux-Tuple{Float64, Float64, Float64, Float64, Float64, ClimateForcingStep, Vararg{Float64, 7}}' href='#GEMB._turbulent_heat_flux-Tuple{Float64, Float64, Float64, Float64, Float64, ClimateForcingStep, Vararg{Float64, 7}}'><span class="jlbinding">GEMB._turbulent_heat_flux</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs, wind_speed, C, pressure_factor, logM, logHT, logHQ, wind_ratio_sq)
```


Core of [`turbulent_heat_flux`](/internals_physics#GEMB.turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}) with the sub-timestep-invariant quantities (`wind_speed`, bulk coefficient `C`, Exner `pressure_factor`, the neutral roughness logs `logM`/`logHT`/`logHQ`, and the squared wind-to-height ratio `wind_ratio_sq`) hoisted to the caller. `calculate_temperature` computes these once per timestep and reuses them across all thermal sub-steps, where only `T_surface` changes. Numerically identical to the inline form.

**Integrated stability functions**

Both branches return `Ψ(ζ) = ∫₀^ζ (1 - φ(z))/z dz`, the integrated correction that enters the transfer coefficients as `coef = log(z/z₀) - Ψ`. Two properties follow from that definition and are what the branches are checked against in `test_turbulent_heat_flux.jl`:
- **`Ψ(0) = 0`.** The integrand vanishes at neutral, so both branches must meet at `Ri = 0`. Without this the fluxes step discontinuously as `T_surface` crosses `T_air` — and since the implicit solver's Newton iteration converges on exactly that crossing, a jump there can leave the surface energy balance with no solution in `T_surface` (Fourteau et al., 2024, Appendix D, who move their own branch point to `Ri_b = 0` for this reason).
  
- **`Ψ > 0` when unstable, `Ψ < 0` when stable.** Convection enhances exchange, so subtracting a positive `Ψ` shrinks `coef` and raises the flux; stable stratification does the reverse.
  

The stable branch is Beljaars & Holtslag (1991) eqs. 28 and 32, integrated in closed form. The unstable branch is Paulson (1970) with Högström's (1988) coefficients. Paulson's argument is the _inverse_ profile function `x = φ⁻¹`, so both exponents are positive: `x_m = (1 - 19ζ)^{1/4}` inverts `φ_m = (1 - 19ζ)^{-1/4}`, and `x_h = (1 - 11.6ζ)^{1/2}` inverts `φ_h = (1 - 11.6ζ)^{-1/2}`. `Ψ_h = 2·ln((1 + x_h)/2)` takes `x_h` to the first power, the square root already being carried by `x_h` itself.

Both closed forms were verified by numerically integrating their own `φ` back through the definition above, agreeing to 6 decimal places at `ζ = ±0.01, ±0.1, ±1, ±5`.

Högström's `φ_h` carries a `0.95` prefactor (his `κ_H/κ_M` ratio) which makes `φ_h(0) = 0.95` rather than 1. It is dropped here: `Ψ` is defined against a profile function that is 1 at neutral, and retaining it would put `Ψ_h(0) = -0.0999` instead of 0 — reintroducing the very discontinuity described above. The ratio belongs to the neutral transfer coefficient, not inside the integral.

`ζ` is bounded below by `ZETA_UNSTABLE_MIN` on the unstable branch. This is a numerical guard, not physics: `Ri` carries `wind_speed^-2`, so at the `min_wind_speed` floor it reaches values far outside the range Monin-Obukhov theory was ever fit over, where `Ψ_h` would grow past `log(z_T/z_Q)` and drive the transfer coefficient through zero. See that constant's comment.

**References**
- Paulson, C. A. (1970). The mathematical representation of wind speed and temperature profiles in the unstable atmospheric surface layer. _J. Appl. Meteorol._ 9, 857-861.
  
- Högström, U. (1988). Non-dimensional wind and temperature profiles in the atmospheric surface layer: a re-evaluation. _Boundary-Layer Meteorol._ 42, 55-78.
  
- Beljaars, A. C. M. & Holtslag, A. A. M. (1991). Flux parameterization over land surfaces for atmospheric models. _J. Appl. Meteorol._ 30, 327-341.
  
- Fourteau, K., Brondex, J., Brun, F. & Dumont, M. (2024). A novel numerical implementation for the surface energy budget of melting snowpacks and glaciers. _Geosci. Model Dev._ 17, 1903-1929. https://doi.org/10.5194/gmd-17-1903-2024
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.turbulent_heat_flux-Tuple{Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}' href='#GEMB.turbulent_heat_flux-Tuple{Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}'><span class="jlbinding">GEMB.turbulent_heat_flux</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
turbulent_heat_flux(T_surface, density_air, z0, zT, zQ, cfs::ClimateForcingStep; min_wind_speed=0.01)
```


Compute sensible and latent heat fluxes using Monin-Obukhov similarity theory. Matches MATLAB's `turbulent_heat_flux.m`.

`min_wind_speed` floors `cfs.wind_speed` (which would otherwise drive the bulk Richardson number to infinity at zero wind); clamping here avoids rebuilding the `ClimateForcingStep` in the caller.

Returns (heat_flux_sensible, heat_flux_latent, latent_heat) [W m-2, W m-2, J kg-1].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

