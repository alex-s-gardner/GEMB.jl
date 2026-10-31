# Model parameters

```@meta
CurrentModule = GEMB
```


Every option is set through [`initialize_parameters`](@ref), which validates each value at
construction. **Defaults are in bold.** Physics defaults follow the two firn-model
intercomparisons — RetMIP (Vandecrux et al., 2020) and FirnMICE (Lundin et al., 2017) — and,
where those are silent, the shipped configurations of the
[Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) and
[IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM). Each option's full citation and
rationale is in the [`ModelParameters`](@ref) docstring. The exception is the snow set: fresh-snow
density, settling and grain growth default to one configuration calibrated together on SNOTEL
seasonal snow (`new_snow_method = :vanKampenhout2017SNOTEL`, `densification_method = :Crocus` with
destructive metamorphism at `destructive_metamorphism_density = 100`, and
`grain_growth_method = :FlannerZender2006`). `:vanKampenhout2017SNOTEL` deposits light snow that
only `:Crocus` settling with destructive metamorphism at that density compacts at the observed
rate, so changing one member alone leaves the others outside the regime they were fitted in:
change them together. Choosing another `densification_method` turns destructive metamorphism
off. The Community Firn Model's configuration is `new_snow_method = :Constant350`,
`densification_method = :Arthern` and `grain_growth_method = :hybrid`, without destructive
metamorphism. The calibration is described under Fresh snow on the Physics Notes page.

## Method options

| Parameter | Options | Notes |
|---|---|---|
| `densification_method` | `:Arthern`, `:ArthernB`, `:HerronLangway`, `:Barnola1991`, **`:Crocus`**, `:CrocusPure`, `:GSFC2020`, `:Simonsen2013`, `:Ligtenberg` | Firn compaction law. FirnMICE finds the raw `:Arthern` (k₀=k₁=1) the most accumulation-sensitive of the eight schemes tested; `:Ligtenberg` and `:GSFC2020` are its recalibrated successors and are worth considering for a specific domain |
| `crocus_handoff_law` | **`:Arthern`**, `:GSFC2020` | The accumulation-driven law `:Crocus` hands firn to |
| `crocus_handoff_rule` | `:density`, **`:rate_match`** | `:density` switches every cell at or above `crocus_handoff_density`; `:rate_match` gives such a cell the faster of the two laws, so the switch falls where their rates cross. Both handoff settings are rejected with any other `densification_method` unless left at their defaults |
| `densification_accumulation` | **`:accumulation`**, `:precipitation` | Whether compaction is driven by mean snowfall or by mean total precipitation |
| `thermal_conductivity_method` | **`:Calonne2019Air`**, `:Calonne2019`, `:Calonne`, `:Sturm`, `:Marchenko2019` | RetMIP §5.1 recommends Calonne et al. (2019) by name; the CFM and IMAU-FDM both ship it. `:Calonne2019Air` is eq. 5 as published, with the air-conductivity ratio on the snow branch, as IMAU-FDM carries it; `:Calonne2019` omits that ratio, as the CFM does |
| `thermal_solver` | **`ExplicitThermal()`**, `ImplicitThermal()` | Thermal integration scheme, selected by *type* rather than by `Symbol`. `ExplicitThermal` sub-steps to the Von Neumann stability limit, so its cost is set by the stiffest cell; `ImplicitThermal` is backward Euler on a tridiagonal system, unconditionally stable, and sub-steps only for accuracy. The explicit default is faster on a well-conditioned column (measured 2.4 s vs 5.8 s over a year of 3-hourly forcing) and is what the implicit path is validated against; the implicit path is the robust choice when thin refrozen lenses collapse the stability limit |
| `thermal_explicit_safety_factor` | **`0.8`**, any value in (0, 1] | Fraction of the diffusive stability limit `ExplicitThermal` may sub-step at; unread by `ImplicitThermal`. Not slack against an imprecise limit: `_max_safe_dt` bounds diffusion only, while the surface cell also carries the surface energy balance's feedback `Λ = dQ/dT₁ ≤ 0` in the same coefficient, and nothing in the model bounds `Λ`. Over a year of 3-hourly forcing at the synthetic site the realized surface amplification factor stayed ≤ 0.955 at 0.8 but exceeded 1 on 0.65% of steps at 0.95, for ~17% fewer sub-steps. Raise it toward 1.0 to trade monotonicity margin for speed on a well-conditioned column; use `ImplicitThermal()` if the step limit should be gone rather than loosened |
| `heat_capacity_method` | **`:CuffeyPaterson`**, `:constant` | `:CuffeyPaterson` is `152.5 + 7.122·T` (Cuffey & Paterson 2010 eq. 9.1, p. 400); `:constant` (`heat_capacity_ice`, 2102 J kg⁻¹ K⁻¹) matches the CFM's melt-enabled path |
| `rain_heat_capacity` | **`:water`**, `:ice` | Heat capacity carrying rain's sensible heat above the melting point |
| `mean_temperature_method` | **`:arithmetic`**, `:arrhenius` | How the mean annual temperature driving densification is averaged |
| `grain_growth_method` | `:hybrid`, `:Arthern`, `:Marbouty`, `:Taillandier2007`, **`:FlannerZender2006`** | Dry non-dendritic grain growth. `:Marbouty` is temperature-gradient metamorphism, inactive above 400 kg m⁻³ and below 16 K m⁻¹; `:hybrid` runs it below 400 kg m⁻³ and `:Arthern` above. The CFM ships `:Arthern`. `:Taillandier2007` replaces Marbouty below 400 kg m⁻³ by measured dry-snow SSA decay (Taillandier et al., 2007, eqs. 9 and 13). `:FlannerZender2006` is the snow aging of SNICAR in CLM5.0 (Flanner and Zender, 2006; Flanner et al., 2007): the SNICAR table's dry law at each cell's temperature, gradient and density, CLM's wet growth and refrozen-water radius, for all snow below 400 kg m⁻³ and all dendritic snow, with radius as state of its own |
| `water_irreducible_method` | **`:ColeouLesaffre`**, `:constant` | Both the CFM and IMAU-FDM use Coléou & Lesaffre. RetMIP §5.4 ties a flat saturation to under-retention in the percolation zone |
| `melt_geometry` | **`:thickness`**, `:density` | What melting does to a cell: shrink `dz` at fixed density (Crocus) or lower density at fixed `dz` (SNOWPACK). Only `:density` makes a melt–refreeze cycle return the geometry as well as the mass |
| `runoff_method` | **`:instantaneous`**, `:ZuoOerlemans`, `:Darcy` | All three RetMIP bucket lineages and both comparison models use instantaneous runoff. The other two give a drainage timescale and permit firn aquifers |
| `new_snow_method` | `:Constant350`, `:Constant315`, `:Constant150`, `:Constant`, `:Fausto`, `:FaustoFit`, `:Pahaut`, `:Kaspers`, `:KuipersMunneke`, `:vanKampenhout2017`, **`:vanKampenhout2017SNOTEL`** | Fresh-snow density. The CFM ships a constant 350. `:Constant*` are bare constants, `:Constant` at `new_snow_density` (50–550 kg m⁻³, default 350); `:Fausto`/`:FaustoFit`/`:Pahaut` also select the Crocus wind-dependent fresh-grain properties. `:Pahaut` is the only alpine-seasonal-snow fit — prefer it for temperate and mid-latitude glaciers, where the polar fits run too dense. `:vanKampenhout2017` is the CLM5 deposition density (temperature plus 10 m wind, van Kampenhout et al., 2017, eqs. 2–4); pair it with `densification_method = :Crocus` (destructive metamorphism is on by default there). `:vanKampenhout2017SNOTEL` is `:vanKampenhout2017` with a temperature term calibrated on SNOTEL snowfall days; see Fresh snow in the Physics Notes |
| `blowing_snow_method` | **`:none`**, `:Crocus`, `:Keenan2021` | Wind rework of snow already on the ground. `:Crocus` is the SURFEX/Crocus `SNOWDRIFT` scheme (Vionnet et al., 2012; Lafaysse et al., 2026 eqs. 59–66). Exactly mass-conserving on its own: it compacts and fragments, it does not erode. Set `blowing_snow_sublimation` for the one term that removes mass, or `drift_rate` for the prescribed path. `:Keenan2021` is the SNOWPACK erosion–redeposition scheme (Keenan et al., 2021, eqs. 3–4) with the Crocus mobility threshold; it replaces the `:Crocus` compaction and is also mass-conserving |
| `albedo_method` | **`:GardnerSharp`**, `:BrunLefebre`, `:GreuellKonzelmann`, `:None` | |
| `emissivity_method` | **`:uniform`**, `:grain_radius_threshold`, `:grain_radius_w_threshold` | |
| `initialize_age` | **`:steady_state`**, `:zero` | |
| `output_frequency` | **`:all`**, `:daily`, `:weekly`, `:monthly`, `:last` | |

## Physical and numerical settings

| Parameter | Default | Units | Notes |
|---|---|---|---|
| `density_ice` | **917.0** | kg m⁻³ | Pure-ice density. Matches the CFM and IMAU-FDM, and the pure-ice density the Calonne (2019) conductivity and Barnola (1991) densification fits were built against |
| `water_irreducible_saturation` | **0.07** | – | Read only under `water_irreducible_method = :constant` |
| `impermeable_density` | **830.0** | kg m⁻³ | Flow-blocking criterion; the same value the CFM uses. RetMIP models span 810–917 |
| `impermeable_thickness` | **0.1** | m | Minimum blocking-lens thickness; the same value the CFM uses |
| `pore_saturation_max` | **1.0** | – | Cap on pore filling; only reachable when runoff is delayed |
| `heat_capacity_ice` | **2102.0** | J kg⁻¹ K⁻¹ | Read only under `heat_capacity_method = :constant` |
| `rain_temperature_threshold` | **273.15** | K | Rain/snow partition temperature |
| `emissivity` | **0.97** | – | |
| `emissivity_grain_radius_large` | **0.97** | – | |
| `emissivity_grain_radius_threshold` | **10.0** | mm | |
| `surface_roughness_effective_ratio` | **0.10** | – | |
| `albedo_snow` / `albedo_ice` / `albedo_fixed` | **0.85** / **0.48** / **0.85** | – | |
| `albedo_density_threshold` | **`Inf`** | kg m⁻³ | `Inf` disables the density-based snow/ice albedo switch |
| `cloud_fraction` | **0.1** | – | |
| `shortwave_subsurface_absorption` | **`false`** | – | Absorb shortwave with depth rather than at the surface only |
| `shortwave_downward_diffuse`, `solar_zenith_angle`, `cloud_optical_thickness` | **0.0** | W m⁻², degrees, – | Used only where the forcing layer of the same name is identically zero |
| `black_carbon_snow`, `black_carbon_ice` | **0.0** | ppmw | Applied to every step |
| `column_ztop` | **10.0** | m | Depth of the uniform near-surface zone |
| `column_dztop` | **0.05** | m | Cell thickness in that zone |
| `column_dzmin` / `column_dzmax` | **0.025** / **0.075** | m | Merge/split bands enforcing the fixed cell count |
| `column_depth_max` | **250.0** | m | *Ceiling* on the constructed column depth, not the depth itself |
| `column_zy` | **1.10** | – | Geometric growth ratio below `column_ztop` |
| `destructive_metamorphism` | **`true` under `:Crocus`/`:CrocusPure`, else `false`** | – | Add Anderson (1976) destructive-metamorphism settling (van Kampenhout et al., 2017, eq. 5; Schneider et al., 2022, eq. 3) to the `:Crocus`/`:CrocusPure` strain rate. Rejected with any other `densification_method` |
| `destructive_metamorphism_density` | **100.0** | kg m⁻³ | `ρ_dm`, above which destructive metamorphism tapers by `exp(−0.046·(ρ − ρ_dm))`. Bounded to [50, 450]. 100 is the CLM4 value; CLM5 uses 175 |
| `crocus_handoff_density` | **300.0** | kg m⁻³ | Density at or above which `:Crocus` hands a cell to `crocus_handoff_law`; under `crocus_handoff_rule = :rate_match`, the floor below which Crocus always governs. Bounded to [200, 830]. The Community Firn Model hands off at 450 under the `:density` rule (`CROCUS_HYBRID_DENSITY`) |
| `crocus_handoff_k0`, `crocus_handoff_k1` | **0.616702**, **0.665705** | – | Multipliers on the `:Crocus` handoff law's rate below / above `crocus_handoff_stage_density`. Bounded to (0, 10] |
| `crocus_handoff_activation_energy_shift` | **−9858.14** | J mol⁻¹ | Added to the handoff law's creep activation energy; the rate is unchanged at 258.15 K (`CROCUS_HANDOFF_REFERENCE_TEMPERATURE`). Bounded to ±60000 |
| `crocus_handoff_accumulation_exponent_scale` | **0.845621** | – | Multiplies the handoff law's accumulation exponent (1 for `:Arthern`, 0.91 / 0.644 for `:GSFC2020`); the rate is unchanged at 600 kg m⁻² yr⁻¹ (`CROCUS_HANDOFF_REFERENCE_ACCUMULATION`). Bounded to (0, 3] |
| `crocus_handoff_stage_density` | **545.319** | kg m⁻³ | The handoff law's stage-1/stage-2 transition. Bounded to [400, 830) |
| `crocus_overburden_factor` | **1.27096** | – | Multiplier on the `:Crocus` viscous strain rate, ramped linearly from 1 at 300 kg m⁻³ to its full value at 350. Bounded to (0, 10]. These six calibration factors are rejected with any other `densification_method` unless left at their defaults. The defaults of the handoff and its factors are `CROCUS_FIRN_CALIBRATION`, the firn-core calibration of the GEMB v2.0 paper; `ModelParameters(; CROCUS_PUBLISHED_HANDOFF...)` runs the published laws |
| `horizontal_strain_rate` | **0.0** | yr⁻¹ | Ice-dynamic layer thinning; 0 disables the term |
| `blowing_snow_sublimation` | **`false`** | – | Whether `blowing_snow_method = :Crocus` also sublimates the suspended fraction (Gordon et al., 2006). Off is Crocus's own default (`OSNOWDRIFT_SUBLIM`), and it is the only part of the computed scheme that removes mass |
| `drift_fetch_length` | **10.0** | m | Fetch length `L` of the `:Keenan2021` saltation flux (Keenan et al., 2021, eq. 3); the reworked mass is inversely proportional to it. Bounded to (0, 1000] |
| `drift_rate` | **0.0** | kg m⁻² yr⁻¹ | Prescribed net drift divergence, positive for erosion; 0 disables the term. A constant-rate shorthand for the `snow_drift` forcing layer, which takes precedence when the forcing carries one. Bounded to ±2000 |
| `surface_slope` | **0.0** | m m⁻¹ | Hydraulic gradient; read only under `:ZuoOerlemans`/`:Darcy` |
| `output_viscosity` | **`false`** | – | Adds a `viscosity` output layer; populated only under `:Crocus`/`:CrocusPure` |

