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
rationale is in the [`ModelParameters`](@ref) docstring.

## Method options

| Parameter | Options | Notes |
|---|---|---|
| `densification_method` | **`:Arthern`**, `:ArthernB`, `:HerronLangway`, `:Barnola1991`, `:Crocus`, `:CrocusPure`, `:GSFC2020`, `:Simonsen2013`, `:Ligtenberg` | Firn compaction law. FirnMICE finds the raw `:Arthern` (k₀=k₁=1) the most accumulation-sensitive of the eight schemes tested; `:Ligtenberg` and `:GSFC2020` are its recalibrated successors and are worth considering for a specific domain |
| `densification_accumulation` | **`:accumulation`**, `:precipitation` | Whether compaction is driven by mean snowfall or by mean total precipitation |
| `thermal_conductivity_method` | **`:Calonne2019`**, `:Calonne2019Air`, `:Calonne`, `:Sturm`, `:Marchenko2019` | RetMIP §5.1 recommends `:Calonne2019` by name; the CFM and IMAU-FDM both ship it. `:Calonne2019Air` carries the air-conductivity ratio of eq. 5 on the snow branch as IMAU-FDM does |
| `thermal_solver` | **`ExplicitThermal()`**, `ImplicitThermal()` | Thermal integration scheme, selected by *type* rather than by `Symbol`. `ExplicitThermal` sub-steps to the Von Neumann stability limit, so its cost is set by the stiffest cell; `ImplicitThermal` is backward Euler on a tridiagonal system, unconditionally stable, and sub-steps only for accuracy. The explicit default is faster on a well-conditioned column (measured 2.4 s vs 5.8 s over a year of 3-hourly forcing) and is what the implicit path is validated against; the implicit path is the robust choice when thin refrozen lenses collapse the stability limit |
| `thermal_explicit_safety_factor` | **`0.8`**, any value in (0, 1] | Fraction of the diffusive stability limit `ExplicitThermal` may sub-step at; unread by `ImplicitThermal`. Not slack against an imprecise limit: `_max_safe_dt` bounds diffusion only, while the surface cell also carries the surface energy balance's feedback `Λ = dQ/dT₁ ≤ 0` in the same coefficient, and nothing in the model bounds `Λ`. Over a year of 3-hourly forcing at the synthetic site the realized surface amplification factor stayed ≤ 0.955 at 0.8 but exceeded 1 on 0.65% of steps at 0.95, for ~17% fewer sub-steps. Raise it toward 1.0 to trade monotonicity margin for speed on a well-conditioned column; use `ImplicitThermal()` if the step limit should be gone rather than loosened |
| `heat_capacity_method` | **`:constant`**, `:CuffeyPaterson` | `:constant` (2102 J kg⁻¹ K⁻¹) matches the CFM's melt-enabled path; `:CuffeyPaterson` is `152.5 + 7.122·T` |
| `rain_heat_capacity` | **`:water`**, `:ice` | Heat capacity carrying rain's sensible heat above the melting point |
| `mean_temperature_method` | **`:arithmetic`**, `:arrhenius` | How the mean annual temperature driving densification is averaged |
| `grain_growth_method` | **`:Arthern`**, `:Marbouty`, `:hybrid` | Dry non-dendritic grain growth. The CFM ships `:Arthern`; `:Marbouty` stops dead at 400 kg m⁻³ |
| `water_irreducible_method` | **`:ColeouLesaffre`**, `:constant` | Both the CFM and IMAU-FDM use Coléou & Lesaffre; RetMIP §5.4 ties the flat 0.07 to under-retention in the percolation zone |
| `melt_geometry` | **`:thickness`**, `:density` | What melting does to a cell: shrink `dz` at fixed density (Crocus) or lower density at fixed `dz` (SNOWPACK). Only `:density` makes a melt–refreeze cycle return the geometry as well as the mass |
| `runoff_method` | **`:instantaneous`**, `:ZuoOerlemans`, `:Darcy` | All three RetMIP bucket lineages and both comparison models run instantaneous. The other two give a drainage timescale and permit firn aquifers |
| `new_snow_method` | **`:Constant350`**, `:Constant315`, `:Constant150`, `:Fausto`, `:FaustoFit`, `:Pahaut`, `:Kaspers`, `:KuipersMunneke` | Fresh-snow density. The CFM ships a constant 350. `:Constant*` are bare constants; `:Fausto`/`:FaustoFit`/`:Pahaut` also select the Crocus wind-dependent fresh-grain properties. `:Pahaut` is the only alpine-seasonal-snow fit — prefer it for temperate and mid-latitude glaciers, where the polar fits run too dense |
| `blowing_snow_method` | **`:none`**, `:Crocus` | Wind rework of snow already on the ground. `:Crocus` is the SURFEX/Crocus `SNOWDRIFT` scheme (Vionnet et al., 2012; Lafaysse et al., 2026 eqs. 59–66) — the only one of the three reference implementations that computes blowing snow at all (the CFM has none; IMAU-FDM imports it from RACMO). Exactly mass-conserving on its own: it compacts and fragments, it does not erode. Set `blowing_snow_sublimation` for the one term that removes mass, or `drift_rate` for the prescribed path |
| `albedo_method` | **`:GardnerSharp`**, `:BrunLefebre`, `:GreuellKonzelmann`, `:None` | |
| `emissivity_method` | **`:uniform`**, `:grain_radius_threshold`, `:grain_radius_w_threshold` | |
| `initialize_age` | **`:steady_state`**, `:zero` | |
| `output_frequency` | **`:all`**, `:daily`, `:weekly`, `:monthly`, `:last` | |

## Physical and numerical settings

| Parameter | Default | Units | Notes |
|---|---|---|---|
| `density_ice` | **917.0** | kg m⁻³ | Pure-ice density. Matches the CFM and IMAU-FDM, and the pure-ice density the Calonne (2019) conductivity and Barnola (1991) densification fits were built against |
| `water_irreducible_saturation` | **0.07** | – | Read only under `water_irreducible_method = :constant` |
| `impermeable_density` | **830.0** | kg m⁻³ | Flow-blocking criterion (CFM `RhoImp`). RetMIP models span 810–917 |
| `impermeable_thickness` | **0.1** | m | Minimum blocking-lens thickness (CFM `ThickImp`) |
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
| `shortwave_downward_diffuse`, `solar_zenith_angle`, `cloud_optical_thickness`, `black_carbon_snow`, `black_carbon_ice` | **0.0** | – | Overridden by the forcing when supplied |
| `column_ztop` | **10.0** | m | Depth of the uniform near-surface zone |
| `column_dztop` | **0.05** | m | Cell thickness in that zone |
| `column_dzmin` / `column_dzmax` | **0.025** / **0.075** | m | Merge/split bands enforcing the fixed cell count |
| `column_depth_max` | **250.0** | m | *Ceiling* on the constructed column depth, not the depth itself |
| `column_zy` | **1.10** | – | Geometric growth ratio below `column_ztop` |
| `horizontal_strain_rate` | **0.0** | yr⁻¹ | Ice-dynamic layer thinning; 0 disables the term |
| `blowing_snow_sublimation` | **`false`** | – | Whether `blowing_snow_method = :Crocus` also sublimates the suspended fraction (Gordon et al., 2006). Off is Crocus's own default (`OSNOWDRIFT_SUBLIM`), and it is the only part of the computed scheme that removes mass |
| `drift_rate` | **0.0** | kg m⁻² yr⁻¹ | Prescribed net drift divergence, positive for erosion; 0 disables the term. A constant-rate shorthand for the `snow_drift` forcing layer, which takes precedence when the forcing carries one. Bounded to ±2000 |
| `surface_slope` | **0.0** | m m⁻¹ | Hydraulic gradient; read only under `:ZuoOerlemans`/`:Darcy` |
| `output_viscosity` | **`false`** | – | Adds a `viscosity` output layer; populated only under `:Crocus`/`:CrocusPure` |

