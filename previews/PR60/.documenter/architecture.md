---
---

# Model Architecture {#Model-Architecture}



After climate forcing, model parameters, and the initial state of the column are defined, [`gemb`](/api#GEMB.gemb) calls `gemb_core` once per forcing timestep. Each call runs the physics modules in sequence, then the two grid controllers restore the column's invariants.

```
initialize_profile → time loop [
    forcing_step = climate_forcing[Ti=At(t)]
    state, flux   = gemb_core(state, forcing_step, mp, verbose)
    accumulate flux → store at output_frequency
] → DimStack
```


## Physics Modules {#Physics-Modules}

Called in this order, once per timestep. Each traces to a published source, cited in the docstring of the function that implements it — see [Internals: Physics](/internals_physics#Internals:-Physics).

|                          Module |                                                         Description |
| -------------------------------:| -------------------------------------------------------------------:|
|          `calculate_grain_size` |     Evolution of effective grain radius, dendricity, and sphericity |
|              `calculate_albedo` | Snow, firn, and ice albedo from grain radius, density, cloud amount |
| `calculate_shortwave_radiation` |               Vertical distribution of absorbed shortwave radiation |
|         `calculate_temperature` |    Temperature profile from energy absorption and thermal diffusion |
|        `calculate_accumulation` |                    Precipitation and deposition added to the column |
|                `calculate_melt` |               Meltwater production, percolation, refreezing, runoff |
|             `calculate_density` |                                             Snow/firn densification |
|      `apply_horizontal_strain!` |   Ice-dynamic layer thinning at constant density (inert by default) |
|        `manage_layer_thickness` |            Layer splitting and merging to maintain grid constraints |


## The vertical grid {#The-vertical-grid}

The column holds a **constant cell count** and a **constant total depth**, both fixed at initialization and enforced every timestep. `manage_layer_thickness` merges and splits cells to keep each within its thickness band and to restore the cell count; total depth is pinned separately at the end of the step. Because the grid is Lagrangian, cells advect with the firn, so `Z` is a **cell index, not a depth** — recover cell-centre depths with [`dz2z`](/api#GEMB.dz2z), or regrid onto fixed depths with [`gemb_interp`](/api#GEMB.gemb_interp).

Profile outputs are sized exactly to the column and top-justified, surface at row 1. There is no padding and no `NaN` scanning.

## Design principles {#Design-principles}

**Cited physics.** Every scheme traces to a published source, cited in the docstring of the function implementing it. Changes are judged against the literature, not against an earlier implementation.

**State as NamedTuple.** The column state is passed between timesteps as a plain `NamedTuple` of `Vector`s. The hot loop carries no `DimStack` overhead; results are packaged back into a `DimStack` on output.

**Symbols for coefficient sets.** Most options that select a coefficient set or a variant use a `Symbol` (`albedo_method=:GardnerSharp`, `output_frequency=:daily`).

**Types for algorithms.** Where an option selects a whole _algorithm_, it is a singleton type under an abstract supertype, dispatched on rather than branched on. [`ModelParameters`](/api#GEMB.ModelParameters) is parameterized on `thermal_solver` ([`ExplicitThermal`](/api#GEMB.ExplicitThermal) / [`ImplicitThermal`](/api#GEMB.ImplicitThermal) `<:`[`AbstractThermalSolver`](/api#GEMB.AbstractThermalSolver)), so selection resolves at compile time and an invalid value is unconstructible rather than caught at validation. A third scheme is a new subtype plus one method, with no central branch to edit.

**Immutable parameters.** [`ModelParameters`](/api#GEMB.ModelParameters) is immutable; construct a new instance to change a setting.

**Conservation as a test.** With `verbose=true`, `gemb_core` validates mass and energy conservation every timestep.

## Surface energy balance numerics {#Surface-energy-balance-numerics}

Three of GEMB's choices in coupling the surface energy balance to the subsurface heat equation are independently endorsed by Fourteau et al. (2024), who set out a finite-volume framework for exactly this coupling and test the alternatives against reference solutions.

**The surface flux is applied as a flux, not as a Dirichlet temperature.** GEMB imposes the net surface energy balance as a flux on the top cell and lets its temperature follow. The alternative — forcing the surface temperature and letting the flux follow — is not conservative, and Fourteau et al. quantify the cost (their Sect. 6.4 and Fig. 14): the spurious energy flux it introduces is −14.5 W m⁻² (σ = 123.5) for a melting glacier surface and +0.34 W m⁻² (σ = 39) for a seasonal snowpack, changing simulated ablation by 40% and 8% respectively. GEMB's error here is structurally zero rather than small: the solve produces a predictor for the implicit face temperatures, which is then applied as one flux per face (`+F` to one cell, `−F` to its neighbour), so pairwise cancellation is exact and the column conserves energy to the last bit independently of the solve residual. This holds on both the explicit and implicit paths.

**Face conductivity is the harmonic mean of the two cells'.** `1/(dz[i+1]/2K[i+1] + dz[i]/2K[i])` — the series resistance of the two half-cells, which is the exact flux for a piecewise-constant conductivity, where an arithmetic mean is not. This is Fourteau et al.'s eq. 6 (after Kadioglu et al., 2008), and a second independent citation for it after Lafaysse et al. (2026) eq. 90. It matters most where GEMB's grid is most graded and where conductivity contrasts are largest, which is the near-surface: a thin fresh-snow cell over dense firn.

**The surface temperature is not an independent unknown.** In Fourteau et al.'s taxonomy GEMB is a _Class 1_ model — the energy balance is applied to the top control volume, and the surface temperature is diagnosed from it (`T_surface = min(CtoK, temperature[1])`) rather than carried as its own degree of freedom. This is the same class as SNTHERM, Crocus, CLM, and CryoGrid, and it is _coupled_: the surface energy balance and the subsurface conduction are solved together, which is what Class 2 skin-layer models (COSIPY, EBFM, SnowModel) give up in exchange for an explicit surface. Their proposed scheme has both, and on their coarser meshes Class 1 diverges from the reference where the coupled-explicit-surface scheme does not. GEMB's `column_dztop = 0.05 m` sits in the regime where that divergence appears in their Figs. 9 and 11, so an explicit surface degree of freedom — a third [`AbstractThermalSolver`](/api#GEMB.AbstractThermalSolver), which would leave existing runs bit-identical — remains a real option rather than a closed question. It is not implemented.

Their Appendix D is the reason the integrated stability functions had to be continuous at neutral stability: they move their own branch point to `Ri_b = 0` so that the surface energy balance stays a well-posed function of `T_surface`, which is precisely the crossing a Newton iteration converges onto.

## Meltwater percolation scheme {#Meltwater-percolation-scheme}

GEMB percolates meltwater with a tipping-bucket scheme: water moves cell by cell from the surface down, refreezing where cold content allows, retained up to the irreducible water content ([`irreducible_saturation`](/internals_physics#GEMB.irreducible_saturation-Tuple{ModelParameters,%20Float64})), and routed to runoff at a contiguous run of cells at or above `impermeable_density` thicker than `impermeable_thickness`. There is no preferential-flow (heterogeneous, "piping") domain and no Richards-equation matrix flow, and none is planned.

That is a deliberate choice, not a missing feature. RetMIP (Vandecrux et al., 2020) intercompared nine firn models at four Greenland sites and found that the three models with explicit deep or preferential percolation (CFM-Cr, CFM-KM, UppsalaUniDeepPerc) performed worse than the bucket schemes at three of the four sites — they infiltrated water too deeply and carried a warm bias in firn temperature at the dry-snow site (Summit) and both percolation sites (Dye-2 mean error +3.6 to +6.2 °C, KAN_U +1.8 to +4.7 °C). At Dye-2 in 2016 the CFM models percolated to 10 m against 2.5 m observed by upward-looking radar, and built multi-metre near-surface ice slabs where none are observed. Their advantage was confined to the firn-aquifer site, where only the deep-percolation schemes recharged the aquifer at all. RetMIP's conclusion (their Sect. 5.2) is that until the physics of preferential flow in firn is better constrained by field and laboratory observation, the more complex schemes do not necessarily give better results than simple bucket schemes.

The levers RetMIP does identify for bucket schemes are _when_ water is blocked and _how fast_ it then leaves, not how deep it goes. Both are exposed:

**The impermeability criterion** — `impermeable_density` and `impermeable_thickness`. See [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) and the physics notes in the README for the range the participating models spanned and how it maps onto their skill at the ice-slab site.

**The runoff timescale** — `runoff_method`, with `surface_slope` as the driving hydraulic gradient. Under the `:instantaneous` default blocked water leaves the column within the timestep and no cell can ever hold more than its irreducible water. The other two let it pond into the pore space above the barrier and drain laterally over a finite timescale:

|  `runoff_method` |                                                Runoff law |                                                                                             Reference |
| ----------------:| ---------------------------------------------------------:| -----------------------------------------------------------------------------------------------------:|
| `:instantaneous` |              All blocked water leaves within the timestep |                     The default; what both comparison models and all three RetMIP bucket lineages use |
|  `:ZuoOerlemans` | `drain = excess · min(1, Δt/τ)`, `τ = c₁ + c₂·exp(−c₃·S)` |                            Zuo and Oerlemans (1996) eqs. 21–22, coefficients via Langen et al. (2017) |
|         `:Darcy` |               `drain = min(excess, ρ_w·Δt·K_sat·K_rel·S)` | Calonne et al. (2012) eq. 6 with van Genuchten (1980) / Yamaguchi et al. (2012) relative permeability |


This is where RetMIP's evidence actually points. The two models with the lowest firn-temperature error at the ice-slab site KAN_U — DMIHH (−1.6 °C) and GEUS (+0.6 °C), against a spread reaching +4.7 °C — were both bucket schemes that _delay_ runoff rather than models that percolate deeper; DMIHH uses the `:ZuoOerlemans` timescale and GEUS a Darcy flux to a virtual downslope neighbour. At the other end, DTU runs water off immediately and produced runoff unrealistic enough that RetMIP excluded it from their multi-model mean.

Delayed runoff is also what makes saturated firn representable at all: with the hard irreducible clamp of `:instantaneous`, a saturated cell cannot exist, and RetMIP (their Sect. 5.4) note that models so constrained "are incapable of modeling actual aquifers" — the firn-aquifer site was dropped from their retention evaluation for exactly this reason. Aquifers form bottom-up under either delayed method, and the `aquifer_thickness` and `aquifer_depth` outputs report the resulting water table. Both are opt-in: they are inert at the `:instantaneous` default.

## Validation {#Validation}

Validation comes from four sources, in rough order of strength:
2. **Closed-form checks** — where a scheme has an analytic solution, the test computes it independently and compares to ~1e-12.
  
3. **Published values** — coefficients and reference points taken from the paper the scheme is cited to.
  
4. **Independent implementations** — agreement with the [Community Firn Model](https://github.com/UWGlaciology/CommunityFirnModel) and [IMAU-FDM](https://github.com/IMAU-ice-and-climate/IMAU-FDM) where the algebra is the same.
  
5. **Invariants** — conservation of mass and energy, monotonicity, bounds, and the physical limits of each scheme.
  
