---
---

# Internals: Grid and Column {#Internals:-Grid-and-Column}

Non-exported grid controllers, layer management, and column initialization. See [Internals](/internals#Internals) for the caveat on depending on these.
<details class='jldocstring custom-block' open>
<summary><a id='GEMB._select_merge_pair-Tuple{Vector{Float64}, Vector{Float64}, Int64, ModelParameters}' href='#GEMB._select_merge_pair-Tuple{Vector{Float64}, Vector{Float64}, Int64, ModelParameters}'><span class="jlbinding">GEMB._select_merge_pair</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_select_merge_pair(dz, dzmax, n, mp) -> Union{Int,Nothing}
```


Choose `i` such that cells `i` and `i+1` are merged. Returns `nothing` when no merge is available (`n < 3`, or only the surface cell would qualify).

Candidates are ranked by two independent preferences, in this order:
2. **In band** — the combined thickness must still fit the lower cell's `dzmax`, so the merge is not immediately undone by next timestep's split test.
  
3. **Deep** — at or below `column_ztop`, leaving the fine near-surface grid (where the radiation and turbulent-flux gradients live) alone where possible.
  

Within each tier the thinnest pair wins. **In-band must dominate deep** — this ordering is load-bearing. The deep region is a finite budget of cells while accumulation demands a merge thousands of times per year, so a deep-first rule eventually has only the bottom pair to offer and returns it every time, folding the deep column into one monolithic cell (measured: 245 m of a 254.6 m column). In-band first self-corrects, since an over-thick bottom cell is by construction out of band.

See also [`_select_split_cell`](/internals_grid#GEMB._select_split_cell-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._select_split_cell-Tuple{Vector{Float64}, Vector{Float64}, Int64, ModelParameters}' href='#GEMB._select_split_cell-Tuple{Vector{Float64}, Vector{Float64}, Int64, ModelParameters}'><span class="jlbinding">GEMB._select_split_cell</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_select_split_cell(dz, dzmin, n, mp) -> Union{Int,Nothing}
```


Choose the cell to split, or `nothing` when none qualifies. Ranked by the same two-tier scheme as [`_select_merge_pair`](/internals_grid#GEMB._select_merge_pair-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters}) (see the note there on why in-band must dominate deep): both halves clearing `dzmin` outranks being below `column_ztop`, and the thickest candidate wins within a tier.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.apply_horizontal_strain!-Tuple{NamedTuple, Float64, ModelParameters}' href='#GEMB.apply_horizontal_strain!-Tuple{NamedTuple, Float64, ModelParameters}'><span class="jlbinding">GEMB.apply_horizontal_strain!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
apply_horizontal_strain!(cols, dt_seconds, mp) -> (mass_lateral, E_lateral)
```


Thin (or thicken) every cell by the horizontal strain rate `mp.horizontal_strain_rate` [yr-1] over `dt_seconds`, at constant density.

Incompressibility ties vertical thinning to horizontal stretching: for a parcel of thickness `h` under a horizontal velocity-gradient trace `D = ε̇_xx + ε̇_yy`, `dh/dt = -h D` at fixed density, whose exact integral over a step of constant `D` is `h *= exp(-D dt)`. The exponential form keeps `dz` positive for any `D dt`, where forward Euler would not — the same reason the `:Crocus` branch of [`calculate_density`](/internals_physics#GEMB.calculate_density-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}) integrates its viscous law in closed form.

`density` and `temperature` are untouched; grain properties are intensive and unchanged. Pore `water` [kg m-2] is a per-unit-area quantity, so it scales by the same factor as the lateral area.

`age` is untouched in both directions. Under divergence every cell loses the same _fraction_ of its mass, which leaves its mean age unchanged; under convergence the material advected in laterally is neighbouring firn at the same depth, so taking the local `age[i]` is the consistent choice.

GEMB's column is per unit area, so the mass this removes leaves **laterally**, not through the base. It is returned signed — negative under divergence (mass leaves), positive under convergence — and is reported to the mass budget separately from [`trim_bottom!`](/internals_grid#GEMB.trim_bottom!-Tuple{NamedTuple,%20Float64,%20ModelParameters})'s basal flux. Folding the two together would make a thinning column read as spurious basal accretion, because `trim_bottom!` refills the column to `z_target` immediately afterwards.

Returns `(0.0, 0.0)` and touches nothing when `horizontal_strain_rate == 0`, which is what keeps the default bit-identical to a run without the term.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.close_slot!-Tuple{NamedTuple, Int64}' href='#GEMB.close_slot!-Tuple{NamedTuple, Int64}'><span class="jlbinding">GEMB.close_slot!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
close_slot!(cols, i)
close_slot!(cols, idx::AbstractVector{Int})
```


Delete cell `i` (or the sorted, unique indices `idx`) across every field. Cells below shift up to close the gap. The vector form deletes in one pass per field, so it is O(n) in total rather than O(n) per index.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.column_bands!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}' href='#GEMB.column_bands!-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters}'><span class="jlbinding">GEMB.column_bands!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
column_bands!(dzmin, dzmax, dz, mp) -> (dzmin, dzmax)
```


Fill the per-cell thickness bands `[dzmin[i], dzmax[i]]` that drive merging and splitting.

Cells within `mp.column_ztop` of the surface get the uniform bands `[mp.column_dzmin, mp.column_dzmax]`. Below that the bands are scaled by `mp.column_zy^k`, where `k` is the geometric grid position **implied by the cell's depth** — that is, the solution of

```julia
z = column_ztop + column_dztop · zy · (zy^k − 1) / (zy − 1)
```


for the cell's lower face depth `z`. This is exactly the stretching [`initialize_grid`](/internals_grid#GEMB.initialize_grid-Tuple{ModelParameters}) uses to build the column, so a freshly-built grid has every cell sitting at the geometric centre of its own band, and a cell is always measured against the thickness appropriate to the depth it has actually reached.

**Why the band must be referenced to depth, not to cell index**

Do not replace this with the obvious version — an accumulator multiplied by `zy` once per cell below `column_ztop` — which indexes the band by _position in the deep sequence_. That is unstable under compaction: as near-surface cells thin, more fit inside `column_ztop`, the deep sequence shortens, and every deep band collapses though no cell has moved. On the default grid `Σdzmax` fell from 382 m to under 50 m within a simulated year; the split pass then shattered the deep column into slivers, the count controller merged them back, and the column degenerated to one cell holding 245.5 m of its 254.6 m by spinup cycle 5. Depth-referencing removes the feedback: a cell's band depends only on where it is, making the geometric grid a fixed point of the merge/split rules.

`dzmin` and `dzmax` must be at least `length(dz)` long; only the first `length(dz)` entries are written.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.column_mass-Tuple{NamedTuple}' href='#GEMB.column_mass-Tuple{NamedTuple}'><span class="jlbinding">GEMB.column_mass</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
column_mass(cols) -> mass [kg m-2]
```


Total column mass: ice/firn `dz*density` plus pore water. Needs no `ModelParameters`, unlike [`column_mass_energy`](/internals_grid#GEMB.column_mass_energy-Tuple{NamedTuple,%20ModelParameters}), so callers that want mass alone can avoid threading one through.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.column_mass_energy-Tuple{NamedTuple, ModelParameters}' href='#GEMB.column_mass_energy-Tuple{NamedTuple, ModelParameters}'><span class="jlbinding">GEMB.column_mass_energy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
column_mass_energy(cols, mp) -> (mass, energy)
```


Total column mass [kg m-2] and enthalpy [J m-2], the reference quantities for every conservation check. Ice/firn contributes `dz*density*h(T)` with `h` the enthalpy integral [`specific_enthalpy`](/internals_physics#GEMB.specific_enthalpy-Tuple{ModelParameters,%20Real}); pore water contributes [`specific_enthalpy_water`](/internals_physics#GEMB.specific_enthalpy_water-Tuple{ModelParameters,%20Real}) — its latent heat plus sensible heat at the melt point.

Absolute energies are measured from a 0 K reference and so are not comparable across `mp.heat_capacity_method`; only differences are meaningful.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.column_state-NTuple{8, Any}' href='#GEMB.column_state-NTuple{8, Any}'><span class="jlbinding">GEMB.column_state</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
column_state(temperature, dz, density, water, grain_radius,
             grain_dendricity, grain_sphericity, age)
```


Bundle the eight per-cell state vectors so the grid primitives can operate on all of them generically. `values(...)` of the result is a homogeneous `NTuple{8,Vector{Float64}}`, so the generic loops below stay type-stable and fully unrolled.

Albedo is deliberately absent: every albedo method is diagnostic in the current column state (see [`calculate_albedo`](/internals_physics#GEMB.calculate_albedo-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ClimateForcingStep,%20ModelParameters})), so there is no per-cell albedo to carry through a merge or a split.

The structural primitives (`open_slot!`, `close_slot!`) need no change when a field is added; `merge_pair!` and `split_cell!` do, since a new field's extensive-vs-intensive merge semantics cannot be inferred. See the module preamble on why the bundle is rebuilt per call site.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.dilute_age-Tuple{Float64, Float64, Float64}' href='#GEMB.dilute_age-Tuple{Float64, Float64, Float64}'><span class="jlbinding">GEMB.dilute_age</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dilute_age(age, mass_old, mass_added) -> age
```


Mass-weighted mean age of a cell that gains `mass_added` [kg m-2] of brand-new (age 0) material, where `mass_old` is the cell's total mass beforehand. The [`mix_age`](/internals_grid#GEMB.mix_age-NTuple{4,%20Float64}) case where one side is exactly zero.

This is the update for every site where mass enters the column from outside: snowfall, rain, and vapour deposition. Returns `age` unchanged for non-positive `mass_added`, so the caller need not branch on the sign of a signed flux (sublimation removes a fraction of the cell and is age-neutral).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.enforce_column_length!-Tuple{NamedTuple, Int64, ModelParameters}' href='#GEMB.enforce_column_length!-Tuple{NamedTuple, Int64, ModelParameters}'><span class="jlbinding">GEMB.enforce_column_length!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
enforce_column_length!(cols, n_target, mp) -> Int
```


Controller 1 (count). Restore the column to exactly `n_target` cells by merging or splitting in the **deep** part of the column, and return the number of merge/split operations performed.

Both operations conserve mass and energy exactly, so this returns no budget terms — it is invisible to the model's mass and energy accounting.

Which cells are touched:
- **Too many cells** → repeatedly merge an adjacent pair, thinnest first ([`_select_merge_pair`](/internals_grid#GEMB._select_merge_pair-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters})).
  
- **Too few cells** → repeatedly split a cell, thickest first ([`_select_split_cell`](/internals_grid#GEMB._select_split_cell-Tuple{Vector{Float64},%20Vector{Float64},%20Int64,%20ModelParameters})).
  

Both selectors rank candidates by staying inside the cell's `dzmin`/`dzmax` band and by lying below `mp.column_ztop`, with **in-band dominating deep** (see `_select_merge_pair` for why), breaking ties toward the deeper candidate. Neither touches cell 1, so the surface cell is never consumed.

Coarsening at depth while fine cells arrive at the top is the history a descending parcel experiences, so this does not fight the grid — [`initialize_grid`](/internals_grid#GEMB.initialize_grid-Tuple{ModelParameters})'s geometric profile is a fixed point of the process. Over a 75-cycle spinup on the default grid the deep structure holds (max cell 22.3 → 30.2 m, cells over 5 m steady at 16) with no cell out of band. Ranking in-band first does let a merge land inside `column_ztop` when no deep pair is admissible: that is the safety valve preventing the deep column from being consumed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.merge_pair!-Tuple{NamedTuple, Int64, Int64, Float64, Float64, ModelParameters}' href='#GEMB.merge_pair!-Tuple{NamedTuple, Int64, Int64, Float64, Float64, ModelParameters}'><span class="jlbinding">GEMB.merge_pair!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
merge_pair!(cols, i, i_target, M_i, M_target, mp) -> M_new
```


Fold the contents of cell `i` into cell `i_target`, conserving mass and energy exactly. `M_i` and `M_target` are the ice/firn masses `dz*density` of the two cells; the combined mass is returned so callers tracking a mass vector can keep it in sync.

Temperature is combined by [`mix_temperature`](/internals_physics#GEMB.mix_temperature-Tuple{ModelParameters,%20Vararg{Real,%204}}), which mixes in enthalpy and inverts. Under a temperature-dependent heat capacity the mass-weighted mean temperature is _not_ energy-conserving — it loses `(b/2)·M·Var_M(T)` joules, thousands of times the conservation tolerance for a realistic deep merge. Thickness and pore water are extensive and are summed, with density recovered from the totals. Grain properties are inherited from cell `i` — that is the historical convention (the _upper_ cell of the pair wins), preserved deliberately.

`age` is mass-weighted on **total** cell mass (`dz*density + water`), matching how the field is defined, so `M_i`/`M_target` alone are not the right weights — the pore water of each cell is added in. Deliberately _not_ the grain-property convention above: inheriting the upper cell's age would discard the older cell's residence time entirely and make deep age drift young with every merge.

This performs only the physics of the merge. The now-redundant cell `i` still occupies a slot; the caller removes it with [`close_slot!`](/internals_grid#GEMB.close_slot!-Tuple{NamedTuple,%20Int64}), either immediately or in a batch.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.mix_age-NTuple{4, Float64}' href='#GEMB.mix_age-NTuple{4, Float64}'><span class="jlbinding">GEMB.mix_age</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
mix_age(age1, M1, age2, M2) -> age
```


Mass-weighted mean age of two masses combined: `(age1·M1 + age2·M2)/(M1+M2)`.

Age is intensive and mixes linearly in mass, so this is the whole rule — the analogue of [`mix_temperature`](/internals_physics#GEMB.mix_temperature-Tuple{ModelParameters,%20Vararg{Real,%204}}) for a quantity with no equation of state. Every site where age changes is one of its three cases: mass arriving from outside at age 0 ([`dilute_age`](/internals_grid#GEMB.dilute_age-Tuple{Float64,%20Float64,%20Float64})), two cells combining ([`merge_pair!`](/internals_grid#GEMB.merge_pair!-Tuple{NamedTuple,%20Int64,%20Int64,%20Float64,%20Float64,%20ModelParameters})), and meltwater joining through-flow in `calculate_melt`. Keeping them one function keeps the zero-mass convention one decision.

Returns `age1` when the total mass is zero, matching [`mix_temperature`](/internals_physics#GEMB.mix_temperature-Tuple{ModelParameters,%20Vararg{Real,%204}}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.open_slot!-Tuple{NamedTuple, Int64}' href='#GEMB.open_slot!-Tuple{NamedTuple, Int64}'><span class="jlbinding">GEMB.open_slot!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
open_slot!(cols, i)
```


Insert a cell at index `i`, duplicating the current contents of cell `i` into the new slot, across every field. Cells `i..n` shift down to `i+1..n+1`.

Used to prepend a fresh-snow cell (`i = 1`) and as the structural half of [`split_cell!`](/internals_grid#GEMB.split_cell!-Tuple{NamedTuple,%20Int64}).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.split_cell!-Tuple{NamedTuple, Int64}' href='#GEMB.split_cell!-Tuple{NamedTuple, Int64}'><span class="jlbinding">GEMB.split_cell!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
split_cell!(cols, i)
```


Split cell `i` into two cells of half its thickness, conserving mass and energy exactly.

Thickness and pore water are halved (both extensive); density, temperature, grain properties and age are intensive and so are simply duplicated. The result occupies indices `i` and `i+1`, and the column grows by one cell.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.trim_bottom!-Tuple{NamedTuple, Float64, ModelParameters}' href='#GEMB.trim_bottom!-Tuple{NamedTuple, Float64, ModelParameters}'><span class="jlbinding">GEMB.trim_bottom!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trim_bottom!(cols, z_target, mp) -> (mass_added, E_added)
```


Controller 2 (mass). Pin the total column depth to `z_target` by adjusting the bottom cell's thickness, and return the mass [kg m-2] and energy [J m-2] added to the column as a result. This is the model's **only** basal mass/energy flux (horizontal strain, see [`apply_horizontal_strain!`](/internals_grid#GEMB.apply_horizontal_strain!-Tuple{NamedTuple,%20Float64,%20ModelParameters}), is the only lateral one).

The adjustment is signed and continuous:
- `Σdz > z_target` (the accumulation regime — mass added at the surface has pushed the column down): the bottom cell shrinks and material passes out through the base. `mass_added` is negative.
  
- `Σdz < z_target` (the ablation regime — mass removed at the surface, or densification has compacted the column): the bottom cell grows and material enters through the base. `mass_added` is positive.
  

Accreted material inherits the bottom cell's properties (the column below is unresolved): `density[n]`, and `temperature[n]`, which is the Dirichlet-pinned basal boundary. Pore water scales with the thickness change; grain properties are intensive and unchanged.

`age[n]` is untouched, so both regimes are handled by doing nothing: removal is proportional (age-neutral), and accreted material inherits `age[n]` like `density` and `temperature`. The consequence is that under sustained ablation the deepest cell's age is a lower bound on its true residence time rather than a measurement of it.

Errors if the adjustment would consume the whole bottom cell — at realistic forcing it is ~1e-3 of the cell, so that signals a bug upstream, not an extreme climate.

The column is an Eulerian window on the firn, not a prognostic ice-thickness model, so `thickness_cumulative` measures basal flux in either regime, not thickness change. For the latter, use the surface terms (`precipitation - runoff + evaporation_condensation`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.manage_layer_thickness-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters, Bool}' href='#GEMB.manage_layer_thickness-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, ModelParameters, Bool}'><span class="jlbinding">GEMB.manage_layer_thickness</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
manage_layer_thickness(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, mp::ModelParameters, verbose::Bool; n_target=length(dz))
```


Return every cell to its per-cell thickness band, then restore the column to exactly `n_target` cells.

Three passes, in order:
2. **Merge** cells thinner than their per-cell `dzmin` band into the cell below.
  
3. **Split** cells thicker than their per-cell `dzmax` band in half.
  
4. **Count control** ([`enforce_column_length!`](/internals_grid#GEMB.enforce_column_length!-Tuple{NamedTuple,%20Int64,%20ModelParameters})) — merge or split in the deep column until the cell count is exactly `n_target`, absorbing whatever net change the two passes above and the surface physics (accumulation, melt-out) produced.
  

All three passes are built from the primitives in `grid_ops.jl` and **conserve mass and energy exactly**, so this function no longer reports a basal mass/energy flux. Column depth is pinned separately, and later in the timestep, by [`trim_bottom!`](/internals_grid#GEMB.trim_bottom!-Tuple{NamedTuple,%20Float64,%20ModelParameters}) — see `grid_ops.jl` for why count and mass are controlled independently.

The Dirichlet temperature boundary condition is re-imposed at the bottom; the energy that requires is the only term returned in `E_added`.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, E_added)`. Arrays are mutated in place and their length on return is `n_target`.

`age` needs no handling here beyond being in the bundle: all three passes go through [`merge_pair!`](/internals_grid#GEMB.merge_pair!-Tuple{NamedTuple,%20Int64,%20Int64,%20Float64,%20Float64,%20ModelParameters}) (which mass-weights it) and [`split_cell!`](/internals_grid#GEMB.split_cell!-Tuple{NamedTuple,%20Int64}) (for which duplication is already correct, age being intensive). The merge pass keeps `M` in sync across a chain of cells folding into one target, so the weighting stays correct over a chain rather than only over a pair.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._derive_column_depth-Tuple{ModelParameters, ClimateSummary}' href='#GEMB._derive_column_depth-Tuple{ModelParameters, ClimateSummary}'><span class="jlbinding">GEMB._derive_column_depth</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_derive_column_depth(mp::ModelParameters, cs::ClimateSummary;
                     margin=1.2, thermal_depths=4)
    -> depth [m]
```


Total column depth this climate needs, replacing the old binary 250-vs-25 m switch with a continuous derivation. Used to build the grid in [`initialize_profile`](/api#GEMB.initialize_profile); it goes no further than that, since the grid is what carries the depth forward.

The configured `mp.column_depth_max` is treated as a **ceiling** rather than a target, so a shallow-firn or ablation site is not forced to carry hundreds of metres of solid ice, while a deep cold site is unaffected. Two requirements set the depth:
2. **Resolve the firn.** Below the depth where the steady-state profile reaches `mp.density_ice` the column carries no further firn information, so `margin × ss.ice_depth` is deep enough — the margin leaves room for the real column to densify more slowly than the guess.
  
3. **Resolve the annual thermal wave.** Even a column of pure ice must be deep enough that its lower boundary does not clamp the seasonal temperature cycle. `thermal_depths` damping depths (`d = sqrt(2K/(ρ·c_p·ω))` at ice density, ≈3.3 m, so ≈13 m) is the floor. This is what keeps an ablation site — where `ice_depth = 0` because nothing is ever buried — from collapsing to a physically useless grid, and it recovers the old 25 m intent from the thermal physics rather than a constant.
  

The result is `clamp(max(firn requirement, thermal floor), mp.column_ztop, mp.column_depth_max)`. Depths are found by marching on a coarse 1 m probe grid, which is all the resolution two scalars need and keeps this independent of the final grid it is used to build.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._uniform_ice_profile-Tuple{ModelParameters, Real}' href='#GEMB._uniform_ice_profile-Tuple{ModelParameters, Real}'><span class="jlbinding">GEMB._uniform_ice_profile</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_uniform_ice_profile(mp::ModelParameters, T_mean) -> DimStack
```


The MATLAB `model_initialize_profile` column: pure ice at a uniform temperature, with non-dendritic faceted ice grains (`grain_radius = 2.5` mm, the non-spherical cap in [`calculate_grain_size`](/internals_physics#GEMB.calculate_grain_size-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters})) and no pore water.

Reached only when both fidelity flags of [`initialize_profile`](/api#GEMB.initialize_profile) are set. Kept as a separate function taking no climate-derived input so no future change to the steady-state scheme can perturb it.

`T_mean` is expected already clamped to `CtoK` by the caller — the clamp lives there so it covers both fidelity paths in one place.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.initialize_grid-Tuple{ModelParameters}' href='#GEMB.initialize_grid-Tuple{ModelParameters}'><span class="jlbinding">GEMB.initialize_grid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
initialize_grid(mp::ModelParameters, depth=mp.column_depth_max)
```


Generate the initial vertical grid layer thicknesses, spanning `depth`. Matches MATLAB's `model_initialize_grid` (local function in model_initialize_profile.m).

`depth` is an argument rather than read from `mp` because [`initialize_profile`](/api#GEMB.initialize_profile) builds the grid on the climate-derived depth ([`_derive_column_depth`](/internals_grid#GEMB._derive_column_depth-Tuple{ModelParameters,%20ClimateSummary})), for which `mp.column_depth_max` is only the ceiling. It is the depth the grid actually spans, hence `depth` and not `column_depth_max`.

Returns a Vector{Float64} of layer thicknesses from surface to depth.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._albedo_residual-Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Float64, ModelParameters}' href='#GEMB._albedo_residual-Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, Float64, Float64, ModelParameters}'><span class="jlbinding">GEMB._albedo_residual</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_albedo_residual(α, ...) -> (g, melt, refreeze)
```


One evaluation of the albedo fixed point: the melt an albedo of `α` produces, the refreeze that implies, and the residual `g = α_target − α` whose root the secant iteration in [`initialize_climate_summary`](/api#GEMB.initialize_climate_summary) seeks.

Returns the melt and refreeze alongside the residual so the caller keeps the pair consistent with the accepted albedo without re-running a forcing sweep.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._cold_content_refreeze-Tuple{Real, Real, Real, Real, ModelParameters}' href='#GEMB._cold_content_refreeze-Tuple{Real, Real, Real, Real, ModelParameters}'><span class="jlbinding">GEMB._cold_content_refreeze</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_cold_content_refreeze(melt, rainfall, accumulation, T_mean, mp) -> R [kg m-2 yr-1]
```


Annual refreeze, limited by both the available liquid water and the cold content of one annual accumulation layer:

```julia
R = min(melt + rainfall, accumulation·(h(273.15) − h(T̄)) / LF)
```


The second term is the standard Pfeffer-style retention capacity: the energy needed to warm one year's accumulation to the melt point, expressed as a mass of refrozen water. Under `:constant` the enthalpy difference is `c_p·(273.15 − T̄)`, the reference's form. A temperate column (`T̄ ≥ 273.15`) has no cold content and refreezes nothing.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._seb_annual_melt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, ModelParameters}' href='#GEMB._seb_annual_melt-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64, ModelParameters}'><span class="jlbinding">GEMB._seb_annual_melt</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_seb_annual_melt(T_air, wind, pressure, sw, lw, vapor, albedo, dt, per_year,
                 zT_obs, zW_obs, mp) -> melt [kg m-2 yr-1]
```


Annual melt from a zero-layer surface energy balance over the forcing series.

For each step, [`_seb_skin_temperature`](/internals_grid#GEMB._seb_skin_temperature-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}) finds the skin temperature closing the flux balance. If that exceeds the melt point the surface is pinned at `273.15 K` and the leftover energy flux `Q` is converted to melt via `Q·Δt/LF`. Sublimation/condensation is deliberately ignored — it is a small term for an initial guess, and including it would need the mass feedback the marcher resolves anyway.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._seb_forcing_step-NTuple{9, Float64}' href='#GEMB._seb_forcing_step-NTuple{9, Float64}'><span class="jlbinding">GEMB._seb_forcing_step</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_seb_forcing_step(dt, T_air, pressure, wind, sw, lw, vapor, zT_obs, zW_obs)
```


Minimal [`ClimateForcingStep`](/api#GEMB.ClimateForcingStep) for the surface-energy-balance estimate. Only the fields [`turbulent_heat_flux`](/internals_physics#GEMB.turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}) reads are meaningful; the climatological-mean and time-varying-parameter slots are unused here and set to zero.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._seb_residual-Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}' href='#GEMB._seb_residual-Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}'><span class="jlbinding">GEMB._seb_residual</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_seb_residual(T_surface, sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs) -> Q [W m-2]
```


Net energy flux into the surface at skin temperature `T_surface`:

```julia
Q = SW_net + LW↓ − εσT⁴ + QH + QE
```


Turbulent terms come from [`turbulent_heat_flux`](/internals_physics#GEMB.turbulent_heat_flux-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}). Positive `Q` warms (or melts) the surface.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._seb_skin_temperature-Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}' href='#GEMB._seb_skin_temperature-Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, ClimateForcingStep}'><span class="jlbinding">GEMB._seb_skin_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_seb_skin_temperature(sw_net, lw_in, ε, density_air, z0, zT, zQ, cfs) -> T [K]
```


Skin temperature closing the surface energy balance, by secant iteration on [`_seb_residual`](/internals_grid#GEMB._seb_residual-Tuple{Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ClimateForcingStep}) bracketed around the air temperature.

The balance is monotone decreasing in `T_surface` (both the `−εσT⁴` radiative term and the turbulent terms oppose warming), so the iteration is well behaved. A few iterations suffice for an initial guess; the result is only used to decide whether the surface reaches melting and by how much energy it overshoots.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._snow_cover_albedo-Tuple{Real, Real, ModelParameters}' href='#GEMB._snow_cover_albedo-Tuple{Real, Real, ModelParameters}'><span class="jlbinding">GEMB._snow_cover_albedo</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_snow_cover_albedo(accumulation_effective, melt, mp) -> α
```


Annual-mean surface albedo implied by an accumulation and a melt rate, blending `mp.albedo_snow` and `mp.albedo_ice` by the fraction of the year the surface is expected to carry snow.

The blend is `f = A_eff / (A_eff + melt)`: with no melt the surface is snow all year; with melt far exceeding accumulation it is bare ice. This is deliberately **continuous**. A binary pick on the sign of the mass balance would reintroduce exactly the threshold this scheme exists to remove — one level up, in the albedo that determines the melt that determines the balance — and would make the whole initialization jump discontinuously as a site crosses `b = 0`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._snow_fraction-Tuple{Real, Real}' href='#GEMB._snow_fraction-Tuple{Real, Real}'><span class="jlbinding">GEMB._snow_fraction</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_snow_fraction(accumulation_effective, melt) -> f
```


Fraction of the year the surface is expected to carry snow, `f = A_eff/(A_eff + melt)`: with no melt the surface is snow all year; with melt far exceeding accumulation it is bare ice.

This is the single definition of the snow/ice mixing weight, used by [`_snow_cover_albedo`](/internals_grid#GEMB._snow_cover_albedo-Tuple{Real,%20Real,%20ModelParameters}) to blend the two albedos.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._snow_fraction_from_albedo-Tuple{Real, ModelParameters}' href='#GEMB._snow_fraction_from_albedo-Tuple{Real, ModelParameters}'><span class="jlbinding">GEMB._snow_fraction_from_albedo</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_snow_fraction_from_albedo(α, mp) -> f
```


Recover the snow fraction from a blended albedo — the exact inverse of [`_snow_cover_albedo`](/internals_grid#GEMB._snow_cover_albedo-Tuple{Real,%20Real,%20ModelParameters}), kept adjacent to it so the pair is maintained together.

[`_seb_annual_melt`](/internals_grid#GEMB._seb_annual_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20Float64,%20Float64,%20Float64,%20Float64,%20ModelParameters}) needs the fraction to interpolate surface roughness, but it is called _with_ a candidate albedo and _before_ the melt that would determine the fraction directly, so inverting is the only order that closes. **If the albedo blend ever stops being linear in `f`, this inverse must change with it** — hence one named function rather than the arithmetic inlined at the use site.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._grain_forcing_step-Tuple{Float64, ClimateSummary}' href='#GEMB._grain_forcing_step-Tuple{Float64, ClimateSummary}'><span class="jlbinding">GEMB._grain_forcing_step</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_grain_forcing_step(dt_seconds, cs::ClimateSummary) -> ClimateForcingStep
```


Minimal `ClimateForcingStep` for the age march's [`calculate_grain_size`](/internals_physics#GEMB.calculate_grain_size-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}) call. Grain growth reads only `dt` from the forcing step (metamorphism depends on the column's own temperature, density and water), so the remaining fields carry the climatological means and are not otherwise used.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._interp_curve-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}' href='#GEMB._interp_curve-Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}'><span class="jlbinding">GEMB._interp_curve</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_interp_curve(values, z_curve, depth; lo, hi) -> Vector{Float64}
```


Linearly interpolate an age-marched curve onto `depth`, with constant extrapolation beyond either end and the result clamped to `[lo, hi]`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._irreducible_water-Tuple{Real, Real, Real, ModelParameters}' href='#GEMB._irreducible_water-Tuple{Real, Real, Real, ModelParameters}'><span class="jlbinding">GEMB._irreducible_water</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_irreducible_water(T, ρ, dz, mp) -> water [kg m-2]
```


Irreducible (capillary-held) water content of a cell at the melt point, zero for a cold cell. Uses the same expression as [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}), `(ρᵢ − ρ)·S·(M/ρ)` with `M = ρ·dz`, so initialization and runtime agree on what "irreducible" means and the first timestep has no water to redistribute.

Two gates are added here, because this seeds a column rather than routing water through an existing one: a cold cell holds nothing, and a cell at or past pore close-off (`DENSITY_PORE_CLOSEOFF`) has no _connected_ pore space to hold water in, so it starts dry rather than with the small amount `(ρᵢ − ρ)` alone would imply. Under `water_irreducible_method = :ColeouLesaffre` [`calculate_melt`](/internals_physics#GEMB.calculate_melt-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Float64,%20ModelParameters,%20Bool}) applies the close-off gate too; under `:constant` it does not.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._march_density_step-Tuple{Any, Vararg{Float64, 6}}' href='#GEMB._march_density_step-Tuple{Any, Vararg{Float64, 6}}'><span class="jlbinding">GEMB._march_density_step</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_march_density_step(p::DensificationCoeffs, ρ, z, b, ρi, dt, T) -> (z, ρ)
```


Advance one step of a steady-state density march by the age step `dt` [yr]:

```julia
ρ  = ρᵢ − (ρᵢ − ρ)·exp(−c·dt)   (Arthern et al. 2010 eq. 1)
z += b·dt/ρ_mid                 (mass balance, trapezoidal in ρ)
```


This is the single definition of the relaxation-and-burial arithmetic, shared by [`steady_state_profile`](/api#GEMB.steady_state_profile) and [`steady_state_density`](/api#GEMB.steady_state_density) so the density law cannot drift between them.

The two differ only in how they choose `dt`, which is deliberate and stays with each caller: `steady_state_profile` derives it from a target _depth_ advance (see [`_march_step_depth`](/internals_grid#GEMB._march_step_depth-Tuple{Real,%20Integer})), while `steady_state_density` uses a fixed age step to reproduce the classic Herron & Langway discretization bit-for-bit.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._march_step_depth-Tuple{Real, Integer}' href='#GEMB._march_step_depth-Tuple{Real, Integer}'><span class="jlbinding">GEMB._march_step_depth</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_march_step_depth(max_depth, n_age) -> dz_step [m]
_march_age_step(dz_step, ρ, b) -> dt [yr]
_march_max_iterations(n_age) -> Int
```


Step-size policy for [`steady_state_profile`](/api#GEMB.steady_state_profile)'s march, which is stepped in **depth** rather than age: `dz_step` is the target advance per step and the age step follows from the local density, `dt = dz_step·ρ/b` (the mass balance `dz = b·dt/ρ`, inverted).

This keeps resolution uniform in the coordinate the result is sampled on, and makes the step count independent of the burial rate — with a fixed age step, a low-accumulation site needed orders of magnitude more steps to reach the same depth, hit the iteration cap partway down, and had the rest of its column filled by extrapolating straight to ice.

Using the current density for `dt` is first-order in `dz_step`, the same order as the density update itself. `_march_max_iterations` is a backstop against a pathological rate law, not the normal exit: the depth bound is what terminates a healthy march.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._steady_state_temperature-Tuple{Real, Real, ClimateSummary, ModelParameters}' href='#GEMB._steady_state_temperature-Tuple{Real, Real, ClimateSummary, ModelParameters}'><span class="jlbinding">GEMB._steady_state_temperature</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_steady_state_temperature(z, ρ, cs::ClimateSummary, mp::ModelParameters) -> T [K]
```


Temperature at depth `z` [m, positive down] for a column in thermal steady state: a damped annual wave superposed on the mean, warmed by the latent heat that refreezing releases.

```julia
T(z) = T̄ + ΔT_latent + A_T·exp(−z/d)·cos(φ − z/d)
```


The damping depth `d = sqrt(2K/(ρ·c_p·ω))` uses the local density's own thermal conductivity from [`thermal_conductivity`](/internals_physics#GEMB.thermal_conductivity-Tuple{AbstractVector,%20AbstractVector,%20ModelParameters}), so the wave decays through a physically consistent column rather than a nominal one. `ω = 2π/yr`.

`ΔT_latent = LF·R/(c_p·A_eff)` is self-limiting: at the cold-content cap it equals `(273.15 − T̄)·accumulation/A_eff`, which is strictly less than the gap to the melt point. The `min(·, CtoK)` clamp is therefore belt-and-braces, and also guards the ablation case where the mean itself is above freezing.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.thermal_damping_depth-Tuple{Real, Real, ModelParameters}' href='#GEMB.thermal_damping_depth-Tuple{Real, Real, ModelParameters}'><span class="jlbinding">GEMB.thermal_damping_depth</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
thermal_damping_depth(T, ρ, mp) -> d [m]
```


E-folding depth of the annual thermal wave in a medium of density `ρ` at temperature `T`:

```julia
d = sqrt(2K/(ρ·c_p·ω)),   ω = 2π/yr
```


with `K` from [`thermal_conductivity`](/internals_physics#GEMB.thermal_conductivity-Tuple{AbstractVector,%20AbstractVector,%20ModelParameters}). About 3.3 m at ice density.

Used both to decay the wave through the initialized column ([`_steady_state_temperature`](/internals_grid#GEMB._steady_state_temperature-Tuple{Real,%20Real,%20ClimateSummary,%20ModelParameters})) and to set the floor on the column depth that resolves it ([`_derive_column_depth`](/internals_grid#GEMB._derive_column_depth-Tuple{ModelParameters,%20ClimateSummary})), so the two agree by construction.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB.gemb_core-Tuple{Any, ClimateForcingStep, ModelParameters, Bool}' href='#GEMB.gemb_core-Tuple{Any, ClimateForcingStep, ModelParameters, Bool}'><span class="jlbinding">GEMB.gemb_core</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
gemb_core(state, cfs::ClimateForcingStep, mp::ModelParameters, verbose::Bool;
          n_target=length(state.dz), z_target=sum(state.dz))
```


Perform a single time-step of the GEMB model. Matches MATLAB's `gemb_core.m`.

The column is returned with exactly `n_target` cells summing to exactly `z_target` metres, so both are invariant across the whole run (see `grid_ops.jl` for the two controllers that enforce them). Cell count is restored by [`manage_layer_thickness`](/internals_grid#GEMB.manage_layer_thickness-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters,%20Bool}) at step 9; column depth is pinned at step 11, _after_ [`calculate_density`](/internals_physics#GEMB.calculate_density-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ClimateForcingStep,%20ModelParameters}), which is the last thing in the timestep to change `dz`.

Returns `(state, flux)` where `state` is the updated column state and `flux` contains energy/mass budget terms for output accumulation.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._assert_grid_feasible-Tuple{Vector{Float64}, Float64, ModelParameters}' href='#GEMB._assert_grid_feasible-Tuple{Vector{Float64}, Float64, ModelParameters}'><span class="jlbinding">GEMB._assert_grid_feasible</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_assert_grid_feasible(dz, z_target, mp)
```


Check at startup that the column the two grid controllers are being handed is one they can actually hold: at least two cells (the thermal solve indexes `temperature[2]`), all thicknesses positive and finite, a positive target depth, and — the substantive check — that `z_target` lies inside `[Σdzmin_i, Σdzmax_i]`, the depth range the `N` cells can span at their band limits.

That last condition is what makes the two controllers compatible. The count controller holds the cell count at `N` while keeping cells inside their bands; the mass controller holds total depth at `z_target`. If `z_target` fell outside the band-limited range, no admissible `N`-cell grid would have the required depth and the two controllers would fight indefinitely — the count controller pushing cells toward their bands, the mass controller pulling total depth away from what those bands permit.

Stable only because [`column_bands!`](/internals_grid#GEMB.column_bands!-Tuple{Vector{Float64},%20Vector{Float64},%20Vector{Float64},%20ModelParameters}) references each cell's band to its **depth** (see there for why). The range then holds steady over a run — measured on the default grid, `Σdzmin` 118–127 m and `Σdzmax` 355–382 m against a 254.6 m target — so the check can be trusted rather than rejecting legitimate configurations after the first cycle.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._compute_output_times-Tuple{Vector{Dates.DateTime}, Symbol}' href='#GEMB._compute_output_times-Tuple{Vector{Dates.DateTime}, Symbol}'><span class="jlbinding">GEMB._compute_output_times</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Compute output times based on output frequency. Returns the last timestep of each day/month, all timesteps, or just the final one.

For `:daily`/`:weekly`/`:monthly`, a timestep is emitted when the next timestep falls in a different day/week/month. The final timestep is compared against a synthetic successor (`times[end] + dt`) rather than emitted unconditionally, so a trailing _partial_ period (e.g. a lone midnight boundary step belonging to a new day) is not saved — only complete days/weeks/months are written. Weeks are grouped by Monday-anchored calendar week (`floor(Date, Week)`). This matches MATLAB's `gemb.m` output indexing (which appends `dates(end) + (dates(end)-dates(end-1))` before diffing).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._gemb_time_loop!-Tuple{Any, Any, Any, Any, Bool, Vector{Dates.DateTime}, Any, Int64, Float64, Float64, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, Vararg{Float64, 7}}' href='#GEMB._gemb_time_loop!-Tuple{Any, Any, Any, Any, Bool, Vector{Dates.DateTime}, Any, Int64, Float64, Float64, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, AbstractVector, Vararg{Float64, 7}}'><span class="jlbinding">GEMB._gemb_time_loop!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_gemb_time_loop!(output, state, model_parameters, mp, verbose, times,
                 output_times, profile_size, z_target, dt_f, <forcing arrays/scalars>)
```


Function-barrier inner loop for [`gemb`](/api#GEMB.gemb). Receives the forcing series as arguments so the compiler specializes on their concrete types instead of re-dispatching every timestep. Mutates `output` in place.

`profile_size` doubles as the fixed cell count and the profile output row count; `z_target` is the fixed total column depth. Both are passed to every `gemb_core` call and hold at every timestep boundary.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='GEMB._extract_profile_at_index-Tuple{DimStack, Int64}' href='#GEMB._extract_profile_at_index-Tuple{DimStack, Int64}'><span class="jlbinding">GEMB._extract_profile_at_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Extract profile at a specific column index from the output DimStack.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/alex-s-gardner/GEMB.jl" target="_blank" rel="noreferrer">source</a></Badge>

</details>

