"""
Shared vertical-grid primitives.

GEMB's column is Lagrangian: cells are created by accumulation and splitting, and
destroyed by melt-out and merging. Historically each of those events edited the nine
per-cell state vectors independently at six separate call sites, so the structural
bookkeeping was duplicated nine ways at each site and adding a tenth state field meant
finding and editing all six.

Everything structural now goes through the primitives below, which take the vectors as a
`ColumnState` NamedTuple and iterate over it generically, so the *structural* bookkeeping
(shifting, merging, splitting) is written once instead of once per field per site.

This removes the per-field duplication within each site, not the number of sites: adding a
field still touches [`column_state`](@ref), the physics that fills it, the
extensive/intensive handling in [`merge_pair!`](@ref) / [`split_cell!`](@ref), and five
plumbing sites that share no source of truth (`initialize_profile`'s two `layers` tuples,
the driver's output allocation and write loop, and `_extract_profile_at_index` — the last
failing *silently* rather than loudly, since a profile field missing there is simply not
carried between spinup cycles). The bundle is
rebuilt per call site because the physics rebinds several vectors mid-timestep, so a
long-lived bundle would go stale. Every field is `Vector{Float64}`, so a mis-ordered
argument list type-checks silently — keep the positional order identical everywhere.

# The two controllers

Column *count* and column *mass* are separate problems with separate controllers:

1. **Count** — [`enforce_column_length!`](@ref) restores the fixed cell count `N` by
   merging adjacent deep cells (when the physics net-created cells) or splitting deep
   cells (when it net-destroyed them). Merging and splitting are **exactly conservative**
   in both mass and energy, so this controller contributes nothing to the model's mass or
   energy budget.
2. **Mass** — [`trim_bottom!`](@ref) pins the total column depth to `Z` by adjusting the
   *bottom cell's thickness*. This is the model's only basal mass/energy flux and is the
   sole grid operation that reports a nonzero `(mass_added, E_added)`.

Together they make `length == N` and `Σdz == Z` true at every timestep boundary, which is
what lets the output arrays be sized exactly `N` rows with no padding.

The mass controller is signed and continuous, so both regimes are the same code path with
opposite sign: the bottom cell shrinks under accumulation (mass leaves the base) and grows
under ablation (basal accretion). Magnitudes are small either way — order 0.1 kg m-2 per step
at Summit, under 2 kg m-2 at 5 m yr-1 ablation — so the flux is smooth rather than the
whole-cell jumps a drop/create scheme would produce.
"""

"""
    column_state(temperature, dz, density, water, grain_radius,
                 grain_dendricity, grain_sphericity, age)

Bundle the eight per-cell state vectors so the grid primitives can operate on all of them
generically. `values(...)` of the result is a homogeneous `NTuple{8,Vector{Float64}}`, so
the generic loops below stay type-stable and fully unrolled.

Albedo is deliberately absent: every albedo method is diagnostic in the current column
state (see [`calculate_albedo`](@ref)), so there is no per-cell albedo to carry through a
merge or a split.

The structural primitives (`open_slot!`, `close_slot!`) need no change when a field is added;
`merge_pair!` and `split_cell!` do, since a new field's extensive-vs-intensive merge semantics
cannot be inferred. See the module preamble on why the bundle is rebuilt per call site.
"""
@inline column_state(temperature, dz, density, water, grain_radius,
    grain_dendricity, grain_sphericity, age) =
    (; temperature, dz, density, water, grain_radius,
       grain_dendricity, grain_sphericity, age)

"""
    mix_age(age1, M1, age2, M2) -> age

Mass-weighted mean age of two masses combined: `(age1·M1 + age2·M2)/(M1+M2)`.

Age is intensive and mixes linearly in mass, so this is the whole rule — the analogue of
[`mix_temperature`](@ref) for a quantity with no equation of state. Every site where age
changes is one of its three cases: mass arriving from outside at age 0
([`dilute_age`](@ref)), two cells combining ([`merge_pair!`](@ref)), and meltwater joining
through-flow in `calculate_melt`. Keeping them one function keeps the zero-mass convention
one decision.

Returns `age1` when the total mass is zero, matching [`mix_temperature`](@ref).
"""
@inline function mix_age(age1::Float64, M1::Float64, age2::Float64, M2::Float64)
    M = M1 + M2
    M <= 0.0 && return age1
    return (age1 * M1 + age2 * M2) / M
end

"""
    dilute_age(age, mass_old, mass_added) -> age

Mass-weighted mean age of a cell that gains `mass_added` [kg m-2] of brand-new (age 0)
material, where `mass_old` is the cell's total mass beforehand. The [`mix_age`](@ref) case
where one side is exactly zero.

This is the update for every site where mass enters the column from outside: snowfall, rain,
and vapour deposition. Returns `age` unchanged for non-positive `mass_added`, so the caller
need not branch on the sign of a signed flux (sublimation removes a fraction of the cell and
is age-neutral).
"""
@inline dilute_age(age::Float64, mass_old::Float64, mass_added::Float64) =
    mass_added <= 0.0 ? age : mix_age(age, mass_old, 0.0, mass_added)

"""
    open_slot!(cols, i)

Insert a cell at index `i`, duplicating the current contents of cell `i` into the new
slot, across every field. Cells `i..n` shift down to `i+1..n+1`.

Used to prepend a fresh-snow cell (`i = 1`) and as the structural half of
[`split_cell!`](@ref).
"""
@inline function open_slot!(cols::NamedTuple, i::Int)
    @inbounds for a in values(cols)
        insert!(a, i, a[i])
    end
    return nothing
end

"""
    close_slot!(cols, i)
    close_slot!(cols, idx::AbstractVector{Int})

Delete cell `i` (or the sorted, unique indices `idx`) across every field. Cells below shift
up to close the gap. The vector form deletes in one pass per field, so it is O(n) in total
rather than O(n) per index.
"""
@inline function close_slot!(cols::NamedTuple, i::Int)
    @inbounds for a in values(cols)
        deleteat!(a, i)
    end
    return nothing
end

@inline function close_slot!(cols::NamedTuple, idx::AbstractVector{Int})
    isempty(idx) && return nothing
    @inbounds for a in values(cols)
        deleteat!(a, idx)
    end
    return nothing
end

"""
    merge_pair!(cols, i, i_target, M_i, M_target, mp) -> M_new

Fold the contents of cell `i` into cell `i_target`, conserving mass and energy exactly.
`M_i` and `M_target` are the ice/firn masses `dz*density` of the two cells; the combined
mass is returned so callers tracking a mass vector can keep it in sync.

Temperature is combined by [`mix_temperature`](@ref), which mixes in enthalpy and inverts.
Under a temperature-dependent heat capacity the mass-weighted mean temperature is *not*
energy-conserving — it loses `(b/2)·M·Var_M(T)` joules, thousands of times the
conservation tolerance for a realistic deep merge. Thickness and pore water are extensive
and are summed, with density recovered from the totals. Grain properties are inherited from
cell `i` — that is the historical convention (the *upper* cell of the pair wins), preserved
deliberately.

`age` is mass-weighted on **total** cell mass (`dz*density + water`), matching how the field
is defined, so `M_i`/`M_target` alone are not the right weights — the pore water of each cell
is added in. Deliberately *not* the grain-property convention above: inheriting the upper
cell's age would discard the older cell's residence time entirely and make deep age drift
young with every merge.

This performs only the physics of the merge. The now-redundant cell `i` still occupies a
slot; the caller removes it with [`close_slot!`](@ref), either immediately or in a batch.
"""
@inline function merge_pair!(cols::NamedTuple, i::Int, i_target::Int,
                             M_i::Float64, M_target::Float64, mp::ModelParameters)
    M_new = M_i + M_target
    @inbounds begin
        cols.temperature[i_target] = mix_temperature(mp,
            cols.temperature[i], M_i, cols.temperature[i_target], M_target)
        # Grain properties come from the upper cell of the merged pair.
        cols.grain_radius[i_target] = cols.grain_radius[i]
        cols.grain_dendricity[i_target] = cols.grain_dendricity[i]
        cols.grain_sphericity[i_target] = cols.grain_sphericity[i]

        # Age, mass-weighted on total cell mass. Computed before `water` is summed below,
        # which needs each cell's pore water separately.
        cols.age[i_target] = mix_age(cols.age[i], M_i + cols.water[i],
            cols.age[i_target], M_target + cols.water[i_target])

        cols.dz[i_target] = cols.dz[i] + cols.dz[i_target]
        cols.density[i_target] = M_new / cols.dz[i_target]
        cols.water[i_target] = cols.water[i] + cols.water[i_target]
    end
    return M_new
end

"""
    split_cell!(cols, i)

Split cell `i` into two cells of half its thickness, conserving mass and energy exactly.

Thickness and pore water are halved (both extensive); density, temperature, grain
properties and age are intensive and so are simply duplicated. The result occupies
indices `i` and `i+1`, and the column grows by one cell.
"""
@inline function split_cell!(cols::NamedTuple, i::Int)
    @inbounds begin
        cols.dz[i] /= 2
        cols.water[i] /= 2
    end
    open_slot!(cols, i)
    return nothing
end

"""
    column_bands!(dzmin, dzmax, dz, mp) -> (dzmin, dzmax)

Fill the per-cell thickness bands `[dzmin[i], dzmax[i]]` that drive merging and splitting.

Cells within `mp.column_ztop` of the surface get the uniform bands
`[mp.column_dzmin, mp.column_dzmax]`. Below that the bands are scaled by `mp.column_zy^k`,
where `k` is the geometric grid position **implied by the cell's depth** — that is, the
solution of

    z = column_ztop + column_dztop · zy · (zy^k − 1) / (zy − 1)

for the cell's lower face depth `z`. This is exactly the stretching
[`initialize_grid`](@ref) uses to build the column, so a freshly-built grid has every cell
sitting at the geometric centre of its own band, and a cell is always measured against the
thickness appropriate to the depth it has actually reached.

# Why the band must be referenced to depth, not to cell index

Do not replace this with the obvious version — an accumulator multiplied by `zy` once per cell
below `column_ztop` — which indexes the band by *position in the deep sequence*. That is
unstable under compaction: as near-surface cells thin, more fit inside `column_ztop`, the deep
sequence shortens, and every deep band collapses though no cell has moved. On the default grid
`Σdzmax` fell from 382 m to under 50 m within a simulated year; the split pass then shattered
the deep column into slivers, the count controller merged them back, and the column degenerated
to one cell holding 245.5 m of its 254.6 m by spinup cycle 5. Depth-referencing removes the
feedback: a cell's band depends only on where it is, making the geometric grid a fixed point of
the merge/split rules.

`dzmin` and `dzmax` must be at least `length(dz)` long; only the first `length(dz)`
entries are written.
"""
function column_bands!(dzmin::Vector{Float64}, dzmax::Vector{Float64},
                       dz::Vector{Float64}, mp::ModelParameters)
    z_cum = 0.0
    ztop = mp.column_ztop + D_TOLERANCE
    # Constant of the geometric stretch, hoisted out of the loop.
    scale = (mp.column_zy - 1.0) / (mp.column_dztop * mp.column_zy)
    @inbounds for i in eachindex(dz)
        z_cum += dz[i]
        if z_cum <= ztop
            dzmin[i] = mp.column_dzmin
            dzmax[i] = mp.column_dzmax
        else
            # `zy^k` for the continuous geometric position `k` implied by this cell's depth.
            # Solving the stretch relation for `k` gives `k = log1p(u)/log(zy)` with
            # `u = (z - ztop)*scale`, so `zy^k` collapses to `1 + u` exactly — the log and the
            # power cancel. Written in the closed form: same value, no transcendentals, and
            # well behaved as `zy → 1` (where `log(zy) → 0` would give `0 * Inf`).
            zy_power = 1.0 + (z_cum - mp.column_ztop) * scale
            dzmin[i] = mp.column_dzmin * zy_power
            dzmax[i] = mp.column_dzmax * zy_power
        end
    end
    return dzmin, dzmax
end

"""
    enforce_column_length!(cols, n_target, mp) -> Int

Controller 1 (count). Restore the column to exactly `n_target` cells by merging or
splitting in the **deep** part of the column, and return the number of merge/split
operations performed.

Both operations conserve mass and energy exactly, so this returns no budget terms — it is
invisible to the model's mass and energy accounting.

Which cells are touched:

- **Too many cells** → repeatedly merge an adjacent pair, thinnest first
  ([`_select_merge_pair`](@ref)).
- **Too few cells** → repeatedly split a cell, thickest first
  ([`_select_split_cell`](@ref)).

Both selectors rank candidates by staying inside the cell's `dzmin`/`dzmax` band and by lying
below `mp.column_ztop`, with **in-band dominating deep** (see `_select_merge_pair` for why),
breaking ties toward the deeper candidate. Neither touches cell 1, so the surface cell is never
consumed.

Coarsening at depth while fine cells arrive at the top is the history a descending parcel
experiences, so this does not fight the grid — [`initialize_grid`](@ref)'s geometric profile is
a fixed point of the process. Over a 75-cycle spinup on the default grid the deep structure
holds (max cell 22.3 → 30.2 m, cells over 5 m steady at 16) with no cell out of band. Ranking
in-band first does let a merge land inside `column_ztop` when no deep pair is admissible: that
is the safety valve preventing the deep column from being consumed.
"""
function enforce_column_length!(cols::NamedTuple, n_target::Int, mp::ModelParameters)
    dz = cols.dz
    n = length(dz)
    n == n_target && return 0

    # Bands are recomputed after each operation: merging and splitting change the
    # cumulative depth of every cell below, and hence which band each falls in.
    dzmin = Vector{Float64}(undef, max(n, n_target) + 1)
    dzmax = Vector{Float64}(undef, max(n, n_target) + 1)

    # Each iteration changes the length by exactly one, so the work is bounded by the initial
    # discrepancy. The cap is a guard against a primitive failing to make progress, not an
    # expected exit. One counter serves both loops: each moves `n` monotonically toward
    # `n_target` and stops on reaching it, so at most one of them ever runs.
    ops = 0
    budget = abs(n - n_target) + 2

    while n > n_target && ops < budget
        column_bands!(dzmin, dzmax, dz, mp)
        i = _select_merge_pair(dz, dzmax, n, mp)
        i === nothing && break
        ii = i::Int
        @inbounds merge_pair!(cols, ii, ii + 1,
            dz[ii] * cols.density[ii], dz[ii+1] * cols.density[ii+1], mp)
        close_slot!(cols, ii)
        n -= 1
        ops += 1
    end

    while n < n_target && ops < budget
        column_bands!(dzmin, dzmax, dz, mp)
        i = _select_split_cell(dz, dzmin, n, mp)
        i === nothing && break
        split_cell!(cols, i::Int)
        n += 1
        ops += 1
    end

    if n != n_target
        error("enforce_column_length!: could not reach the target cell count " *
              "(n = $n, target = $n_target). No admissible deep merge/split candidate " *
              "was available, which means the column geometry and the " *
              "column_dzmin/column_dzmax/column_ztop settings are inconsistent.")
    end

    return ops
end

# Index of the first cell whose *upper* face lies below `column_ztop`, i.e. the start of the
# deep (geometrically stretched) region. Returns `n + 1` when the whole column is shallower
# than `column_ztop`.
#
# Deliberately offset by one cell from `column_bands!`, which classifies by the *lower* face:
# the straddling cell gets a stretched band there but is not counted deep here. This is only a
# preference for coarsening away from the surface (band admissibility comes from
# `column_bands!`), so excluding the boundary cell is the conservative direction. Keep the
# offset in this direction if either rule changes.
@inline function _deep_start(dz::Vector{Float64}, n::Int, mp::ModelParameters)
    z_cum = 0.0
    ztop = mp.column_ztop + D_TOLERANCE
    @inbounds for i in 1:n
        if z_cum > ztop
            return i
        end
        z_cum += dz[i]
    end
    return n + 1
end

"""
    _select_merge_pair(dz, dzmax, n, mp) -> Union{Int,Nothing}

Choose `i` such that cells `i` and `i+1` are merged. Returns `nothing` when no merge is
available (`n < 3`, or only the surface cell would qualify).

Candidates are ranked by two independent preferences, in this order:

 1. **In band** — the combined thickness must still fit the lower cell's `dzmax`, so the
    merge is not immediately undone by next timestep's split test.
 2. **Deep** — at or below `column_ztop`, leaving the fine near-surface grid (where the
    radiation and turbulent-flux gradients live) alone where possible.

Within each tier the thinnest pair wins. **In-band must dominate deep** — this ordering is
load-bearing. The deep region is a finite budget of cells while accumulation demands a merge
thousands of times per year, so a deep-first rule eventually has only the bottom pair to
offer and returns it every time, folding the deep column into one monolithic cell (measured:
245 m of a 254.6 m column). In-band first self-corrects, since an over-thick bottom cell is
by construction out of band.

See also [`_select_split_cell`](@ref).
"""
function _select_merge_pair(dz::Vector{Float64}, dzmax::Vector{Float64},
                            n::Int, mp::ModelParameters)
    n < 3 && return nothing
    deep = _deep_start(dz, n, mp)

    best = (0, 0)                        # (tier, index); higher tier wins
    best_dz = Inf
    @inbounds for i in 2:(n-1)           # never merge away the surface cell
        combined = dz[i] + dz[i+1]
        in_band = combined <= dzmax[i+1] + D_TOLERANCE
        tier = (in_band ? 2 : 0) + (i >= deep ? 1 : 0)
        # `<=` not `<` so equal-thickness ties go to the *deeper* pair: on a uniform column
        # every pair scores identically, and coarsening belongs away from the surface.
        if tier > best[1] || (tier == best[1] && combined <= best_dz)
            best = (tier, i)
            best_dz = combined
        end
    end
    best[2] == 0 && return nothing
    return best[2]
end

"""
    _select_split_cell(dz, dzmin, n, mp) -> Union{Int,Nothing}

Choose the cell to split, or `nothing` when none qualifies. Ranked by the same two-tier
scheme as [`_select_merge_pair`](@ref) (see the note there on why in-band must dominate
deep): both halves clearing `dzmin` outranks being below `column_ztop`, and the thickest
candidate wins within a tier.
"""
function _select_split_cell(dz::Vector{Float64}, dzmin::Vector{Float64},
                            n::Int, mp::ModelParameters)
    n < 2 && return nothing
    deep = _deep_start(dz, n, mp)

    best = (0, 0)
    best_dz = -Inf
    @inbounds for i in 2:n               # leave the surface cell intact
        d = dz[i]
        in_band = d / 2 >= dzmin[i] - D_TOLERANCE
        tier = (in_band ? 2 : 0) + (i >= deep ? 1 : 0)
        # `>=` breaks equal-thickness ties toward the deeper cell, as in `_select_merge_pair`.
        if tier > best[1] || (tier == best[1] && d >= best_dz)
            best = (tier, i)
            best_dz = d
        end
    end
    best[2] == 0 && return nothing
    return best[2]
end

"""
    trim_bottom!(cols, z_target, mp) -> (mass_added, E_added)

Controller 2 (mass). Pin the total column depth to `z_target` by adjusting the bottom
cell's thickness, and return the mass [kg m-2] and energy [J m-2] added to the column as a
result. This is the model's **only** basal mass/energy flux (horizontal strain, see
[`apply_horizontal_strain!`](@ref), is the only lateral one).

The adjustment is signed and continuous:

- `Σdz > z_target` (the accumulation regime — mass added at the surface has pushed the
  column down): the bottom cell shrinks and material passes out through the base.
  `mass_added` is negative.
- `Σdz < z_target` (the ablation regime — mass removed at the surface, or densification
  has compacted the column): the bottom cell grows and material enters through the base.
  `mass_added` is positive.

Accreted material inherits the bottom cell's properties (the column below is unresolved):
`density[n]`, and `temperature[n]`, which is the Dirichlet-pinned basal boundary. Pore water
scales with the thickness change; grain properties are intensive and unchanged.

`age[n]` is untouched, so both regimes are handled by doing nothing: removal is proportional
(age-neutral), and accreted material inherits `age[n]` like `density` and `temperature`. The
consequence is that under sustained ablation the deepest cell's age is a lower bound on its
true residence time rather than a measurement of it.

Errors if the adjustment would consume the whole bottom cell — at realistic forcing it is
~1e-3 of the cell, so that signals a bug upstream, not an extreme climate.

The column is an Eulerian window on the firn, not a prognostic ice-thickness model, so
`thickness_cumulative` measures basal flux in either regime, not thickness change. For the
latter, use the surface terms (`precipitation - runoff + evaporation_condensation`).
"""
function trim_bottom!(cols::NamedTuple, z_target::Float64, mp::ModelParameters)
    dz = cols.dz
    n = length(dz)
    n == 0 && error("trim_bottom!: empty column")

    z_total = 0.0
    @inbounds for i in 1:n
        z_total += dz[i]
    end

    delta = z_total - z_target          # > 0 ⇒ too deep, shrink the bottom cell
    if abs(delta) <= D_TOLERANCE
        return 0.0, 0.0
    end

    @inbounds begin
        dz_old = dz[n]
        dz_new = dz_old - delta
        if dz_new <= 0.0
            error("trim_bottom!: basal adjustment of $(delta) m would consume the " *
                  "entire bottom cell ($(dz_old) m). Column depth is $(z_total) m " *
                  "against a target of $(z_target) m.")
        end

        # Fraction of the bottom cell leaving (positive) or being accreted (negative).
        frac = delta / dz_old
        water_delta = frac * cols.water[n]
        mass_delta = delta * cols.density[n] + water_delta

        mass_added = -mass_delta
        E_added = -(delta * cols.density[n] * specific_enthalpy(mp, cols.temperature[n]) +
                    water_delta * specific_enthalpy_water(mp))

        dz[n] = dz_new
        cols.water[n] -= water_delta
    end

    return mass_added, E_added
end

"""
    apply_horizontal_strain!(cols, dt_seconds, mp) -> (mass_lateral, E_lateral)

Thin (or thicken) every cell by the horizontal strain rate `mp.horizontal_strain_rate`
[yr-1] over `dt_seconds`, at constant density.

Incompressibility ties vertical thinning to horizontal stretching: for a parcel of
thickness `h` under a horizontal velocity-gradient trace `D = ε̇_xx + ε̇_yy`, `dh/dt = -h D`
at fixed density, whose exact integral over a step of constant `D` is `h *= exp(-D dt)`.
The exponential form keeps `dz` positive for any `D dt`, where forward Euler would not —
the same reason the `:Crocus` branch of [`calculate_density`](@ref) integrates its viscous
law in closed form.

`density` and `temperature` are untouched; grain properties are intensive and unchanged.
Pore `water` [kg m-2] is a per-unit-area quantity, so it scales by the same factor as the
lateral area.

`age` is untouched in both directions. Under divergence every cell loses the same *fraction*
of its mass, which leaves its mean age unchanged; under convergence the material advected in
laterally is neighbouring firn at the same depth, so taking the local `age[i]` is the
consistent choice.

GEMB's column is per unit area, so the mass this removes leaves **laterally**, not through
the base. It is returned signed — negative under divergence (mass leaves), positive under
convergence — and is reported to the mass budget separately from [`trim_bottom!`](@ref)'s
basal flux. Folding the two together would make a thinning column read as spurious basal
accretion, because `trim_bottom!` refills the column to `z_target` immediately afterwards.

Returns `(0.0, 0.0)` and touches nothing when `horizontal_strain_rate == 0`, which is what
keeps the default bit-identical to a run without the term.
"""
function apply_horizontal_strain!(cols::NamedTuple, dt_seconds::Float64,
    mp::ModelParameters)
    strain_rate = mp.horizontal_strain_rate
    strain_rate == 0.0 && return 0.0, 0.0

    # yr-1 -> per-step factor. `SECONDS_PER_YEAR` (Julian, 365.25 d) rather than
    # `DAYS_PER_YEAR_DENSIFICATION`: this is an observed ice-flow strain rate, not a
    # densification coefficient fitted against a 365-day year.
    shrink = expm1(-strain_rate * dt_seconds / SECONDS_PER_YEAR)   # < 0 under divergence

    mass_lateral = 0.0
    E_lateral = 0.0
    enthalpy_water = specific_enthalpy_water(mp)
    @inbounds for i in eachindex(cols.dz)
        dz_delta = cols.dz[i] * shrink
        water_delta = cols.water[i] * shrink
        mass_ice = dz_delta * cols.density[i]
        mass_lateral += mass_ice + water_delta
        E_lateral += mass_ice * specific_enthalpy(mp, cols.temperature[i]) +
                     water_delta * enthalpy_water
        cols.dz[i] += dz_delta
        cols.water[i] += water_delta
    end
    return mass_lateral, E_lateral
end

"""
    column_mass(cols) -> mass [kg m-2]

Total column mass: ice/firn `dz*density` plus pore water. Needs no `ModelParameters`,
unlike [`column_mass_energy`](@ref), so callers that want mass alone can avoid threading
one through.
"""
function column_mass(cols::NamedTuple)
    mass = 0.0
    @inbounds for i in eachindex(cols.dz)
        mass += cols.dz[i] * cols.density[i] + cols.water[i]
    end
    return mass
end

"""
    column_mass_energy(cols, mp) -> (mass, energy)

Total column mass [kg m-2] and enthalpy [J m-2], the reference quantities for every
conservation check. Ice/firn contributes `dz*density*h(T)` with `h` the enthalpy integral
[`specific_enthalpy`](@ref); pore water contributes [`specific_enthalpy_water`](@ref) —
its latent heat plus sensible heat at the melt point.

Absolute energies are measured from a 0 K reference and so are not comparable across
`mp.heat_capacity_method`; only differences are meaningful.
"""
function column_mass_energy(cols::NamedTuple, mp::ModelParameters)
    mass = 0.0
    energy = 0.0
    h_water = specific_enthalpy_water(mp)
    @inbounds for i in eachindex(cols.dz)
        mi = cols.dz[i] * cols.density[i]
        mass += mi + cols.water[i]
        energy += mi * specific_enthalpy(mp, cols.temperature[i]) + cols.water[i] * h_water
    end
    return mass, energy
end
