"""
    manage_layer_thickness(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, mp::ModelParameters, verbose::Bool; dt, n_target=length(dz))

Return every cell to its per-cell thickness band, then restore the column to exactly
`n_target` cells.

Three passes, in order:

1. **Merge** cells thinner than their per-cell `dzmin` band into the cell below (the bottom
   cell into the deepest surviving cell above it). A surface snow cell is the exception: it
   is kept, however thin, while it holds at least [`snow_layer_min_mass`](@ref)`(dt)`, `dt`
   [s] being the forcing step, and is distinct
   from cell 2 ([`snow_layers_distinct`](@ref)), which is how fresh snow that fell on
   dissimilar snow, firn or ice stays a layer of its own until it has aged or compacted
   toward what lies below. Snow merged into snow mixes the grain state by mass
   ([`merge_snow_pair!`](@ref)); surface snow merged into firn or ice is absorbed at the
   firn or ice density ([`absorb_snow!`](@ref)); every other merge takes the grain state of
   the thin cell ([`merge_pair!`](@ref)).
2. **Split** cells thicker than their per-cell `dzmax` band in half.
3. **Count control** ([`enforce_column_length!`](@ref)) — merge or split in the deep
   column until the cell count is exactly `n_target`, absorbing whatever net change the
   two passes above and the surface physics (accumulation, melt-out) produced.

All three passes are built from the primitives in `grid_ops.jl` and **conserve mass and
energy exactly**, so this function no longer reports a basal mass/energy flux. Column depth
is pinned separately, and later in the timestep, by [`trim_bottom!`](@ref) — see
`grid_ops.jl` for why count and mass are controlled independently.

The surface cell can therefore enter the next thermal solve with as little as
`snow_layer_min_mass(dt)`. Both solvers are stable there; the explicit one shortens its sub-step
through the surface relaxation bound (`THERMAL_SURFACE_STEP_FRACTION`).

The Dirichlet temperature boundary condition is re-imposed at the bottom; the energy that
requires is the only term returned in `E_added`.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity,
grain_sphericity, age, E_added)`. Arrays are mutated in place and their
length on return is `n_target`.

Two keywords serve the snow over the ground of `lower_boundary = :ground`, whose cell count
varies: `n_max > 0` replaces the fixed count by a ceiling, the count controller merging only
while there are more than `n_max` cells, and `dirichlet_bottom = false` leaves the bottom cell's
temperature alone, the snow's base not being a fixed-temperature boundary. A column of one cell
is not merged.

`age` needs no handling here beyond being in the bundle: all three passes go through
[`merge_pair!`](@ref), [`merge_snow_pair!`](@ref) or [`absorb_snow!`](@ref) (which mass-weight
it) and [`split_cell!`](@ref) (for which duplication is already correct, age being
intensive). The merge pass keeps `M` in sync across a chain of cells folding into one target,
so the weighting stays correct over a chain rather than only over a pair.
"""
function manage_layer_thickness(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, age::Vector{Float64},
    mp::ModelParameters, verbose::Bool; dt::Float64, n_target::Int=length(dz),
    workspace::ColumnWorkspace=ColumnWorkspace(), n_max::Int=0,
    dirichlet_bottom::Bool=true)

    cols = column_state(temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, age)

    if verbose
        M_total_initial, E_total_initial = column_mass_energy(cols, mp)
    end

    T_bottom = temperature[end]
    m = length(temperature)

    # Per-cell thickness bands; dzmax is carried through the merge pass because merging
    # shifts cells up and their band assignment must shift with them. Both come from the
    # workspace, so they are usually longer than the column: `column_bands!` writes
    # `eachindex(dz)` and every read below is bounded by the cell count, never by the
    # buffer length.
    _resize_workspace!(workspace, max(m, n_target) + 1)
    column_dzmin2 = workspace.dzmin
    column_dzmax2 = workspace.dzmax
    column_bands!(column_dzmin2, column_dzmax2, dz, mp)

    # Cell masses, kept in sync through the merge pass so a chain of merges into the same
    # target weights correctly. Filled cell by cell rather than broadcast, because the
    # buffer is longer than `dz` and `.=` would be a length mismatch.
    M = workspace.mass
    @inbounds for i in 1:m
        M[i] = dz[i] * density[i]
    end

    ## MERGE CELLS BELOW THEIR MINIMUM THICKNESS
    # A view over the live cells, so `findlast` cannot run off into the surplus tail — where an
    # unwritten `false` would read as the last surviving cell and pick the wrong merge target.
    delete_cell = view(workspace.delete_cell, 1:m)
    fill!(delete_cell, false)
    for i in 1:m
        if m > 1 && dz[i] < (column_dzmin2[i] - D_TOLERANCE)
            # A distinct snow layer at the surface is kept however thin, down to
            # `snow_layer_min_mass(dt)`.
            (i == 1 && m > 1 && _surface_layer_kept(cols, mp, dt)) && continue
            delete_cell[i] = true

            # Merge downward, except for the bottom cell, which merges into the
            # deepest surviving cell above it.
            i_target = i == m ? findlast(!, delete_cell) : i + 1
            M[i_target] = _merge_thin_cell!(cols, i, i_target, M[i], M[i_target], mp)
        end
    end

    to_delete = workspace.to_delete
    empty!(to_delete)
    for i in 1:m
        delete_cell[i] && push!(to_delete, i)
    end
    if !isempty(to_delete)
        close_slot!(cols, to_delete)
        # Every index is <= m, so deleting them from the longer buffer shifts the live prefix
        # exactly as it would a buffer sized to the column; only the inert tail differs.
        deleteat!(column_dzmax2, to_delete)
    end

    m = length(temperature)

    ## SPLIT CELLS ABOVE THEIR MAXIMUM THICKNESS
    # Collected first, then applied back-to-front so earlier indices stay valid.
    f = workspace.to_split
    empty!(f)
    @inbounds for i in 1:m
        if dz[i] > column_dzmax2[i] + D_TOLERANCE
            push!(f, i)
        end
    end

    @inbounds for k in length(f):-1:1
        split_cell!(cols, f[k])
    end

    ## COUNT CONTROL — restore the fixed cell count (exactly conservative)
    n_max > 0 && (n_target = min(length(dz), n_max))
    enforce_column_length!(cols, n_target, mp; workspace)

    # The count controller conserves mass and energy exactly, and depth is pinned by
    # `trim_bottom!` later in the timestep, so the only budget term arising here is the
    # Dirichlet bottom boundary condition.
    E_added = 0.0
    if dirichlet_bottom
        E_added = (dz[end] * density[end]) *
                  (specific_enthalpy(mp, T_bottom) - specific_enthalpy(mp, temperature[end]))
        temperature[end] = T_bottom
    end

    ## CHECK FOR MASS AND ENERGY CONSERVATION
    if verbose
        M_total_final, E_total_final = column_mass_energy(cols, mp)

        M_delta = M_total_initial - M_total_final
        E_delta = E_total_initial - E_total_final + E_added

        E_tol = energy_tolerance(E_total_initial)

        if (abs(M_delta) > M_TOLERANCE) || (abs(E_delta) > E_tol)
            error("Mass and/or energy are not conserved in manage_layer_thickness:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
    age, E_added
end

# Whether a surface cell thinner than its band is kept as a layer of its own: it is snow, holds
# at least `snow_layer_min_mass(dt)`, and is distinct from cell 2 (`snow_layers_distinct`, always
# so over firn or ice).
@inline function _surface_layer_kept(cols::NamedTuple, mp::ModelParameters, dt::Float64)
    density = cols.density
    density[1] < DENSITY_PORE_CLOSEOFF - D_TOLERANCE || return false
    cols.dz[1] * density[1] >= snow_layer_min_mass(dt) || return false
    D = snow_layer_distance(cols.grain_dendricity[1], cols.grain_sphericity[1],
        cols.grain_radius[1], cols.water[1] > WATER_TOLERANCE, cols.age[1], density[1],
        cols.grain_dendricity[2], cols.grain_sphericity[2], cols.grain_radius[2],
        cols.water[2] > WATER_TOLERANCE, cols.age[2], density[2])
    return snow_layers_distinct(cols.dz[1], D, mp)
end

# Fold a cell thinner than its band into `i_target`. Snow into snow mixes the grain state by mass
# (`merge_snow_pair!`); surface snow into firn or ice is absorbed at the firn or ice density
# (`absorb_snow!`); every other pair keeps the grain state of cell `i` (`merge_pair!`).
@inline function _merge_thin_cell!(cols::NamedTuple, i::Int, i_target::Int,
    M_i::Float64, M_target::Float64, mp::ModelParameters)
    snow_i = cols.density[i] < DENSITY_PORE_CLOSEOFF - D_TOLERANCE
    snow_target = cols.density[i_target] < DENSITY_PORE_CLOSEOFF - D_TOLERANCE
    if snow_i && snow_target
        return merge_snow_pair!(cols, i, i_target, M_i, M_target, mp)
    elseif snow_i && i == 1
        return absorb_snow!(cols, i, i_target, M_i, M_target, mp)
    end
    return merge_pair!(cols, i, i_target, M_i, M_target, mp)
end
