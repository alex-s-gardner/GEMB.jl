"""
    manage_layer_thickness(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, mp::ModelParameters, verbose::Bool; n_target=length(dz))

Return every cell to its per-cell thickness band, then restore the column to exactly
`n_target` cells.

Three passes, in order:

1. **Merge** cells thinner than their per-cell `dzmin` band into the cell below.
2. **Split** cells thicker than their per-cell `dzmax` band in half.
3. **Count control** ([`enforce_column_length!`](@ref)) — merge or split in the deep
   column until the cell count is exactly `n_target`, absorbing whatever net change the
   two passes above and the surface physics (accumulation, melt-out) produced.

All three passes are built from the primitives in `grid_ops.jl` and **conserve mass and
energy exactly**, so this function no longer reports a basal mass/energy flux. Column depth
is pinned separately, and later in the timestep, by [`trim_bottom!`](@ref) — see
`grid_ops.jl` for why count and mass are controlled independently.

The Dirichlet temperature boundary condition is re-imposed at the bottom; the energy that
requires is the only term returned in `E_added`.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity,
grain_sphericity, albedo, albedo_diffuse, E_added)`. Arrays are mutated in place and their
length on return is `n_target`.
"""
function manage_layer_thickness(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, albedo::Vector{Float64},
    albedo_diffuse::Vector{Float64},
    mp::ModelParameters, verbose::Bool; n_target::Int=length(dz))

    cols = column_state(temperature, dz, density, water, grain_radius,
        grain_dendricity, grain_sphericity, albedo, albedo_diffuse)

    if verbose
        M_total_initial, E_total_initial = column_mass_energy(cols)
    end

    T_bottom = temperature[end]
    m = length(temperature)

    # Per-cell thickness bands; dzmax is carried through the merge pass because merging
    # shifts cells up and their band assignment must shift with them.
    column_dzmin2 = Vector{Float64}(undef, m)
    column_dzmax2 = Vector{Float64}(undef, m)
    column_bands!(column_dzmin2, column_dzmax2, dz, mp)

    # Cell masses, kept in sync through the merge pass so a chain of merges into the same
    # target weights correctly.
    M = dz .* density

    ## MERGE CELLS BELOW THEIR MINIMUM THICKNESS
    delete_cell = falses(m)
    for i in 1:m
        if dz[i] < (column_dzmin2[i] - D_TOLERANCE)
            delete_cell[i] = true

            # Merge downward, except for the bottom cell, which merges into the
            # deepest surviving cell above it.
            i_target = i == m ? findlast(.!delete_cell) : i + 1
            M[i_target] = merge_pair!(cols, i, i_target, M[i], M[i_target])
        end
    end

    to_delete = findall(delete_cell)
    if !isempty(to_delete)
        close_slot!(cols, to_delete)
        deleteat!(column_dzmax2, to_delete)
    end

    m = length(temperature)

    ## SPLIT CELLS ABOVE THEIR MAXIMUM THICKNESS
    # Collected first, then applied back-to-front so earlier indices stay valid.
    f = Int[]
    @inbounds for i in 1:m
        if dz[i] > column_dzmax2[i] + D_TOLERANCE
            push!(f, i)
        end
    end

    @inbounds for k in length(f):-1:1
        split_cell!(cols, f[k])
    end

    ## COUNT CONTROL — restore the fixed cell count (exactly conservative)
    enforce_column_length!(cols, n_target, mp)

    # The count controller conserves mass and energy exactly, and depth is pinned by
    # `trim_bottom!` later in the timestep, so the only budget term arising here is the
    # Dirichlet bottom boundary condition.
    E_added = ((T_bottom - temperature[end]) * (dz[end] * density[end]) * C_ICE)
    temperature[end] = T_bottom

    ## CHECK FOR MASS AND ENERGY CONSERVATION
    if verbose
        M_total_final, E_total_final = column_mass_energy(cols)

        M_delta = M_total_initial - M_total_final
        E_delta = E_total_initial - E_total_final + E_added

        if (abs(M_delta) > 1e-3) || (abs(E_delta) > 1e-3)
            error("Mass and/or energy are not conserved in manage_layer_thickness:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, albedo, albedo_diffuse, E_added
end
