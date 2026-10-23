# Layout of a `lower_boundary = :ground` column.
#
# The profile and the `gemb` output hold a fixed `N = ground_snow_cells_max + ng` rows, surface
# first: the `ns` snow cells (`material = MATERIAL_SNOW`, possibly none), then the `ng` ground
# cells (`MATERIAL_GROUND`), then `N − ns − ng` empty slots (`MATERIAL_EMPTY`: zero thickness,
# density, pore water, grains and age, temperature `NaN`). Row 1 is therefore always the surface,
# and `dz2z(dz)` gives true depths, the empty slots adding none.
#
# A ground cell's `density` is its solids plus pore ice per unit volume, and its `water` its
# unfrozen pore water [kg m-2], so `dz·density + water` is the cell's constant mass through
# freezing and thawing.

"""
    ground_snow_slots(mp) -> Int

Number of snow slots of a `lower_boundary = :ground` column, `mp.ground_snow_cells_max`.
"""
ground_snow_slots(mp::ModelParameters) = mp.ground_snow_cells_max

"""
    initialize_profile_ground(mp) -> profile::DimStack

The initial column of `lower_boundary = :ground`: no snow, the soil column of `mp.ground_depth`
on [`ground_grid`](@ref), every ground cell at `mp.ground_temperature_deep` (the mean annual
ground temperature, which the bottom cell keeps for the run) with its pore water liquid, or
frozen where that temperature is below 273.15 K. The profile carries the `RESTART_LAYERS` and
`material`.
"""
function initialize_profile_ground(mp::ModelParameters)
    validate_parameters(mp)
    g = ground_column(mp, ground_grid(mp.ground_depth))
    ng = length(g.dz)
    T0 = mp.ground_temperature_deep
    liquid = T0 >= CtoK ? g.water_total : zeros(ng)
    state = _ground_state_vectors(mp, g, Float64[], Float64[], Float64[], Float64[], Float64[],
        Float64[], Float64[], Float64[], fill(T0, ng), liquid)
    zdim = Z(1:length(state.dz); metadata=cf_layer_index_attributes())
    layers = NamedTuple(l => DimArray(getfield(state, l), (zdim,))
                        for l in (RESTART_LAYERS..., :material))
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end

"""
    _ground_state_vectors(mp, g, temperature, dz, density, water, grain_radius,
                          grain_dendricity, grain_sphericity, age, temperature_ground,
                          water_ground) -> NamedTuple

Assemble the fixed-length column of the module preamble from the `ns` snow cells and the ground
state (temperature [K] and unfrozen pore water [kg m-2] per ground cell). Throws when the snow
exceeds its slots.
"""
function _ground_state_vectors(mp::ModelParameters, g::GroundColumn, temperature, dz, density,
    water, grain_radius, grain_dendricity, grain_sphericity, age, temperature_ground,
    water_ground)
    ns = length(dz)
    ng = length(g.dz)
    S = ground_snow_slots(mp)
    ns <= S || error("ground column: $(ns) snow cells exceed the $(S) snow slots " *
                     "(ground_snow_cells_max)")
    N = S + ng
    out = (temperature=fill(NaN, N), dz=zeros(N), density=zeros(N), water=zeros(N),
        grain_radius=zeros(N), grain_dendricity=zeros(N), grain_sphericity=zeros(N),
        age=zeros(N), material=fill(MATERIAL_EMPTY, N))
    for i in 1:ns
        out.temperature[i] = temperature[i]; out.dz[i] = dz[i]; out.density[i] = density[i]
        out.water[i] = water[i]; out.grain_radius[i] = grain_radius[i]
        out.grain_dendricity[i] = grain_dendricity[i]
        out.grain_sphericity[i] = grain_sphericity[i]; out.age[i] = age[i]
        out.material[i] = MATERIAL_SNOW
    end
    solids = (1 - g.porosity) * mp.ground_solid_density
    for j in 1:ng
        i = ns + j
        out.temperature[i] = temperature_ground[j]
        out.dz[i] = g.dz[j]
        out.density[i] = solids + (g.water_total[j] - water_ground[j]) / g.dz[j]
        out.water[i] = water_ground[j]
        out.material[i] = MATERIAL_GROUND
    end
    return out
end

"""
    _ground_state_split(state, g) -> (snow::NamedTuple, temperature_ground, water_ground)

Inverse of [`_ground_state_vectors`](@ref): fresh vectors of the snow cells, as
[`column_state`](@ref) bundles them, and the ground temperatures and unfrozen pore water. Throws
unless the column has the layout of the module preamble and its ground cells are `g`'s.
"""
function _ground_state_split(state, g::GroundColumn)
    mat = state.material
    ng = length(g.dz)
    ns = 0
    while ns < length(mat) && mat[ns+1] == MATERIAL_SNOW
        ns += 1
    end
    length(mat) >= ns + ng || error("ground column: $(length(mat)) rows cannot hold $(ns) snow " *
                                    "and $(ng) ground cells")
    for j in 1:ng
        mat[ns+j] == MATERIAL_GROUND || error("ground column: row $(ns + j) is not ground")
        state.dz[ns+j] == g.dz[j] || error("ground column: ground cell $(j) thickness " *
            "$(state.dz[ns+j]) m differs from the ground grid's $(g.dz[j]) m")
    end
    for i in ns+ng+1:length(mat)
        (mat[i] == MATERIAL_EMPTY && state.dz[i] == 0.0 && state.density[i] == 0.0 &&
         state.water[i] == 0.0) || error("ground column: row $(i) below the ground is not an " *
                                         "empty slot")
    end
    r = 1:ns
    snow = column_state(state.temperature[r], state.dz[r], state.density[r], state.water[r],
        state.grain_radius[r], state.grain_dendricity[r], state.grain_sphericity[r],
        state.age[r])
    rg = ns+1:ns+ng
    return snow, state.temperature[rg], state.water[rg]
end

"""
    ground_soil_temperature(temperature_ground, g, depth) -> T [K]

Temperature at `depth` [m] below the ground surface, linear between the centres of the ground
cells, and the top (bottom) cell's temperature above (below) its centre.
"""
function ground_soil_temperature(temperature_ground::AbstractVector{Float64}, g::GroundColumn,
    depth::Real)
    z_prev = g.dz[1] / 2
    depth <= z_prev && return temperature_ground[1]
    z = 0.0
    for j in 1:length(g.dz)-1
        zc = z + g.dz[j] / 2
        zn = z + g.dz[j] + g.dz[j+1] / 2
        if depth <= zn
            w = (depth - zc) / (zn - zc)
            return temperature_ground[j] + w * (temperature_ground[j+1] - temperature_ground[j])
        end
        z += g.dz[j]
    end
    return temperature_ground[end]
end
