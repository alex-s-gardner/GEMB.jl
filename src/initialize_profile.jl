"""
    initialize_profile(mp::ModelParameters, cf::ClimateForcing)

Initialize a GEMB firn column profile as a DimStack.
Matches MATLAB's `model_initialize_profile.m`.

Returns a DimStack with Z dimension containing:
- z_center, dz, temperature, density, water, grain_radius,
  grain_dendricity, grain_sphericity, albedo, albedo_diffuse
"""
function initialize_profile(mp::ModelParameters, cf::ClimateForcing)
    T_mean = cf.temperature_air_mean

    @assert T_mean > 0 "temperature_air_mean must exceed 0 K."
    if T_mean < 100
        @warn "temperature_air_mean should be in kelvin, but is below 100, suggesting an error."
    end

    # Initialize grid
    dz = initialize_grid(mp)
    z_center = dz2z(dz)
    m = length(dz)

    # Create Z dimension
    zdim = Z(1:m)

    return DimStack((
        z_center=DimArray(z_center, (zdim,)),
        dz=DimArray(dz, (zdim,)),
        temperature=DimArray(fill(T_mean, m), (zdim,)),
        density=DimArray(fill(mp.density_ice, m), (zdim,)),
        water=DimArray(zeros(m), (zdim,)),
        grain_radius=DimArray(fill(2.5, m), (zdim,)),
        grain_dendricity=DimArray(ones(m), (zdim,)),
        grain_sphericity=DimArray(fill(0.5, m), (zdim,)),
        albedo=DimArray(fill(mp.albedo_snow, m), (zdim,)),
        albedo_diffuse=DimArray(fill(mp.albedo_snow, m), (zdim,)),
    ))
end

"""
    initialize_grid(mp::ModelParameters)

Generate the initial vertical grid layer thicknesses.
Matches MATLAB's `model_initialize_grid` (local function in model_initialize_profile.m).

Returns a Vector{Float64} of layer thicknesses from surface to depth.
"""
function initialize_grid(mp::ModelParameters)
    d_tolerance = 1e-11

    # Calculate number of top grid points
    n_top = mp.column_ztop / mp.column_dztop
    @assert mod(n_top, 1) == 0 "Top grid cell structure length does not go evenly into specified top structure depth."

    n_top = Int(n_top)

    if mp.column_dztop < 0.05 - d_tolerance
        @warn "Initial top grid cell length (column_dztop) is < 0.05 m."
    end

    # Initialize top grid (constant spacing)
    dzT = fill(mp.column_dztop, n_top)

    # Build bottom grid (geometrically stretched)
    dzB = Float64[]
    gp0 = mp.column_dztop
    z0 = mp.column_ztop

    while mp.column_zmax > (z0 + d_tolerance)
        dz_new = gp0 * mp.column_zy
        push!(dzB, dz_new)
        gp0 = dz_new
        z0 += gp0
    end

    # Combine top and bottom
    return vcat(dzT, dzB)
end
