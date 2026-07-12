"""
    gemb_profile(out::DimStack)
    gemb_profile(out::DimStack, time_extract::DateTime)

Extract a column state from GEMB output as a Profile DimStack.

If `time_extract` is not provided, the last time step is used.
If `time_extract` does not exactly match any output time, the nearest time step is used.

Matches MATLAB's `gemb_profile.m`.
"""
function gemb_profile(out::DimStack, time_extract::DateTime)
    # Get output times
    out_times = collect(dims(out[:temperature], Ti))
    n_times = length(out_times)

    @assert time_extract >= out_times[1] "time_extract cannot be before the first time step of output."
    @assert time_extract <= out_times[end] "time_extract cannot be after the last time step of output."

    # Find nearest time index
    if time_extract == out_times[end]
        col_idx = n_times
    else
        _, col_idx = findmin(abs.(Dates.value.(out_times .- time_extract)))
    end

    return _extract_profile_at_index(out, col_idx)
end

function gemb_profile(out::DimStack)
    n_times = size(out[:temperature], 2)
    return _extract_profile_at_index(out, n_times)
end

"""
Extract profile at a specific column index from the output DimStack.
"""
function _extract_profile_at_index(out::DimStack, col_idx::Int)
    # Extract column, removing NaN padding
    temp_col = parent(out[:temperature])[:, col_idx]
    valid = .!isnan.(temp_col)
    m = sum(valid)

    zdim = Z(1:m)

    # Extract albedo from profile output if available, otherwise use defaults
    if haskey(out, :albedo)
        albedo_col = parent(out[:albedo])[valid, col_idx]
        albedo_diffuse_col = parent(out[:albedo_diffuse])[valid, col_idx]
    else
        # Fallback for output without albedo profiles
        albedo_col = fill(0.85, m)
        albedo_diffuse_col = fill(0.85, m)
    end

    return DimStack((
        z_center=DimArray(dz2z(parent(out[:dz])[valid, col_idx]), (zdim,)),
        dz=DimArray(parent(out[:dz])[valid, col_idx], (zdim,)),
        temperature=DimArray(parent(out[:temperature])[valid, col_idx], (zdim,)),
        density=DimArray(parent(out[:density])[valid, col_idx], (zdim,)),
        water=DimArray(parent(out[:water])[valid, col_idx], (zdim,)),
        grain_radius=DimArray(parent(out[:grain_radius])[valid, col_idx], (zdim,)),
        grain_dendricity=DimArray(parent(out[:grain_dendricity])[valid, col_idx], (zdim,)),
        grain_sphericity=DimArray(parent(out[:grain_sphericity])[valid, col_idx], (zdim,)),
        albedo=DimArray(albedo_col, (zdim,)),
        albedo_diffuse=DimArray(albedo_diffuse_col, (zdim,)),
    ))
end
