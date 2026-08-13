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
    # Profile output is top-justified with a fixed row count, so every row is a real cell
    # and the whole column is taken as-is.
    m = size(parent(out[:temperature]), 1)

    zdim = Z(1:m)

    # Extract albedo from profile output if available, otherwise use defaults
    if haskey(out, :albedo)
        albedo_col = parent(out[:albedo])[:, col_idx]
        albedo_diffuse_col = parent(out[:albedo_diffuse])[:, col_idx]
    else
        # Fallback for output without albedo profiles
        albedo_col = fill(0.85, m)
        albedo_diffuse_col = fill(0.85, m)
    end

    # `z_center` is intentionally omitted so the extracted profile matches the
    # output layout (a slice of it); recompute via `dz2z(profile[:dz])` if needed.
    return DimStack((
        dz=DimArray(parent(out[:dz])[:, col_idx], (zdim,)),
        temperature=DimArray(parent(out[:temperature])[:, col_idx], (zdim,)),
        density=DimArray(parent(out[:density])[:, col_idx], (zdim,)),
        water=DimArray(parent(out[:water])[:, col_idx], (zdim,)),
        grain_radius=DimArray(parent(out[:grain_radius])[:, col_idx], (zdim,)),
        grain_dendricity=DimArray(parent(out[:grain_dendricity])[:, col_idx], (zdim,)),
        grain_sphericity=DimArray(parent(out[:grain_sphericity])[:, col_idx], (zdim,)),
        albedo=DimArray(albedo_col, (zdim,)),
        albedo_diffuse=DimArray(albedo_diffuse_col, (zdim,)),
    ))
end
