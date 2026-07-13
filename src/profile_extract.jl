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

    @assert time_extract >= out_times[1] "time_extract cannot be before the first time step of output."
    @assert time_extract <= out_times[end] "time_extract cannot be after the last time step of output."

    # Find nearest time (use At for semantic time-based indexing)
    nearest_time = if time_extract == out_times[end]
        out_times[end]
    else
        _, idx = findmin(abs.(Dates.value.(out_times .- time_extract)))
        out_times[idx]
    end

    return _extract_profile_at_time(out, nearest_time)
end

function gemb_profile(out::DimStack)
    # Extract profile at last timestep
    out_times = collect(dims(out[:temperature], Ti))
    return _extract_profile_at_time(out, out_times[end])
end

"""
Extract profile at a specific time from the output DimStack using At() indexing.
"""
function _extract_profile_at_time(out::DimStack, time::DateTime)
    # Extract column using At() for semantic time indexing, removing NaN padding
    temp_col = out[:temperature][:, At(time)]
    valid = .!isnan.(temp_col)
    m = sum(valid)

    zdim = Z(1:m)

    # Extract albedo from profile output if available, otherwise use defaults
    if haskey(out, :albedo)
        albedo_col = out[:albedo][valid, At(time)]
        albedo_diffuse_col = out[:albedo_diffuse][valid, At(time)]
    else
        # Fallback for output without albedo profiles
        albedo_col = fill(0.85, m)
        albedo_diffuse_col = fill(0.85, m)
    end

    return DimStack((
        z_center=DimArray(dz2z(out[:dz][valid, At(time)]), (zdim,)),
        dz=DimArray(out[:dz][valid, At(time)], (zdim,)),
        temperature=DimArray(out[:temperature][valid, At(time)], (zdim,)),
        density=DimArray(out[:density][valid, At(time)], (zdim,)),
        water=DimArray(out[:water][valid, At(time)], (zdim,)),
        grain_radius=DimArray(out[:grain_radius][valid, At(time)], (zdim,)),
        grain_dendricity=DimArray(out[:grain_dendricity][valid, At(time)], (zdim,)),
        grain_sphericity=DimArray(out[:grain_sphericity][valid, At(time)], (zdim,)),
        albedo=DimArray(albedo_col, (zdim,)),
        albedo_diffuse=DimArray(albedo_diffuse_col, (zdim,)),
    ))
end
