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

    zdim = Z(1:m; metadata=cf_layer_index_attributes())

    # Exactly `RESTART_LAYERS`, taken from that constant rather than listed again, so what `gemb`
    # requires of a profile and what this writes into one cannot drift apart. `z_center` is
    # absent from it by design, so the extracted profile matches the output layout (a slice of
    # it); recompute via `dz2z(profile[:dz])` if needed. Carrying `age` is what makes it
    # accumulate across `gemb_spinup` cycles, which round-trip the column through here.
    layers = NamedTuple(l => DimArray(parent(out[l])[:, col_idx], (zdim,)) for l in RESTART_LAYERS)

    # Carry the CF attributes through, minus `cell_methods`: the extracted column has no
    # `Ti` dimension, so an interval reduction referencing a `time` coordinate that isn't
    # there would be wrong. Stack-level metadata is deliberately left unset — `gemb_spinup`
    # owns that slot (`rebuild(profile; metadata=prov)` in `spinup.jl`), and `rebuild`
    # preserves the layer metadata set here.
    return DimStack(layers; layermetadata=cf_layermetadata(layers; time_axis=false))
end
