"""
    gemb_interp(z_center::AbstractMatrix, A::AbstractMatrix, z_target::AbstractVector; interp_method=:linear)

Regularize GEMB Lagrangian output onto a consistent vertical grid.

This function is necessary because the vertical spacing of GEMB output
evolves with every timestep.

# Arguments
- `z_center`: MxN matrix of grid cell center heights (from `dz2z(OutData.dz)`)
- `A`: MxN matrix (or `DimArray` with `Z` and `Ti` dimensions) of data to regrid
- `z_target`: Vector of target depth coordinates defining the output grid
- `interp_method`: Interpolation method (`:linear` default, `:nearest`)

# Returns
A `DimArray` of size `length(z_target) × N` with `Z(z_target)` and `Ti` dimensions.
When `A` is a `DimArray` with a `Ti` dimension, the time coordinates are preserved.

Target heights outside the column — above the top cell centre or below the bottom cell
centre — are **not extrapolated**; they are left as the missing value `NaN`. The column
depth and the surface height both migrate over a run (the surface cell is centimetres
thick and the base is pinned only to `z_target`'s own tolerance), so a fixed `z_target`
range like `(-20, 0)` inevitably reaches past the column. Extrapolating there continued
the near-basal temperature gradient into rock and produced firn temperatures above the
melt point — values the model itself never holds. `NaN` says "no column here" instead,
and Makie renders it as blank.

"""
function gemb_interp(z_center::AbstractMatrix, A::AbstractMatrix, z_target::AbstractVector; interp_method::Symbol=:linear)
    @assert size(z_center) == size(A) "Dimensions of z_center and A must agree."
    @assert size(z_center, 1) > 1 "Inputs z_center and A must contain multiple rows representing profile depth."

    n_times = size(A, 2)

    # Extract time coordinates if A carries them
    ti_dim = A isa DimArray ? dims(A, Ti) : Ti(1:n_times)

    # Preallocate output. Unlike the `Z` of `gemb` output — a bare cell index — `z_target`
    # is a genuine vertical coordinate, so it gets real CF coordinate attributes. The
    # regridded values are the same quantity as the input, so carry the input's own
    # attributes across; a plain matrix input simply has none to carry.
    zdim = Z(z_target; metadata=cf_height_attributes())
    A_regularized = DimArray(fill(NaN, length(z_target), n_times), (zdim, ti_dim);
        metadata=A isa DimArray ? DD.metadata(A) : DD.NoMetadata())

    # Loop through each timestep. GEMB profile output is top-justified with a fixed row
    # count, so every row is a real cell; the finite check only guards against genuinely
    # missing data (e.g. an output array sliced before it was fully written).
    for k in 1:n_times
        isf = isfinite.(A[:, k]) .& isfinite.(view(z_center, :, k))

        if count(isf) < 2
            continue
        end

        z_src = z_center[isf, k]
        a_src = A[isf, k]

        if interp_method == :linear
            A_regularized[:, k] = _interp1(z_src, a_src, z_target)
        elseif interp_method == :nearest
            A_regularized[:, k] = _interp1_nearest(z_src, a_src, z_target)
        else
            error("interp_method must be :linear or :nearest")
        end
    end

    return A_regularized
end

"""
    _interp1(x, y, xi)

1D linear interpolation, **without** extrapolation: query points outside `[min(x), max(x)]`
return `NaN` rather than a continuation of the end interval. `x` must be sorted (ascending
or descending). See [`gemb_interp`](@ref) for why the ends are missing rather than
extrapolated.
"""
function _interp1(x::AbstractVector, y::AbstractVector, xi::AbstractVector)
    n = length(x)
    result = similar(xi, Float64)

    # Determine sort order - x may be ascending or descending
    ascending = x[end] > x[1]

    for i in eachindex(xi)
        xq = xi[i]

        if ascending
            # Find bracketing interval
            if xq < x[1] || xq > x[end]
                # Outside the data: no extrapolation.
                result[i] = NaN
                continue
            elseif xq == x[1]
                idx_lo, idx_hi = 1, 2
            elseif xq == x[end]
                idx_lo, idx_hi = n - 1, n
            else
                # Binary search for interval
                idx_hi = searchsortedfirst(x, xq)
                idx_lo = idx_hi - 1
            end
        else
            # Descending x
            if xq > x[1] || xq < x[end]
                # Outside the data: no extrapolation.
                result[i] = NaN
                continue
            elseif xq == x[1]
                idx_lo, idx_hi = 1, 2
            elseif xq == x[end]
                idx_lo, idx_hi = n - 1, n
            else
                # Linear search (descending)
                idx_lo = 1
                for j in 1:(n-1)
                    if x[j] >= xq >= x[j+1]
                        idx_lo = j
                        break
                    end
                end
                idx_hi = idx_lo + 1
            end
        end

        # Linear interpolation/extrapolation
        dx = x[idx_hi] - x[idx_lo]
        if abs(dx) < eps(Float64)
            result[i] = y[idx_lo]
        else
            t = (xq - x[idx_lo]) / dx
            result[i] = y[idx_lo] + t * (y[idx_hi] - y[idx_lo])
        end
    end

    return result
end

"""
    _interp1_nearest(x, y, xi)

1D nearest-neighbor interpolation, **without** extrapolation: query points outside
`[min(x), max(x)]` return `NaN`, matching [`_interp1`](@ref).
"""
function _interp1_nearest(x::AbstractVector, y::AbstractVector, xi::AbstractVector)
    result = similar(xi, Float64)
    x_lo, x_hi = extrema(x)
    for i in eachindex(xi)
        if xi[i] < x_lo || xi[i] > x_hi
            result[i] = NaN
            continue
        end
        _, idx = findmin(abs.(x .- xi[i]))
        result[i] = y[idx]
    end
    return result
end
