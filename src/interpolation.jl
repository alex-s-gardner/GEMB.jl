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

Matches MATLAB's `gemb_interp.m`.
"""
function gemb_interp(z_center::AbstractMatrix, A::AbstractMatrix, z_target::AbstractVector; interp_method::Symbol=:linear)
    @assert size(z_center) == size(A) "Dimensions of z_center and A must agree."
    @assert size(z_center, 1) > 1 "Inputs z_center and A must contain multiple rows representing profile depth."

    n_times = size(A, 2)

    # Extract time coordinates if A carries them
    ti_dim = A isa DimArray ? dims(A, Ti) : Ti(1:n_times)

    # Preallocate output
    A_regularized = fill(NaN, Z(z_target), ti_dim)

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
            A_regularized[:, k] = _interp1_extrap(z_src, a_src, z_target)
        elseif interp_method == :nearest
            A_regularized[:, k] = _interp1_nearest(z_src, a_src, z_target)
        else
            error("interp_method must be :linear or :nearest")
        end
    end

    return A_regularized
end

"""
    _interp1_extrap(x, y, xi)

1D linear interpolation with extrapolation.
x must be sorted (ascending or descending). Extrapolates beyond bounds.
"""
function _interp1_extrap(x::AbstractVector, y::AbstractVector, xi::AbstractVector)
    n = length(x)
    result = similar(xi, Float64)

    # Determine sort order - x may be ascending or descending
    ascending = x[end] > x[1]

    for i in eachindex(xi)
        xq = xi[i]

        if ascending
            # Find bracketing interval
            if xq <= x[1]
                # Extrapolate below
                idx_lo, idx_hi = 1, 2
            elseif xq >= x[end]
                # Extrapolate above
                idx_lo, idx_hi = n - 1, n
            else
                # Binary search for interval
                idx_hi = searchsortedfirst(x, xq)
                idx_lo = idx_hi - 1
            end
        else
            # Descending x
            if xq >= x[1]
                idx_lo, idx_hi = 1, 2
            elseif xq <= x[end]
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

1D nearest-neighbor interpolation with extrapolation.
"""
function _interp1_nearest(x::AbstractVector, y::AbstractVector, xi::AbstractVector)
    result = similar(xi, Float64)
    for i in eachindex(xi)
        _, idx = findmin(abs.(x .- xi[i]))
        result[i] = y[idx]
    end
    return result
end
