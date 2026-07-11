"""
    gemb_interp(Z_center::AbstractMatrix, A::AbstractMatrix, profile::DimStack; interp_method=:linear)

Regularize GEMB Lagrangian output onto a consistent vertical grid.

This function is necessary because the vertical spacing of GEMB output
evolves with every timestep.

# Arguments
- `Z_center`: MxN matrix of grid cell center heights (from `dz2z(OutData.dz)`)
- `A`: MxN matrix of data to regrid (e.g., `OutData.temperature`)
- `profile`: Profile DimStack containing `z_center` field defining the target grid
- `interp_method`: Interpolation method (`:linear` default, `:nearest`)

# Returns
An `M_regularized x N` matrix where the vertical postings correspond to `profile[:z_center]`.

Matches MATLAB's `gemb_interp.m`.
"""
function gemb_interp(Z_center::AbstractMatrix, A::AbstractMatrix, profile::DimStack; interp_method::Symbol=:linear)
    @assert size(Z_center) == size(A) "Dimensions of Z_center and A must agree."
    @assert size(Z_center, 1) > 1 "Inputs Z_center and A must contain multiple rows representing profile depth."

    # Target z values from the profile
    z_target = collect(Float64, parent(profile[:z_center]))
    n_target = length(z_target)
    n_times = size(A, 2)

    # Preallocate output
    A_regularized = fill(NaN, n_target, n_times)

    # Loop through each timestep
    for k in 1:n_times
        # Get finite mask for this column
        isf = isfinite.(A[:, k])

        if sum(isf) < 2
            continue
        end

        # Source data for this column (only finite values)
        z_src = Z_center[isf, k]
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
