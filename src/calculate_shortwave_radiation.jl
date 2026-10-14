"""
    calculate_shortwave_radiation(dz, density, grain_radius, albedo_broadband, albedo_diffuse, cfs::ClimateForcingStep, mp::ModelParameters)

Distribute absorbed shortwave radiation vertically within snow/ice.

Depending on model configuration:
1. Surface Absorption: All net shortwave energy is absorbed by the top grid cell
   (`shortwave_subsurface_absorption = false`).
2. Subsurface Penetration: Shortwave energy penetrates and is absorbed by deeper layers
   (`shortwave_subsurface_absorption = true`), using either:
   - Density-dependent extinction (Bassford, 2002)
   - Spectral-dependent extinction for "BrunLefebre" (Lefebre et al., 2003)

Returns `shortwave_flux` vector [W m-2] of absorbed shortwave radiation per grid cell.

Every branch absorbs the same total, the net surface shortwave of
[`shortwave_net_surface`](@ref): `albedo_broadband` applied to the direct beam
`SW↓ − SW↓_diffuse` and `albedo_diffuse` to `SW↓_diffuse`. The branches differ only in how that
total is distributed with depth, and on the subsurface branches the part transmitted through
the base of the column is not absorbed:

- surface absorption: all of it in the top cell;
- density-dependent extinction: 36% in the top cell (wavelengths > 0.8 µm), the remaining
  64% attenuated as `exp(−B·z)` with `B` linear in density from 10 m⁻¹ at 300 kg m⁻³ to
  1.3 m⁻¹ at `mp.density_ice`;
- `:BrunLefebre` spectral extinction: the net is split over three bands in proportion to
  `sF_k·(1 − a_k)`, the spectral fractions and spectral albedos of
  [`_albedo_brun_lefebre`](@ref); the 1.5–2.8 µm band is absorbed in the top cell and the
  other two are attenuated with their own grain-size-dependent coefficients.

# References
- Lefebre, F., Gallée, H., van Ypersele, J.-P., and Greuell, W. (2003). Modeling of snow and
  ice melt at ETH Camp (West Greenland): a study of surface albedo. *J. Geophys. Res.* 108,
  4231.
- Bassford, R. P. (2002). *Geophysical and numerical modelling investigations of the ice caps
  on Severnaya Zemlya*. PhD thesis, University of Bristol. Source of the `Bs`/`Bi`
  extinction coefficients of the density-dependent branch. Inherited from the MATLAB
  implementation, which cites this work as "Bassford, 2002" for the formulation and
  "Bassford, 2006" for the same two coefficients; the thesis is the only work of that name
  carrying them, so the date is unified here. The 2006 *Arct. Antarct. Alp. Res.* papers of
  Bassford et al. are the published mass-balance results of the same modelling and do not
  tabulate these coefficients.
"""
calculate_shortwave_radiation(dz::Vector{Float64}, density::Vector{Float64},
    grain_radius::Vector{Float64}, albedo_broadband::Float64, albedo_diffuse::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters) =
    calculate_shortwave_radiation!(zeros(length(density)), dz, density, grain_radius,
        albedo_broadband, albedo_diffuse, cfs, mp)

"""
    calculate_shortwave_radiation!(shortwave_flux, dz, density, grain_radius,
                                   albedo_broadband, albedo_diffuse, cfs, mp) -> shortwave_flux

Write the absorbed-shortwave profile of [`calculate_shortwave_radiation`](@ref) into
`shortwave_flux`, which must be at least the column length. It is zeroed in full on entry, because
the surface-absorption branches set only the top cell and leave the rest at zero — which also means a
buffer longer than the column contributes nothing to the caller's `sum`.
"""
function calculate_shortwave_radiation!(shortwave_flux::AbstractVector, dz::Vector{Float64},
    density::Vector{Float64},
    grain_radius::Vector{Float64}, albedo_broadband::Float64,
    albedo_diffuse::Float64,
    cfs::ClimateForcingStep, mp::ModelParameters)

    m = length(density)
    fill!(shortwave_flux, 0.0)

    shortwave_net = shortwave_net_surface(albedo_broadband, albedo_diffuse, cfs)

    if (!mp.shortwave_subsurface_absorption) ||
       ((mp.density_ice - density[1]) < D_TOLERANCE)
        # all sw radiation is absorbed by the top grid cell
        shortwave_flux[1] = shortwave_net

    else  # sw radiation is absorbed at depth within the glacier

        if mp.albedo_method == :BrunLefebre
            # convert effective radius [mm] to grain size [m]
            gsz = (grain_radius .* 2) ./ 1000

            # Spectral fractions and surface band albedos, which set how the net is shared
            # between the bands.
            sF = BRUN_LEFEBRE_SPECTRAL_FRACTIONS
            a1, a2, a3 = _brun_lefebre_band_albedos(grain_radius[1])
            w = (sF[1] * (1 - a1), sF[2] * (1 - a2), sF[3] * (1 - a3))
            w_total = w[1] + w[2] + w[3]
            (all(>=(0.0), w) && w_total > 0.0) ||
                error("BrunLefebre spectral absorption weights must be non-negative with a " *
                      "positive sum: $(w) at grain radius $(grain_radius[1]) mm")
            swfS = (shortwave_net * w[1] / w_total, shortwave_net * w[2] / w_total,
                shortwave_net * w[3] / w_total)

            # absorption coefficient for spectral range
            h = density ./ (gsz .^ 0.5)
            B1 = 0.0192 .* h                 # 0.3 - 0.8um
            B2 = 0.1098 .* h                 # 0.8 - 1.5um

            # cumulative extinction factors
            B1_cum = ones(m + 1)
            B2_cum = ones(m + 1)
            _extinction_profile!(B1_cum, B1, dz)
            _extinction_profile!(B2_cum, B2, dz)

            # flux across grid cell boundaries
            Qs1 = swfS[1] .* B1_cum
            Qs2 = swfS[2] .* B2_cum

            # net energy flux to each grid cell
            @inbounds for i in 1:m
                shortwave_flux[i] = (Qs1[i] - Qs1[i+1]) + (Qs2[i] - Qs2[i+1])
            end

            # add flux absorbed at surface
            shortwave_flux[1] = shortwave_flux[1] + swfS[3]

        else  # function of grid cell density
            # fraction of sw radiation absorbed in top grid cell (wavelength > 0.8um)
            SWs = 0.36

            # calculate surface shortwave radiation fluxes [W m-2]
            swf_s = SWs * shortwave_net

            # calculate subsurface shortwave radiation fluxes [W m-2]
            swf_ss = (1 - SWs) * shortwave_net

            # SW extinction coefficients
            Bs = 10.0    # snow SW extinction coefficient [m-1] (Bassford, 2002)
            Bi = 1.3     # ice SW extinction coefficient [m-1] (Bassford, 2002)

            # calculate extinction coefficient B [m-1] vector
            B = Bs .+ (300 .- density) .* ((Bs - Bi) / (mp.density_ice - 300))

            # cumulative extinction factor
            B_cum = ones(m + 1)
            _extinction_profile!(B_cum, B, dz)

            # flux across grid cell boundaries
            Qs = swf_ss .* B_cum

            # net energy flux to each grid cell
            @inbounds for i in 1:m
                shortwave_flux[i] = Qs[i] - Qs[i+1]
            end

            # add flux absorbed at surface
            shortwave_flux[1] = shortwave_flux[1] + swf_s
        end
    end

    return shortwave_flux
end

"""
    shortwave_net_surface(albedo, albedo_diffuse, cfs) -> SW_net [W m-2]

Shortwave absorbed by the surface, `(1 − albedo)·(SW↓ − SW↓_diffuse) +
(1 − albedo_diffuse)·SW↓_diffuse`. Throws when `SW↓_diffuse > SW↓`.
"""
@inline function shortwave_net_surface(albedo::Float64, albedo_diffuse::Float64,
    cfs::ClimateForcingStep)
    direct = cfs.shortwave_downward - cfs.shortwave_downward_diffuse
    direct >= 0.0 || error("shortwave_downward_diffuse ($(cfs.shortwave_downward_diffuse) " *
                           "W m-2) exceeds shortwave_downward ($(cfs.shortwave_downward) W m-2)")
    return (1.0 - albedo) * direct + (1.0 - albedo_diffuse) * cfs.shortwave_downward_diffuse
end

"""
    effective_albedo(albedo, albedo_diffuse, cfs) -> α

Shortwave-weighted surface albedo, `1 − SW_net/SW↓` with `SW_net` from
[`shortwave_net_surface`](@ref): the reflected fraction of the incident shortwave. When
`SW↓ = 0` there is nothing to weight by, and the direct-beam `albedo` is returned.
"""
@inline function effective_albedo(albedo::Float64, albedo_diffuse::Float64,
    cfs::ClimateForcingStep)
    sw = cfs.shortwave_downward
    sw > 0.0 || return albedo
    return 1.0 - shortwave_net_surface(albedo, albedo_diffuse, cfs) / sw
end

"""
    _extinction_profile!(out, B, dz) -> out

Write the cumulative transmission down the column into `out[2:end]`, leaving `out[1]` as the
caller set it (1.0 — no attenuation above the surface).

`out[i+1] = prod(exp(-B[j] * dz[j]) for j in 1:i)`, accumulated as a running product in the
same order, and so to the same last bit, as `cumprod`. `Base.cumprod` is avoided deliberately:
it reaches `Base._accumulate!`'s `dims` machinery, whose `CartesianIndices` construction the
`--trim` verifier cannot resolve, which would take the whole column-step call tree out of the
statically compiled C API (see `capi/README.md`).
"""
function _extinction_profile!(out::AbstractVector{Float64}, B::AbstractVector{Float64},
    dz::AbstractVector{Float64})
    axes(B) == axes(dz) ||
        throw(DimensionMismatch("B and dz must match: $(axes(B)) vs $(axes(dz))"))
    length(out) == length(B) + 1 ||
        throw(DimensionMismatch("out must be one longer than B: $(length(out)) vs $(length(B))"))

    running = 1.0
    o = firstindex(out)
    for i in eachindex(B, dz)
        running *= exp(-B[i] * dz[i])
        out[o+1] = running
        o += 1
    end
    return out
end
