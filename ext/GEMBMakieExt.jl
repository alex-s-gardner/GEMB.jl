"""
Extension for GEMB.jl providing diagnostic plotting via Makie.

This extension is automatically loaded when both GEMB and a Makie backend
(e.g. CairoMakie or GLMakie) are loaded. It implements [`gemb_plot_output`](@ref).
"""
module GEMBMakieExt

using GEMB
using DimensionalData
using Dates
using Statistics
using Makie

#=============================================================================
# Metadata: labels, colormaps, and the grouping of scalar variables.
=============================================================================#

# Short, unit-bearing labels for scalar variables (used as y-axis labels or
# legend entries). Falls back to the raw name when a key is missing.
const _LABELS = Dict{Symbol,String}(
    :melt => "melt",
    :runoff => "runoff",
    :refreeze => "refreeze",
    :evaporation_condensation => "evap/cond",
    :precipitation => "precip",
    :shortwave_net => "SW net",
    :longwave_net => "LW net",
    :heat_flux_sensible => "sensible",
    :heat_flux_latent => "latent",
    :albedo_surface => "surface albedo",
    :temperature_air => "air temperature",
    :firn_air_content => "firn air content",
    :thickness_cumulative => "cumulative Δ thickness",
    :densification_from_compaction => "compaction",
    :densification_from_melt => "melt",
    :valid_profile_length => "valid layers",
    # 2-D (profile) variables
    :temperature => "temperature [K]",
    :density => "density [kg m⁻³]",
    :water => "water [kg m⁻²]",
    :grain_radius => "grain radius [mm]",
    :grain_dendricity => "dendricity",
    :grain_sphericity => "sphericity",
    :albedo => "albedo",
    :albedo_diffuse => "diffuse albedo",
)
_label(v::Symbol) = get(_LABELS, v, string(v))

# Perceptually-uniform colormap per 2-D variable, chosen so each field reads at
# maximum contrast (cmocean maps: thermal/dense/haline/matter/ice).
const _CMAP = Dict{Symbol,Symbol}(
    :temperature => :thermal,
    :density => :dense,
    :water => :haline,
    :grain_radius => :matter,
    :grain_dendricity => :speed,
    :grain_sphericity => :speed,
    :albedo => :ice,
    :albedo_diffuse => :ice,
)
_cmap(v::Symbol) = get(_CMAP, v, :viridis)

# Scalar variables grouped by physical theme + shared unit so multiple series
# share one axis (increasing data density, removing near-duplicate panels).
# Each entry: (panel title with units, [member variables]).
const _SCALAR_GROUPS = [
    ("Energy fluxes [W m⁻²]",
        [:shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent]),
    ("Mass fluxes [kg m⁻²]",
        [:melt, :runoff, :refreeze, :evaporation_condensation, :precipitation]),
    ("Surface albedo [–]", [:albedo_surface]),
    ("Air temperature [K]", [:temperature_air]),
    ("Firn air content [m]", [:firn_air_content]),
    ("Densification [m]", [:densification_from_compaction, :densification_from_melt]),
]

#=============================================================================
# Helpers
=============================================================================#

# Convert a date-limit endpoint (DateTime or decimal year) to a decimal year.
_to_decyear(x::DateTime) = GEMB.datetime2decyear([x])[1]
_to_decyear(x::Real) = Float64(x)

# Robust color limits: clip to the (plo, phi) percentiles of finite data so a
# handful of outliers don't wash out the map. For a (near-)constant field, widen
# minimally around the value so the colorbar reflects its true magnitude rather
# than an arbitrary span.
function _robust_limits(A::AbstractArray; plo=2, phi=98)
    v = filter(isfinite, vec(A))
    isempty(v) && return (0.0, 1.0)
    lo, hi = quantile(v, plo / 100), quantile(v, phi / 100)
    # A (near-)constant field gives lo ≈ hi. Widen the range so it survives
    # Makie's Float32 color pipeline, which errors when cmin == cmax (and can
    # collapse a range narrower than Float32 precision, e.g. 0.85 ± 1e-16).
    if hi - lo <= abs(hi) * 1e-4 + 1e-9
        pad = max(abs(lo) * 1e-3, 1e-6)
        return (lo - pad, hi + pad)
    end
    return (lo, hi)
end

# Even index stride so a heatmap never exceeds `maxcols` columns — a raster
# figure a few thousand px wide cannot resolve more, and full width would make
# rendering (and memory) explode for multi-decade, sub-daily runs.
function _time_stride(ntime::Int, maxcols::Int)
    ntime <= maxcols && return 1:ntime
    step = cld(ntime, maxcols)
    return 1:step:ntime
end

#=============================================================================
# Main entry point
=============================================================================#

function GEMB.gemb_plot_output(output::DimStack;
    datelims=nothing, depthlims=(-10, 0), variables=nothing, title="")

    # ---- Select variables -------------------------------------------------
    present = collect(keys(output))
    wanted = variables === nothing ? present : intersect(present, collect(variables))
    profile_vars = [v for v in wanted if ndims(output[v]) == 2]
    # Drop from the default profile column unless the caller asks for a variable
    # explicitly: `dz` is a grid-management diagnostic (not a physical field), and
    # the depth-resolved albedo fields duplicate the surface-albedo time series.
    if variables === nothing
        default_drop = (:dz, :albedo, :albedo_diffuse)
        profile_vars = filter(v -> v ∉ default_drop, profile_vars)
    end
    scalar_groups = [(t, [v for v in g if v in wanted]) for (t, g) in _SCALAR_GROUPS]
    scalar_groups = [(t, g) for (t, g) in scalar_groups if !isempty(g)]

    # ---- Shared time axis (decimal year) ----------------------------------
    times = collect(dims(output, Ti))
    decyear = GEMB.datetime2decyear(times)
    xlims = datelims === nothing ? nothing :
            (_to_decyear(datelims[1]), _to_decyear(datelims[2]))

    # ---- Fixed vertical grid for regridding 2-D variables -----------------
    # Range from the evolving cell centres (or `depthlims`); resolution is a
    # bounded render target, not the (padded) layer count.
    z_center = dz2z(parent(output[:dz]))
    if depthlims === nothing
        finite = filter(isfinite, z_center)
        ztop, zbot = maximum(finite), minimum(finite)
    else
        ztop, zbot = maximum(depthlims), minimum(depthlims)
    end
    nz = 240
    z_target = collect(range(ztop, zbot; length=nz))  # descending: surface first

    # Heatmaps: stride time so the raster stays a sane width.
    tcols = _time_stride(length(times), 2400)
    x_hm = decyear[tcols]

    # ---- Figure scaffold --------------------------------------------------
    nprof, ngrp = length(profile_vars), length(scalar_groups)
    nrows = max(nprof, ngrp, 1)
    fig = Figure(size=(1500, 200 + 175 * nrows), fontsize=14)

    # Header: title + traceability metadata banner.
    _header!(fig[1, 1:2], output, times, z_center, tcols, title)

    body = fig[2, 1:2] = GridLayout()
    left = body[1, 1] = GridLayout()   # profile heatmaps
    right = body[1, 2] = GridLayout()  # grouped scalar series
    colsize!(body, 1, Relative(0.56))
    colgap!(body, 28)

    prof_axes = Axis[]
    scal_axes = Axis[]

    # ---- Left column: profile heatmaps ------------------------------------
    for (i, v) in enumerate(profile_vars)
        gridded = gemb_interp(z_center[:, tcols], parent(output[v])[:, tcols], z_target)
        ax = Axis(left[i, 1]; ylabel="depth [m]")
        push!(prof_axes, ax)
        lo, hi = _robust_limits(parent(gridded))
        hm = heatmap!(ax, x_hm, collect(dims(gridded, Z)), parent(gridded');
            colormap=_cmap(v), colorrange=(lo, hi))
        Colorbar(left[i, 2], hm; label=_label(v), width=12)
        ylims!(ax, zbot, ztop)
    end

    # ---- Right column: grouped scalar time series -------------------------
    for (i, (gtitle, members)) in enumerate(scalar_groups)
        ax = Axis(right[i, 1]; ylabel=gtitle)
        push!(scal_axes, ax)
        # With ~10⁵ sub-daily points, overplotted opaque lines saturate to a
        # solid block; thin, semi-transparent strokes let density show through.
        alpha = length(members) > 1 ? 0.75 : 0.9
        for v in members
            lines!(ax, decyear, Float64.(parent(output[v]));
                label=_label(v), linewidth=0.6, alpha=alpha)
        end
        length(members) > 1 && axislegend(ax; position=:rt, framevisible=false,
            padding=(4, 4, 2, 2), rowgap=0, labelsize=11, patchsize=(12, 8))
    end

    # ---- Link x-axes; strip redundant decorations -------------------------
    all_axes = vcat(prof_axes, scal_axes)
    !isempty(all_axes) && linkxaxes!(all_axes...)
    if xlims !== nothing
        for ax in all_axes
            xlims!(ax, xlims...)
        end
    end
    # Only the bottom axis of each column keeps its x ticks + label.
    _label_bottom_only!(prof_axes)
    _label_bottom_only!(scal_axes)

    return fig
end

# Title + one-line metadata banner for traceability.
function _header!(pos, output, times, z_center, tcols, title)
    ver = try
        string(pkgversion(GEMB))
    catch
        "?"
    end
    t0, t1 = first(times), last(times)
    # Median sampling cadence in hours.
    dt_h = length(times) > 1 ?
           median(diff(Dates.value.(times))) / 3.6e6 : NaN
    cad = isnan(dt_h) ? "?" :
          dt_h >= 24 ? "$(round(dt_h/24, digits=1)) d" : "$(round(dt_h, digits=1)) h"
    finite = filter(isfinite, z_center)
    zmax = isempty(finite) ? NaN : -minimum(finite)
    strided = length(tcols) < length(times)

    meta = string(
        "GEMB v$ver  │  ",
        Dates.format(t0, "yyyy-mm-dd"), " → ", Dates.format(t1, "yyyy-mm-dd"),
        "  │  ", length(times), " steps (~", cad, ")",
        "  │  max depth ", round(zmax, digits=1), " m",
        strided ? "  │  heatmaps strided to $(length(tcols)) cols" : "",
        "  │  generated ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"),
    )
    head = title == "" ? "GEMB diagnostic output" : title
    Label(pos, rich(head * "\n", fontsize=20, font=:bold,
            rich(meta, fontsize=11, color=:gray40, font=:regular));
        justification=:left, halign=:left, padding=(4, 0, 6, 0))
    return nothing
end

# Keep x ticks/label only on the bottom-most axis of a linked column.
function _label_bottom_only!(axes::Vector{Axis})
    isempty(axes) && return
    for ax in axes[1:end-1]
        hidexdecorations!(ax; grid=false)
    end
    axes[end].xlabel = "year"
    return nothing
end

end # module
