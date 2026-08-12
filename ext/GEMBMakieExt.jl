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
    :temperature => "temperature [°C]",
    :density => "density [kg m⁻³]",
    :water => "water [kg m⁻²]",
    :grain_radius => "grain radius [mm]",
    :grain_dendricity => "dendricity",
    :grain_sphericity => "sphericity",
    :albedo => "albedo",
    :albedo_diffuse => "diffuse albedo",
)
_label(v::Symbol) = get(_LABELS, v, string(v))

# Colormap per 2-D variable. Most are perceptually-uniform named maps; `water`
# and `density` use hand-built blue ramps so their endpoints carry meaning:
#   water   — white (no liquid water) → deep blue (maximum water content)
#   density — pale blue (low-density fresh snow) → medium blue (glacier ice)
# Dendricity and sphericity share the grain-radius map so the grain-metamorphism
# panels read as one family.
const _CMAP = Dict{Symbol,Any}(
    :temperature => :thermal,
    :density => cgrad([colorant"#eef4fb", colorant"#3b7dc4"]),
    :water => cgrad([colorant"#ffffff", colorant"#08306b"]),
    :grain_radius => :matter,
    :grain_dendricity => :matter,
    :grain_sphericity => :matter,
    :albedo => :ice,
    :albedo_diffuse => :ice,
)
_cmap(v::Symbol) = get(_CMAP, v, :viridis)

# Scalar variables grouped by physical theme + shared unit so multiple series
# share one axis (increasing data density, removing near-duplicate panels).
# Each entry: (panel title with units, [member variables]).
# Air temperature is placed first so it sits at the top of the right column,
# above (aligned with) the temperature profile heatmap in the left column.
const _SCALAR_GROUPS = [
    ("Air temperature [°C]", [:temperature_air]),
    ("Energy fluxes [W m⁻²]",
        [:shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent]),
    ("Mass fluxes [kg m⁻²]",
        [:melt, :runoff, :refreeze, :evaporation_condensation, :precipitation]),
    ("Surface albedo [–]", [:albedo_surface]),
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

# Map a median sampling cadence (hours) to an input run-frequency word. Falls
# back to a rounded interval when it doesn't match a common cadence.
function _freq_word(dt_h::Real)
    isnan(dt_h) && return "?"
    approx(target) = abs(dt_h - target) <= 0.05 * target
    approx(1) && return "hourly"
    approx(3) && return "3-hourly"
    approx(6) && return "6-hourly"
    approx(12) && return "12-hourly"
    approx(24) && return "daily"
    approx(24 * 7) && return "weekly"
    approx(24 * 30) && return "monthly"
    approx(24 * 365) && return "annual"
    return dt_h >= 24 ? "$(round(dt_h/24, digits=1))-daily" :
           "$(round(dt_h, digits=1))-hourly"
end

# Human-readable forcing provenance from the output stack metadata, e.g.
# "ERA5-Land @ 72.58°N, 38.46°W, +100 m". Every piece is optional: a stack with
# no provenance metadata yields "" and the banner renders as it did before.
function _provenance_str(md)
    md === nothing && return ""
    getmd(k) = md isa AbstractDict ? get(md, k, nothing) :
               (haskey(md, k) ? md[k] : nothing)

    parts = String[]
    ds = getmd("dataset")
    ds isa AbstractString && !isempty(ds) && push!(parts, ds)

    lat, lon = getmd("latitude"), getmd("longitude")
    if lat isa Real && lon isa Real && isfinite(lat) && isfinite(lon)
        lats = string(round(abs(lat), digits=2), "°", lat >= 0 ? "N" : "S")
        lons = string(round(abs(lon), digits=2), "°", lon >= 0 ? "E" : "W")
        coord = "$lats, $lons"
        isempty(parts) ? push!(parts, coord) : (parts[end] *= " @ " * coord)
    end

    off = getmd("elevation_offset")
    if off isa Real && isfinite(off) && abs(off) > 0
        push!(parts, string(off > 0 ? "+" : "", round(off, digits=1), " m"))
    end

    return join(parts, ", ")
end

# Spinup provenance from the output stack metadata: the climatology averaging
# window and whether the spinup converged, e.g.
# "spinup 1950-01-01→1980-12-31, converged". If the run was not spun up, reports
# "no spinup"; if no spinup metadata is present at all, yields "".
function _spinup_str(md)
    md === nothing && return ""
    getmd(k) = md isa AbstractDict ? get(md, k, nothing) :
               (haskey(md, k) ? md[k] : nothing)

    performed = getmd("spinup_performed")
    performed === false && return "no spinup"
    performed === nothing && return ""   # metadata predates provenance support

    parts = String[]
    t0, t1 = getmd("climatology_window_start"), getmd("climatology_window_stop")
    if t0 isa DateTime && t1 isa DateTime
        push!(parts, string("spinup ", Dates.format(t0, "yyyy-mm-dd"),
                            "→", Dates.format(t1, "yyyy-mm-dd")))
    else
        push!(parts, "spinup")
    end
    conv = getmd("spinup_converged")
    conv isa Bool && push!(parts, conv ? "converged" : "not converged")
    return join(parts, ", ")
end

# Average-annual surface mass-budget summary for the metadata banner, e.g.
# "↓ snow 420  ↓ rain 35  → runoff 180  ↺ refreeze 60  ↑ evap 12  [kg m⁻² yr⁻¹]".
# The flux outputs are per-output-step mass accumulations (kg m⁻²); summing the
# record and dividing by its length in years gives the annual mean. Leading
# arrows mark direction: ↓ mass input (falls onto/into the column), ↑ mass loss
# to the atmosphere (evaporation/sublimation), → runoff (liquid leaving laterally),
# ↺ internal phase change (refreeze). Evaporation/condensation is signed, so
# condensation reads as a ↓ input and evaporation/sublimation as a ↑ loss.
function _mass_budget_str(output, times)
    length(times) < 2 && return ""
    years = Dates.value(last(times) - first(times)) / (365.25 * 86400 * 1000)
    years <= 0 && return ""

    present = keys(output)
    total(v) = v in present ? sum(filter(isfinite, output[v])) : NaN
    ann(v) = total(v) / years

    parts = String[]
    # Snowfall is the model's precipitation minus its liquid (rain) fraction.
    if :precipitation in present && :rain in present
        snow = ann(:precipitation) - ann(:rain)
        isfinite(snow) && push!(parts, string("↓ snow ", round(Int, snow)))
    elseif :precipitation in present
        p = ann(:precipitation)
        isfinite(p) && push!(parts, string("↓ precip ", round(Int, p)))
    end
    if :rain in present
        r = ann(:rain)
        isfinite(r) && push!(parts, string("↓ rain ", round(Int, r)))
    end
    if :runoff in present
        ro = ann(:runoff)
        isfinite(ro) && push!(parts, string("→ runoff ", round(Int, ro)))
    end
    if :refreeze in present
        rf = ann(:refreeze)
        isfinite(rf) && push!(parts, string("↺ refreeze ", round(Int, rf)))
    end
    if :evaporation_condensation in present
        ec = ann(:evaporation_condensation)
        if isfinite(ec)
            # Positive = condensation/deposition (mass gain, ↓ input); negative =
            # evaporation/sublimation (mass loss to atmosphere, ↑).
            arrow = ec >= 0 ? "↓" : "↑"
            label = ec >= 0 ? "cond" : "evap"
            push!(parts, string(arrow, " ", label, " ", round(Int, abs(ec))))
        end
    end

    isempty(parts) && return ""
    return join(parts, "  ") * "  [kg m⁻² yr⁻¹]"
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
        # `dz` is a grid-management diagnostic and the depth-resolved albedo fields
        # duplicate the surface-albedo series; dendricity/sphericity are dropped by
        # default (grain radius already represents the grain-metamorphism family).
        default_drop = (:dz, :albedo, :albedo_diffuse, :grain_dendricity, :grain_sphericity)
        profile_vars = filter(v -> v ∉ default_drop, profile_vars)
    end
    scalar_groups = [(t, [v for v in g if v in wanted]) for (t, g) in _SCALAR_GROUPS]
    scalar_groups = [(t, g) for (t, g) in scalar_groups if !isempty(g)]
    # Densification is dropped from the default panel set (kept available when the
    # caller lists its variables explicitly via `variables`).
    if variables === nothing
        scalar_groups = [(t, g) for (t, g) in scalar_groups if t != "Densification [m]"]
    end

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
        vals = parent(output[v])[:, tcols]
        # Temperature is stored in K; show it in °C so the profile shares the
        # air-temperature panel's units.
        v === :temperature && (vals = vals .- 273.15)
        gridded = gemb_interp(z_center[:, tcols], vals, z_target)
        ax = Axis(left[i, 1]; ylabel="depth [m]")
        push!(prof_axes, ax)
        lo, hi = _robust_limits(parent(gridded))
        # Anchor water's ramp at zero so "no liquid water" reads as white.
        v === :water && (lo = 0.0)
        hm = heatmap!(ax, x_hm, collect(dims(gridded, Z)), parent(gridded');
            colormap=_cmap(v), colorrange=(lo, hi))
        Colorbar(left[i, 2], hm; label=_label(v), width=12)
        ylims!(ax, zbot, ztop)
    end

    # ---- Right column: grouped scalar time series -------------------------
    for (i, (gtitle, members)) in enumerate(scalar_groups)
        ax = Axis(right[i, 1]; ylabel=gtitle)
        push!(scal_axes, ax)
        # Air temperature: show in °C with the curve filled to zero — warm
        # (> 0 °C) in red, cold (≤ 0 °C) in blue — so melt-permissive periods
        # read at a glance.
        if members == [:temperature_air]
            tc = Float64.(parent(output[:temperature_air])) .- 273.15
            band!(ax, decyear, min.(tc, 0.0), 0.0;
                color=(:blue, 0.5))
            band!(ax, decyear, 0.0, max.(tc, 0.0);
                color=(:red, 0.5))
            continue
        end
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
    # Default to the exact data span (Makie otherwise auto-pads ~5% past the
    # ends); an explicit `datelims` overrides.
    xl = xlims === nothing ? (first(decyear), last(decyear)) : xlims
    for ax in all_axes
        xlims!(ax, xl...)
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
    # Median sampling cadence in hours, expressed as an input run frequency word
    # (e.g. "hourly", "daily") prefixed to the time span.
    dt_h = length(times) > 1 ?
           median(diff(Dates.value.(times))) / 3.6e6 : NaN
    freq = _freq_word(dt_h)
    strided = length(tcols) < length(times)
    md = DimensionalData.metadata(output)
    prov = _provenance_str(md)
    spin = _spinup_str(md)
    budget = _mass_budget_str(output, times)

    meta = string(
        "GEMB v$ver  │  ",
        prov == "" ? "" : prov * "  │  ",
        freq, " ",
        Dates.format(t0, "yyyy-mm-dd"), " → ", Dates.format(t1, "yyyy-mm-dd"),
        spin == "" ? "" : "  │  " * spin,
        budget == "" ? "" : "  │  " * budget,
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
