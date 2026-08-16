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
# Metadata: labels, units, colormaps, and the grouping of scalar variables.
#
# Units are *not* duplicated here. They are read from each layer's CF metadata
# (`GEMB.GEMB_CF_ATTRIBUTES`, attached by `gemb`), so this file holds only what
# is genuinely a plotting concern: short display names sized for cramped axes and
# legends, and the display units a reader expects to see (°C rather than K).
=============================================================================#

# Short display names. A CF `long_name` is a sentence-length description; a legend
# entry has room for two words, so these are separate by design. Falls back to the
# raw variable name when a key is missing.
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
    :albedo_broadband => "broadband albedo",
    :temperature_air => "air temperature",
    :firn_air_content => "firn air content",
    :firn_air_content_10m => "top 10 m",
    :firn_air_content_20m => "top 20 m",
    :percolation_depth => "percolation depth",
    :ice_slab_thickness => "slab thickness",
    :ice_slab_depth => "depth to slab",
    :aquifer_thickness => "saturated thickness",
    :aquifer_depth => "water table",
    :thickness_cumulative => "cumulative Δ thickness",
    :densification_from_compaction => "compaction",
    :densification_from_melt => "melt",
    :strain_thinning => "horizontal strain",
    :valid_profile_length => "valid layers",
    # 2-D (profile) variables
    :temperature => "temperature",
    :density => "density",
    :water => "water",
    :grain_radius => "grain radius",
    :grain_dendricity => "dendricity",
    :grain_sphericity => "sphericity",
)
_short(v::Symbol) = get(_LABELS, v, string(v))

# CF units → (display unit, value conversion). GEMB stores temperature in kelvin,
# which is right for the model and for CF, but glaciological plots are read in °C —
# so the conversion lives here, keyed by unit rather than by variable name, and
# applies to every temperature field automatically. Units absent from this table
# are displayed as stored.
const _DISPLAY_UNITS = Dict{String,Tuple{String,Any}}(
    "K" => ("°C", x -> x - 273.15),
)

# CF units are ASCII/UDUNITS ("kg m-3"); render the exponents as superscripts for
# display, and show the dimensionless unit "1" as an en dash.
const _SUPERSCRIPTS = Dict{Char,Char}(
    '-' => '⁻', '0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
)
function _pretty_units(u::AbstractString)
    isempty(u) && return ""
    u == "1" && return "–"
    # Superscript a trailing exponent on each whitespace-separated factor, e.g.
    # "m-3" → "m⁻³". A bare factor with no exponent ("kg", "W") passes through.
    factors = map(split(u)) do f
        i = findfirst(c -> c == '-' || isdigit(c), f)
        i === nothing && return f
        # A factor that is all digits (rare) is not an exponent — leave it alone.
        i == 1 && return f
        String(f[1:i-1]) * map(c -> get(_SUPERSCRIPTS, c, c), f[i:end])
    end
    return join(factors, " ")
end

# CF `units` for a layer, from the stack's own metadata when it has any, else from
# the package table, else "". The fallback keeps the extension working on output
# stacks produced before layer metadata existed.
function _units(output, v::Symbol)
    md = DimensionalData.metadata(output[v])
    if !(md isa DimensionalData.Dimensions.Lookups.NoMetadata) && haskey(md, "units")
        return String(md["units"])
    end
    haskey(GEMB.GEMB_CF_ATTRIBUTES, v) && return GEMB.GEMB_CF_ATTRIBUTES[v].units
    return ""
end

# Display units for a layer and the conversion that produces them.
function _display_units(output, v::Symbol)
    u = _units(output, v)
    haskey(_DISPLAY_UNITS, u) && return _DISPLAY_UNITS[u]
    return (_pretty_units(u), identity)
end

# Short display name with its display units appended, e.g. "density [kg m⁻³]".
function _label(output, v::Symbol)
    du, _ = _display_units(output, v)
    return isempty(du) ? _short(v) : string(_short(v), " [", du, "]")
end

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
    # Age reads as depth-of-residence: young mass near the surface, old at the base.
    :age => :acton,
)
_cmap(v::Symbol) = get(_CMAP, v, :viridis)

# Display order for the profile (left) column — the counterpart of `_SCALAR_GROUPS` below,
# which orders the right column. Storage order would put `age` last (it was added last) and
# `dz` second.
#
# `water` sits below `grain_radius`: grain radius is set by the dry-snow metamorphism that
# follows from `temperature` and `density` above it, so those three read as one continuous
# microstructure block, and `water` then opens the melt-and-refreeze story that `age` closes.
# Variables absent from this list sort after the listed ones, keeping their stack order.
const _PROFILE_ORDER = [:temperature, :density, :grain_radius, :water, :age,
    :grain_dendricity, :grain_sphericity, :dz]
_prof_rank(v::Symbol) =
    something(findfirst(==(v), _PROFILE_ORDER), length(_PROFILE_ORDER) + 1)

# Scalar variables grouped by physical theme + shared unit so multiple series
# share one axis (increasing data density, removing near-duplicate panels).
# Each entry: (theme name, [member variables]). The unit suffix is *not* stored —
# it is derived from the members' own CF units at plot time by `_group_title`,
# so the panel can never disagree with the data it draws.
# Air temperature is placed first so it sits at the top of the right column,
# above (aligned with) the temperature profile heatmap in the left column.
const _SCALAR_GROUPS = [
    ("Air temperature", [:temperature_air]),
    ("Energy fluxes",
        [:shortwave_net, :longwave_net, :heat_flux_sensible, :heat_flux_latent]),
    ("Mass fluxes",
        [:melt, :runoff, :refreeze, :evaporation_condensation, :precipitation]),
    ("Broadband albedo", [:albedo_broadband]),
    ("Firn air content",
        [:firn_air_content, :firn_air_content_20m, :firn_air_content_10m]),
    # All depths/thicknesses below the surface, so they share an axis. `ice_slab_depth` and
    # `aquifer_depth` are NaN whenever nothing qualifies, which Makie draws as a gap — the
    # honest rendering of "there is nothing to plot here", and the reason neither layer is
    # interval-averaged. The aquifer series are flat zero / all-NaN under the
    # `runoff_method = :instantaneous` default, where no water can stand.
    ("Firn properties",
        [:percolation_depth, :ice_slab_depth, :ice_slab_thickness,
            :aquifer_depth, :aquifer_thickness]),
    ("Densification", [:densification_from_compaction, :densification_from_melt]),
    ("Strain thinning", [:strain_thinning]),
]

# Panel title for a scalar group: the theme name plus the members' shared display
# units. A group whose members disagree on units gets no suffix rather than a
# misleading one — the per-series legend then carries the distinction.
function _group_title(output, theme::AbstractString, members::Vector{Symbol})
    us = unique(first(_display_units(output, v)) for v in members)
    (length(us) == 1 && !isempty(us[1])) || return theme
    return string(theme, " [", us[1], "]")
end

# Height in points the figure allots each panel row. A rotated y-label is bounded by this,
# not by the figure's width, so every group title has to stay short enough to read on one
# line here — which is why the themes above are two or three words, not sentences.
const _ROW_HEIGHT = 175

#=============================================================================
# Isochrones: annual accumulation surfaces drawn over the age heatmap.
=============================================================================#

# The month/day whose annual recurrence is drawn as an isochrone: 1 May in the northern
# hemisphere, 1 November in the southern, both marking the end of the accumulation season, so
# each line is that year's snow surface. Chosen from the run's latitude metadata; defaults to
# the northern date when latitude is absent or unusable, since a wrong-by-half-a-year line is
# still a valid annual horizon.
function _isochrone_monthday(md)
    md === nothing && return (5, 1)
    lat = md isa AbstractDict ? get(md, "latitude", nothing) :
          (haskey(md, "latitude") ? md["latitude"] : nothing)
    return (lat isa Real && isfinite(lat) && lat < 0) ? (11, 1) : (5, 1)
end

# Legend text naming the date the isochrones mark. Without it the overlaid lines are unexplained
# — and which date they are is hemisphere-dependent, so it cannot be a fixed string.
function _isochrone_label(md)
    mo, dy = _isochrone_monthday(md)
    return string(dy, " ", monthabbr(mo))
end

# The isochrone levels, as decimal years, for every such date inside the record.
function _isochrone_dates(times, md)
    isempty(times) && return Float64[]
    mo, dy = _isochrone_monthday(md)
    t0, t1 = first(times), last(times)
    dates = [DateTime(y, mo, dy) for y in year(t0):year(t1)]
    filter!(t -> t0 <= t <= t1, dates)
    return GEMB.datetime2decyear(dates)
end

# Depth of one isochrone at each timestep: where the column's *deposition year* crosses
# `level`.
#
# The traced field is deposition year (`t - age`), not age. Age is a running clock, so a given
# age contour slides downward through the firn with time and traces no physical surface.
# Deposition year is fixed for a parcel of snow and travels with it, so its level set is
# exactly a buried annual surface.
#
# The crossing is solved on the model's own layers rather than on the regridded raster the
# heatmap draws. Regridding onto the 240-level display grid spreads each crossing over a grid
# cell, so a contour of it comes out as a band whose width is the grid spacing, not a horizon;
# and where several horizons fall inside the top metre — any column carrying spinup age — the
# bands merge into one blob. Solving here is also exact: deposition year is monotonic down the
# column, so the crossing is a single linear inverse per timestep.
#
# `NaN` marks "this horizon is not resolved here", which breaks the line — either the level
# postdates the record's surface at that time, it lies below the deepest modelled layer, or the
# bracketing interface is an unconformity (see `max_gap`).
#
# `max_gap`, in years, is what keeps the ablation zone honest. There, annual snow rests
# directly on ice tens of years older, so one layer interface spans that entire age jump and
# *every* level in the record brackets across it. Without the test each horizon is placed at
# the same snow/ice contact, and the lines stack into the solid band that gives the
# ablation zone away. But such an interface is a hiatus: the intervening years were never
# deposited, or have already melted away, so there is no surface to draw. Interpolating across
# it would invent a depth for snow that does not exist. Terminating is the honest rendering,
# and it leaves exactly the horizons the column really preserves.
function _isochrone_depth(z_col, deposition, level, max_gap)
    nz, nt = size(deposition)
    depth = fill(NaN, nt)
    for j in 1:nt
        for i in 1:(nz - 1)
            d1, d2 = deposition[i, j], deposition[i + 1, j]
            (isnan(d1) || isnan(d2)) && continue
            # Deposition year decreases with depth (deeper is older), so the crossing is
            # bracketed when the level sits between a layer and the one beneath it. `>=` on
            # the upper side catches a level landing exactly on a layer centre.
            if d1 >= level >= d2
                # Bracketed, but across a gap too large to be a real annual surface: leave
                # this timestep blank rather than pinning the horizon to the unconformity.
                d1 - d2 <= max_gap || break
                z1, z2 = z_col[i, j], z_col[i + 1, j]
                w = d1 == d2 ? 0.0 : (d1 - level) / (d1 - d2)
                depth[j] = z1 + w * (z2 - z1)
                break
            end
        end
    end
    return depth
end

# Isochrone stroke style: a hairline dark core over a wider translucent white casing.
# The age colormap runs from near-black at the surface to near-white at depth, so no
# single flat colour reads across it — the yellow this used to be disappeared over the
# pale, oldest firn. A cased line reads on both ends, for the same reason a road on a
# map is drawn with a halo. Widths stay small so a dozen horizons annotate the field
# rather than dominate it.
const _ISO_CORE = RGBAf(0.10, 0.10, 0.12, 0.90)
const _ISO_CASE = RGBAf(1.00, 1.00, 1.00, 0.55)
const _ISO_LW = 0.7
const _ISO_CASE_LW = 2.2

# Draw one thin line per annual horizon over an already-plotted heatmap.
#
# A layer interface may legitimately span a couple of years of deposition age once the model
# has merged thin layers deep in the column, so the gap tolerance is a few times the annual
# spacing rather than exactly one year — loose enough not to chop up well-resolved firn,
# tight enough to reject the decades-wide jumps of an ablation surface.
function _isochrones!(ax, x, z_col, deposition, levels, label; max_gap=3.0)
    isempty(levels) && return nothing
    for level in levels
        depth = _isochrone_depth(z_col, deposition, level, max_gap)
        # Casing first, core second, so every core sits above every casing — otherwise a
        # closely-spaced pair of horizons has the upper one's casing wash out the lower's line.
        lines!(ax, x, depth; color=_ISO_CASE, linewidth=_ISO_CASE_LW)
        lines!(ax, x, depth; color=_ISO_CORE, linewidth=_ISO_LW)
    end
    # The legend entry is built explicitly rather than harvested from a `label=` on the first
    # line: at legend scale a hairline is invisible (the reason the old yellow swatch could not
    # be read), so the swatch is drawn as one thicker cased line — a schematic of the horizon,
    # not a literal copy of its width.
    axislegend(ax,
        [[LineElement(color=_ISO_CASE, linewidth=4.0),
          LineElement(color=_ISO_CORE, linewidth=1.6)]],
        [label];
        position=:rb, framevisible=false, labelsize=11, labelcolor=:gray20,
        patchsize=(20, 10), padding=(5, 5, 3, 3), backgroundcolor=(:white, 0.85))
    return nothing
end

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
# "ERA5-Land @ 72.58°N 38.46°W 3216 m, +100 m" — the elevation is the source dataset's own
# surface and the trailing Δ is the offset applied to reach the target. Every piece is optional: a stack with
# no provenance metadata yields "" and the banner renders as it did before.
function _provenance_str(md)
    md === nothing && return ""
    getmd(k) = md isa AbstractDict ? get(md, k, nothing) :
               (haskey(md, k) ? md[k] : nothing)

    parts = String[]
    ds = getmd("dataset")
    ds isa AbstractString && !isempty(ds) && push!(parts, ds)

    off = getmd("elevation_offset")
    lat, lon = getmd("latitude"), getmd("longitude")
    if lat isa Real && lon isa Real && isfinite(lat) && isfinite(lon)
        lats = string(round(abs(lat), digits=2), "°", lat >= 0 ? "N" : "S")
        lons = string(round(abs(lon), digits=2), "°", lon >= 0 ? "E" : "W")
        coord = "$lats $lons"
        # Surface elevation the source dataset itself represents — what the target was
        # lapse-rate adjusted *from*. Carried explicitly as `elevation_native` rather than
        # derived from target-minus-offset, so a missing target can't silently poison it.
        # Sits with the dataset/coordinates; the Δ that follows gives the target.
        # Older outputs predate the key: fall back to `elevation` (equal to the native
        # elevation when no adjustment was applied) and omit the field entirely otherwise.
        zc = getmd("elevation_native")
        zc isa Real && isfinite(zc) || (zc = getmd("elevation"))
        zc isa Real && isfinite(zc) && (coord *= " " * string(round(Int, zc), " m"))
        isempty(parts) ? push!(parts, coord) : (parts[end] *= " @ " * coord)
    end

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
    # The terms are per-interval mass accumulations divided by the record length, so the
    # unit is a mass layer's own unit per year. Any of them will do; `parts` being
    # non-empty guarantees at least one is present.
    mass_var = findfirst(in(present),
        [:precipitation, :rain, :runoff, :refreeze, :evaporation_condensation])
    base = mass_var === nothing ? "" :
           _units(output, [:precipitation, :rain, :runoff, :refreeze,
                           :evaporation_condensation][mass_var])
    suffix = isempty(base) ? "" : "  [" * _pretty_units(base * " yr-1") * "]"
    return join(parts, "  ") * suffix
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
    # dendricity/sphericity are redundant with grain radius, which already represents
    # the grain-metamorphism family.
    if variables === nothing
        default_drop = (:dz, :grain_dendricity, :grain_sphericity)
        profile_vars = filter(v -> v ∉ default_drop, profile_vars)
    end
    # Physical order rather than storage order (see `_PROFILE_ORDER`). Stable, so any
    # variable the list doesn't name keeps its position relative to the other unnamed ones.
    sort!(profile_vars; by=_prof_rank, alg=MergeSort)
    scalar_groups = [(t, [v for v in g if v in wanted]) for (t, g) in _SCALAR_GROUPS]
    scalar_groups = [(t, g) for (t, g) in scalar_groups if !isempty(g)]
    # Densification and strain thinning are dropped from the default panel set (kept
    # available when the caller lists their variables explicitly via `variables`). Strain
    # thinning is identically zero unless `horizontal_strain_rate` is set, so it would
    # otherwise render an empty panel for every run.
    if variables === nothing
        default_drop_groups = ("Densification", "Strain thinning")
        scalar_groups = [(t, g) for (t, g) in scalar_groups if t ∉ default_drop_groups]
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
    fig = Figure(size=(1500, 200 + _ROW_HEIGHT * nrows), fontsize=14)

    # Header: title + traceability metadata banner.
    _header!(fig[1, 1:2], output, times, z_center, title)

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
        # Convert to display units (e.g. K → °C, so the profile shares the
        # air-temperature panel's units). Driven by the layer's CF `units`, so it
        # applies to any temperature field rather than to `:temperature` by name.
        _, convert_units = _display_units(output, v)
        convert_units === identity || (vals = convert_units.(vals))
        gridded = gemb_interp(z_center[:, tcols], vals, z_target)
        ax = Axis(left[i, 1]; ylabel="depth [m]")
        push!(prof_axes, ax)
        lo, hi = _robust_limits(parent(gridded))
        # Anchor water's ramp at zero so "no liquid water" reads as white.
        v === :water && (lo = 0.0)
        hm = heatmap!(ax, x_hm, collect(dims(gridded, Z)), parent(gridded');
            colormap=_cmap(v), colorrange=(lo, hi))
        Colorbar(left[i, 2], hm; label=_label(output, v), width=12)
        # Age is the one profile whose heatmap has annual structure worth marking: overlay the
        # buried accumulation surfaces. Traced on the native layer centres, not on `gridded` —
        # see `_isochrone_depth` for why regridding first turns each horizon into a band.
        if v === :age
            md = DimensionalData.metadata(output)
            age_yr = parent(output[:age])[:, tcols] ./ 365.25
            _isochrones!(ax, x_hm, z_center[:, tcols], reshape(x_hm, 1, :) .- age_yr,
                _isochrone_dates(times, md), _isochrone_label(md))
        end
        ylims!(ax, zbot, ztop)
    end

    # ---- Right column: grouped scalar time series -------------------------
    for (i, (theme, members)) in enumerate(scalar_groups)
        ax = Axis(right[i, 1]; ylabel=_group_title(output, theme, members))
        push!(scal_axes, ax)
        # Air temperature: shown in its display units (°C) with the curve filled to
        # zero — warm (> 0 °C) in red, cold (≤ 0 °C) in blue — so melt-permissive
        # periods read at a glance. The zero crossing is meaningful only in °C, which
        # is why this panel is special-cased rather than drawn as a plain series.
        if members == [:temperature_air]
            _, convert_units = _display_units(output, :temperature_air)
            tc = convert_units.(Float64.(parent(output[:temperature_air])))
            band!(ax, decyear, min.(tc, 0.0), 0.0;
                color=(:blue, 0.5))
            band!(ax, decyear, 0.0, max.(tc, 0.0);
                color=(:red, 0.5))
            continue
        end
        # With ~10⁵ sub-daily points, overplotted opaque lines saturate to a
        # solid block; thin, semi-transparent strokes let density show through.
        alpha = length(members) > 1 ? 0.75 : 0.9
        # The axis title carries the shared units, so legend entries are the short
        # names alone. If the members' units disagree the title drops its suffix, and
        # each legend entry has to carry its own units instead.
        shared_units = _group_title(output, theme, members) != theme
        for v in members
            _, convert_units = _display_units(output, v)
            vals = Float64.(parent(output[v]))
            convert_units === identity || (vals = convert_units.(vals))
            lines!(ax, decyear, vals;
                label=shared_units ? _short(v) : _label(output, v),
                linewidth=0.6, alpha=alpha)
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
function _header!(pos, output, times, z_center, title)
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
        "  │  generated ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM"),
    )
    head = title == "" ? "GEMB.jl glacier firn model output" : title
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
