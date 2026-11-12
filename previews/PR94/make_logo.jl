# Generates docs/src/assets/logo.svg — the GEMB.jl mark.
#
# The mark is a firn column drawn as stacked rows of grains. Grain size grows with depth,
# which is the way snow actually evolves — small fresh crystals at the surface (Julia purple)
# coarsen by sintering and vapour transport into firn (green) and then into the large crystals
# of close-off ice (red).
#
# Rows rather than a jammed packing: at logo size the row structure reads as stratigraphy —
# the layering a firn core actually shows — and stays legible at favicon scale, where an
# irregular packing turns to noise. A thin rule on each row boundary makes those layers
# explicit.
#
# The background is transparent, so the mark sits on whatever it is placed on (VitePress
# light and dark themes, README, GitHub avatar) without carrying a pale square with it.
#
# Row radii are fitted to the available height, so the top and bottom rows always land flush
# against the margins for any ROWS / R_TOP / R_BOT.
#
# Run:  julia --project=docs docs/src/assets/make_logo.jl

const W = 128.0
const H = 128.0

const PURPLE = "#9558b2"   # Julia purple — surface snow
const GREEN = "#389826"    # Julia green  — firn
const RED = "#cb3c33"      # Julia red    — close-off ice

const MARGIN = 4.0         # keeps whole grains clear of the rounded column edge
const ROWS = 7             # number of grain rows, surface to base
const R_TOP = 1.0          # top-row grain radius, as a *relative* weight...
const R_BOT = 2.2          # ...and bottom-row weight; absolute sizes are fitted to H

const LINE = "#3d4451"     # row-boundary rule: near-black, slightly cool
const LINE_W = 0.6
const LINE_OPACITY = 0.4
const LINE_INSET = 3.0     # rules stop short of the edge, so they read as layers not ruling

"""
    rows() -> Vector{NamedTuple}

Return the grain rows as `(; y, r)` in SVG coordinates, surface first.

Radii grow linearly with row index from `R_TOP` to `R_BOT`, then the whole set is scaled so
the stack of diameters exactly fills the height between the margins. That fitting is why the
top and bottom rows are always flush: the constants set the *ratio* of surface to basal grain
size, and the geometry follows.
"""
function rows()
    weights = [R_TOP + (R_BOT - R_TOP) * (i - 1) / (ROWS - 1) for i in 1:ROWS]
    avail = H - 2 * MARGIN
    scale = avail / (2 * sum(weights))
    radii = weights .* scale

    out = NamedTuple{(:y, :r),Tuple{Float64,Float64}}[]
    y = MARGIN
    for r in radii
        push!(out, (; y = y + r, r))
        y += 2 * r
    end
    return out
end

# Snow, firn and ice by row. Splitting on row index rather than depth keeps each regime a
# whole number of rows, so no row is bicoloured. The bands are uneven in row count but even
# in *height*: the deep rows are so much thicker that two of them span as much of the mark as
# the three shallow ones.
row_color(i) = i <= 3 ? PURPLE : (i <= 5 ? GREEN : RED)

r2(v) = round(v, digits = 2)

function write_logo(path)
    rs = rows()
    n_grains = 0
    open(path, "w") do io
        print(io, """
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128" width="128" height="128"
             role="img" aria-labelledby="gemb-logo-title gemb-logo-desc">
          <title id="gemb-logo-title">GEMB.jl</title>
          <desc id="gemb-logo-desc">A firn column in the Julia colors: rows of grains that coarsen
          with depth, from small fresh snow crystals at the surface through firn to the large crystals
          of close-off ice at the base, with each layer boundary ruled.</desc>
          <defs><clipPath id="gemb-col"><rect x="0" y="0" width="128" height="128" rx="12"/></clipPath></defs>
          <g clip-path="url(#gemb-col)">
        """)

        for (i, row) in enumerate(rs)
            fill = row_color(i)
            # Lay the row out with whole grains only, centred, so none is clipped by the
            # column edge: take the count that fits inside the margins, then centre them.
            # Centring every row (rather than staggering alternate ones) keeps the column
            # symmetric — with grains this large, offsetting a row by a fraction of its
            # pitch visibly throws the whole mark off axis.
            pitch = 2 * row.r
            n = max(1, floor(Int, (W - 2 * MARGIN) / pitch))
            span = (n - 1) * pitch
            x0 = (W - span) / 2
            for k in 0:(n - 1)
                print(io, """    <circle cx="$(r2(x0 + k * pitch))" cy="$(r2(row.y))" """,
                      """r="$(r2(row.r))" fill="$fill"/>\n""")
                n_grains += 1
            end
        end

        # Row boundaries, drawn last so they read as stratigraphy over the grains. Interior
        # boundaries only — the outer two would just trace the column edge.
        print(io, """    <g stroke="$LINE" stroke-width="$LINE_W" stroke-opacity="$LINE_OPACITY">\n""")
        for i in 1:(length(rs) - 1)
            y = rs[i].y + rs[i].r        # bottom of row i == top of row i+1
            print(io, """      <line x1="$LINE_INSET" y1="$(r2(y))" """,
                  """x2="$(r2(W - LINE_INSET))" y2="$(r2(y))"/>\n""")
        end
        print(io, "    </g>\n")

        print(io, "  </g>\n</svg>\n")
    end
    return n_grains
end

n = write_logo(joinpath(@__DIR__, "logo.svg"))
println("wrote logo.svg with $n grains in $ROWS rows")
