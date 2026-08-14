#=============================================================================
# CF-style attributes for the `gemb` output layers.
#
# Concise by design: `units` (UDUNITS-parsable), `long_name`, `cell_methods`
# (how the value was reduced over the output interval), and `standard_name`
# ONLY where a CF standard name matches GEMB's quantity exactly. Most GEMB
# variables have none. Several near-misses exist — `temperature_in_surface_snow`,
# `surface_snow_density`, `liquid_water_content_of_surface_snow` are all real CF
# names — but each denotes a bulk whole-pack value, while GEMB reports per-cell
# values down the column. Omitting `standard_name` is CF-legal; misapplying one
# is not, so those are deliberately left off.
#
# Names here were checked against the CF Standard Name Table XML v94
# (2026-06-09) and cross-checked on the NERC vocabulary server (collection P07).
#
# The table is a load-time constant, read only while `gemb` allocates its output
# stack — never inside the time loop.
=============================================================================#

"""
    CFAttrs(units, long_name, cell_methods; standard_name="", comment="")

CF attributes for one output variable. `cell_methods`, `standard_name`, and
`comment` are empty strings when they do not apply, and [`cf_attributes`](@ref)
omits empty fields from the emitted attribute dictionary rather than writing
blanks.

`cell_methods` is held separately from the rest so it can be dropped for products
with no time coordinate — a `cell_methods` of `"time: point"` must not appear on
a variable whose `time` coordinate does not exist.
"""
struct CFAttrs
    units::String
    long_name::String
    cell_methods::String   # "" when there is no time reduction
    standard_name::String  # "" when no CF standard name applies exactly
    comment::String        # "" when nothing needs saying
end

CFAttrs(units, long_name, cell_methods; standard_name="", comment="") =
    CFAttrs(units, long_name, cell_methods, standard_name, comment)

"""
    GEMB_CF_ATTRIBUTES::Dict{Symbol,CFAttrs}

CF attributes for every layer of the [`gemb`](@ref) output stack, keyed by layer
name. `cell_methods` records how each layer was reduced over the output interval:
`"time: sum"` for mass/thickness accumulations, `"time: mean"` for interval
averages, `"time: point"` for instantaneous snapshots taken at the output time.
"""
const GEMB_CF_ATTRIBUTES = Dict{Symbol,CFAttrs}(
    # ---- monolevel, interval sums -------------------------------------------
    :melt => CFAttrs("kg m-2", "meltwater produced", "time: sum";
        standard_name="surface_snow_melt_amount",
        comment="Column total, including bare-ice melt in the ablation regime."),
    :runoff => CFAttrs("kg m-2", "meltwater runoff leaving the column",
        "time: sum"; standard_name="runoff_amount"),
    :refreeze => CFAttrs("kg m-2", "meltwater refrozen within the column",
        "time: sum";
        comment="CF has only surface_snow_and_ice_refreezing_flux (kg m-2 s-1); " *
                "there is no amount form, so standard_name is omitted."),
    :evaporation_condensation => CFAttrs("kg m-2",
        "net surface water vapour mass exchange", "time: sum";
        comment="Positive = mass gain by the column (condensation/deposition); " *
                "negative = evaporation/sublimation. This is the opposite sign " *
                "to CF's surface_upward_water_vapor_flux_in_air, which is why " *
                "no standard_name is given."),
    :precipitation => CFAttrs("kg m-2", "total precipitation", "time: sum";
        standard_name="precipitation_amount"),
    :rain => CFAttrs("kg m-2", "rainfall", "time: sum";
        standard_name="rainfall_amount"),
    :densification_from_compaction => CFAttrs("m",
        "column thickness change from dry compaction", "time: sum"),
    :densification_from_melt => CFAttrs("m",
        "column thickness change from melt and wet compaction", "time: sum"),
    :strain_thinning => CFAttrs("m",
        "column thickness change from horizontal ice-dynamic strain", "time: sum";
        comment="Positive = thinning under horizontal divergence. Zero unless " *
                "horizontal_strain_rate is set. No CF standard name exists for " *
                "this quantity."),

    # ---- monolevel, interval means ------------------------------------------
    :shortwave_net => CFAttrs("W m-2", "net shortwave radiation at the surface",
        "time: mean"; standard_name="surface_net_downward_shortwave_flux",
        comment="Positive downward (toward the surface)."),
    :longwave_net => CFAttrs("W m-2", "net longwave radiation at the surface",
        "time: mean"; standard_name="surface_net_downward_longwave_flux",
        comment="Positive downward (toward the surface)."),
    :heat_flux_sensible => CFAttrs("W m-2", "sensible heat flux at the surface",
        "time: mean"; standard_name="surface_downward_sensible_heat_flux",
        comment="Positive downward (toward the surface)."),
    :heat_flux_latent => CFAttrs("W m-2", "latent heat flux at the surface",
        "time: mean"; standard_name="surface_downward_latent_heat_flux",
        comment="Positive downward (toward the surface)."),
    :albedo_broadband => CFAttrs("1", "broadband surface albedo", "time: mean";
        standard_name="surface_albedo"),
    :temperature_air => CFAttrs("K", "near-surface air temperature",
        "time: mean"; standard_name="air_temperature"),
    :firn_air_content => CFAttrs("m", "firn air content", "time: mean";
        comment="No CF standard name exists for firn air content."),
    :thickness_cumulative => CFAttrs("m", "cumulative column thickness change",
        "time: mean";
        comment="Interval mean of a running total since the start of the run."),

    # ---- monolevel, instantaneous -------------------------------------------
    :valid_profile_length => CFAttrs("1", "number of valid column layers",
        "time: point"),

    # ---- profile (Z x Ti), instantaneous snapshots --------------------------
    :temperature => CFAttrs("K", "layer temperature", "time: point";
        comment="Per-cell value; CF's temperature_in_surface_snow and " *
                "land_ice_temperature denote bulk quantities, so neither applies."),
    :dz => CFAttrs("m", "layer thickness", "time: point"),
    :density => CFAttrs("kg m-3", "layer bulk density", "time: point";
        comment="Per-cell value; CF's surface_snow_density is a whole-pack quantity."),
    :water => CFAttrs("kg m-2", "layer liquid water content", "time: point";
        comment="Per-cell value; CF's liquid_water_content_of_surface_snow is " *
                "a whole-pack quantity."),
    :grain_radius => CFAttrs("mm", "snow grain effective radius", "time: point"),
    :grain_dendricity => CFAttrs("1", "snow grain dendricity", "time: point"),
    :grain_sphericity => CFAttrs("1", "snow grain sphericity", "time: point"),
)

"""
    cf_attributes(key::Symbol; time_axis=true)

CF attributes for output layer `key` as a `Dict{String,String}` ready to attach
as `DimArray` metadata (or to write as NetCDF variable attributes). Attributes
that do not apply to the variable are absent rather than empty.

Pass `time_axis=false` for a product with no time coordinate (e.g. a single
column from [`gemb_profile`](@ref)) to drop `cell_methods`, which would otherwise
reference a `time` coordinate that is not there.

A fresh `Dict` is built on each call, so [`GEMB_CF_ATTRIBUTES`](@ref) cannot be
mutated through the returned value.
"""
function cf_attributes(key::Symbol; time_axis::Bool=true)
    a = GEMB_CF_ATTRIBUTES[key]
    d = Dict{String,String}("units" => a.units, "long_name" => a.long_name)
    time_axis && !isempty(a.cell_methods) && (d["cell_methods"] = a.cell_methods)
    isempty(a.standard_name) || (d["standard_name"] = a.standard_name)
    isempty(a.comment) || (d["comment"] = a.comment)
    return d
end

"""
    cf_layermetadata(layers::NamedTuple; time_axis=true)

Build the `layermetadata` NamedTuple for a `DimStack` whose layers are `layers`.

Mapping over `keys(layers)` rather than a fixed order is deliberate:
DimensionalData requires the `layermetadata` keys to match the layer keys both in
membership and in order, and the extracted-profile stack orders its layers
differently from the `gemb` output stack.
"""
cf_layermetadata(layers::NamedTuple; time_axis::Bool=true) =
    NamedTuple{keys(layers)}(map(k -> cf_attributes(k; time_axis), keys(layers)))

"""
    cf_time_attributes()

CF attributes for the `Ti` dimension of GEMB output.
"""
cf_time_attributes() = Dict{String,String}(
    "standard_name" => "time", "long_name" => "time", "axis" => "T")

"""
    cf_layer_index_attributes()

CF attributes for the `Z` dimension of GEMB output, which is a bare 1-based cell
index rather than a coordinate. It therefore carries no `standard_name` and no
`positive`: claiming either would present an index as a physical height. Depths
come from the `dz` layer via [`dz2z`](@ref).
"""
cf_layer_index_attributes() = Dict{String,String}(
    "units" => "1",
    "long_name" => "vertical layer index (1 = surface layer)",
    "axis" => "Z",
    "comment" => "Layer index, surface first. Layer thickness is the `dz` " *
                 "variable; cell-centre heights via dz2z(dz).")

"""
    cf_height_attributes()

CF attributes for a genuine vertical coordinate, as produced by
[`gemb_interp`](@ref) when it regrids onto fixed heights. [`dz2z`](@ref) returns
heights measured negative-downward from the surface, hence `positive = "up"`.
"""
cf_height_attributes() = Dict{String,String}(
    "units" => "m",
    "long_name" => "height relative to the snow/firn surface",
    "axis" => "Z",
    "positive" => "up")

"""
    GEMB_CF_GLOBAL_ATTRIBUTES::Dict{String,Any}

CF global attributes stamped onto the [`gemb`](@ref) output stack alongside the
run's own provenance. A `history` attribute is deliberately omitted: a generation
timestamp would make the output non-deterministic and break the numerical
regression fingerprint.
"""
const GEMB_CF_GLOBAL_ATTRIBUTES = Dict{String,Any}(
    "Conventions" => "CF-1.11",
    "title" => "GEMB (Glacier Energy and Mass Balance) point simulation output",
    "source" => "GEMB.jl",
)
