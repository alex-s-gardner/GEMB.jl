# Compilation entry point for the GEMB C API shared library.
#
# Built by `capi/build.jl`, which drives JuliaC. Nothing here is loaded by a normal
# `using GEMB` — the `Base.@ccallable` entry points, the handle registry and the trim
# constraints they impose live in this project, not in `src/`.
#
# Scope of the compiled surface: one column, one timestep. `GEMB.gemb_core` and the leaf
# physics it calls. See `capi/README.md`.

module GEMBCApi

using GEMB

# ---------------------------------------------------------------------------------------------
# Status codes. Mirrored in `capi/gemb.h`.
# ---------------------------------------------------------------------------------------------

const GEMB_OK = Cint(0)
const GEMB_ERR_BAD_ARGUMENT = Cint(-1)
const GEMB_ERR_BAD_HANDLE = Cint(-2)
const GEMB_ERR_BAD_OPTION = Cint(-3)
const GEMB_ERR_PHYSICS = Cint(-4)
const GEMB_ERR_UNEXPECTED = Cint(-5)

# ---------------------------------------------------------------------------------------------
# Handle registries.
#
# `ModelParameters` carries `Symbol` fields and a `Vector{Float64}`, and `ThermalWorkspace`
# owns growable buffers; neither has a C layout. Both stay on the Julia side and are named
# across the boundary by an opaque `Cint` handle.
#
# These tables are mutable module-level state, which the model in `src/` deliberately has none
# of. Confining it here is the point: `src/` stays free of it, so a `ModelParameters` remains a
# pure description that any number of threads can read.
# ---------------------------------------------------------------------------------------------

const PARAMS = Dict{Cint,GEMB.ModelParameters{GEMB.ExplicitThermal}}()
const WORKSPACES = Dict{Cint,GEMB.ThermalWorkspace}()
const NEXT_HANDLE = Ref(Cint(1))
const LAST_ERROR = Ref{String}("")

function _take_handle()
    h = NEXT_HANDLE[]
    NEXT_HANDLE[] = h + Cint(1)
    return h
end

"""
    _record_error(context, e)

Record a message for `gemb_last_error` describing a caught exception.

Only the exception's type name is reported, not its message: rendering an arbitrary exception
needs `showerror` on an `Any`, which the `--trim` verifier cannot resolve. `context` is what
makes the result useful — it names the operation that failed.
"""
function _record_error(context::String, @nospecialize(e))
    LAST_ERROR[] = context * ": " * string(nameof(typeof(e)))
    return nothing
end

# ---------------------------------------------------------------------------------------------
# Option code tables.
#
# One tuple per `Symbol`-valued `ModelParameters` field a single column step reads, in the
# order the `options` array declares them. A code indexes into the field's tuple, so the
# integer a host passes is checked here rather than deep in the physics.
#
# `output_frequency` and `initialize_age` are absent: they select output cadence and profile
# initialization, neither of which a single step consults.
# ---------------------------------------------------------------------------------------------

const OPTION_FIELDS = (
    :densification_method,
    :densification_coeffs_M01,
    :densification_accumulation,
    :mean_temperature_method,
    :new_snow_method,
    :emissivity_method,
    :thermal_conductivity_method,
    :heat_capacity_method,
    :rain_heat_capacity,
    :grain_growth_method,
    :water_irreducible_method,
    :melt_geometry,
    :runoff_method,
    :albedo_method,
    :blowing_snow_method,
)

const OPTION_VALUES = (
    (:HerronLangway, :Arthern, :ArthernB, :Barnola1991, :Crocus, :CrocusPure, :GSFC2020,
     :Simonsen2013, :Ligtenberg),
    (:Ant_ERA5_GS_SW0, :Ant_ERA5v4_Paolo23, :Ant_ERA5_BF_SW1, :Ant_RACMO_GS_SW0,
     :Ant_Ligtenberg, :Gre_ERA5_GS_SW0, :Gre_RACMO_GS_SW0, :Gre_RACMO_GB_SW1,
     :Gre_KuipersMunneke),
    (:accumulation, :precipitation),
    (:arithmetic, :arrhenius),
    (:Constant150, :Constant315, :Constant350, :Fausto, :FaustoFit, :Pahaut, :Kaspers,
     :KuipersMunneke),
    (:uniform, :grain_radius_threshold, :grain_radius_w_threshold),
    (:Sturm, :Calonne, :Calonne2019, :Calonne2019Air, :Marchenko2019),
    (:constant, :CuffeyPaterson),
    (:water, :ice),
    (:Marbouty, :Arthern, :hybrid),
    (:constant, :ColeouLesaffre),
    (:thickness, :density),
    (:instantaneous, :ZuoOerlemans, :Darcy),
    (:None, :GardnerSharp, :BrunLefebre, :GreuellKonzelmann),
    (:none, :Crocus),
)

const N_OPTIONS = length(OPTION_FIELDS)

# Numeric `ModelParameters` fields a single column step reads, in the order the `numeric`
# array declares them.
const NUMERIC_FIELDS = (
    :density_ice,
    :rain_temperature_threshold,
    :emissivity,
    :emissivity_grain_radius_large,
    :emissivity_grain_radius_threshold,
    :surface_roughness_effective_ratio,
    :heat_capacity_ice,
    :water_irreducible_saturation,
    :impermeable_density,
    :impermeable_thickness,
    :pore_saturation_max,
    :albedo_density_threshold,
    :albedo_snow,
    :albedo_ice,
    :albedo_fixed,
    :column_ztop,
    :column_dztop,
    :column_dzmin,
    :column_dzmax,
    :column_depth_max,
    :column_zy,
    :horizontal_strain_rate,
    :surface_slope,
    :drift_rate,
    :thermal_explicit_safety_factor,
)

const N_NUMERIC = length(NUMERIC_FIELDS)

# Boolean fields, in the order the `flags` array declares them.
const FLAG_FIELDS = (
    :shortwave_subsurface_absorption,
    :output_viscosity,
    :blowing_snow_sublimation,
)

const N_FLAGS = length(FLAG_FIELDS)

const N_FORCING = fieldcount(GEMB.ClimateForcingStep)
const N_FLUX = 21  # the flux scalars; `viscosity` comes back through its own buffer

"""
    _resolve_option(i, code) -> Symbol

The `Symbol` that `code` selects for the `i`th entry of `OPTION_FIELDS`.

Throws when the code is out of range for that field, so an invalid option is rejected at the
boundary instead of reaching a physics branch that would silently fall through.
"""
function _resolve_option(i::Int, code::Cint)
    values = OPTION_VALUES[i]
    1 <= code <= length(values) ||
        throw(ArgumentError("option $(OPTION_FIELDS[i]): code $code is out of range 1:$(length(values))"))
    return values[code]
end

"""
    _build_parameters(dt, options, numeric, flags) -> ModelParameters{ExplicitThermal}

Assemble a `ModelParameters` from the three C-side spec arrays and derive `dt_divisors` for a
forcing step of `dt` seconds.

The derivation is not optional: `dt_divisors` is empty in a default `ModelParameters` (it is
listed in `GEMB.DERIVED_PARAMETERS`) and the thermal sub-step search indexes into it, so a
column stepped with the empty vector fails on a bounds error. Doing it here is why parameters
are created through a call rather than filled in as a struct by the host.

Every field is named explicitly. Reconstructing by iterating `fieldnames` yields the `UnionAll`
`ModelParameters` rather than a concrete `ModelParameters{ExplicitThermal}`, which leaves the
whole downstream call tree unresolvable to the `--trim` verifier.
"""
function _build_parameters(dt::Float64, options::Vector{Cint}, numeric::Vector{Float64},
    flags::Vector{Cint})

    o(i) = _resolve_option(i, options[i])
    f(i) = flags[i] != 0

    return GEMB.ModelParameters(;
        densification_method=o(1),
        densification_coeffs_M01=o(2),
        densification_accumulation=o(3),
        mean_temperature_method=o(4),
        new_snow_method=o(5),
        emissivity_method=o(6),
        thermal_conductivity_method=o(7),
        heat_capacity_method=o(8),
        rain_heat_capacity=o(9),
        grain_growth_method=o(10),
        water_irreducible_method=o(11),
        melt_geometry=o(12),
        runoff_method=o(13),
        albedo_method=o(14),
        blowing_snow_method=o(15),
        density_ice=numeric[1],
        rain_temperature_threshold=numeric[2],
        emissivity=numeric[3],
        emissivity_grain_radius_large=numeric[4],
        emissivity_grain_radius_threshold=numeric[5],
        surface_roughness_effective_ratio=numeric[6],
        heat_capacity_ice=numeric[7],
        water_irreducible_saturation=numeric[8],
        impermeable_density=numeric[9],
        impermeable_thickness=numeric[10],
        pore_saturation_max=numeric[11],
        albedo_density_threshold=numeric[12],
        albedo_snow=numeric[13],
        albedo_ice=numeric[14],
        albedo_fixed=numeric[15],
        column_ztop=numeric[16],
        column_dztop=numeric[17],
        column_dzmin=numeric[18],
        column_dzmax=numeric[19],
        column_depth_max=numeric[20],
        column_zy=numeric[21],
        horizontal_strain_rate=numeric[22],
        surface_slope=numeric[23],
        drift_rate=numeric[24],
        thermal_explicit_safety_factor=numeric[25],
        shortwave_subsurface_absorption=f(1),
        output_viscosity=f(2),
        blowing_snow_sublimation=f(3),
        dt_divisors=GEMB.fast_divisors(round(Int, dt * 10000)) ./ 10000,
    )
end

# ---------------------------------------------------------------------------------------------
# Entry points.
#
# Every one returns a status code: a Julia exception must not unwind into C. The message the
# code cannot carry is held in `LAST_ERROR` for `gemb_last_error`.
# ---------------------------------------------------------------------------------------------

"""
    gemb_defaults(options, numeric, flags) -> Cint

Fill the three spec arrays with the package defaults, so a host can take them and override only
the entries it cares about.

`options` holds `GEMB_N_OPTIONS` ints, `numeric` `GEMB_N_NUMERIC` doubles, `flags`
`GEMB_N_FLAGS` ints.
"""
Base.@ccallable function gemb_defaults(options::Ptr{Cint}, numeric::Ptr{Cdouble},
    flags::Ptr{Cint})::Cint
    (options == C_NULL || numeric == C_NULL || flags == C_NULL) && return GEMB_ERR_BAD_ARGUMENT
    try
        defaults = GEMB.ModelParameters()

        for i in 1:N_OPTIONS
            # An explicit loop, not `findfirst`: `OPTION_VALUES` is a heterogeneous tuple, and
            # searching it with a closure leaves the call unresolvable to the `--trim` verifier.
            code = 0
            values = OPTION_VALUES[i]
            value = getfield(defaults, OPTION_FIELDS[i])
            for j in eachindex(values)
                values[j] === value && (code = j; break)
            end
            code == 0 && return GEMB_ERR_BAD_OPTION
            unsafe_store!(options, Cint(code), i)
        end
        for i in 1:N_NUMERIC
            unsafe_store!(numeric, getfield(defaults, NUMERIC_FIELDS[i]), i)
        end
        for i in 1:N_FLAGS
            unsafe_store!(flags, Cint(getfield(defaults, FLAG_FIELDS[i])), i)
        end
        return GEMB_OK
    catch e
        _record_error("gemb_defaults", e)
        return GEMB_ERR_UNEXPECTED
    end
end

"""
    gemb_params_create(dt, options, numeric, flags, handle_out) -> Cint

Build a parameter set for a forcing step of `dt` seconds and write its handle to `handle_out`.

The handle stays valid until `gemb_params_destroy`. The parameters are immutable, so one handle
may be shared by any number of threads.
"""
Base.@ccallable function gemb_params_create(dt::Cdouble, options::Ptr{Cint},
    numeric::Ptr{Cdouble}, flags::Ptr{Cint}, handle_out::Ptr{Cint})::Cint
    (options == C_NULL || numeric == C_NULL || flags == C_NULL || handle_out == C_NULL) &&
        return GEMB_ERR_BAD_ARGUMENT
    dt > 0.0 || return GEMB_ERR_BAD_ARGUMENT
    try
        opts = [unsafe_load(options, i) for i in 1:N_OPTIONS]
        nums = [unsafe_load(numeric, i) for i in 1:N_NUMERIC]
        flgs = [unsafe_load(flags, i) for i in 1:N_FLAGS]

        mp = _build_parameters(dt, opts, nums, flgs)
        h = _take_handle()
        PARAMS[h] = mp
        unsafe_store!(handle_out, h)
        return GEMB_OK
    catch e
        _record_error("gemb_params_create", e)
        return e isa ArgumentError ? GEMB_ERR_BAD_OPTION : GEMB_ERR_UNEXPECTED
    end
end

"""
    gemb_params_destroy(handle) -> Cint

Release a parameter set.
"""
Base.@ccallable function gemb_params_destroy(handle::Cint)::Cint
    haskey(PARAMS, handle) || return GEMB_ERR_BAD_HANDLE
    delete!(PARAMS, handle)
    return GEMB_OK
end

"""
    gemb_workspace_create(handle_out) -> Cint

Create thermal-solver scratch and write its handle to `handle_out`.

One workspace belongs to one column being stepped: it is mutated in place, so concurrent
threads must not share a handle. Its buffers grow to the column on first use and are reused
across steps, which is the reason to create one per column rather than one per call.
"""
Base.@ccallable function gemb_workspace_create(handle_out::Ptr{Cint})::Cint
    handle_out == C_NULL && return GEMB_ERR_BAD_ARGUMENT
    try
        h = _take_handle()
        WORKSPACES[h] = GEMB.ThermalWorkspace()
        unsafe_store!(handle_out, h)
        return GEMB_OK
    catch e
        _record_error("gemb_workspace_create", e)
        return GEMB_ERR_UNEXPECTED
    end
end

"""
    gemb_workspace_destroy(handle) -> Cint

Release thermal-solver scratch.
"""
Base.@ccallable function gemb_workspace_destroy(handle::Cint)::Cint
    haskey(WORKSPACES, handle) || return GEMB_ERR_BAD_HANDLE
    delete!(WORKSPACES, handle)
    return GEMB_OK
end

"""
    gemb_step(n, temperature, dz, density, water, grain_radius, grain_dendricity,
              grain_sphericity, age, scalars, forcing, params_handle, workspace_handle,
              flux_out, viscosity_out) -> Cint

Advance one `n`-cell column by one forcing step.

The eight profile arrays and `scalars` are read and then overwritten with the new column state.
The column keeps its cell count and its total depth, so `n` and `sum(dz)` are the same on
return as on entry.

`scalars` is `[evaporation_condensation, melt_surface]`. `forcing` is `GEMB_N_FORCING` doubles
in `ClimateForcingStep` field order, `flux_out` receives `GEMB_N_FLUX` doubles. `viscosity_out`
takes `n` doubles when the parameter set requested viscosity, and may be `NULL` otherwise.
"""
Base.@ccallable function gemb_step(
    n::Cint,
    temperature::Ptr{Cdouble},
    dz::Ptr{Cdouble},
    density::Ptr{Cdouble},
    water::Ptr{Cdouble},
    grain_radius::Ptr{Cdouble},
    grain_dendricity::Ptr{Cdouble},
    grain_sphericity::Ptr{Cdouble},
    age::Ptr{Cdouble},
    scalars::Ptr{Cdouble},
    forcing::Ptr{Cdouble},
    params_handle::Cint,
    workspace_handle::Cint,
    flux_out::Ptr{Cdouble},
    viscosity_out::Ptr{Cdouble},
)::Cint

    n >= 2 || return GEMB_ERR_BAD_ARGUMENT
    for p in (temperature, dz, density, water, grain_radius, grain_dendricity,
              grain_sphericity, age, scalars, forcing, flux_out)
        p == C_NULL && return GEMB_ERR_BAD_ARGUMENT
    end

    mp = get(PARAMS, params_handle, nothing)
    mp === nothing && return GEMB_ERR_BAD_HANDLE
    ws = get(WORKSPACES, workspace_handle, nothing)
    ws === nothing && return GEMB_ERR_BAD_HANDLE

    try
        m = Int(n)

        # Copies, not `unsafe_wrap` views: the physics modules rebind these names to freshly
        # allocated arrays as cells merge and split, so the results are written back below
        # rather than accumulating in the caller's memory.
        state = (
            temperature=[unsafe_load(temperature, i) for i in 1:m],
            dz=[unsafe_load(dz, i) for i in 1:m],
            density=[unsafe_load(density, i) for i in 1:m],
            water=[unsafe_load(water, i) for i in 1:m],
            grain_radius=[unsafe_load(grain_radius, i) for i in 1:m],
            grain_dendricity=[unsafe_load(grain_dendricity, i) for i in 1:m],
            grain_sphericity=[unsafe_load(grain_sphericity, i) for i in 1:m],
            age=[unsafe_load(age, i) for i in 1:m],
            evaporation_condensation=unsafe_load(scalars, 1),
            melt_surface=unsafe_load(scalars, 2),
        )

        cfs = GEMB.ClimateForcingStep(
            unsafe_load(forcing, 1), unsafe_load(forcing, 2), unsafe_load(forcing, 3),
            unsafe_load(forcing, 4), unsafe_load(forcing, 5), unsafe_load(forcing, 6),
            unsafe_load(forcing, 7), unsafe_load(forcing, 8), unsafe_load(forcing, 9),
            unsafe_load(forcing, 10), unsafe_load(forcing, 11), unsafe_load(forcing, 12),
            unsafe_load(forcing, 13), unsafe_load(forcing, 14), unsafe_load(forcing, 15),
            unsafe_load(forcing, 16), unsafe_load(forcing, 17), unsafe_load(forcing, 18),
            unsafe_load(forcing, 19), unsafe_load(forcing, 20), unsafe_load(forcing, 21),
            unsafe_load(forcing, 22),
        )

        new_state, flux = GEMB.gemb_core(state, cfs, mp, false; thermal_workspace=ws)

        for i in 1:m
            unsafe_store!(temperature, new_state.temperature[i], i)
            unsafe_store!(dz, new_state.dz[i], i)
            unsafe_store!(density, new_state.density[i], i)
            unsafe_store!(water, new_state.water[i], i)
            unsafe_store!(grain_radius, new_state.grain_radius[i], i)
            unsafe_store!(grain_dendricity, new_state.grain_dendricity[i], i)
            unsafe_store!(grain_sphericity, new_state.grain_sphericity[i], i)
            unsafe_store!(age, new_state.age[i], i)
        end
        unsafe_store!(scalars, new_state.evaporation_condensation, 1)
        unsafe_store!(scalars, new_state.melt_surface, 2)

        unsafe_store!(flux_out, flux.albedo_broadband, 1)
        unsafe_store!(flux_out, flux.shortwave_net, 2)
        unsafe_store!(flux_out, flux.heat_flux_sensible, 3)
        unsafe_store!(flux_out, flux.heat_flux_latent, 4)
        unsafe_store!(flux_out, flux.heat_flux_basal, 5)
        unsafe_store!(flux_out, flux.longwave_upward, 6)
        unsafe_store!(flux_out, flux.rain, 7)
        unsafe_store!(flux_out, flux.melt, 8)
        unsafe_store!(flux_out, flux.runoff, 9)
        unsafe_store!(flux_out, flux.refreeze, 10)
        unsafe_store!(flux_out, flux.mass_added, 11)
        unsafe_store!(flux_out, flux.mass_lateral, 12)
        unsafe_store!(flux_out, flux.mass_blowing_snow, 13)
        unsafe_store!(flux_out, flux.E_added, 14)
        unsafe_store!(flux_out, flux.densification_from_compaction, 15)
        unsafe_store!(flux_out, flux.densification_from_melt, 16)
        unsafe_store!(flux_out, flux.percolation_depth, 17)
        unsafe_store!(flux_out, flux.ice_slab_thickness, 18)
        unsafe_store!(flux_out, flux.ice_slab_depth, 19)
        unsafe_store!(flux_out, flux.aquifer_thickness, 20)
        unsafe_store!(flux_out, flux.aquifer_depth, 21)

        if viscosity_out != C_NULL
            for i in eachindex(flux.viscosity)
                unsafe_store!(viscosity_out, flux.viscosity[i], i)
            end
        end

        return GEMB_OK
    catch e
        _record_error("gemb_step", e)
        return GEMB_ERR_PHYSICS
    end
end

"""
    gemb_last_error(buf, len) -> Cint

Copy the most recent error message into `buf`, NUL-terminated, and return its length in bytes
excluding the terminator. A message longer than `len - 1` bytes is truncated.
"""
Base.@ccallable function gemb_last_error(buf::Ptr{Cchar}, len::Cint)::Cint
    (buf == C_NULL || len < 1) && return GEMB_ERR_BAD_ARGUMENT
    bytes = codeunits(LAST_ERROR[])
    k = min(length(bytes), Int(len) - 1)
    for i in 1:k
        unsafe_store!(buf, reinterpret(Cchar, bytes[i]), i)
    end
    unsafe_store!(buf, Cchar(0), k + 1)
    return Cint(k)
end

end # module
