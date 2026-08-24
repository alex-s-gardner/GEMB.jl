# The column, forcing and parameters both the Julia reference run and the C driver use.
#
# Kept in one place so the two cannot drift apart. `capi/test/test_capi.c` restates these values
# as literals — it cannot include Julia — and `capi/test/runtests.jl` checks that the C file's
# literals still match what is here.

using Printf

const N_CELLS = 20
const N_STEPS = 240
const DT = 10800.0

# Flux fields accumulated over the run, in the order `capi/gemb.h` declares them.
const FLUX_NAMES = (
    :albedo_broadband, :shortwave_net, :heat_flux_sensible, :heat_flux_latent,
    :heat_flux_basal, :longwave_upward, :rain, :melt, :runoff, :refreeze,
    :mass_added, :mass_lateral, :mass_blowing_snow, :E_added,
    :densification_from_compaction, :densification_from_melt, :percolation_depth,
    :ice_slab_thickness, :ice_slab_depth, :aquifer_thickness, :aquifer_depth,
)

"""
    reference_setup() -> (state, forcing, mp)

A melting-season column and one forcing step held constant across the run.

The forcing is warm enough and bright enough to melt, so the run exercises percolation,
refreezing and wet compaction rather than only cold diffusion — the paths where a boundary fault
is most likely to show up.
"""
function reference_setup()
    state = (
        temperature=fill(263.0, N_CELLS),
        dz=fill(0.5, N_CELLS),
        # A fixed increment per cell rather than `range(350, 800; length=N_CELLS)`: the C driver
        # has to reproduce this column exactly, and `range` computes its interior points to
        # better than naive interpolation, so the two would differ in the last bits.
        density=[350.0 + 20.0 * (i - 1) for i in 1:N_CELLS],
        water=zeros(N_CELLS),
        grain_radius=fill(0.5, N_CELLS),
        grain_dendricity=zeros(N_CELLS),
        grain_sphericity=fill(0.5, N_CELLS),
        age=zeros(N_CELLS),
        evaporation_condensation=0.0,
        melt_surface=0.0,
    )

    forcing = GEMB.ClimateForcingStep(
        DT,       # dt
        272.0,    # temperature_air
        90000.0,  # pressure_air
        0.5,      # precipitation
        5.0,      # wind_speed
        350.0,    # shortwave_downward
        280.0,    # longwave_downward
        400.0,    # vapor_pressure
        258.0,    # temperature_air_mean
        5.0,      # wind_speed_mean
        300.0,    # precipitation_mean
        2.0,      # temperature_observation_height
        10.0,     # wind_observation_height
        0.0,      # black_carbon_snow
        0.0,      # black_carbon_ice
        0.0,      # cloud_optical_thickness
        0.0,      # solar_zenith_angle
        0.0,      # shortwave_downward_diffuse
        0.1,      # cloud_fraction
    )

    mp = GEMB.ModelParameters(; dt_divisors=GEMB.fast_divisors(round(Int, DT * 10000)) ./ 10000)

    return state, forcing, mp
end
