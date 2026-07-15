"""
    gemb_spinup(profile::DimStack, cf::ClimateForcing, mp::ModelParameters, n_cycles::Int; verbose=false)

Run GEMB for multiple spinup cycles to reach quasi-steady state.

Forces `output_frequency=:last` internally to minimize memory usage during spinup.
Returns the spun-up profile DimStack.

Matches MATLAB's `gemb_spinup.m`.
"""
function gemb_spinup(profile::DimStack, cf::ClimateForcing, mp::ModelParameters, n_cycles::Int; verbose::Bool=false)
    # Force output_frequency to :last for spinup efficiency
    mp_spinup = ModelParameters(;
        (field => getfield(mp, field) for field in fieldnames(ModelParameters) if field != :output_frequency)...,
        output_frequency=:last
    )

    current_profile = profile
    for cycle in 1:n_cycles
        out = gemb(current_profile, cf, mp_spinup; verbose=verbose)
        current_profile = gemb_profile(out)
    end
    return current_profile
end
