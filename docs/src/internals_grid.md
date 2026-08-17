# Internals: Grid and Column

Non-exported grid controllers, layer management, and column initialization. See
[Internals](@ref) for the caveat on depending on these.

```@autodocs
Modules = [GEMB]
Public = false
Pages = [
    "grid_ops.jl",
    "manage_layer_thickness.jl",
    "initialize_profile.jl",
    "initialize_climate_summary.jl",
    "initialize_forcing.jl",
    "initialize_parameters.jl",
    "steady_state_profile.jl",
    "gemb_core.jl",
    "gemb_driver.jl",
    "spinup.jl",
    "profile_extract.jl",
    "forcing_climatology.jl",
]
Order = [:module, :type, :constant, :function, :macro]
```
