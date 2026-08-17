# Internals: Physics

Non-exported physics kernels, called once per timestep from `gemb_core`. See
[Internals](@ref) for the caveat on depending on these.

```@autodocs
Modules = [GEMB]
Public = false
Pages = [
    "calculate_accumulation.jl",
    "calculate_albedo.jl",
    "calculate_density.jl",
    "calculate_grain_size.jl",
    "calculate_melt.jl",
    "calculate_shortwave_radiation.jl",
    "calculate_temperature.jl",
    "densification_lookup.jl",
    "heat_capacity.jl",
    "hydrology.jl",
    "thermal_conductivity.jl",
    "turbulent_heat_flux.jl",
]
Order = [:module, :type, :constant, :function, :macro]
```
