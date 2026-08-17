# Internals

These pages document GEMB's non-exported internals: the physics kernels called by
[`gemb`](@ref), the grid controllers, and the thermodynamic and hydrologic helpers.

They are **not** part of the public API. They are documented because the physics lives
here — each scheme's docstring carries the citation for the law it implements — and
because the narrative pages link to them. Their names and signatures may change in a
minor release; depend on the public API in [API Reference](@ref) instead.

The internals are split across three pages by subsystem:

- [Internals: Physics](@ref) — the per-timestep physics kernels and their schemes
- [Internals: Grid and Column](@ref) — grid controllers, layer management, profile setup
- [Internals: Support](@ref) — constants, metadata, interpolation, and other helpers
