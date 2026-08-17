---
---

# Internals {#Internals}

These pages document GEMB's non-exported internals: the physics kernels called by [`gemb`](/api#GEMB.gemb), the grid controllers, and the thermodynamic and hydrologic helpers.

They are **not** part of the public API. They are documented because the physics lives here — each scheme's docstring carries the citation for the law it implements — and because the narrative pages link to them. Their names and signatures may change in a minor release; depend on the public API in [API Reference](/api#API-Reference) instead.

The internals are split across three pages by subsystem:
- [Internals: Physics](/internals_physics#Internals:-Physics) — the per-timestep physics kernels and their schemes
  
- [Internals: Grid and Column](/internals_grid#Internals:-Grid-and-Column) — grid controllers, layer management, profile setup
  
- [Internals: Support](/internals_support#Internals:-Support) — constants, metadata, interpolation, and other helpers
  
