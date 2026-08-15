# Specific heat capacity and the enthalpy it integrates to.
#
# Every energy budget in GEMB goes through `specific_enthalpy` rather than a bare
# `T * c_p` product. That is not stylistic: for a temperature-dependent `c_p` the two
# differ by a term comparable to the latent heat of fusion (see the `specific_enthalpy`
# docstring), so mixing the two forms breaks conservation outright.

"""
    heat_capacity(mp::ModelParameters, T) -> c_p [J kg-1 K-1]

Specific heat capacity of the ice matrix at temperature `T` [K], per
`mp.heat_capacity_method`:

- `:constant` — `mp.heat_capacity_ice` (default 2102, the melting-point value). MATLAB
  behaviour.
- `:CuffeyPaterson` — Cuffey and Paterson (2010) eq. 9.1, `c_p = 152.5 + 7.122·T`. Gives
  2097.9 at the melting point but only 1862 at 240 K and 1648 at 210 K, so the constant
  value over-permits refreezing in cold firn by +12.6% and +27.5% respectively.

Applies to snow, firn, and ice alike: the load-bearing matrix is ice, and pore air and
pore water contribute nothing to this term.

For anything that accumulates or balances joules, use [`specific_enthalpy`](@ref) — not
`T * heat_capacity(mp, T)`.
"""
@inline function heat_capacity(mp::ModelParameters, T::Real)
    if mp.heat_capacity_method === :constant
        return mp.heat_capacity_ice
    end
    return HEAT_CAPACITY_CUFFEY_A + HEAT_CAPACITY_CUFFEY_B * Float64(T)
end

"""
    specific_enthalpy(mp::ModelParameters, T) -> h [J kg-1]

Specific enthalpy of the ice matrix at temperature `T` [K], measured from a 0 K
reference: `h(T) = ∫₀ᵀ c_p dT′`.

- `:constant` — `c·T` exactly, so the default path is bit-identical to the MATLAB
  `M·T·C_ICE` accounting.
- `:CuffeyPaterson` — `a·T + (b/2)T²`.

!!! warning
    Never use `T * heat_capacity(mp, T)` in an energy budget. For `c_p = a + bT` that
    product is `a·T + b·T²`, which overstates the enthalpy by `(b/2)T² = 2.66e5 J kg-1` at
    the melting point — 0.79 × the latent heat of fusion. It is a dominant error, not a
    correction. `test/test_heat_capacity.jl` guards against reintroducing it.

The 0 K reference makes absolute enthalpies incomparable *between* methods (574 061 vs
307 385 J kg-1 at the melting point), so a reported `E_thermal` jumps when
`heat_capacity_method` changes. Only differences are physically meaningful.

See [`temperature_from_specific_enthalpy`](@ref) for the inverse.
"""
@inline function specific_enthalpy(mp::ModelParameters, T::Real)
    if mp.heat_capacity_method === :constant
        return mp.heat_capacity_ice * Float64(T)
    end
    Tf = Float64(T)
    return HEAT_CAPACITY_CUFFEY_A * Tf + (HEAT_CAPACITY_CUFFEY_B / 2) * Tf * Tf
end

"""
    temperature_from_specific_enthalpy(mp::ModelParameters, h) -> T [K]

Inverse of [`specific_enthalpy`](@ref). `h` is in J kg-1 from a 0 K reference.

For `:CuffeyPaterson` this inverts the monotone quadratic `h = aT + (b/2)T²` using the
cancellation-free root

    T = 2h / (a + √(a² + 2bh))

which is well-conditioned for small `h`, unlike the textbook quadratic formula. `c_p > 0`
over the physical range, so `h` is strictly increasing and the root is unique.
"""
@inline function temperature_from_specific_enthalpy(mp::ModelParameters, h::Real)
    if mp.heat_capacity_method === :constant
        return Float64(h) / mp.heat_capacity_ice
    end
    hf = Float64(h)
    a = HEAT_CAPACITY_CUFFEY_A
    b = HEAT_CAPACITY_CUFFEY_B
    return 2 * hf / (a + sqrt(a * a + 2 * b * hf))
end

"""
    enthalpy_temperature_scale(mp::ModelParameters) -> s

Divisor that [`temperature_from_scaled_enthalpy`](@ref) expects to have been applied to a
specific enthalpy already: `mp.heat_capacity_ice` under `:constant`, `1.0` otherwise.

Exists so a caller with a per-cell mass reciprocal can fold `1/c` into it once and turn the
inverse map into a multiply. See [`temperature_from_scaled_enthalpy`](@ref).
"""
@inline enthalpy_temperature_scale(mp::ModelParameters) =
    mp.heat_capacity_method === :constant ? mp.heat_capacity_ice : 1.0

"""
    temperature_from_scaled_enthalpy(mp::ModelParameters, hs) -> T [K]

[`temperature_from_specific_enthalpy`](@ref) for a caller that has already divided the
specific enthalpy by [`enthalpy_temperature_scale`](@ref). Under `:constant` that division
*is* the whole inverse, so this is the identity; under `:CuffeyPaterson` the scale is 1 and
this is the plain inverse.

The point is to get the division off an inner loop. The thermal solver converts enthalpy to
temperature twice per cell per sub-step; folding `1/c` into the mass reciprocal it already
precomputes replaces both divisions with multiplies, which is worth ~40% of the solver's
runtime — division is the one arithmetic op on that loop that does not pipeline.

Folding the two constants into one reciprocal is a reassociation, so results differ from
`h/M/c` in the last bit or two. That is not gated by the bench fingerprint: the 75-year
non-converged spinup it drives amplifies any 1-ULP perturbation to O(1e-2) in the melt
fields, as a control reassociation in the *unmodified* solver confirms.
"""
@inline function temperature_from_scaled_enthalpy(mp::ModelParameters, hs::Real)
    mp.heat_capacity_method === :constant && return Float64(hs)
    return temperature_from_specific_enthalpy(mp, hs)
end

"""
    specific_enthalpy_water(mp::ModelParameters) -> h [J kg-1]

Specific enthalpy of liquid water at the melting point: the enthalpy of ice at `CtoK`
plus the latent heat of fusion. This is the `LF + CtoK * C_ICE` of the MATLAB budgets,
generalized to a temperature-dependent `c_p`.

GEMB carries pore water only at the melting point, so no liquid heat capacity is needed.

See the two-argument form for liquid *above* the melting point, which only rain is.
"""
@inline specific_enthalpy_water(mp::ModelParameters) = LF + specific_enthalpy(mp, CtoK)

"""
    specific_enthalpy_water(mp::ModelParameters, T) -> h [J kg-1]

Specific enthalpy of liquid water at temperature `T` [K]: melting-point water plus its
sensible heat above the melting point, `h_water + c_water·(T − CtoK)`.

Rain is the only liquid GEMB carries above `CtoK`, and it is the only caller. Everything
else — pore water, refreezing, runoff — is isothermal at the melting point and uses the
one-argument form, which this reduces to exactly at `T == CtoK`.

`mp.rain_heat_capacity` selects the heat capacity of that sensible term:

- `:water` (default) — [`HEAT_CAPACITY_WATER`](@ref) ≈ 4220 J kg-1 K-1, the physical value.
- `:ice` — `heat_capacity(mp, CtoK)` ≈ 2102, which is what MATLAB and GEMB.jl before this
  option used. It understates the rain's sensible heat by about a factor two: ~3.9 kJ kg-1
  for rain at 275 K, against 334.5 kJ kg-1 of `LF`. Small (~1% of the rain's energy) but
  systematic, and growing linearly with `T_air − CtoK`.

Below `CtoK` this returns the melting-point value unchanged rather than extrapolating: liquid
colder than the melting point is not a state the column represents, and the callers reach
this only on the rain branch, which is gated on `T_air > mp.rain_temperature_threshold`.
"""
@inline function specific_enthalpy_water(mp::ModelParameters, T::Real)
    dT = Float64(T) - CtoK
    dT <= 0.0 && return specific_enthalpy_water(mp)
    c = mp.rain_heat_capacity === :water ? HEAT_CAPACITY_WATER : heat_capacity(mp, CtoK)
    return specific_enthalpy_water(mp) + c * dT
end

"""
    mix_temperature(mp::ModelParameters, T1, M1, T2, M2) -> T [K]
    mix_temperature(mp::ModelParameters, h1, M1, h2, M2, ::Val{:enthalpy}) -> T [K]

Temperature of the mixture of two masses, conserving enthalpy: solve
`M₁h(T₁) + M₂h(T₂) = (M₁+M₂)·h(T)` for `T`.

This is the only correct way to combine cells of different temperature under a
temperature-dependent `c_p`. The mass-weighted mean temperature
`(M₁T₁ + M₂T₂)/(M₁+M₂)` is *not* equivalent: `h` is convex, so the mean loses
`(b/2)·M·Var_M(T)` joules — about 2 225 J for an equal-mass pair 5 K apart, against an
energy-conservation tolerance of `E_TOLERANCE = 1e-3` J.

The second form takes specific enthalpies directly, for callers that already hold them
(the mixing partner's enthalpy is often `specific_enthalpy_water` or `h(T_air) + LF`
rather than a temperature).

Under `:constant` the first form mixes in temperature space. That is not an approximation:
for constant `c_p` the `c` cancels out of the enthalpy balance, so the mass-weighted mean
*is* the exact enthalpy-conserving mixture. Doing it directly rather than via a
`h`/`h⁻¹` round-trip keeps the default path bit-identical to the MATLAB arithmetic, which
would otherwise drift at ~1e-13 per merge.

Returns `T1` when the total mass is zero.
"""
@inline function mix_temperature(mp::ModelParameters, T1::Real, M1::Real, T2::Real, M2::Real)
    if mp.heat_capacity_method === :constant
        M = Float64(M1) + Float64(M2)
        M <= 0.0 && return Float64(T1)
        return (Float64(T1) * Float64(M1) + Float64(T2) * Float64(M2)) / M
    end
    return mix_temperature(mp, specific_enthalpy(mp, T1), M1,
        specific_enthalpy(mp, T2), M2, Val(:enthalpy))
end

@inline function mix_temperature(mp::ModelParameters, h1::Real, M1::Real, h2::Real,
    M2::Real, ::Val{:enthalpy})
    M = Float64(M1) + Float64(M2)
    M <= 0.0 && return temperature_from_specific_enthalpy(mp, h1)
    h = (Float64(M1) * Float64(h1) + Float64(M2) * Float64(h2)) / M
    return temperature_from_specific_enthalpy(mp, h)
end

"""
    mix_temperature_liquid(mp::ModelParameters, T_liquid, M_liquid, T_solid, M_solid) -> T [K]

Temperature after mixing `M_liquid` [kg] of liquid water at `T_liquid` [K] into `M_solid`
[kg] of ice matrix at `T_solid` [K], conserving enthalpy. The liquid carries the latent heat
of fusion, so its specific enthalpy is [`specific_enthalpy_water`](@ref)`(mp, T_liquid)`.

Two callers, and they sit on opposite sides of the melting point:

- **Rain on snow** (`calculate_accumulation`), with `T_liquid = T_air > CtoK`. The rain's
  sensible heat above the melting point is carried at `mp.rain_heat_capacity` — see
  [`specific_enthalpy_water`](@ref).
- **Refreezing** ([`refreeze_temperature`](@ref)), with `T_liquid = CtoK` exactly. The
  sensible term vanishes and this reduces to the melting-point form.

Under `:constant` the `c` cancels out of the enthalpy balance, so both are evaluated as a
weighted mean of temperatures with the liquid at an effective temperature `h_liquid/c` —
`T_liquid + LF/c` when the liquid's heat capacity is the matrix's, which is the grouping the
MATLAB reference uses. Keeping that grouping is what makes the refreeze path (and
`rain_heat_capacity = :ice`) bit-identical to the reference rather than drifting through an
`h`/`h⁻¹` round-trip.

Returns `T_solid` when the total mass is zero.
"""
@inline function mix_temperature_liquid(mp::ModelParameters, T_liquid::Real, M_liquid::Real,
    T_solid::Real, M_solid::Real)
    # Melting-point liquid, or liquid whose sensible heat is deliberately carried at the ice
    # heat capacity: the reference's arithmetic, preserved bit-for-bit.
    if mp.rain_heat_capacity === :ice || Float64(T_liquid) <= CtoK
        if mp.heat_capacity_method === :constant
            M = Float64(M_liquid) + Float64(M_solid)
            M <= 0.0 && return Float64(T_solid)
            return (Float64(M_liquid) * (Float64(T_liquid) + LF / mp.heat_capacity_ice) +
                    Float64(T_solid) * Float64(M_solid)) / M
        end
        return mix_temperature(mp, specific_enthalpy(mp, T_liquid) + LF, M_liquid,
            specific_enthalpy(mp, T_solid), M_solid, Val(:enthalpy))
    end
    # Above-freezing liquid (rain) with the liquid-water heat capacity.
    if mp.heat_capacity_method === :constant
        M = Float64(M_liquid) + Float64(M_solid)
        M <= 0.0 && return Float64(T_solid)
        c = mp.heat_capacity_ice
        T_liquid_effective = CtoK + LF / c + (HEAT_CAPACITY_WATER / c) * (Float64(T_liquid) - CtoK)
        return (Float64(M_liquid) * T_liquid_effective +
                Float64(T_solid) * Float64(M_solid)) / M
    end
    return mix_temperature(mp, specific_enthalpy_water(mp, T_liquid), M_liquid,
        specific_enthalpy(mp, T_solid), M_solid, Val(:enthalpy))
end

"""
    cold_content_mass(mp::ModelParameters, T, M) -> [kg]

Mass of melting-point water that mass `M` [kg] of ice matrix at temperature `T` [K] can
refreeze before reaching the melting point: `M·(h(CtoK) − h(T))/LF`. Zero for `T >= CtoK`.

This is the `freeze_max` of the melt equations. Making it the enthalpy difference rather
than `M·c·(CtoK−T)` is what keeps it consistent with
[`refreeze_temperature`](@ref) — if the two used different energy maps, refreezing the full
cold content would not land the cell exactly at the melting point.
"""
@inline function cold_content_mass(mp::ModelParameters, T::Real, M::Real)
    if mp.heat_capacity_method === :constant
        # Grouped as the MATLAB reference groups it, so the default path is bit-identical.
        return max(0.0, -((Float64(T) - CtoK) * Float64(M) * mp.heat_capacity_ice) / LF)
    end
    return max(0.0, Float64(M) * (specific_enthalpy(mp, CtoK) - specific_enthalpy(mp, T)) / LF)
end

"""
    cold_content_mass(mp::ModelParameters, T, density, dz) -> [kg]

As above, for a cell given by its density [kg m-3] and thickness [m] rather than its mass.

Not merely `cold_content_mass(mp, T, density*dz)`: the reference multiplies `ρ` and `dz` into
the product separately at this site and pre-multiplied at the other, and floating-point
multiplication is not associative. Keeping both groupings is what makes the `:constant` path
bit-identical at *both* call sites.
"""
@inline function cold_content_mass(mp::ModelParameters, T::Real, density::Real, dz::Real)
    if mp.heat_capacity_method === :constant
        return max(0.0, -((Float64(T) - CtoK) * Float64(density) * Float64(dz) *
                          mp.heat_capacity_ice) / LF)
    end
    return cold_content_mass(mp, T, Float64(density) * Float64(dz))
end

"""
    melt_mass_from_excess(mp::ModelParameters, T_excess, density, dz) -> [kg]

Mass of ice [kg] that the enthalpy held above the melting point can melt, for a cell of
`density` [kg m-3] and thickness `dz` [m] whose temperature exceeded `CtoK` by `T_excess`
[K]: `h_excess·ρ·dz/LF`.

Grouped to match the reference's `T_excess·ρ·dz·c/LF` on the `:constant` path.
"""
@inline function melt_mass_from_excess(mp::ModelParameters, T_excess::Real, density::Real,
    dz::Real)
    if mp.heat_capacity_method === :constant
        return Float64(T_excess) * Float64(density) * Float64(dz) * mp.heat_capacity_ice / LF
    end
    return excess_specific_enthalpy(mp, T_excess) * Float64(density) * Float64(dz) / LF
end

"""
    refreeze_temperature(mp::ModelParameters, T, M_new, freeze) -> T [K]

Temperature after `freeze` [kg] of melting-point pore water refreezes into ice matrix that
was at temperature `T` [K]. `M_new` is the matrix mass *after* the refrozen mass is added,
matching the melt equations' bookkeeping order.

Releases both the latent heat and the sensible heat of cooling the new ice from `CtoK` to
the mixture temperature — i.e. it is [`mix_temperature_liquid`](@ref) with the liquid at the
melting point, written in the reference's incremental form.
"""
@inline function refreeze_temperature(mp::ModelParameters, T::Real, M_new::Real, freeze::Real)
    Mn = Float64(M_new)
    Mn <= 0.0 && return Float64(T)
    if mp.heat_capacity_method === :constant
        c = mp.heat_capacity_ice
        return Float64(T) + (Float64(freeze) * (LF + (CtoK - Float64(T)) * c) / (Mn * c))
    end
    return mix_temperature_liquid(mp, CtoK, freeze, T, Mn - Float64(freeze))
end

"""
    excess_specific_enthalpy(mp::ModelParameters, T_excess) -> [J kg-1]

Specific enthalpy held above the melting point by ice that was `T_excess` [K] warmer than
`CtoK`: `h(CtoK + T_excess) − h(CtoK)`. This is the energy available to melt the cell, so
comparing it directly against `LF` replaces the reference's `LF/c_p` temperature threshold —
a quantity that only exists because `h` was linear.

Evaluated in the cancellation-free form `a·ΔT + (b/2)·ΔT·(ΔT + 2·CtoK)` rather than as a
difference of two large enthalpies.
"""
@inline function excess_specific_enthalpy(mp::ModelParameters, T_excess::Real)
    dT = Float64(T_excess)
    dT <= 0.0 && return 0.0
    if mp.heat_capacity_method === :constant
        return mp.heat_capacity_ice * dT
    end
    return HEAT_CAPACITY_CUFFEY_A * dT +
           (HEAT_CAPACITY_CUFFEY_B / 2) * dT * (dT + 2 * CtoK)
end

"""
    column_enthalpy(mp::ModelParameters, M, temperature, water) -> E [J]

Total enthalpy [J] of a column given per-cell matrix mass `M` [kg], `temperature` [K], and
pore water `water` [kg m-2]: `Σ Mᵢh(Tᵢ) + (Σwaterᵢ)·h_water`.

Written as an explicit loop rather than a generator. The verbose energy budgets sit inside
large functions (`calculate_melt` especially), and a generator's closure there was enough to
push inference over its limit and infer *every* local in the enclosing function as `Any` —
boxing the whole hot path for an 18x allocation regression, even with `verbose=false`.
"""
@inline function column_enthalpy(mp::ModelParameters, M::AbstractVector{Float64},
    temperature::AbstractVector{Float64})
    E = 0.0
    @inbounds for i in eachindex(M)
        E += M[i] * specific_enthalpy(mp, temperature[i])
    end
    return E
end

@inline function column_enthalpy(mp::ModelParameters, M::AbstractVector{Float64},
    temperature::AbstractVector{Float64}, water::AbstractVector{Float64})
    return column_enthalpy(mp, M, temperature) + sum(water) * specific_enthalpy_water(mp)
end

"""
    surplus_energy(mp::ModelParameters, T_excess, T_surplus, M) -> Q [J]

Energy held by mass `M` [kg] *beyond* what melting its entire mass would consume, given a
temperature excess `T_excess` [K] above the melting point and the corresponding surplus
`T_surplus = max(0, T_excess − full_melt_excess_temperature(mp))`.

Zero when `T_surplus` is zero, so this agrees exactly with the kelvin-space mask the melt
equations branch on. Under `:constant` it evaluates as `T_surplus·c·M` (the reference's
grouping); otherwise as `(h_excess(T_excess) − LF)·M`.
"""
@inline function surplus_energy(mp::ModelParameters, T_excess::Real, T_surplus::Real, M::Real)
    Float64(T_surplus) <= 0.0 && return 0.0
    if mp.heat_capacity_method === :constant
        return Float64(T_surplus) * mp.heat_capacity_ice * Float64(M)
    end
    return (excess_specific_enthalpy(mp, T_excess) - LF) * Float64(M)
end

"""
    excess_temperature_from_specific_enthalpy(mp::ModelParameters, e) -> ΔT [K]

Inverse of [`excess_specific_enthalpy`](@ref): the temperature excess above `CtoK` that
holds `e` [J kg-1] of excess specific enthalpy. Uses the cancellation-free root of
`(b/2)ΔT² + (a + b·CtoK)ΔT − e = 0`.

Needed because the melt equations branch on kelvin thresholds (`T_TOLERANCE` is in kelvin)
while the arithmetic is in joules, so the two representations must be interconvertible.
"""
@inline function excess_temperature_from_specific_enthalpy(mp::ModelParameters, e::Real)
    ef = Float64(e)
    ef <= 0.0 && return 0.0
    if mp.heat_capacity_method === :constant
        return ef / mp.heat_capacity_ice
    end
    a = HEAT_CAPACITY_CUFFEY_A + HEAT_CAPACITY_CUFFEY_B * CtoK
    return 2 * ef / (a + sqrt(a * a + 2 * HEAT_CAPACITY_CUFFEY_B * ef))
end

"""
    add_energy_temperature(mp::ModelParameters, T, M, Q) -> T [K]

Temperature of mass `M` [kg] at temperature `T` [K] after absorbing `Q` joules:
`h⁻¹(h(T) + Q/M)`.

For constant `c_p` this is evaluated as `Q/M/c + T` rather than via an `h`/`h⁻¹` round-trip,
so the default path is bit-identical to the reference arithmetic. Returns `T` unchanged when
`M <= 0`.
"""
@inline function add_energy_temperature(mp::ModelParameters, T::Real, M::Real, Q::Real)
    Mf = Float64(M)
    Mf <= 0.0 && return Float64(T)
    if mp.heat_capacity_method === :constant
        return Float64(Q) / Mf / mp.heat_capacity_ice + Float64(T)
    end
    return temperature_from_specific_enthalpy(mp, specific_enthalpy(mp, T) + Float64(Q) / Mf)
end

"""
    full_melt_excess_temperature(mp::ModelParameters) -> ΔT [K]

Temperature excess above the melting point at which a cell holds exactly enough enthalpy to
melt its entire mass — the `LF/c_p` of the reference, generalized. A cell hotter than this
has surplus energy to pass downward.

Kept in kelvin because the melt equations' branch tests use a kelvin tolerance
(`T_tolerance = 1e-10`); reformulating those tests in joules would silently rescale them by
~2100.
"""
@inline function full_melt_excess_temperature(mp::ModelParameters)
    if mp.heat_capacity_method === :constant
        return LF / mp.heat_capacity_ice
    end
    return excess_temperature_from_specific_enthalpy(mp, LF)
end
