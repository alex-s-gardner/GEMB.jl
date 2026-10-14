# Tabulated Monin-Obukhov stability for the thermal solvers. `obukhov_stability` solves
# z_w/L = Ri·F_m²/F_h by four Newton steps on every flux evaluation; within a run that relation
# depends only on the measurement heights and the roughness lengths, so the solvers instead read
# its converged solution from a table built once per (heights, roughness) and cached in the
# `ThermalWorkspace`.

# Transfer functions and their logarithmic derivatives at stability `ζ`:
# q = ln|Ri(ζ)| with Ri = ζ·F_h/F_m², G = ln(F_m·F_h), and dq/ds, dG/ds with s = ln|ζ|.
@inline function _stability_point(ζ::Float64, logM::Float64, logHT::Float64, r0::Float64,
    rT::Float64, rTs::Float64)
    ψm, φm = _stability_m(ζ)
    ψm0, φm0 = _stability_m(ζ * r0)
    ψhT, φhT = _stability_h(ζ * rT)
    ψh0, φh0 = _stability_h(ζ * rTs)
    Fm = logM - ψm + ψm0
    Fh = logHT - ψhT + ψh0
    dlnFm = (φm - φm0) / Fm         # ζ·F_m'/F_m = d ln F_m / ds
    dlnFh = (φhT - φh0) / Fh
    q = log(abs(ζ)) + log(Fh / (Fm * Fm))
    return q, 1 + dlnFh - 2 * dlnFm, log(Fm * Fh), dlnFm + dlnFh
end

# Converged `s = ln|ζ|` with `ln|Ri(ζ)| = q` for the sign `sgn` of Ri, by Newton from `s0`.
function _stability_solve(q::Float64, sgn::Float64, s0::Float64, logM::Float64,
    logHT::Float64, r0::Float64, rT::Float64, rTs::Float64)
    s = s0
    for _ in 1:100
        q_s, dq, _, _ = _stability_point(sgn * exp(s), logM, logHT, r0, rT, rTs)
        Δ = (q_s - q) / dq
        s -= Δ
        abs(Δ) <= 1e-14 * max(1.0, abs(s)) && return s
    end
    error("_stability_solve: Newton on ln|ζ| did not converge for ln|Ri| = $(q)")
end

"""
    _stability_table(z_w, z_T, z0, zT) -> StabilityTable

Build the [`StabilityTable`](@ref) for wind height `z_w`, temperature and humidity height
`z_T`, momentum roughness `z0` and scalar roughness `zT` [m]. Nodes are uniform in
`q = ln|Ri|` from `ln(STABILITY_TABLE_RI_MIN)`, spaced at most `STABILITY_TABLE_DQ`: up to
`ln(STABILITY_TABLE_RI_MAX)` for stable `Ri`, and up to the `Ri` at which `ζ` reaches
`ZETA_UNSTABLE_MIN` for unstable `Ri`. At each node `ζ` is the converged solution of
`z_w/L = Ri·F_m²/F_h` (the same Beljaars-Holtslag / Paulson-Högström functions as
[`obukhov_stability`](@ref)), and `ds/dq` and `dG/dq` are exact.
"""
function _stability_table(z_w::Float64, z_T::Float64, z0::Float64, zT::Float64)
    logM = log(z_w / z0)
    logHT = log(z_T / zT)
    r0, rT, rTs = z0 / z_w, z_T / z_w, zT / z_w
    q_min = log(STABILITY_TABLE_RI_MIN)
    s_neutral = 2 * log(logM) - log(logHT)          # ln|ζ| − ln|Ri| at neutral

    function nodes(sgn, q_max)
        n = ceil(Int, (q_max - q_min) / STABILITY_TABLE_DQ) + 1
        h = (q_max - q_min) / (n - 1)
        s = Vector{Float64}(undef, n); ds = similar(s); G = similar(s); dG = similar(s)
        s_prev = q_min + s_neutral
        for k in 1:n
            q = q_min + (k - 1) * h
            s[k] = _stability_solve(q, sgn, s_prev, logM, logHT, r0, rT, rTs)
            _, dq, G[k], dGds = _stability_point(sgn * exp(s[k]), logM, logHT, r0, rT, rTs)
            ds[k] = 1 / dq
            dG[k] = dGds / dq
            s_prev = s[k]
        end
        return h, s, ds, G, dG
    end

    stable = nodes(1.0, log(STABILITY_TABLE_RI_MAX))
    q_clamp, = _stability_point(ZETA_UNSTABLE_MIN, logM, logHT, r0, rT, rTs)
    unstable = nodes(-1.0, q_clamp)
    return StabilityTable(z_w, z_T, z0, zT, logM, logHT, log(logM * logHT), q_min,
        stable..., unstable...)
end

"""
    _stability_table!(workspace, z_w, z_T, z0, zT) -> StabilityTable

The cached [`StabilityTable`](@ref) for these heights and roughness lengths, built and stored in
`workspace` on first use. A workspace belongs to one run, so the cache needs no locking.
"""
function _stability_table!(ws::ThermalWorkspace, z_w::Float64, z_T::Float64, z0::Float64,
    zT::Float64)
    for tab in ws.stability_tables
        tab.z_w == z_w && tab.z_T == z_T && tab.z0 == z0 && tab.zT == zT && return tab
    end
    tab = _stability_table(z_w, z_T, z0, zT)
    push!(ws.stability_tables, tab)
    return tab
end

# Cubic Hermite interpolation of node values `y` with node derivatives `d` (per unit q) on a
# grid of spacing `h` from `q_min`; returns the value and its derivative with respect to q.
@inline function _hermite(y::Vector{Float64}, d::Vector{Float64}, q_min::Float64, h::Float64,
    q::Float64)
    t = (q - q_min) / h
    k = clamp(floor(Int, t) + 1, 1, length(y) - 1)
    x = t - (k - 1)
    x2 = x * x
    x3 = x2 * x
    y0, y1 = y[k], y[k+1]
    m0, m1 = h * d[k], h * d[k+1]
    v = (2x3 - 3x2 + 1) * y0 + (x3 - 2x2 + x) * m0 + (-2x3 + 3x2) * y1 + (x3 - x2) * m1
    dv = ((6x2 - 6x) * y0 + (3x2 - 4x + 1) * m0 + (-6x2 + 6x) * y1 + (3x2 - 2x) * m1) / h
    return v, dv
end

"""
    _stability_lookup(tab, Ri) -> (G, dG_dRi)

`G = ln(F_m·F_h)` and `dG/dRi` for bulk Richardson number `Ri`, from the table:

- `|Ri| < STABILITY_TABLE_RI_MIN`: linear in `Ri` between the neutral value and the first node;
- `Ri > STABILITY_TABLE_RI_MAX`: the converged solution, computed for this `Ri`;
- `Ri` at or beyond the unstable end, where `ζ` reaches `ZETA_UNSTABLE_MIN`: the value there,
  constant in `Ri` (the bound of [`obukhov_stability`](@ref));
- otherwise cubic Hermite interpolation in `ln|Ri|`.
"""
@inline function _stability_lookup(tab::StabilityTable, Ri::Float64)
    isfinite(Ri) || error("_stability_lookup: bulk Richardson number is $(Ri)")
    a = abs(Ri)
    if a < STABILITY_TABLE_RI_MIN
        G1 = Ri >= 0.0 ? tab.stable_G[1] : tab.unstable_G[1]
        slope = (G1 - tab.G_neutral) / STABILITY_TABLE_RI_MIN
        return tab.G_neutral + slope * a, Ri >= 0.0 ? slope : -slope
    end
    q = log(a)
    if Ri > 0.0
        q > tab.q_min + (length(tab.stable_G) - 1) * tab.stable_h &&
            return _stability_exact_G(tab, q)
        G, dG = _hermite(tab.stable_G, tab.stable_dG, tab.q_min, tab.stable_h, q)
    else
        q >= tab.q_min + (length(tab.unstable_G) - 1) * tab.unstable_h &&
            return tab.unstable_G[end], 0.0
        G, dG = _hermite(tab.unstable_G, tab.unstable_dG, tab.q_min, tab.unstable_h, q)
    end
    return G, dG / Ri                     # dq/dRi = 1/Ri on either side
end

# Converged `(G, dG/dRi)` for a stable `Ri = exp(q)` beyond the table.
@noinline function _stability_exact_G(tab::StabilityTable, q::Float64)
    r0, rT, rTs = tab.z0 / tab.z_w, tab.z_T / tab.z_w, tab.zT / tab.z_w
    s0 = tab.stable_s[end] + (q - (tab.q_min + (length(tab.stable_s) - 1) * tab.stable_h)) *
         tab.stable_ds[end]
    s = _stability_solve(q, 1.0, s0, tab.logM, tab.logHT, r0, rT, rTs)
    _, dq, G, dGds = _stability_point(exp(s), tab.logM, tab.logHT, r0, rT, rTs)
    return G, (dGds / dq) / exp(q)
end

"""
    _stability_zeta(tab, Ri) -> ζ

Stability parameter `ζ = z_w/L` interpolated from the table, with the same regimes as
[`_stability_lookup`](@ref): `ζ` linear in `Ri` below `STABILITY_TABLE_RI_MIN`, exact above
`STABILITY_TABLE_RI_MAX`, and `ZETA_UNSTABLE_MIN` beyond the unstable end.
"""
function _stability_zeta(tab::StabilityTable, Ri::Float64)
    Ri == 0.0 && return 0.0
    a = abs(Ri)
    q = log(a)
    if a < STABILITY_TABLE_RI_MIN
        s1 = Ri > 0.0 ? tab.stable_s[1] : tab.unstable_s[1]
        return sign(Ri) * exp(s1) * a / STABILITY_TABLE_RI_MIN
    end
    if Ri > 0.0
        if q > tab.q_min + (length(tab.stable_s) - 1) * tab.stable_h
            r0, rT, rTs = tab.z0 / tab.z_w, tab.z_T / tab.z_w, tab.zT / tab.z_w
            return exp(_stability_solve(q, 1.0, tab.stable_s[end], tab.logM, tab.logHT, r0,
                rT, rTs))
        end
        return exp(_hermite(tab.stable_s, tab.stable_ds, tab.q_min, tab.stable_h, q)[1])
    end
    q >= tab.q_min + (length(tab.unstable_s) - 1) * tab.unstable_h && return ZETA_UNSTABLE_MIN
    return -exp(_hermite(tab.unstable_s, tab.unstable_ds, tab.q_min, tab.unstable_h, q)[1])
end

"""
    _table_turbulent_heat_flux(tab, T_surface, density_air, C, ri_scale, cfs)
        -> (heat_flux_sensible, heat_flux_latent, latent_heat, dshf_dT, dlhf_dT)

The fluxes of [`_turbulent_heat_flux`](@ref) with `ln(F_m·F_h)` read from `tab` (heat and
moisture share their roughness length), and their exact derivatives with respect to
`T_surface` through the Richardson number, the mean air temperature and the saturation vapor
pressure. The roughness logs and heights are those the table was built for.
"""
@inline function _table_turbulent_heat_flux(tab::StabilityTable, T_surface::Float64,
    density_air::Float64, C::Float64, ri_scale::Float64, cfs::ClimateForcingStep)

    T_air = cfs.temperature_air
    T_mean = 0.5 * (T_air + T_surface)
    ΔT = T_air - T_surface
    Ri = ri_scale * ΔT / T_mean
    G, dG_dRi = _stability_lookup(tab, Ri)
    inv_FF = exp(-G)
    # dRi/dT_s = −ri_scale·T_air/T_mean²
    dG_dT = dG_dRi * (-ri_scale * T_air / (T_mean * T_mean))

    k_sh = density_air * C * HEAT_CAPACITY_AIR
    heat_flux_sensible = k_sh * ΔT * inv_FF
    dshf = k_sh * inv_FF * (-1.0 - ΔT * dG_dT)

    if T_surface >= CtoK - T_TOLERANCE
        latent_heat = LV
        eS = saturation_vapor_pressure_water(T_surface)
        deS = saturation_vapor_pressure_water_derivative(T_surface)
    else
        latent_heat = LS
        eS = saturation_vapor_pressure_ice(T_surface)
        deS = saturation_vapor_pressure_ice_derivative(T_surface)
    end
    Δe = cfs.vapor_pressure - eS
    k_lh = C * latent_heat / 461.9
    heat_flux_latent = k_lh * Δe / T_mean * inv_FF
    dlhf = k_lh * inv_FF * (-deS / T_mean - Δe / (2 * T_mean * T_mean) - Δe / T_mean * dG_dT)

    return heat_flux_sensible, heat_flux_latent, latent_heat, dshf, dlhf
end
