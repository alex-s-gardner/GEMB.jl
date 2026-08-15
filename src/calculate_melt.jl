"""
    irreducible_saturation(mp::ModelParameters, density) -> S_wi [-]

Irreducible (capillary-held) water saturation of the pore space, per
`mp.water_irreducible_method`. The retention of a cell is `(ρᵢ − ρ) · S_wi · dz`
[kg m-2] at every site that uses it, so only `S_wi` varies between methods.

- `:constant` — `mp.water_irreducible_saturation` at every density (Colbeck, 1974).
  This is the MATLAB behaviour and the default.
- `:ColeouLesaffre` — Coléou and Lesaffre (1998) eq. 3 via Langen et al. (2017) eq. 4,
  where irreducible water mass per unit total mass is
  `wmi = 0.057(ρᵢ − ρ)/ρ + 0.017` and

      S_wi = wmi/(1 − wmi) · ρᵢ·ρ / (ρ_w(ρᵢ − ρ))

  `mp.water_irreducible_saturation` is ignored. Retention rises with density —
  ~0.069 at ρ = 300 against ~0.163 at ρ = 800 (ρᵢ = 917) — which is where the constant
  value under-retains most, in the percolation zone.

Both methods are gated at `DENSITY_PORE_CLOSEOFF`, where there is no connected pore
space left to hold water: `:ColeouLesaffre` returns zero there, as the Community Firn
Model and [`_irreducible_water`](@ref) do. `:constant` is deliberately *not* gated, so
the default path stays bit-identical to MATLAB. The gate also removes the `ρ → ρᵢ`
singularity in `S_wi`; the retention product `(ρᵢ − ρ)·S_wi` is finite in that limit,
but the saturation alone is not.

The gate is the *constant* `DENSITY_PORE_CLOSEOFF`, deliberately not the tunable
`mp.impermeable_density`. The two thresholds answer different questions: this one is where
capillary retention ceases for want of connected pore space, while `impermeable_density` is
where a lens stops conducting flow at the scale the model represents. Tying them together
would make lowering the flow criterion silently change retention too.
"""
@inline function irreducible_saturation(mp::ModelParameters, density::Float64)
    if mp.water_irreducible_method === :constant
        return mp.water_irreducible_saturation
    end
    # :ColeouLesaffre
    density >= DENSITY_PORE_CLOSEOFF - D_TOLERANCE && return 0.0
    ρi = mp.density_ice
    wmi = 0.057 * (ρi - density) / density + 0.017
    return wmi / (1 - wmi) * ρi * density / (DENSITY_WATER * (ρi - density))
end

"""
    calculate_melt(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity, age, rain, mp::ModelParameters, verbose::Bool)

Compute meltwater generation, percolation, refreezing, and runoff using a tipping bucket approach.

Processes:
1. Initial Refreeze: Existing pore water in cold layers is refrozen.
2. Melt Generation: Excess energy above 0 C is converted to liquid meltwater.
3. Percolation: Liquid water percolates downward, refreezing in cold layers,
   being retained as pore water, or running off at impermeable ice lenses.

Water is routed to runoff at a contiguous run of cells at or above `mp.impermeable_density`
thicker than `mp.impermeable_thickness` (defaults 830 kg m-3 and 0.1 m, as in MATLAB; across
the RetMIP models the density criterion spans 810–917 kg m-3 — see `ModelParameters`).

# Ponding above a barrier

What happens to water that reaches a barrier depends on `mp.runoff_method`:

- `:instantaneous` (the default, and MATLAB's behaviour) — it leaves the column in the same
  timestep. Every cell is capped at its irreducible saturation, so the column can never hold
  standing water.
- `:ZuoOerlemans` / `:Darcy` — it *backs up*. After the downward percolation loop, an upward
  pass fills the cells above each barrier to their pore capacity
  `mp.pore_saturation_max·(ρᵢ − ρ)·dz` [kg m-2], deepest first; only what is left after cell 1
  is full genuinely leaves as runoff, the column then being saturated to the surface. Water in
  excess of irreducible instead drains laterally, over a finite timescale, in
  [`apply_lateral_drainage!`](@ref). Barrier cells are excluded from the fill: at or above
  `mp.impermeable_density` there is no connected pore space to pond in, which is the premise
  of their being impermeable in the first place.

  The capacity uses the same `(ρᵢ − ρ)·S·dz` form as [`irreducible_saturation`](@ref) with
  `S = mp.pore_saturation_max`, so `water/capacity` is exactly the saturation `S_wi` is
  defined against and a cell at `S_wi = 1` reads as full rather than over-full.

  The upward pass does not refreeze what it deposits, even into a cell that is still cold
  (a cell can pass water down while retaining cold content if its refreeze was capped by
  `d_max` rather than by cold content). Any such cell is resolved at the top of the *next*
  timestep by the REFREEZE PORE WATER block above, so the only cost is a one-step delay; mass
  and energy are conserved either way, because pore water and runoff carry the same
  melt-point enthalpy.

  Firn aquifers then form bottom-up with no further machinery: water still drains to
  irreducible everywhere it passes, piles up at the deepest barrier, and backs up from there.
  RetMIP (Vandecrux et al., 2020, Sect. 5.4) excluded the firn-aquifer site from its retention
  evaluation for exactly this reason — models that cannot hold water above irreducible
  saturation cannot represent an aquifer at all.

  This deviates from MATLAB, which has only the `:instantaneous` behaviour.

# No preferential-flow domain

This is a single-domain bucket scheme: water descends cell by cell through the matrix, and
there is no second, fast, heterogeneous ("piping") flow path, nor a Richards-equation matrix
solver. That is deliberate. RetMIP (Vandecrux et al., 2020) found that the three of nine
models carrying an explicit deep- or preferential-percolation scheme (CFM-Cr, CFM-KM,
UppsalaUniDeepPerc) did worse than the bucket schemes at three of the four Greenland sites:
they infiltrated too deeply and ran warm — firn-temperature mean error +3.6 to +6.2 °C at
Dye-2 and +1.8 to +4.7 °C at KAN_U, plus a warm bias at the near-melt-free Summit site — and
at Dye-2 in 2016 they percolated to 10 m against the 2.5 m the upward-looking radar observed,
growing near-surface ice slabs several metres thick where none exist. Their one advantage was
the firn-aquifer site, which only they recharged. RetMIP's conclusion (their Sect. 5.2) is
that until preferential flow in firn is better constrained observationally, the more complex
schemes do not necessarily give better results than simple bucket schemes; the productive
levers for a bucket scheme are instead the impermeability criterion and the runoff timescale
above — the two best-performing models at the ice-slab site KAN_U were bucket schemes that
delay runoff rather than models that percolate deeper.

Returns `(temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
age, melt_total, melt_surface, runoff_total, freeze_total, percolation_depth)`.

# Age transport

`age` [d] is the mass-weighted mean age of all mass in a cell, matrix and pore water
together, so a phase change that stays inside one cell (the initial pore-water refreeze) does
not move it at all. Melt removal is likewise age-neutral: it takes a *fraction* of the cell,
which leaves the mean age of the remainder unchanged.

What does move age is percolation. Meltwater carries the age of the firn it came from, so the
percolation loop transports the age *moment* `mass × age` alongside `flux_dn`, in a companion
`flux_age` vector holding the mass-weighted mean age of the water in transit. At each cell the
water generated there and any water the cell expels (both at `age[i]`) mix with the through-flow
from above (at `flux_age[i]`) into one pool; whatever refreezes enters the cell at the pool's
age, and whatever continues down or runs off leaves at it. Refreezing therefore makes a cell
*younger* only to the extent that genuinely younger water reached it, rather than resetting it.

`percolation_depth` [m] is a diagnostic that takes no part in the mass or energy budget: the
base of the deepest cell water reached this timestep, 0 when no water moved. Note that the
percolation loop can walk past cells water never entered (its exit also waits on the deepest
cell holding pre-existing excess pore water), so this tracks water arrival per cell rather
than the loop index. The companion slab diagnostics are taken in
[`gemb_core`](@ref) via [`ice_slab_diagnostics`](@ref), after the grid controllers run.

Arrays may shrink (cells deleted when mass=0).
"""
function calculate_melt(temperature::Vector{Float64}, dz::Vector{Float64},
    density::Vector{Float64}, water::Vector{Float64},
    grain_radius::Vector{Float64}, grain_dendricity::Vector{Float64},
    grain_sphericity::Vector{Float64}, age::Vector{Float64}, rain::Float64,
    mp::ModelParameters, verbose::Bool)

    # Note: arrays are modified in-place. May shrink via deleteat! when cells lose all mass.

    T_tolerance = 1e-10
    d_tolerance = 1e-11
    water_tolerance = 1e-13

    # The density-based impermeability criterion. Both default to the MATLAB values
    # (`DENSITY_PORE_CLOSEOFF` = 830 and 0.1 m) — see `ModelParameters` for why they are
    # tunable and what range the literature supports.
    d_phc = mp.impermeable_density       # pore hole close off density [kg m-3]
    ice_layer_dzmin = mp.impermeable_thickness   # minimum ice layer thickness for runoff [m]

    m = length(temperature)
    water_delta = zeros(m)

    # store initial mass [kg]
    M = dz .* density

    if verbose
        M_total_initial = sum(water) + sum(M)
        E_total_initial = column_enthalpy(mp, M, temperature, water)
    end

    # initialize melt and runoff scalars
    runoff_total = 0.0
    melt_total = 0.0
    melt_surface = 0.0

    # Diagnostics (read-only; they take no part in the mass or energy budget). Initialized
    # here because the percolation block below is conditional — a timestep with no water
    # anywhere skips it entirely and must still report a wetting front of zero depth.
    percolation_depth = 0.0

    # calculate temperature excess above 0 degC
    T_excess = max.(0.0, temperature .- CtoK)

    # new grid point center temperature. Rebind to a fresh array (do not mutate
    # the caller's temperature vector) to preserve the previous behavior.
    temperature = min.(temperature, CtoK)

    ## REFREEZE PORE WATER

    if sum(water) > water_tolerance
        # Fused per-cell refreeze: freeze pore water and update snow/ice
        # properties. Numerically identical (same per-element ops and
        # left-to-right ordering) to the previous broadcast chain, but computes
        # everything in a single pass. density, dz, and water are rebound to
        # fresh arrays (matching the previous `density = M ./ dz` etc.) so the
        # caller's inputs are not mutated; M and temperature are mutated in
        # place since they are already function-local fresh arrays.
        #
        # `age` needs no update here: this moves mass from `water[i]` to `M[i]` within one
        # cell, and `age` is defined on their sum, so the cell's total mass and hence its
        # mean age are both unchanged.
        dz_orig = dz
        density = similar(density)
        dz = similar(dz_orig)
        water = copy(water)
        @inbounds for i in 1:m
            # maximum freeze amount [kg]
            freeze_max_i = cold_content_mass(mp, temperature[i], M[i])

            # freeze pore water and change snow/ice properties
            wd = min(freeze_max_i, water[i])
            water_delta[i] = wd
            water[i] = water[i] - wd
            M[i] = M[i] + wd
            density[i] = M[i] / dz_orig[i]
            if M[i] > water_tolerance
                temperature[i] = refreeze_temperature(mp, temperature[i], M[i], wd)
            end

            # if pore water froze in ice then adjust density and dz thickness
            if density[i] > mp.density_ice - d_tolerance
                density[i] = mp.density_ice
            end
            dz[i] = M[i] / density[i]
        end
    end

    # squeeze water from snow pack (compute water_excess without materializing
    # the water_irreducible temporary)
    water_excess = Vector{Float64}(undef, m)
    @inbounds for i in 1:m
        water_irreducible = (mp.density_ice - density[i]) * irreducible_saturation(mp, density[i]) * (M[i] / density[i])
        water_excess[i] = max(0.0, water[i] - water_irreducible)
    end

    ## MELT, PERCOLATION AND REFREEZE

    # Seed freeze with the pore-water refreeze accumulated above, then reset
    # water_delta for reuse in the percolation loop.
    freeze = copy(water_delta)
    fill!(water_delta, 0.0)

    # run melt algorithm if there is melt water or excess pore water
    if (sum(T_excess) > T_tolerance) || (sum(water_excess) > water_tolerance)

        # Check to see if thermal energy exceeds energy to melt entire cell.
        # `T_excess`/`T_surplus` stay in kelvin: the branch tests below compare against
        # `T_tolerance`, a kelvin tolerance, so moving them to joules would rescale those
        # thresholds by ~c_p. The surplus *energy* is computed from them per cell.
        T_full_melt = full_melt_excess_temperature(mp)
        T_surplus = max.(0.0, T_excess .- T_full_melt)

        if sum(T_surplus) > T_tolerance
            # calculate surplus energy. Built into a concretely-typed vector rather than a
            # comprehension: the comprehension infers as `Vector{Any}`, which then makes
            # `T_surplus` itself infer as `Any` through the shared loop below and boxes every
            # arithmetic op in the hot path (18x the allocations in `calculate_melt`).
            E_surplus = Vector{Float64}(undef, m)
            @inbounds for i in 1:m
                E_surplus[i] = surplus_energy(mp, T_excess[i], T_surplus[i], M[i])
            end
            i = 1

            while (sum(E_surplus) > T_tolerance) && (i < (m + 1))
                if i < m
                    # use surplus energy to increase the temperature of lower cell
                    temperature[i+1] = add_energy_temperature(mp, temperature[i+1], M[i+1],
                        E_surplus[i])

                    T_excess[i+1] = max(0.0, temperature[i+1] - CtoK) + T_excess[i+1]
                    temperature[i+1] = min(CtoK, temperature[i+1])

                    T_surplus[i+1] = max(0.0, T_excess[i+1] - T_full_melt)
                    E_surplus[i+1] = surplus_energy(mp, T_excess[i+1], T_surplus[i+1], M[i+1])
                else
                    error("surplus energy reached the base of gemb column (i.e. entire column melted out in a single time step)")
                end

                # adjust current cell properties: the cell keeps exactly the excess that
                # melts it out, having handed the rest downward.
                T_excess[i] = T_full_melt
                E_surplus[i] = 0.0
                i += 1
            end
        end

        # Convert temperature excess to melt [kg] and compute the max refreeze amount in a
        # single fused pass. The melt mass is the excess enthalpy divided by the latent heat;
        # `freeze_max` is the cold content, the same quantity with the sign reversed. Also
        # tracks the running melt sum and the deepest cell with melt/excess pore water,
        # avoiding the melt_maximum, freeze_max, and findlast BitVector temporaries.
        melt = Vector{Float64}(undef, m)
        freeze_max = Vector{Float64}(undef, m)
        melt_sum = 0.0
        X = 1
        @inbounds for i in 1:m
            melt_max_i = melt_mass_from_excess(mp, T_excess[i], density[i], dz[i])
            mi = min(melt_max_i, M[i])
            melt[i] = mi
            melt_sum += mi
            freeze_max[i] = cold_content_mass(mp, temperature[i], density[i], dz[i])
            if mi > water_tolerance || water_excess[i] > water_tolerance
                X = i
            end
        end
        melt_surface = melt[1]
        melt_total = max(0.0, melt_sum - rain)

        # initialize refreeze, runoff, flux_dn and water_delta vectors
        runoff = zeros(m)
        flux_dn = zeros(m + 1)
        # Mass-weighted mean age [d] of the water in `flux_dn`, same indexing. Only entries
        # with `flux_dn > 0` are ever read, so the zero fill needs no sentinel.
        flux_age = zeros(m + 1)

        Xi = 1
        m = length(temperature)

        # Deepest cell water actually entered, for `percolation_depth`. Tracked separately
        # from `Xi` (the cell the loop stopped at) because the loop's break condition also
        # waits for `i > X`, the deepest cell holding *pre-existing* excess pore water: a
        # single wet cell at depth keeps the loop running through every cell above it, none
        # of which the water ever crossed. 0 when no cell receives water.
        i_wet = 0

        # meltwater percolation
        for i in 1:m
            # calculate total melt water entering cell
            melt_input = melt[i] + flux_dn[i]
            if abs(melt_input) > water_tolerance
                i_wet = i
            end

            # Age moment of the cell on entry. `M[i]` still includes `melt[i]` (the branches
            # below subtract it) and `water_delta[i]` has not been applied, so this is the
            # cell's full mass before any of this timestep's transfers.
            A_entry = (M[i] + water[i]) * age[i]

            ice_depth = 0.0
            # If this grid cell's density exceeds the pore closeoff density:
            if density[i] >= d_phc - d_tolerance
                for l in i:m
                    if density[l] >= d_phc - d_tolerance
                        ice_depth += dz[l]
                        if ice_depth > ice_layer_dzmin + d_tolerance
                            break
                        end
                    else
                        break
                    end
                end
            end

            # break loop if there is no meltwater and if depth is > mw_depth
            if abs(melt_input) < water_tolerance && i > X
                break

            # if reaches impermeable ice layer all liquid water runs off
            elseif (density[i] >= (mp.density_ice - d_tolerance)) ||
                   ((density[i] >= d_phc - d_tolerance) && (ice_depth > ice_layer_dzmin + d_tolerance))

                M[i] = M[i] - melt[i]
                water_irr = (mp.density_ice - density[i]) * irreducible_saturation(mp, density[i]) * (M[i] / density[i])
                water_delta[i] = max(min(melt_input, water_irr - water[i]), -water[i])
                runoff[i] = max(0.0, melt_input - water_delta[i])

            # check if no energy to refreeze meltwater
            elseif abs(freeze_max[i]) < d_tolerance

                M[i] = M[i] - melt[i]
                water_irr = (mp.density_ice - density[i]) * irreducible_saturation(mp, density[i]) * (M[i] / density[i])
                water_delta[i] = max(min(melt_input, water_irr - water[i]), -1 * water[i])
                flux_dn[i+1] = max(0.0, melt_input - water_delta[i])
                runoff[i] = 0.0

            # some or all meltwater refreezes
            else
                M[i] = M[i] - melt[i]
                dz_0 = M[i] / density[i]
                d_max = (mp.density_ice - density[i]) * dz_0
                freeze1 = min(min(melt_input, d_max), freeze_max[i])
                M[i] = M[i] + freeze1
                density[i] = M[i] / dz_0

                # pore water. `density[i]` is the post-refreeze value, so the saturation
                # is evaluated on it too under `:ColeouLesaffre`.
                water_irr = (mp.density_ice - density[i]) * irreducible_saturation(mp, density[i]) * dz_0
                water_delta[i] = max(min(melt_input - freeze1, water_irr - water[i]), -1 * water[i])
                freeze2 = 0.0

                if water_delta[i] < 0.0 - water_tolerance
                    d_max = (mp.density_ice - density[i]) * dz_0
                    freeze2_max = min(d_max, freeze_max[i] - freeze1)
                    freeze2 = min(-1.0 * water_delta[i], freeze2_max)
                    M[i] = M[i] + freeze2
                    density[i] = M[i] / dz_0
                end

                freeze[i] = freeze[i] + freeze1 + freeze2
                flux_dn[i+1] = max(0.0, melt_input - freeze1 - water_delta[i])

                if M[i] > water_tolerance
                    temperature[i] = refreeze_temperature(mp, temperature[i], M[i],
                        freeze1 + freeze2)
                end

                # check if an ice layer forms
                if abs(density[i] - mp.density_ice) < d_tolerance
                    runoff[i] = flux_dn[i+1]
                    flux_dn[i+1] = 0.0
                end
            end

            # --- age transport ---------------------------------------------------------
            # Applied once for all three branches above, which have by now set
            # `water_delta[i]`, `runoff[i]` and `flux_dn[i+1]` and adjusted `M[i]`.
            #
            # The mobile pool at this cell is the water generated here plus any the cell is
            # expelling (`water_delta < 0`), both at `age[i]`, mixed with the through-flow
            # from above at `flux_age[i]`. Everything leaving the cell — downward flux and
            # runoff — leaves at the pool's age; everything retained or refrozen enters at it.
            own = melt[i] + max(0.0, -water_delta[i])
            A_in = flux_dn[i] * flux_age[i]
            pool_age = mix_age(age[i], own, flux_age[i], flux_dn[i])
            flux_age[i+1] = pool_age

            # Post-transfer cell mass, read back from the arrays the branches just wrote
            # rather than derived from an assumed balance, so the mean stays bounded.
            S_new = M[i] + water[i] + water_delta[i]
            A_new = A_entry + A_in - (runoff[i] + flux_dn[i+1]) * pool_age
            age[i] = S_new > water_tolerance ? max(0.0, A_new / S_new) : 0.0

            Xi = Xi + 1
        end

        # Check for negative pore water
        if verbose
            if any(water .< 0.0 - water_tolerance)
                error("Negative pore water generated in melt equations.")
            end
        end

        # Let blocked water back up into the pore space above the barrier instead of leaving
        # in this timestep, unless `runoff_method` is `:instantaneous` (the MATLAB default,
        # under which this returns the blocked mass unchanged and touches nothing).
        runoff_blocked = pond_blocked_water!(M, density, water, water_delta, age,
            runoff, flux_dn, flux_age, Xi, mp)

        # adjust pore water
        water = water .+ water_delta

        # calculate runoff_total
        runoff_total = sum(runoff) + runoff_blocked

        # Wetting-front depth [m], for comparison against upward-looking radar (e.g. Heilig et
        # al., 2018, as used by RetMIP): the base of the deepest cell water reached. Taken
        # before cells are deleted and `dz` rebuilt, so the indices still line up with the
        # column the loop walked.
        @inbounds for i in 1:i_wet
            percolation_depth += dz[i]
        end

        # Delete all cells that melted out entirely. `dz` is rebuilt from the conserved
        # cell mass immediately below, so it is excluded from the shift here and the
        # column state is assembled with `M` standing in for it.
        to_delete = findall(M .<= water_tolerance)
        if !isempty(to_delete)
            close_slot!(column_state(temperature, M, density, water, grain_radius,
                grain_dendricity, grain_sphericity, age), to_delete)
        end

        # calculate new grid lengths
        dz = M ./ density
    end

    freeze_total = sum(freeze)

    ## CHECK FOR MASS AND ENERGY CONSERVATION
    if verbose
        E_total_runoff = runoff_total * specific_enthalpy_water(mp)

        M_total_final = sum(water) + sum(M) + runoff_total
        E_total_final = column_enthalpy(mp, M, temperature, water)

        M_delta = M_total_initial - M_total_final
        E_delta = E_total_initial - E_total_final - E_total_runoff

        E_tol = energy_tolerance(E_total_initial)

        if (abs(M_delta) > 1e-3) || (abs(E_delta) > E_tol)
            error("Mass and/or energy are not conserved in melt equations:\n M_delta: $(M_delta) E_delta: $(E_delta)\n")
        end

        if any(water .< 0.0 - water_tolerance)
            error("Negative pore water generated in melt equations.")
        end
    end

    return temperature, dz, density, water, grain_radius, grain_dendricity, grain_sphericity,
    age, melt_total, melt_surface, runoff_total, freeze_total, percolation_depth
end

"""
    pore_capacity(mp, density, dz) -> capacity [kg m-2]

Maximum water a cell can hold, `mp.pore_saturation_max·(ρᵢ − ρ)·dz`.

Deliberately the same `(ρᵢ − ρ)·S·dz` form as the retention in
[`irreducible_saturation`](@ref), with `S = mp.pore_saturation_max` in place of `S_wi`, so
that `water/capacity` is exactly the saturation `S_wi` is defined against: a cell holding its
irreducible water sits at `S_wi/pore_saturation_max` of capacity, and one at `S_wi = 1` reads
as full rather than over-full. The physically exact pore mass `ρ_w·dz·(1 − ρ/ρᵢ)` would
disagree with the retention expression unless `ρ_w = ρᵢ`, and mixing the two conventions
would let a cell exceed capacity while still under saturation 1.

Returns 0 for a cell at or above `mp.density_ice`, which has no pore space at all.
"""
@inline pore_capacity(mp::ModelParameters, density::Float64, dz::Float64) =
    max(0.0, mp.pore_saturation_max * (mp.density_ice - density) * dz)

"""
    pond_blocked_water!(M, density, water, water_delta, age, runoff, flux_dn, flux_age, Xi, mp)
        -> runoff_blocked [kg m-2]

Back water that percolation could not pass — the runoff each cell generated, plus the basal
outflow `flux_dn[Xi]` — up into the pore space above, and return only the mass that still
leaves the column.

Called from [`calculate_melt`](@ref) after the percolation loop and before `water_delta` is
applied, so it works on `water .+ water_delta` as the current pore water. A no-op under
`mp.runoff_method === :instantaneous`, where it returns `flux_dn[Xi]` with `runoff`
untouched, leaving the total identical to MATLAB's `sum(runoff) + flux_dn[Xi]`.

Otherwise the blocked water is pooled (mass-weighted by [`mix_age`](@ref), since each
contribution leaves its cell at that cell's outflow age `flux_age[i+1]`), the contributing
`runoff` entries are zeroed, and the pool is spent filling cells to
[`pore_capacity`](@ref) from the deepest contributing cell upward. Cells at or above
`mp.impermeable_density` are skipped: at pore close-off the pore space is disconnected, which
is the premise of their blocking flow. Whatever the column cannot hold is returned as runoff —
at that point every cell up to the surface is at capacity, so it genuinely leaves.

Mass is conserved by construction (every kilogram is either placed in a cell or returned),
and so is energy: pore water and runoff carry the same melt-point enthalpy
[`specific_enthalpy_water`](@ref), so moving mass between them changes neither budget term's
per-kilogram value. Nothing refreezes here; see the ponding section of
[`calculate_melt`](@ref) on why the next timestep handles that.
"""
function pond_blocked_water!(M::Vector{Float64}, density::Vector{Float64},
    water::Vector{Float64}, water_delta::Vector{Float64}, age::Vector{Float64},
    runoff::Vector{Float64}, flux_dn::Vector{Float64}, flux_age::Vector{Float64},
    Xi::Int, mp::ModelParameters)

    # `flux_dn[Xi]` is basal outflow only when the loop ran to completion; when it broke
    # early `flux_dn[Xi]` is the sub-tolerance inflow to the cell it stopped at, which is
    # why this is summed rather than treated as a separate case.
    mp.runoff_method === :instantaneous && return flux_dn[Xi]

    water_tolerance = 1e-13
    d_tolerance = 1e-11
    d_phc = mp.impermeable_density

    # Pool the blocked water. Each contribution left its cell at that cell's outflow age,
    # `flux_age[i+1]`; the basal term left the last cell visited, at `flux_age[Xi]`.
    blocked = 0.0
    blocked_age = 0.0
    i_start = 0
    @inbounds for i in eachindex(runoff)
        if runoff[i] > water_tolerance
            blocked_age = mix_age(blocked_age, blocked, flux_age[i+1], runoff[i])
            blocked += runoff[i]
            runoff[i] = 0.0
            i_start = i
        end
    end
    if flux_dn[Xi] > water_tolerance
        blocked_age = mix_age(blocked_age, blocked, flux_age[Xi], flux_dn[Xi])
        blocked += flux_dn[Xi]
        i_start = max(i_start, Xi - 1)
    end

    blocked <= water_tolerance && return blocked
    i_start = min(i_start, length(M))

    # Fill upward from the deepest contributing cell. Water blocked at a shallow barrier never
    # reached a deeper one, so a single upward sweep from the deepest source covers every
    # barrier that received water this timestep.
    @inbounds for i in i_start:-1:1
        blocked <= water_tolerance && break
        density[i] >= d_phc - d_tolerance && continue

        dz_i = M[i] / density[i]
        room = pore_capacity(mp, density[i], dz_i) - (water[i] + water_delta[i])
        room <= water_tolerance && continue

        add = min(blocked, room)
        # Mass already in the cell, matrix and pore water together, which is what `age` is
        # the mean over.
        S = M[i] + water[i] + water_delta[i]
        age[i] = mix_age(age[i], S, blocked_age, add)
        water_delta[i] += add
        blocked -= add
    end

    return blocked
end

"""
    ice_slab_diagnostics(dz, density, mp::ModelParameters) -> (thickness, depth)

Summarize the ice slabs in a column under the same criterion the percolation scheme uses.

`thickness` [m] is the total thickness of all cells at or above `mp.impermeable_density`,
whether or not they form a flow-blocking run — the quantity a firn core measures.

`depth` [m] is the depth to the top of the shallowest run that would actually block
percolation, matching both clauses of the impermeable branch of [`calculate_melt`](@ref): a
contiguous run thicker than `mp.impermeable_thickness`, or any single cell at
`mp.density_ice`, which blocks unconditionally however thin it is. It is `NaN` when nothing
qualifies, so "no slab" stays distinguishable from "slab at the surface".

Both thresholds are capped at `mp.density_ice`, since `mp.impermeable_density` may be set
above it (see [`ice_slab_diagnostics`](@ref) source) and solid ice must always count.

Diagnostic only: nothing in the physics reads either value. Called from
[`gemb_core`](@ref) after the grid controllers, not from [`calculate_melt`](@ref), so the
column scanned is the one the profile output records at that timestamp.
"""
function ice_slab_diagnostics(dz::Vector{Float64}, density::Vector{Float64},
    mp::ModelParameters)

    d_tolerance = 1e-11
    # Solid ice is dense enough to count whatever `impermeable_density` is set to: the
    # criterion may legitimately exceed `density_ice` (the validator allows up to 917 so that
    # lowering `density_ice` does not reject the default), and in that configuration every
    # pure-ice cell still blocks flow via the unconditional `density_ice` clause in
    # `calculate_melt`. Taking the smaller of the two keeps this scan and that branch in
    # agreement instead of reporting "no slab" for a column of solid ice.
    d_phc = min(mp.impermeable_density, mp.density_ice)
    dzmin = mp.impermeable_thickness

    thickness = 0.0
    depth = NaN
    z_top = 0.0      # depth to the top of the current run
    run_dz = 0.0     # thickness of the current run
    z = 0.0          # depth to the top of cell i

    @inbounds for i in eachindex(dz)
        if density[i] >= d_phc - d_tolerance
            thickness += dz[i]
            if run_dz == 0.0
                z_top = z
            end
            run_dz += dz[i]
            # First qualifying run wins; keep scanning for `thickness` but never overwrite.
            # A cell at `density_ice` blocks on its own, however thin, matching the
            # unconditional clause in `calculate_melt`.
            if isnan(depth) && (run_dz > dzmin + d_tolerance ||
                                density[i] >= mp.density_ice - d_tolerance)
                depth = z_top
            end
        else
            run_dz = 0.0
        end
        z += dz[i]
    end

    return thickness, depth
end

"""
    aquifer_diagnostics(dz, density, water, mp) -> (thickness, depth)

Summarize the standing (super-irreducible) water in a column.

`thickness` [m] is the total thickness of cells holding more water than capillary forces can
retain — the saturated thickness a borehole or a firn core would find wet. `depth` [m] is the
depth to the top of the shallowest such cell, `NaN` when there are none, so "dry column" stays
distinguishable from "water table at the surface".

Both are zero and `NaN` under `mp.runoff_method === :instantaneous`, where no cell can hold
more than its irreducible water. They are what makes the ponding in [`calculate_melt`](@ref)
observable in output, and the quantity RetMIP (Vandecrux et al., 2020, Sect. 5.4) evaluates
aquifer representation against.

A cell counts as saturated only once it exceeds irreducible by [`AQUIFER_TOLERANCE`](@ref) of
its pore space, rather than by the arithmetic `WATER_TOLERANCE` the physics uses. `calculate_melt`
leaves a retaining cell at exactly irreducible, but `calculate_density` then compacts it within
the same timestep, shrinking the pore space that retention was computed against and leaving the
cell genuinely — if negligibly — over-saturated. That residue accumulates between melt events
and would otherwise read as a metres-thick water table in a run that cannot have one. It is
separable from real ponding by saturation but not by mass: the residue sits ~1e-7 of pore space
above irreducible against 1e-1 or more for ponded water.

Diagnostic only: nothing in the physics reads either value. Called from [`gemb_core`](@ref)
after the grid controllers, alongside [`ice_slab_diagnostics`](@ref), so the column scanned is
the one the profile output records at that timestamp.
"""
function aquifer_diagnostics(dz::Vector{Float64}, density::Vector{Float64},
    water::Vector{Float64}, mp::ModelParameters)

    thickness = 0.0
    depth = NaN
    z = 0.0

    @inbounds for i in eachindex(dz)
        pore = (mp.density_ice - density[i]) * dz[i]
        irreducible = pore * irreducible_saturation(mp, density[i])
        if water[i] > irreducible + AQUIFER_TOLERANCE * pore
            thickness += dz[i]
            isnan(depth) && (depth = z)
        end
        z += dz[i]
    end

    return thickness, depth
end
