# Thermal solvers

```@meta
CurrentModule = GEMB
```

GEMB carries two schemes for advancing the subsurface temperature profile over one forcing
timestep, selected by the type of `mp.thermal_solver` rather than by a symbol (see
[`AbstractThermalSolver`](@ref)):

| | [`ExplicitThermal`](@ref) | [`ImplicitThermal`](@ref) |
|---|---|---|
| Scheme | Explicit finite volume | Backward Euler, Thomas algorithm |
| Stability | Conditional — Von Neumann limit | Unconditional |
| Sub-steps set by | The stiffest cell in the column | [`THERMAL_IMPLICIT_DT_TARGET`](@ref) alone |
| Sub-steps buy | Stability | Accuracy |
| Order in time | First | First |
| Default | Yes | |

`ExplicitThermal` is the default and the only scheme the synthetic regression fingerprints are
pinned to. This page is the measured record behind that choice: which scheme to use when, what
the implicit path costs, and which optimizations of it were tried and rejected. The schemes'
algebra and invariants are documented on the two `GEMB._thermal_solve!` methods.

## Which to use

Use the default `ExplicitThermal` unless the column develops a cell thin enough to dominate the
stability limit. Because that limit scales with the *smallest* cell, one thin layer sets the cost
of the whole column: a single 1e-4 m refrozen lens drops the limit from 903 s to 3.98 s on the
264-cell benchmark column — a 227x sub-step increase for one cell (measured in
`test/test_calculate_temperature.jl`). The implicit count does not move at all.

That is the case `ImplicitThermal` exists for, and the only one in which it is also *faster*. On a
column with no stiff cell it is 1.5-2.4x slower (below), so unconditional stability is the reason
to select it and throughput is not.

## Measured runtime

A year of 3-hourly synthetic forcing over a 264-cell column, whole-model runtime:

| Solver | `DT_TARGET` | Runtime |
|---|---|---|
| `ExplicitThermal` | — | 2.44 s |
| `ImplicitThermal` | 1800 s (6 sub-steps) | 3.59 s |
| `ImplicitThermal` | 900 s (12 sub-steps, default) | 5.75 s |
| `ImplicitThermal` | 450 s (24 sub-steps) | 9.62 s |

So 1.5x at the coarsest usable accuracy and 2.4x at the calibrated default.

### Static condensation: 2.67x

The implicit path's nonlinearity is confined to row 1, the surface energy balance; under the
default `:constant` heat capacity rows `2 … n` are linear in the iterate. Re-eliminating them once
per Newton iteration was wasted work. Condensing the interior **bottom-up, once per sub-step**
reduces each row to an affine function of the cell above it, so Newton iterates on a single scalar
equation in `T_1` at O(1) per iteration and one forward substitution propagates the answer back
down the column.

Measured: **15.20 s to 5.70 s (2.67x)**, with whole-model output agreeing to 1.2e-6 K and annual
melt to 3.2e-8 kg m-2 — the same converged answer to round-off, as the `Λ`-independence property
of the scheme requires. This is what the exponential-integrator idea in the design notes reduces
to in cheap form: the operator is constant across all Newton iterations of a sub-step, so it is
factorized once and reused.

### What remains

The Newton iteration itself: **4.19 iterations per sub-step solve**, over 1.12M solves. Capping the
iteration at 1 takes the run to 4.98 s, so the iteration is most of the remaining cost. It is also
irreducible — surface-only and whole-profile convergence criteria give identical counts (4.190 vs
4.190), i.e. cell 1 is always the last cell to converge.

## Rejected optimizations

Recorded so they are not retried.

**Analytic turbulent-flux derivatives.** The `Λ` slope entering the Newton diagonal takes
[`THERMAL_BC_DERIVATIVE_STEP`](@ref)-sized finite differences through the Beljaars-Holtslag
stability branches. Replacing them analytically requires holding the transfer coefficients
`coefM`/`coefHT`/`coefHQ` fixed, and they cannot be: they depend on `T_surface` through the bulk
Richardson number, and the term dropped by freezing them substantially cancels the rest. Over 675
sampled surface states the frozen-coefficient slope overshot the true derivative by up to 10x,
with 9% of states outside a factor of two. Overshoot under-relaxes Newton: iteration count more
than doubled, hit `THERMAL_IMPLICIT_MAX_ITERATIONS`, and whole-model runtime *rose* 43% (15.2 s to
21.8 s) while the unconverged output drifted 15 K. Differentiating the stability branches too
would be sound, and remains open.

**The turbulent-flux calls are not the cost.** Direct timing puts `_surface_energy_balance` at
65 ns and `_surface_energy_balance_slope` at 69 ns, so all flux evaluations together were 0.70 s of
the pre-condensation 15.20 s. An earlier attribution of the cost to doubled flux calls was wrong;
the cost was the O(n) sweep per iteration, which condensation removed.

**Relaxing the convergence tolerance.** Taking `THERMAL_IMPLICIT_T_TOLERANCE` from 1e-10 to
1e-3 — seven orders of magnitude, far looser than the model's other branch tolerances — cut the
pre-condensation runtime only from 15.2 s to 10.8 s.

**Sub-step halving for non-convergence.** The remedy Fourteau et al. (2026) Sect. 3.1.3 pairs with
their own Newton solve. Rejected here in favour of step damping; see
`THERMAL_IMPLICIT_DAMPING_FLOOR` for the comparison.

## Convergence

Linearizing the surface flux into the diagonal is not sufficient on its own — a fraction of solves
end in small limit cycles across the melting-point switches rather than converging, which is why
the iteration is damped. That measurement (39% non-convergence before damping, 0.9% after, and why
sub-step halving was rejected in favour of it) is recorded on
`GEMB.THERMAL_IMPLICIT_DAMPING_FLOOR`, the constant it justifies.

## References

- Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing.
- Versteeg, H. K. & Malalasekera, W. (2007). *An Introduction to Computational Fluid Dynamics:
  The Finite Volume Method*, 2nd ed. Pearson, Ch. 8.
- Thomas, L. H. (1949). *Elliptic Problems in Linear Difference Equations over a Network*.
  Watson Scientific Computing Laboratory, Columbia University.
- Fourteau, K., Brondex, J., Cancès, C., and Dumont, M. (2026). Numerical strategies for
  representing Richards' equation and its couplings in snowpack models. *Geosci. Model Dev.*
  19, 3193-3212.
