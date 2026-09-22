# Experimental direct SVD fit for 2D Algoim surface moments

This describes the current `Mesh2` triangle implementation in
`cpp/problem/AlgoimCutFEM.tpp`, controlled by `ProblemOption` in
`cpp/solver/solver.hpp`. It is an implementation description, not an
accuracy theorem. The direct SVD fit is opt-in through
`algoim_ibp_surface_svd_`; the default is the older Gram-matrix fit.
The SVD option compares both candidates and may return the legacy one.
The volume-weight correction still uses the older solver.

## Why a surface fit is needed

On a cut triangle `K`, the volume rule integrates `{phi_B < 0}` and the
surface rule integrates its implicit boundary. For a polynomial vector
field `F`, the two rules and the element-edge rule should obey

```text
sum_volume w_q div F(x_q)
  = sum_surface w_q F(y_q)·n_q
    + sum_inside_edges v_q·F(z_q).
```

The edge rule integrates the portions of the three triangle edges inside
`{phi_B < 0}`. It uses the same Bernstein level set `phi_B` as the cut rules.
The surface stage first enforces the identity for divergence-free fields,
where the volume term vanishes. The volume stage then fits the remaining
scalar moments against that corrected boundary rule. The Stokes driver
currently requests degree `D = k+m` to cover products of a potential in
`P_m` and a velocity test in `[P_k]^2`; this choice does not ensure accurate
integration of every curved-interface term.

## Surface constraints

Algoim supplies fixed surface nodes `y_q`, positive scalar weights `w_q`,
and unit normals `n_q`, for `q=1,...,N`. The unknowns are the **two
components** of each vector weight

```text
omega_q = w_q n_q,       omega in R^(2N).
```

The fit does not move nodes. It may change both the lengths and directions
of vector weights. For every monomial
`psi_ab = xi^a eta^b`, with `1 <= a+b <= D+1`, the code constructs

```text
xi = (x-c_x)/h_K,   eta = (y-c_y)/h_K,
F_ab = (b xi^a eta^(b-1), -a xi^(a-1) eta^b).
```

These fields are divergence-free and span the intended polynomial class.
There are `(D+1)(D+4)/2` scalar constraints. The missing physical `1/h_K`
factor in the curl would multiply both sides of each corresponding identity
by the same constant, so it does not change that identity. Define
`C_i,2q = F_i(y_q)_x`, `C_i,2q+1 = F_i(y_q)_y`, and let `g_i` be the sum of
`F_i·(n ds)` over the inside-edge rule. The target is

```text
C omega + g = 0.
```

The implementation sets `rho_i = ||C_i,:||_2`, stores `A_i,: = C_i,:/rho_i`
when `rho_i > 0`, and evaluates the row-scaled residual
`r_i(omega) = -(C_i omega + g_i)/rho_i`. Acceptance uses
`max_i |r_i|`. A zero surface row is currently assigned zero residual;
the edge target for such a row is not independently checked. That is an
important degeneracy for review.

## Direct SVD update

For a raw-rule surface mass `W = sum_q w_q`, the code assigns both components
at node `q` the scale

```text
s_q = max(|w_q|, 0.1 W/N, h_K*1e-14),
S = diag(s_1,s_1,...,s_N,s_N).
```

It seeks `delta_omega = S z`, so the norm of `z` measures relative changes
to the quadrature weights with a floor at very small nodes. Starting from
`A delta_omega = r`, it forms `A S`. Each row of `A S` and the matching
right-hand side are divided by that row's Euclidean norm. Call the resulting
matrix `B` and right-hand side `b`. `LAPACKE_dgesvd` computes the thin
factorization `B = U diag(sigma) V^T` directly, without forming `B B^T`.

Only singular values satisfying

```text
sigma_j > sigma_max * max(1e-12, 64*epsilon*max(N_constraints,2N))
```

are retained. On that truncated subspace, the update is

```text
z = sum_retained v_j (u_j^T b)/sigma_j,
delta_omega = S z.
```

This is the minimum `||S^-1 delta_omega||_2` update for the retained
constraints. Discarded components of `b` are not enforced. The diagnostic
`ibp_svd_discarded_rhs` reports their relative Euclidean norm; a value near
one can arise when the attempted right-hand side is already at roundoff,
so it is not by itself a failed-rule certificate. Direct SVD avoids the
squared condition number of the older `A A^T` eigensolve. The older
`1e-12` eigenvalue cutoff corresponds to a relative singular-value cutoff
around `1e-6` on `A`; the new cutoff admits smaller singular directions,
subject to the update bound and final moment check.

The SVD update is rejected if its Euclidean norm is nonfinite or exceeds
`0.25*(W + h_K*1e-8)`. Up to four correction passes are attempted. A pass is
committed only when it strictly decreases the maximum row-scaled residual.
Failure of the factorization, the correction cap, or monotonicity leaves the
last accepted vector weights in place.

## Packing, acceptance, and candidate choice

After fitting, each vector weight becomes a nonnegative scalar weight
`||omega_q||_2` and its unit direction. For the SVD candidate, a vector
shorter than `1e-14*h_K` is set to zero weight with its previous normal;
the legacy candidate treats such an entry as invalid. The code then
recomputes every moment residual from the **actual packed weights and
normals**. It accepts a candidate only if its weights are valid and

```text
max_i |r_i| <= 100 * 3e-16 * (1 + W).
```

With `algoim_ibp_surface_svd_ = true`, the code copies the same input rule,
fits one copy by SVD and the other by the legacy solver, and selects in this
order:

1. Keep the legacy rule if it passes.
2. Otherwise keep the SVD rule if it passes.
3. If both fail, keep the candidate with smaller measured final residual.

The selection checks polynomial moments. It does not compare general
surface integrals or penalize changes in the normals. Both solvers can rotate
or reverse a normal; a small moment residual does not certify geometric
fidelity. The whole-rule diagnostics include the relative vector-weight
change, relative surface-mass change, and minimum old/new normal alignment.
The current cap bounds each attempted update's **aggregate norm over all
nodes**, not each node's rotation or weight change.

## Subdivision, volume correction, and strict gate

If the chosen surface rule fails, `quadGenSurf` can subdivide the triangle
into four midpoint children, recursively up to
`algoim_subdivision_depth_`. It concatenates the child rules and fits them
again against the **original parent** constraints. Children alone cannot
certify the parent-cell divergence theorem. `quadGenVol` regenerates its
matching surface rule, fits volume weights through degree `D-1` using its
existing Gram-matrix solver, and likewise refits after subdivision.

The optional `algoim_ibp_require_valid_parent_` throws when the final
surface or volume parent fit fails. It is disabled in recursive child calls
so that a parent can repair failed children. The gate is off by default;
without it, a failed final rule is returned to assembly with
`ibp_correction_ok = false`. The gate only checks the measured moment
tolerance. It does not reject normal reversals or certify quadrature
accuracy for a nonpolynomial curved-interface integrand.

## Diagnostic and generalization limits

`AlgoimQuadratureRule` carries counts of subdivision, capped failures,
SVD attempts/selections, final correction status, rank, discarded right-hand
side, update norm, and change measures. Some SVD quantities describe the
selected candidate's last update, rather than a complete history of all
attempts and children. For reproducible per-cell comparisons, use the
workfiles `algoim_refinement_survey` and an independent flux or form test.

The implementation is restricted to the `Mesh2` triangle surface path.
Its cost includes dense SVD and, when enabled, both candidate fits. It has
no automatic choice of the Algoim one-dimensional Gauss count `q`, no node
enrichment, no positivity or angular constraint on the corrected normal,
and no guarantee of uniform numerical rank for short arcs or slivers.
Extending it to other geometries or higher dimensions requires defining
the moment space, trusted boundary targets, scaling, node placement,
acceptance tolerance, and normal-fidelity criterion anew.
