# Research-facing library invariants

These rules summarize implementation constraints in the current drivers. For historical
evidence, locate the relevant private project through the machine's path map.

- Distinct trial/test active meshes must be mapped through the background element before
  selecting basis functions and DOFs. A fix in one assembly path does not establish all
  contribution overloads are safe for cross-mesh forms.
- Algoim-backed geometry must use its matching stored volume, surface and cut-face rules.
  Do not silently route a geometry-dependent contribution through legacy cut partitioning.
- When an IBP-corrected triangle rule is subdivided, refit the concatenated surface and
  volume rules against the parent cell's moments before assembly. Child rules can each
  be rank-deficient on short arcs, so child corrections alone do not guarantee the
  divergence theorem on the original cell. Check the returned rule's low-degree moments.
- The opt-in direct-SVD surface fit and final-parent moment gate are described in
  [ALGOIM_DIRECT_SVD_FIT.md](ALGOIM_DIRECT_SVD_FIT.md). Moment acceptance does not
  constrain corrected-normal orientation or general surface-integral accuracy.
- Parametric Mesh2 cut rules (`ProblemOption::parametric_cut_rule_`,
  `problem/ParametricCutQuadrature.hpp`): build `AlgoimInterface` from the same
  `ProblemOption`; it then marks an element cut exactly when the rule has an arc. Do not add
  thresholds (arc length, root distance, narrow band) to the classification or the root
  finding: neighbours must agree on shared edges, or the divergence theorem summed over the
  mesh fails. Both sides of the interface use one arc; the 4-argument `quadGenVol` does not
  negate the level set.
- The parametric visibility test samples arc nodes: it guarantees positive weights, not that
  nodes stay inside the triangle. `workfiles:src/stokes/parametric_cut_probe.cpp` checks both.
- IBP-corrected (moment-fitted) Algoim triangle rules can carry non-positive weights; do not
  pass them to helpers that assert positive weights (`L2_norm_cut`, `integral_algoim`).
- The integration helpers in `problem/AlgoimIntegration.hpp` take cut rules from the global
  `defaultProblemOption` (a `static`, one copy per translation unit); set it explicitly where
  norms are computed. Stationary `L2_norm_cut` takes `order_space` for uncut elements (default
  4, i.e. the degree-5 rule); `QF_Simplex<R2>` stops at degree 9.
- The aggregate extension extrapolates root polynomials over chains of cells and so amplifies
  noise in its input. Check polynomial reproduction with analytic data, not with derivatives
  of an interpolant.
- Patch stabilization integrates its intended full background patch, not physical cut volume.
- Preserve all CutFEMParameter coefficients in unified assembly.
- `ListItemVF::reduce()` sums only items that `ItemVF::operator==` deems identical up to the
  scalar `c`; that comparison must include every factor of the term (parameter lists `coefu`/
  `coefv` and `pfunU`/`pfunV` included). Only exactly cancelled items are dropped, and a fully
  cancelled list keeps one zero item because assembly takes the FE space from `VF[0]`.
  `cpp/example/tests/test_itemvf_reduce.cpp` checks this.
- Structured GridPhi evaluation must address the explicitly resolved FE element rather than
  relying on mutable cached-element state in FunFEM.
- Maintain MPI-consistent replicated geometry/fallback data when every rank evaluates a
  complete background field.
- A passing build does not establish integration-by-parts compatibility, conservation,
  cut-uniform stability, or a mathematical inf-sup condition.

Detailed private project findings belong in the hub; keep this document implementation-focused.
