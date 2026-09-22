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
