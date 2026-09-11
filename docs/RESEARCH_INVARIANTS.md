# Research-facing library invariants

These rules summarize implementation constraints in the current drivers. For historical
evidence, locate the relevant private project through the machine's path map.

- Distinct trial/test active meshes must be mapped through the background element before
  selecting basis functions and DOFs. A fix in one assembly path does not establish all
  contribution overloads are safe for cross-mesh forms.
- Algoim-backed geometry must use its matching stored volume, surface and cut-face rules.
  Do not silently route a geometry-dependent contribution through legacy cut partitioning.
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
