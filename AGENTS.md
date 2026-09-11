# CutFEM library instructions

This is the shared CutFEM library, currently on the `development` branch. Its origin is
the public CutFEM/CutFEM-Library repository. Keep unpublished scientific conclusions,
private manuscripts, and personal planning in the private research hub.

## Shared scientific context

Read `~/.config/phd/paths.json` to locate `research-hub`, related code, and manuscripts.
Read the hub's `WORKFLOW.md` and relevant `projects/<project>/STATUS.md` before substantive
work, then follow links to the needed knowledge, proof register, and investigations.
`alias:path` references are resolved through the local path map; do not assume Mac/Linux paths.
Automatically update the appropriate project status/knowledge after substantive work.
Record evidence and unresolved issues; never call an agent-written proof human-checked.
Keep investigations out of AGENTS.md/CLAUDE.md. Git actions follow the user's current
authorization; never stage or commit pre-existing changes as part of an unrelated task.

## Library orientation

- `cpp/common/`: meshes, geometry, quadrature, active interfaces, common data.
- `cpp/FESpace/`: finite elements, FE functions, operators and expressions.
- `cpp/problem/`: weak forms, assembly and CutFEM contributions.
- `cpp/solver/`, `cpp/parallel/`: solvers and MPI/OpenMP integration.
- `docs/BUILD.md`: existing library build conventions and local presets.
- `docs/RESEARCH_INVARIANTS.md`: compact reusable invariants from current research drivers.

Changes to this library can invalidate results in several projects. Rebuild affected library
targets and the dependent workfiles executable, and run the relevant small verification.
Record the library commit or uncommitted patch used by any important experiment.

Do not infer push/publication permission from the presence of an origin remote. Follow the
user's current task instructions. Maintain implementation documentation automatically.
