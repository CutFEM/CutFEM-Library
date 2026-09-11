# Library build configuration

Shared CMakePresets.json provides orientation and disables optional examples/documentation.
Ignored CMakeUserPresets.json captures this machine's compiler, MPI, MUMPS and other dependencies.
Create a local preset from the existing working CMakeCache.txt with the hub's
`scripts/local_presets.py --role development --repo <development> --cache <cache> --apply`.

Use `cmake --preset local`, then `cmake --build --preset local --target solver -j 4`
when the task calls for configuring/rebuilding the library. Build all required library targets
for a fresh checkout. Do not silently reconfigure a running production build.

## Deliberately retained legacy constraint

The current library creates cpp/cutFEMConfig.h using a path relative to build/, and cfmpi.hpp
includes that source-tree header. Consequently these presets keep the library binary directory
at <development>/build. Use a separate source checkout for an independent library configuration.
Moving generated headers entirely into the binary tree needs a separate, tested code change;
the workspace restructuring did not modify the user's uncommitted library implementation.

The library's existing FindMUMPS accepts explicit cached include/library paths and the MUMPS_DIR
environment variable. Preserve the local workstation's MPI/MUMPS/BLAS linkage configuration.
The restructuring does not claim a tested Linux or Dardel library rebuild.
