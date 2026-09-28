# Reproducibility and feasibility record

Test host: Ubuntu/OpenFOAM v2412, AMD Ryzen 5 7500F (6 physical cores),
15 GiB RAM, no swap. Tests were run on 2026-09-27 and 2026-09-28.

## Mesh verification

- volume cells: 35,840;
- finite-area bed faces: 1,024;
- water-mesh bounding box: 0.05 x 0.05 x 0.05 m;
- minimum/maximum cell volume: 1e-9 / 6.25e-6 m3;
- maximum aspect ratio: 4.075;
- maximum/average non-orthogonality: 0 / 0 deg;
- maximum skewness: 5.93e-14;
- `checkMesh`: Mesh OK;
- `checkFaMesh`: completed successfully.

## Solver checks

The basic smoke test advanced two parallel steps and exited with status 0.
The coupled smoke test then forced both release times to zero and verified, in
the same run:

- suspension equation: `Max(Cs)=8.91e-4` after two steps;
- moving Exner bed: maximum vertex motion `1.29e-6 m` in the latest regression;
- six-DoF body: centre changed from `(0.025 0.025 0.010)` to
  `(0.025003456 0.025 0.0099998236)`;
- six-face contact force was finite and included tangential friction
  `(-3.31e-4, 0, 0.02132) N` after the second regression step;
- buoyancy was integrated from the actually submerged body-mask volume
  (about 2e-6 m3) and remains pose dependent;
- parallel run exited with status 0 and no swap activity;
- latest 5 cm regression elapsed time was 2.21 s for two coupled start-up
  steps; `/usr/bin/time` returned exit status 0.

The production release times (`3 s`) and production end time (`63 s`) were
automatically restored after the test.

The orientation-aware contact and partial-submergence upgrades were compiled
with OpenFOAM v2412 and passed a fresh two-step, six-rank fully coupled
regression: suspension reached `Max(Cs)=8.91e-4`, the Exner bed moved, the body translated, and the solver
exited with status 0. The quadrature verifier independently confirms 54 points,
0.001 m2 total cuboid surface area, 0.0001 m2 top/bottom faces and 0.0002 m2 long
side faces.

The contact gap and support/friction directions were subsequently upgraded
from a fixed vertical direction to the nearest evolving finite-area bed tangent
plane. The solver compiled and repeated the same coupled regression with
identical flat-bed reference forces and exit status 0, as expected before the
bed develops a slope.

The authoritative remote evidence for the resized case is preserved as
`log.couplingSmoke.5cmRegression`,
`log.couplingSmoke.resources.5cmRegression` and `5cmRegression.sha256` in
the Ubuntu case directory; it is not committed because runtime logs are
deliberately excluded from Git.

Smoke/restart dictionary handling was tested separately with a zero-duration
coupling run. SHA-256 checks confirmed that `system/controlDict` and
`constant/immersedBodyProperties` were restored byte-for-byte, so validation
runs no longer expand includes or silently reformat production inputs.

The quality-control function objects (`yPlus`, bed `wallShearStress` and
`fieldMinMax`) were then instantiated in a one-step six-rank solver run. The
solver recognized the bed patch, advanced normally and exited with status 0;
the production end time was again restored to 63 s.

## Parallel decomposition benchmark

The original generic Scotch split created upper-water partitions with zero
finite-area bed faces at 12 ranks. It was replaced by a generated hierarchical
x/y split with `nz=1`, so every process spans the depth and owns bed faces.

| configuration | decomposition | two-step wall time | result |
|---|---:|---:|---|
| 6 physical cores | 3 x 2 x 1 | 2.21 s | selected |

The earlier, larger 0.50 m case showed that 12 hardware threads were slower
than six physical cores. Six physical cores are therefore retained as the
desktop default; any HPC layout still requires a fresh strong-scaling test.

## Full-run feasibility on this host

At a 1 mm minimum cell and 0.50 m/s flow, `maxDeltaT=0.0005 s` implies at
least 126,000 steps for 63 s. Startup timing is pessimistic, but it projects
roughly 1--2 days of uninterrupted wall time on this six-core CPU. This is a
start-up extrapolation rather than a guarantee; the new 35,840-cell mesh easily
fits on the current desktop host.

The local Ubuntu host is now a practical production target. HPC is optional;
if used, strong-scaling should be measured first because 64 MPI ranks would be
excessive for only 35,840 cells.
