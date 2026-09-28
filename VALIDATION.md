# Reproducibility and feasibility record

Test host: Ubuntu/OpenFOAM v2412, AMD Ryzen 5 7500F (6 physical cores),
15 GiB RAM, no swap. Tests were run on 2026-09-27 and 2026-09-28.

## Mesh verification

- volume cells: 3,072,832;
- finite-area bed faces: 23,104;
- bounding box: 0.50 x 0.50 x 0.25 m;
- minimum/maximum cell volume: 1e-9 / 6.25e-6 m3;
- maximum aspect ratio: 25;
- maximum/average non-orthogonality: 0 / 0 deg;
- maximum skewness: 5.54e-13;
- `checkMesh`: Mesh OK;
- `checkFaMesh`: completed successfully.

## Solver checks

The basic smoke test advanced two parallel steps and exited with status 0.
The coupled smoke test then forced both release times to zero and verified, in
the same run:

- suspension equation: `Max(Cs)=8.80e-4` after two steps;
- moving Exner bed: maximum vertex motion `2.37e-6 m`;
- six-DoF body: centre changed from `(0.25 0.25 0.10)` to
  `(0.25000281 0.25 0.099999807)`;
- contact force was finite and included tangential friction;
- parallel run exited with status 0 and no swap activity;
- elapsed time was 92.43 s for two start-up steps at about 586% CPU.

The production release times (`3 s`) and production end time (`63 s`) were
automatically restored after the test.

## Parallel decomposition benchmark

The original generic Scotch split created upper-water partitions with zero
finite-area bed faces at 12 ranks. It was replaced by a generated hierarchical
x/y split with `nz=1`, so every process spans the depth and owns bed faces.

| configuration | decomposition | two-step wall time | result |
|---|---:|---:|---|
| 6 physical cores | 3 x 2 x 1 | 82.32 s | selected |
| 12 hardware threads | 4 x 3 x 1 | 93.15 s | slower |

At 12 ranks the maximum bed-face count was only 0.66% above the mean, proving
that the finite-area work was balanced. Six physical cores are retained as the
desktop default because simultaneous multithreading added overhead.

## Full-run feasibility on this host

At a 1 mm minimum cell and 0.50 m/s flow, `maxDeltaT=0.0005 s` implies at
least 126,000 steps for 63 s. Startup timing is pessimistic, but it projects
roughly 3--8 weeks of uninterrupted wall time on this six-core CPU. The mesh
itself fits in memory (observed total system use about 8 GiB), but the current
machine is a validation/pilot host, not an efficient production host for this
resolution.

A practical production target is a Linux HPC node or allocation with at least
64 CPU cores, 128 GiB RAM and 100 GiB scratch space. Actual strong-scaling and
storage behavior must be measured with a 0.1 s benchmark before submitting the
63 s job.
