# 可复现性与可行性验证记录

本文件保存网格、短时耦合、数值稳定性、并行效率和资源需求的验证证据。

Test host: Ubuntu/OpenFOAM v2412, AMD Ryzen 5 7500F (6 physical cores),
15 GiB RAM, no swap. Tests were run on 2026-09-27 and 2026-09-28.

## Mesh verification

- target volume cells: 125,000 (uniform 50 x 50 x 50; pending fresh OpenFOAM check);
- target finite-area bed faces: 2,500 (uniform 50 x 50; pending fresh OpenFOAM check);
- water-mesh bounding box: 0.05 x 0.05 x 0.05 m;
- target minimum/maximum cell volume: 1e-9 / 1e-9 m3;
- target maximum aspect ratio: 1.0;
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

The orientation-aware contact and dynamic-buoyancy upgrades were compiled
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

### Release-instability regression

The first 5 cm production attempt failed after release because the 15.74 g
body left the short domain and the unconstrained moving-bed/RANS feedback later
produced unbounded fields. The original crash evidence is preserved remotely
as `log.crash.5cm.unbounded` plus its SHA-256 record.

The corrected solver adds a body-only horizontal retaining frame, resolved
linear/angular damping, bounded `Cs/k/omega/nut`, a 1e-6 m per-step Exner guard,
and a pressure reference point inside the resized domain. A deliberately harsh
test released both the body and sediment at time zero and reached 2.5000831 s
on six ranks with exit status 0:

- body centre `(0.0400961 0.0243947 0.0047722) m`, inside the test section;
- body velocity below 0.001 m/s and angular speed below 0.36 rad/s at exit;
- hydrodynamic force remained O(0.03 N), without large-number growth;
- `Cs` remained non-negative, with maximum 0.18495;
- sediment thickness remained 0.04718--0.05244 m;
- turbulence remained finite (`k` maximum about 0.034 m2/s2);
- elapsed wall time 310.75 s, maximum per-rank RSS 106288 KiB, no swap;
- no floating-point exception, fatal error, segmentation fault or OOM.

The evidence is preserved remotely as `log.couplingSmoke.5cmStability2p5s`,
`log.couplingSmoke.resources.5cmStability2p5s` and
`5cmStability2p5s.sha256`. The Exner limiter was active, so this run establishes
numerical robustness, not site-calibrated morphological accuracy.

That short regression was insufficient: the subsequent production run
reproduced a floating-point failure at 6.2786213 s. Once the adaptive CFD time
step collapsed, the fixed per-step Exner cap implied an unbounded bed-mesh
velocity and amplified the pressure/immersed-boundary feedback. The current
candidate therefore adds a second `2.5e-3 m/s` bed-change-rate cap, keeps the
diffuse solid mask 0.003 m inside the horizontal CFD boundaries, resolves the
0.0002 s Brinkman time with `maxDeltaT=0.0002 s`, and uses two PIMPLE outer
loops plus three pressure correctors. Stability beyond 6.2786213 s remains a
required acceptance test, not an assumed result.

The rebuilt uniform 1 mm candidate passed a harsher immediate-release test to
0.2 s on six physical cores. It completed 1,000+ coupled steps with exit status
0 in 357.85 s wall time and 142,696 KiB peak RSS per rank. At the final state:

- maximum cell Courant number was 0.18413;
- body centre was `(0.0370289 0.0246834 0.0103146) m`, after a stable encounter
  with the inset retainer;
- hydrodynamic force remained O(0.06 N), contact support O(0.14 N);
- maximum velocity was 0.71355 m/s and maximum `Cs` was 0.14047;
- sediment thickness remained 0.04951--0.05030 m;
- no NaN, fatal error, segmentation fault or floating-point signal occurred.

The checksummed evidence is stored on the simulation disk under
`validation_evidence/1mmFix2_0p2s`. This establishes short-horizon coupled
stability; the running 63 s production case remains the acceptance test beyond
the previous 6.2786213 s failure time.

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

At a 1 mm minimum cell and 0.50 m/s flow, `maxDeltaT=0.0002 s` implies at
least 315,000 steps for 63 s. The 0.2 s uniform-grid benchmark projects about
31.3 h on this six-core CPU; reserve a 36--48 h uninterrupted window because
bed motion, output and reconstruction add workload later. The measured memory
use is small relative to 15 GiB, but the final wall time remains provisional
until the full run completes.

The local Ubuntu host is now a practical production target. HPC is optional;
if used, strong-scaling should be measured first; six physical cores remain the
desktop baseline and 8--16 ranks are the first sensible HPC measurements.
