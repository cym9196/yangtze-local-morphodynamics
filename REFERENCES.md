# 参数来源与模型参考文献

本文件记录各物理参数、经验公式和模型实现的出处，便于复核与论文引用。

The parameter table deliberately distinguishes observations from assumptions.
No site-independent value should be presented as a calibrated Yangtze value.

## Yangtze observations

1. Hu et al. (2024), *Dune Development Dominates Flow Resistance Increase in
   a Large Dammed River*, Water Resources Research,
   https://doi.org/10.1029/2023WR036660. The alluvial middle-Yangtze subreaches
   are described as medium-fine sand with surface-bed `D50 approximately
   0.23 mm`. This is the direct basis for the baseline bed diameter.
2. Yuan et al. (2021), *Hydrodynamics, Sediment Transport and Morphological
   Features at the Confluence Between the Yangtze River and the Poyang Lake*,
   Water Resources Research, https://doi.org/10.1029/2020WR028284. Reported
   Yangtze bed-material D50 values include 162, 186 and 279 micrometres,
   supporting a sensitivity range around, but not a universal replacement for,
   0.23 mm.
3. The 2021 estuary review, *Declines in suspended sediment concentration and
   their geomorphological and biological impacts in the Yangtze River Estuary
   and adjacent sea*, https://doi.org/10.1016/j.ecss.2021.107724, reports an
   annual Datong SSC decrease from about 0.36 to 0.16 g/L between 1990--1999
   and 2000--2020. The 0.16 kg/m3 value is recorded only as a future
   multi-fraction inlet scenario.
4. Chinese environmental monitoring documentation reports Yangtze-estuary
   navigation-channel averages of about 0.32 kg/m3 in 2008--2014 and annual
   values of 0.10 and 0.16 kg/m3 in 2015 and 2016:
   https://www.mee.gov.cn/ywdt/gsgg/gongshi/wqgs_1/201809/W020210309392981491095.pdf.

The suspended load is much finer than the sandy bed at many sites. Therefore
the current single-fraction run uses zero imposed inlet concentration instead
of incorrectly assigning fine wash load the 0.23 mm bed-material diameter.

## Numerical and transport models

1. Renaud et al. (2026), *sedExnerFoam 2412: a 3D Exner-based sediment
   transport and morphodynamics model*, Geoscientific Model Development 19,
   2299--2333, https://doi.org/10.5194/gmd-19-2299-2026.
2. Upstream solver repository and validation cases:
   https://github.com/SedFoam/sedExnerFoam.
3. Meyer-Peter and Muller (1948), *Formulas for Bed-Load Transport*, original
   IAHR report record: https://repository.tudelft.nl/record/uuid:4fda9b61-be28-4703-ab06-43cdc2a21bd7.
4. Soulsby (1997), *Dynamics of Marine Sands*, threshold-of-motion chapter,
   https://doi.org/10.1680/doms.25844.0006.
5. OpenFOAM `nutkRoughWallFunction` documentation defines `Ks` as sand-grain
   roughness and gives `Cs=0.5--1.0`:
   https://doc.openfoam.com/2312/tools/processing/boundary-conditions/rtm/derived/wall/nutkRoughWallFunction/.

## Explicit engineering priors

The iron density, dimensions, mass and inertia are material or geometric
properties. In contrast, friction coefficients, Winkler stiffness, damping and
bearing capacity are engineering priors, not Yangtze measurements. They must be
replaced by submerged drag, plate-bearing and sinkage tests using the actual
bed sample. The present values are intended to make uncertainty visible and to
provide a reproducible starting point, not to claim field calibration.
