#!/bin/sh
# Export every reconstructed bed surface and render the 0.05 m-deep solid
# sand layer. Run this after reconstructing the desired MPI time directories.

cd "${0%/*}" || exit 1

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    . /usr/lib/openfoam/openfoam2412/etc/bashrc
fi

set -e

# -noZero omits the undeformed initial field; each exported frame retains its
# physical OpenFOAM time in the generated .vtm.series index.
foamToVTK -overwrite -noZero -no-internal -no-fields -patches '(bed)' \
    -name VTK_bed_animation > log.foamToVTK 2>&1

rm -rf frames_solid
python3 render_solid_sand.py \
    --series VTK_bed_animation/yangtzeIronBlock3D_1to25mm.vtm.series \
    --motion processor0/postProcessing/immersedBodyMotion/motion.dat \
    --output frames_solid \
    --exaggeration 8

echo "Solid-sand frames and metrics.csv completed"
