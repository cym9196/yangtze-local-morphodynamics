#!/bin/sh
# Print a compact, read-only production status report.

cd "${0%/*}" || exit 1

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    . /usr/lib/openfoam/openfoam2412/etc/bashrc
fi

latest=$(foamListTimes -processor -latestTime 2>/dev/null || true)
[ -n "$latest" ] || latest=0

progress=$(awk -v t="$latest" 'BEGIN {
    p=(t<=3 ? 0 : (t-3)/60*100);
    if (p>100) p=100;
    printf "%.3f", p
}')

if pgrep -f '(^|/)immersedSedExnerFoam( |$)' >/dev/null 2>&1; then
    state=RUNNING
else
    state=STOPPED
fi

echo "state=$state"
echo "latestPhysicalTime=$latest s"
echo "coupledIntervalProgress=$progress %"
df -h . | awk 'NR==2 {print "diskAvailable=" $4}'

if [ -f log.immersedSedExnerFoam.63s ]; then
    echo "lastSolverMessages:"
    tail -n 8 log.immersedSedExnerFoam.63s
fi
