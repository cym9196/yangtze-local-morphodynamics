#!/bin/sh
# 中文说明：只读显示当前物理时刻、耦合进度、运行进程和剩余磁盘空间。
# Print a compact, read-only production status report.

cd "${0%/*}" || exit 1

if [ -z "${WM_PROJECT_DIR:-}" ]; then
    . /usr/lib/openfoam/openfoam2412/etc/bashrc
fi

saved=$(foamListTimes -processor -latestTime 2>/dev/null || true)
[ -n "$saved" ] || saved=0

# Saved directories advance only every 2 s. Read the most recent solver time
# as well so status remains useful between output writes without scanning the
# complete (potentially very large) log.
live=0
if [ -f log.immersedSedExnerFoam.63s ]; then
    live=$(tail -n 5000 log.immersedSedExnerFoam.63s \
        | awk '/^Time =/{t=$3} END{if (t=="") t=0; print t}')
fi
latest=$(awk -v a="$saved" -v b="$live" \
    'BEGIN {print (a>b ? a : b)}')

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
echo "latestSavedTime=$saved s"
echo "coupledIntervalProgress=$progress %"
df -h . | awk 'NR==2 {print "diskAvailable=" $4}'

if [ -f log.immersedSedExnerFoam.63s ]; then
    echo "lastSolverMessages:"
    tail -n 8 log.immersedSedExnerFoam.63s
fi
