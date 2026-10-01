#!/bin/bash
# 中文说明：固定使用 /dev/sdc2 上的用户程序、动态库和当前算例，前台执行正式计算。

case_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
user_dir=$(CDPATH= cd -- "$case_dir/../.." && pwd)

# 先载入 OpenFOAM 系统环境，再覆盖用户目录；否则 bashrc 会默认回到
# /home/cym/OpenFOAM/cym-v2412。
. /usr/lib/openfoam/openfoam2412/etc/bashrc
set -e
export WM_PROJECT_USER_DIR="$user_dir"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"

cd "$case_dir"
NPROCS="${NPROCS:-6}"
export NPROCS

# Reuse a valid decomposed checkpoint when one exists.  This keeps an SSH
# reconnect or machine restart from silently deleting an in-progress run.
latest=$(foamListTimes -processor -latestTime 2>/dev/null || true)
if [ -n "$latest" ]; then
    echo "Found decomposed checkpoint at physical time $latest; resuming."
    exec ./Allresume63s
fi

echo "No decomposed checkpoint found; starting a new production run."
exec ./Allrun63s
