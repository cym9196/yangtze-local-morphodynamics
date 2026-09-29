#!/bin/sh
# 中文说明：在 SSH 断开后继续运行全域 1 mm、3+60 s 正式计算。

set -e
case_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$case_dir"

if pgrep -af "immersedSedExnerFoam -parallel" | grep -v grep >/dev/null 2>&1
then
    echo "已有 immersedSedExnerFoam 并行任务，拒绝重复启动。" >&2
    pgrep -af "immersedSedExnerFoam|mpirun" >&2 || true
    exit 1
fi

nohup "$case_dir/run_from_sdc2.sh" \
    > "$case_dir/log.driver.1mm60s" 2>&1 < /dev/null &
pid=$!
printf '%s\n' "$pid" > "$case_dir/production.pid"
echo "后台正式计算已启动，包装进程 PID=$pid"
