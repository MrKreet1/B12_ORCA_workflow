#!/usr/bin/env bash
set -u

ROOT_DIR="calculations/B12"
ORCA_BIN="/opt/orca/orca"
LOG="run_all_timeout.log"
TIME_LIMIT="4h"

echo "============================================================" | tee -a "$LOG"
echo "B12 ORCA timeout run started: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG"
echo "TIME_LIMIT=$TIME_LIMIT per calculation" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

mapfile -t JOBS < <(find "$ROOT_DIR" -type f -name input.inp | grep -v "final_refinement" | sort)
TOTAL="${#JOBS[@]}"
COUNT=0

for inp in "${JOBS[@]}"; do
    COUNT=$((COUNT + 1))
    dir="$(dirname "$inp")"
    out="$dir/output.out"

    if [ -f "$out" ] && grep -q "ORCA TERMINATED NORMALLY" "$out"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP completed $COUNT/$TOTAL: $dir" | tee -a "$LOG"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $COUNT/$TOTAL: $dir" | tee -a "$LOG"

    (
        cd "$dir" || exit 100
        timeout "$TIME_LIMIT" "$ORCA_BIN" input.inp > output.out
    )
    rc=$?

    if [ "$rc" -eq 0 ] && grep -q "ORCA TERMINATED NORMALLY" "$out"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK $dir" | tee -a "$LOG"
    elif [ "$rc" -eq 124 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] TIMEOUT $dir" | tee -a "$LOG"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAIL rc=$rc $dir" | tee -a "$LOG"
    fi
done
