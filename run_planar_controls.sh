#!/usr/bin/env bash
# Run representative planar/quasi-planar B12 controls without running the full tree.

set -uo pipefail

LOG_FILE="${LOG_FILE:-planar_controls.log}"
CONTROLS=(
  "calculations/B12/planar_double_ring"
  "calculations/B12/quasi_planar_buckled_double_ring"
)

rc_total=0
for control_dir in "${CONTROLS[@]}"; do
  if [[ ! -d "${control_dir}" ]]; then
    echo "ERROR: ${control_dir} not found. Run python3 generate_inputs.py first." | tee -a "${LOG_FILE}"
    rc_total=1
    continue
  fi

  ROOT_DIR="${control_dir}" LOG_FILE="${LOG_FILE}" ./run_all.sh
  rc=$?
  if [[ "${rc}" -ne 0 ]]; then
    rc_total="${rc}"
  fi
done

exit "${rc_total}"
