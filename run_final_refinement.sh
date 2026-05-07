#!/usr/bin/env bash
# Run the final PBE0-D3BJ/def2-TZVP Opt Freq refinement separately.

set -uo pipefail

FINAL_DIR="${FINAL_DIR:-calculations/B12/final_refinement}"
LOG_FILE="${LOG_FILE:-final_refinement.log}"

if [[ -z "${ORCA_BIN:-}" ]]; then
  ORCA_BIN="$(command -v orca || true)"
fi

if [[ -z "${ORCA_BIN}" ]]; then
  echo "ERROR: ORCA executable not found in PATH. Set ORCA_BIN=/full/path/to/orca" | tee -a "${LOG_FILE}"
  exit 1
fi

if [[ ! -f "${FINAL_DIR}/input.inp" ]]; then
  echo "ERROR: ${FINAL_DIR}/input.inp not found. Run python3 analyze_results.py first." | tee -a "${LOG_FILE}"
  exit 1
fi

if grep -qi "Grid5" "${FINAL_DIR}/input.inp"; then
  echo "ERROR: ${FINAL_DIR}/input.inp still contains Grid5; fix the input before running." | tee -a "${LOG_FILE}"
  exit 1
fi

if ! grep -Eq '(^|[[:space:]])Freq([[:space:]]|$)' "${FINAL_DIR}/input.inp"; then
  echo "ERROR: ${FINAL_DIR}/input.inp does not request Freq." | tee -a "${LOG_FILE}"
  exit 1
fi

{
  echo "============================================================"
  echo "B12 final refinement started: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "FINAL_DIR=${FINAL_DIR}"
  echo "ORCA_BIN=${ORCA_BIN}"
  echo "============================================================"
} | tee -a "${LOG_FILE}"

(
  cd "${FINAL_DIR}" || exit 2
  "${ORCA_BIN}" input.inp > output.out 2> output.err
)
rc=$?

if [[ "${rc}" -eq 0 ]] && grep -q "ORCA TERMINATED NORMALLY" "${FINAL_DIR}/output.out"; then
  echo "B12 final refinement finished normally: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${LOG_FILE}"
  exit 0
fi

echo "B12 final refinement failed rc=${rc}: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${LOG_FILE}"
exit "${rc}"
