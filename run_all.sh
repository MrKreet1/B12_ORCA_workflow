#!/usr/bin/env bash
# run_all.sh
# Sequential ORCA runner for B12-only workflow.
# It never starts multiple ORCA jobs in parallel.

set -uo pipefail

ROOT_DIR="${ROOT_DIR:-calculations/B12}"
LOG_FILE="${LOG_FILE:-run_all.log}"
RESUME_MODE="${RESUME_MODE:-1}"
SLEEP_BETWEEN_JOBS="${SLEEP_BETWEEN_JOBS:-2}"

# For parallel ORCA jobs using %pal, a full path is safer than a bare "orca".
# Override when needed:
#   ORCA_BIN=/full/path/to/orca ./run_all.sh
if [[ -z "${ORCA_BIN:-}" ]]; then
  ORCA_BIN="$(command -v orca || true)"
fi

if [[ -z "${ORCA_BIN}" ]]; then
  echo "ERROR: ORCA executable not found in PATH. Set ORCA_BIN=/full/path/to/orca" | tee -a "${LOG_FILE}"
  exit 1
fi

if [[ ! -d "${ROOT_DIR}" ]]; then
  echo "ERROR: ${ROOT_DIR} does not exist. Run: python3 generate_inputs.py" | tee -a "${LOG_FILE}"
  exit 1
fi

start_ts="$(date '+%Y-%m-%d %H:%M:%S')"
{
  echo "============================================================"
  echo "B12 ORCA sequential run started: ${start_ts}"
  echo "ROOT_DIR=${ROOT_DIR}"
  echo "ORCA_BIN=${ORCA_BIN}"
  echo "RESUME_MODE=${RESUME_MODE}"
  echo "No parallel job dispatch is used; one input.inp is run at a time."
  echo "============================================================"
} | tee -a "${LOG_FILE}"

mapfile -t INPUTS < <(find "${ROOT_DIR}" -type f -name input.inp | sort)

if [[ "${#INPUTS[@]}" -eq 0 ]]; then
  echo "ERROR: no input.inp files found under ${ROOT_DIR}" | tee -a "${LOG_FILE}"
  exit 1
fi

TOTAL="${#INPUTS[@]}"
COUNT=0
OK=0
FAIL=0
SKIP=0

for input_path in "${INPUTS[@]}"; do
  COUNT=$((COUNT + 1))
  calc_dir="$(dirname "${input_path}")"

  # Skip template folders that intentionally contain placeholders.
  if grep -q "REPLACE_WITH_BEST_MULTIPLICITY" "${input_path}"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP template ${calc_dir}" | tee -a "${LOG_FILE}"
    SKIP=$((SKIP + 1))
    continue
  fi

  if [[ "${RESUME_MODE}" == "1" && -f "${calc_dir}/output.out" ]] && grep -q "ORCA TERMINATED NORMALLY" "${calc_dir}/output.out"; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP completed ${COUNT}/${TOTAL}: ${calc_dir}" | tee -a "${LOG_FILE}"
    SKIP=$((SKIP + 1))
    continue
  fi

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] START ${COUNT}/${TOTAL}: ${calc_dir}" | tee -a "${LOG_FILE}"

  (
    cd "${calc_dir}" || exit 2

    # Remove stale temporary files from interrupted previous jobs, but keep useful outputs.
    rm -f ./*.tmp ./*.tmp.* ./input.tmp* ./input.*.tmp 2>/dev/null || true

    job_start="$(date +%s)"
    "${ORCA_BIN}" input.inp > output.out 2> output.err
    rc=$?
    job_end="$(date +%s)"
    elapsed=$((job_end - job_start))

    if [[ "${rc}" -eq 0 ]] && grep -q "ORCA TERMINATED NORMALLY" output.out; then
      echo "SUCCESS elapsed_seconds=${elapsed}"
      exit 0
    fi

    echo "FAILED rc=${rc} elapsed_seconds=${elapsed}"
    if [[ -s output.err ]]; then
      echo "--- output.err tail ---"
      tail -n 30 output.err || true
    fi
    if [[ -s output.out ]]; then
      echo "--- output.out tail ---"
      tail -n 50 output.out || true
    fi
    exit 1
  ) >> "${LOG_FILE}" 2>&1

  job_rc=$?
  if [[ "${job_rc}" -eq 0 ]]; then
    OK=$((OK + 1))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK ${calc_dir}" | tee -a "${LOG_FILE}"
  else
    FAIL=$((FAIL + 1))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAIL ${calc_dir}" | tee -a "${LOG_FILE}"
  fi

  sleep "${SLEEP_BETWEEN_JOBS}"
done

end_ts="$(date '+%Y-%m-%d %H:%M:%S')"
{
  echo "============================================================"
  echo "B12 ORCA sequential run finished: ${end_ts}"
  echo "Total input files: ${TOTAL}"
  echo "Successful: ${OK}"
  echo "Failed: ${FAIL}"
  echo "Skipped: ${SKIP}"
  echo "Log file: ${LOG_FILE}"
  echo "============================================================"
} | tee -a "${LOG_FILE}"

if [[ "${FAIL}" -gt 0 ]]; then
  exit 1
fi

exit 0
