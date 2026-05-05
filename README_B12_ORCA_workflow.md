# B12 ORCA 6.1 3D workflow

This package prepares a reproducible ORCA 6.1 workflow for B12 only.

## Files

- `generate_inputs.py` — creates the calculation tree, 3D B12 XYZ files, and ORCA `input.inp` files.
- `run_all.sh` — runs all ORCA calculations sequentially.
- `analyze_results.py` — parses ORCA outputs, creates `summary.csv`, `best_structure_B12.txt`, and final refinement input when a true minimum is found.
- `starter_xyz_examples/` — ready example XYZ files for d = 2.0 Å.
- `final_refinement_template.inp` — template for PBE0-D3BJ/def2-TZVP Opt Freq.
- `report_template_B12.md` — final report template.
- `summary_template.csv` — requested summary columns.

## Quick start on Linux VPS

```bash
cd B12_ORCA_workflow
python3 generate_inputs.py
chmod +x run_all.sh

# Recommended: run inside tmux
tmux new -s b12_orca
./run_all.sh
# detach: Ctrl-b then d

python3 analyze_results.py
```

## ORCA command

The script resolves the full ORCA path via `command -v orca` and then runs:

```bash
/path/to/orca input.inp > output.out
```

This is safer for parallel ORCA jobs using `%pal nprocs 8 end`.

## Selection criterion

The best B12 geometry is the lowest-energy structure among calculations that:

1. terminated normally;
2. have a converged SCF;
3. have no imaginary vibrational frequencies;
4. have the expected 30 vibrational modes for nonlinear B12, when the output contains the full mode list.
