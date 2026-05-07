# B12 ORCA 6.1 workflow

This package prepares a reproducible ORCA 6.1 workflow for B12 only.
It includes 3D, planar, and quasi-planar starting geometries so that a
preliminary 3D minimum is not treated as the final scientific conclusion before
the planar controls have been processed.

## Files

- `generate_inputs.py` — creates the calculation tree, B12 XYZ files, and ORCA `input.inp` files.
- `run_all.sh` — runs preliminary ORCA calculations sequentially and skips `final_refinement`.
- `run_planar_controls.sh` — runs representative planar/quasi-planar controls only.
- `run_final_refinement.sh` — runs `calculations/B12/final_refinement/input.inp` separately.
- `analyze_results.py` — parses ORCA outputs, creates `summary.csv`, `best_structure_B12.txt`, and final refinement input when a true minimum is found.
- `starter_xyz_examples/` — ready example XYZ files for d = 2.0 Å.
- `final_refinement_template.inp` — template for PBE0-D3BJ/def2-TZVP Opt Freq.
- `report_template_B12.md` — final report template.
- `summary_template.csv` — requested summary columns.

## Quick start on Linux VPS

```bash
cd B12_ORCA_workflow
python3 generate_inputs.py
chmod +x run_all.sh run_planar_controls.sh run_final_refinement.sh

# Recommended: run inside tmux
tmux new -s b12_orca
./run_all.sh
# detach: Ctrl-b then d

python3 analyze_results.py
```

## Required follow-up before final reporting

1. Run the planar/quasi-planar controls or at least two representative controls,
   for example `planar_double_ring` and `quasi_planar_buckled_double_ring`.
   The helper command is `./run_planar_controls.sh`.
2. Run `calculations/B12/final_refinement/input.inp` separately after
   `analyze_results.py` creates it for the current best preliminary minimum.
3. Base the final scientific conclusion only on completed, frequency-checked
   PBE0-D3BJ/def2-TZVP final-refinement results.

## ORCA command

The script resolves the full ORCA path via `command -v orca` and then runs:

```bash
/path/to/orca input.inp > output.out
```

This is safer for parallel ORCA jobs using `%pal nprocs 8 end`.

## Selection criterion

The preliminary best B12 geometry is the lowest-energy structure among
calculations that:

1. terminated normally;
2. have a converged SCF;
3. have no imaginary vibrational frequencies;
4. have the expected 30 vibrational modes for nonlinear B12.
