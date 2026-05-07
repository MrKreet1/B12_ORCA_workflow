# Computational Methods

All calculations were performed with ORCA 6.1.

The neutral B12 cluster was studied with charge 0 and multiplicities 1, 3, and 5. Initial geometries include 3D, planar, and quasi-planar motifs.

The preliminary level was PBE0-D3BJ/def2-SVP. The current report keeps the preliminary screening inputs as `Opt Freq`, because this directly validates each completed candidate as a true minimum or excludes it.

True minima were identified by:

1. normal ORCA termination;
2. SCF convergence;
3. absence of imaginary frequencies;
4. 30 vibrational modes for nonlinear B12;
5. lowest total energy among calculations satisfying the criteria above.

## Recommended Staged Workflow

For faster completion, future production runs should use a staged workflow:

1. Stage 1: PBE0-D3BJ/def2-SVP `Opt` for all starting geometries and multiplicities.
2. Stage 2: `Freq` only for the 10-20 lowest-energy optimized candidates.
3. Stage 3: PBE0-D3BJ/def2-TZVP `Opt Freq` for the final candidate set.

This keeps the current report scientifically conservative while avoiding full frequency calculations for every high-energy or failed preliminary structure in future runs.
