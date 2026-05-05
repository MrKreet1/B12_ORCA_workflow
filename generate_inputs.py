#!/usr/bin/env python3
"""
generate_inputs.py

Fully reproducible generator for 3D B12 ORCA 6.1 input folders.

Scope:
- B12 only
- neutral cluster only: charge = 0
- 3D starting geometries only
- multiplicities: 1, 3, 5
- distances d: 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0 Angstrom
- preliminary method: PBE0-D3BJ/def2-SVP, Opt Freq
- alternative method is documented in comments: B3LYP-D3BJ/def2-SVP
- ORCA resources: %pal nprocs 8 end, %maxcore 2500

This script DOES NOT invent energies or optimized structures. It only creates
initial XYZ geometries and ORCA input files.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

Atom = Tuple[str, float, float, float]
Coord = Tuple[float, float, float]

CLUSTER = "B12"
N_ATOMS = 12
CHARGE = 0
MULTIPLICITIES = [1, 3, 5]
DISTANCES_A = [5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0]
ROOT = Path("calculations") / CLUSTER

NPROCS = 8
MAXCORE_MB = 2500

# Preliminary search level requested in the workflow.
METHOD = "PBE0"
DISPERSION = "D3BJ"
BASIS = "def2-SVP"

# Alternative preliminary level, not used by default. To use it, change METHOD.
ALTERNATIVE_METHOD_COMMENT = "Alternative screening level: B3LYP-D3BJ/def2-SVP"

RANDOM_SEED = 12012026


@dataclass(frozen=True)
class StructureSpec:
    name: str
    builder: Callable[[float], List[Atom]]
    description: str


def center_atoms(atoms: Sequence[Atom]) -> List[Atom]:
    """Translate atoms so that the centroid is at the origin."""
    cx = sum(a[1] for a in atoms) / len(atoms)
    cy = sum(a[2] for a in atoms) / len(atoms)
    cz = sum(a[3] for a in atoms) / len(atoms)
    return [(el, x - cx, y - cy, z - cz) for el, x, y, z in atoms]


def scale_coords(coords: Sequence[Coord], factor: float) -> List[Coord]:
    return [(x * factor, y * factor, z * factor) for x, y, z in coords]


def atoms_from_coords(coords: Sequence[Coord]) -> List[Atom]:
    return [("B", x, y, z) for x, y, z in coords]


def distance(a: Atom, b: Atom) -> float:
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def min_distance(atoms: Sequence[Atom]) -> float:
    md = float("inf")
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            md = min(md, distance(atoms[i], atoms[j]))
    return md


def pair_distance_stats(atoms: Sequence[Atom]) -> Dict[str, float]:
    ds: List[float] = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            ds.append(distance(atoms[i], atoms[j]))
    return {
        "min_distance_A": min(ds),
        "max_distance_A": max(ds),
        "mean_distance_A": sum(ds) / len(ds),
    }


def validate_b12_3d(name: str, atoms: Sequence[Atom]) -> None:
    if len(atoms) != N_ATOMS:
        raise ValueError(f"{name}: expected {N_ATOMS} atoms, got {len(atoms)}")
    if any(atom[0] != "B" for atom in atoms):
        raise ValueError(f"{name}: only boron atoms are allowed")
    xs = [a[1] for a in atoms]
    ys = [a[2] for a in atoms]
    zs = [a[3] for a in atoms]
    if max(zs) - min(zs) < 1e-6:
        raise ValueError(f"{name}: structure is planar in z; 3D geometry required")
    if max(xs) - min(xs) < 1e-6 or max(ys) - min(ys) < 1e-6:
        raise ValueError(f"{name}: degenerate coordinate range; 3D geometry required")
    if min_distance(atoms) < 0.45:
        raise ValueError(f"{name}: atoms are unrealistically close")


def icosahedron(d: float) -> List[Atom]:
    """12 vertices of a regular icosahedron; nearest-neighbor edge length ~= d."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    base: List[Coord] = []
    for y in (-1.0, 1.0):
        for z in (-phi, phi):
            base.append((0.0, y, z))
    for x in (-1.0, 1.0):
        for y in (-phi, phi):
            base.append((x, y, 0.0))
    for x in (-phi, phi):
        for z in (-1.0, 1.0):
            base.append((x, 0.0, z))
    # In this coordinate convention the edge length is 2.
    coords = scale_coords(base, d / 2.0)
    return center_atoms(atoms_from_coords(coords))


def distorted_icosahedron(d: float, seed: int, amplitude: float) -> List[Atom]:
    """Slightly distorted icosahedron to remove excessive symmetry."""
    rng = random.Random(seed)
    atoms = icosahedron(d)
    distorted: List[Atom] = []
    for idx, (el, x, y, z) in enumerate(atoms):
        # Deterministic radial + Cartesian perturbation. The amplitude is a fraction of d.
        r = math.sqrt(x * x + y * y + z * z) or 1.0
        radial = 1.0 + rng.uniform(-0.035, 0.035)
        dx = rng.uniform(-amplitude, amplitude) * d
        dy = rng.uniform(-amplitude, amplitude) * d
        dz = rng.uniform(-amplitude, amplitude) * d
        distorted.append((el, x * radial + dx, y * radial + dy, z * radial + dz))
    return center_atoms(distorted)


def compact_3d(d: float) -> List[Atom]:
    """
    Compact non-planar 12-atom cluster: distorted tetrahedral/cage-like shell.
    The coordinates are scaled so that the shortest B-B distance is close to d.
    """
    raw: List[Coord] = [
        (0.00, 0.00, 1.28),
        (1.18, 0.00, 0.43),
        (-0.59, 1.02, 0.43),
        (-0.59, -1.02, 0.43),
        (0.72, 1.18, -0.32),
        (-0.72, 1.18, -0.32),
        (0.72, -1.18, -0.32),
        (-0.72, -1.18, -0.32),
        (1.12, 0.00, -1.08),
        (-1.12, 0.00, -1.08),
        (0.00, 0.98, -1.28),
        (0.00, -0.98, -1.28),
    ]
    atoms0 = atoms_from_coords(raw)
    scale = d / min_distance(atoms0)
    return center_atoms(atoms_from_coords(scale_coords(raw, scale)))


def bilayer_3d(d: float) -> List[Atom]:
    """Two staggered 6-atom rings separated along z; explicitly non-planar."""
    coords: List[Coord] = []
    radius = d
    zsep = 0.62 * d
    for k in range(6):
        ang = 2.0 * math.pi * k / 6.0
        coords.append((radius * math.cos(ang), radius * math.sin(ang), +zsep / 2.0))
    for k in range(6):
        ang = 2.0 * math.pi * (k + 0.5) / 6.0
        # Slightly different radius avoids artificial high symmetry.
        coords.append((0.94 * radius * math.cos(ang), 0.94 * radius * math.sin(ang), -zsep / 2.0))
    return center_atoms(atoms_from_coords(coords))


def cuboctahedral_fragment(d: float) -> List[Atom]:
    """12 vertices of a cuboctahedron-like high-symmetry 3D structure."""
    raw: List[Coord] = []
    for x in (-1.0, 1.0):
        for y in (-1.0, 1.0):
            raw.append((x, y, 0.0))
    for x in (-1.0, 1.0):
        for z in (-1.0, 1.0):
            raw.append((x, 0.0, z))
    for y in (-1.0, 1.0):
        for z in (-1.0, 1.0):
            raw.append((0.0, y, z))
    # Cuboctahedron edge length in these coordinates is sqrt(2).
    coords = scale_coords(raw, d / math.sqrt(2.0))
    return center_atoms(atoms_from_coords(coords))


def distorted_cage(d: float) -> List[Atom]:
    """Deformed 3D cage made from an anisotropically distorted icosahedral shell."""
    base = icosahedron(d)
    coords: List[Coord] = []
    for i, (_, x, y, z) in enumerate(base):
        sx = 1.10 + 0.03 * math.sin(i)
        sy = 0.88 + 0.04 * math.cos(2 * i)
        sz = 1.18 + 0.05 * math.sin(3 * i)
        twist = 0.10 * math.sin(i + 1)
        x2 = sx * x + twist * y
        y2 = sy * y - twist * x
        z2 = sz * z + 0.04 * d * math.cos(i)
        coords.append((x2, y2, z2))
    return center_atoms(atoms_from_coords(coords))


def random_3d(d: float, seed: int) -> List[Atom]:
    """Random 3D structure with controlled minimum B-B distance."""
    rng = random.Random(seed)
    target_min = 0.72 * d
    # Large box for large d; compact enough for smaller d. Rejection sampling.
    box_half = 1.75 * d
    points: List[Coord] = []
    attempts = 0
    max_attempts = 200000
    while len(points) < N_ATOMS and attempts < max_attempts:
        attempts += 1
        x = rng.uniform(-box_half, box_half)
        y = rng.uniform(-box_half, box_half)
        z = rng.uniform(-box_half, box_half)
        # Avoid accidentally flat clouds.
        if abs(z) < 0.05 * d:
            z += math.copysign(0.12 * d, z if z != 0 else rng.uniform(-1, 1))
        candidate = (x, y, z)
        ok = True
        for px, py, pz in points:
            dd = math.sqrt((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2)
            if dd < target_min:
                ok = False
                break
        if ok:
            points.append(candidate)
    if len(points) != N_ATOMS:
        raise RuntimeError(f"Could not generate random_3d with seed={seed}, d={d}")
    atoms = center_atoms(atoms_from_coords(points))
    return atoms


def build_structure_specs() -> List[StructureSpec]:
    return [
        StructureSpec("icosahedron", icosahedron, "Regular icosahedral B12 vertex shell."),
        StructureSpec("compact_3d", compact_3d, "Compact non-planar B12 cage-like cluster."),
        StructureSpec("bilayer_3d", bilayer_3d, "Two staggered non-planar B6 layers."),
        StructureSpec("cuboctahedral_fragment", cuboctahedral_fragment, "Cuboctahedron-like 12-vertex 3D structure."),
        StructureSpec("distorted_cage", distorted_cage, "Anisotropically deformed 3D cage."),
        StructureSpec("distorted_icosahedron_01", lambda d: distorted_icosahedron(d, seed=101, amplitude=0.025), "Low-symmetry distorted icosahedron variant 01."),
        StructureSpec("distorted_icosahedron_02", lambda d: distorted_icosahedron(d, seed=202, amplitude=0.040), "Low-symmetry distorted icosahedron variant 02."),
        StructureSpec("random_3d_01", lambda d: random_3d(d, seed=RANDOM_SEED + 1 + int(round(d * 10))), "Random 3D B12 cloud with controlled minimum distance 01."),
        StructureSpec("random_3d_02", lambda d: random_3d(d, seed=RANDOM_SEED + 2 + int(round(d * 10))), "Random 3D B12 cloud with controlled minimum distance 02."),
        StructureSpec("random_3d_03", lambda d: random_3d(d, seed=RANDOM_SEED + 3 + int(round(d * 10))), "Random 3D B12 cloud with controlled minimum distance 03."),
    ]


def write_xyz(path: Path, atoms: Sequence[Atom], comment: str) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment.strip() + "\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def orca_input_text(multiplicity: int) -> str:
    """Return a complete ORCA input file for preliminary Opt Freq."""
    return f"""# B12 only. Neutral 3D boron cluster. Preliminary geometry optimization + frequency check.
# {ALTERNATIVE_METHOD_COMMENT}
# Run with: orca input.inp > output.out

! {METHOD} {DISPERSION} {BASIS} Opt Freq TightSCF Grid5 NoAutoStart XYZFile

%pal
nprocs {NPROCS}
end

%maxcore {MAXCORE_MB}

%geom
MaxIter 300
end

%scf
MaxIter 500
end

%output
PrintLevel Normal
Print[P_Cartesian] 1
end

* xyzfile {CHARGE} {multiplicity} start.xyz
"""


def final_refinement_template_text() -> str:
    return f"""# Final refinement template for the best preliminary B12 structure.
# Before running, copy the selected optimized.xyz to this folder as best_start.xyz.
# Use the multiplicity of the selected best preliminary minimum.
# Run with: orca input.inp > output.out

! PBE0 D3BJ def2-TZVP Opt Freq TightSCF Grid5 NoAutoStart XYZFile

%pal
nprocs {NPROCS}
end

%maxcore {MAXCORE_MB}

%geom
MaxIter 500
end

%scf
MaxIter 700
end

%output
PrintLevel Normal
Print[P_Cartesian] 1
end

* xyzfile {CHARGE} REPLACE_WITH_BEST_MULTIPLICITY best_start.xyz
"""


def write_metadata(path: Path, metadata: Dict[str, object]) -> None:
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_summary_template(root: Path) -> None:
    headers = [
        "cluster",
        "structure_name",
        "initial_distance_A",
        "charge",
        "multiplicity",
        "method",
        "basis",
        "total_energy_Eh",
        "relative_energy_eV",
        "terminated_normally",
        "scf_converged",
        "has_imaginary_frequencies",
        "number_of_imaginary_frequencies",
        "min_real_frequency_cm-1",
        "number_of_vibrational_modes",
        "optimized_xyz_path",
        "output_path",
    ]
    (root / "summary_template.csv").write_text(",".join(headers) + "\n", encoding="utf-8")


def write_report_template(root: Path) -> None:
    text = """# Итоговый отчёт: B12 ORCA 6.1 3D workflow

## 1. Цель расчёта
Цель — найти наиболее стабильную 3D-геометрию нейтрального кластера B12 по минимальной полной энергии среди истинных минимумов, подтверждённых частотным расчётом.

## 2. Условия расчёта
- ПО: ORCA 6.1
- ОС: Linux VPS
- CPU: 8 cores
- RAM: 24 GB
- Disk: 400 GB SSD
- Запуск: `orca input.inp > output.out`
- Ресурсы ORCA: `%pal nprocs 8 end`, `%maxcore 2500`
- Режим запуска: последовательный, без параллельного запуска нескольких ORCA jobs.

## 3. Методика
Стартовые геометрии B12 генерируются только в 3D. Для каждой геометрии перебираются расстояния d = 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0 Å и мультиплетности 1, 3, 5. Предварительный уровень: PBE0-D3BJ/def2-SVP Opt Freq. Альтернатива для скрининга: B3LYP-D3BJ/def2-SVP. Финальное уточнение лучшей структуры: PBE0-D3BJ/def2-TZVP Opt Freq.

## 4. Результаты B12
После выполнения `python3 analyze_results.py` вставить таблицу `summary.csv`.

Обязательные поля:
- энергия Eh;
- относительная энергия eV;
- наличие/отсутствие мнимых частот;
- число колебательных мод;
- путь к optimized.xyz;
- путь к output.out.

## 5. Сравнение 3D-структур
Сравнить только реально оптимизированные структуры. Лучшая структура не выбирается по стартовой геометрии: стартовая конфигурация — только начальная точка оптимизации. Сравнение выполняется по финальной полной энергии после оптимизации и по результатам частотного расчёта. Структуры с мнимыми частотами не считаются истинными минимумами.

## 6. Итоговый вывод
Указать структуру с минимальной полной энергией среди истинных минимумов, её мультиплетность, энергию, относительную энергию, 30 нормальных колебательных мод B12, отсутствие мнимых частот и рекомендации для дополнительных расчётов.

## Пояснения
- Нельзя выбирать лучшую структуру только по стартовой геометрии, потому что разные начальные конфигурации могут оптимизироваться в один и тот же минимум или в разные локальные минимумы.
- Оптимизация обязательна, потому что энергия стартовой геометрии не соответствует стационарной точке потенциальной поверхности.
- Частотная проверка обязательна, потому что она отличает минимум от седловой точки.
- Мнимая частота означает отрицательную кривизну вдоль нормальной координаты; такая структура не является устойчивым минимумом.
- Энергии разных мультиплетностей сравниваются только при одинаковом заряде, уровне теории, базисе и после одинаковой процедуры Opt Freq.
- Для нелинейного B12 ожидается 3N − 6 = 30 нормальных колебаний.
- Последовательный запуск нужен, потому что один ORCA job уже использует все 8 CPU cores.
- `%maxcore 2500` означает примерно 2.5 GB на процесс; при 8 процессах это около 20 GB, что безопаснее, чем использовать все 24 GB RAM.
"""
    (root / "report_template_B12.md").write_text(text, encoding="utf-8")


def write_readme(project_root: Path) -> None:
    text = f"""# B12 ORCA 6.1 3D workflow

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
"""
    (project_root / "README_B12_ORCA_workflow.md").write_text(text, encoding="utf-8")


def main() -> None:
    specs = build_structure_specs()
    ROOT.mkdir(parents=True, exist_ok=True)
    write_summary_template(ROOT)
    write_report_template(ROOT)

    # Example XYZ files for direct inspection/copying.
    examples_dir = Path("starter_xyz_examples")
    examples_dir.mkdir(parents=True, exist_ok=True)

    manifest: List[Dict[str, object]] = []

    for spec in specs:
        example_atoms = spec.builder(2.0)
        validate_b12_3d(spec.name, example_atoms)
        write_xyz(
            examples_dir / f"{spec.name}_d_2.0.xyz",
            example_atoms,
            f"B12 {spec.name}; 3D initial geometry; d=2.0 A; generated by generate_inputs.py",
        )

        for d in DISTANCES_A:
            atoms = spec.builder(d)
            validate_b12_3d(spec.name, atoms)
            stats = pair_distance_stats(atoms)
            for mult in MULTIPLICITIES:
                calc_dir = ROOT / spec.name / f"d_{d:.1f}" / f"mult_{mult}"
                calc_dir.mkdir(parents=True, exist_ok=True)
                xyz_comment = (
                    f"B12 {spec.name}; neutral charge={CHARGE}; multiplicity={mult}; "
                    f"3D initial geometry; nominal d={d:.1f} A; {spec.description}"
                )
                write_xyz(calc_dir / "start.xyz", atoms, xyz_comment)
                (calc_dir / "input.inp").write_text(orca_input_text(mult), encoding="utf-8")
                metadata = {
                    "cluster": CLUSTER,
                    "structure_name": spec.name,
                    "description": spec.description,
                    "initial_distance_A": d,
                    "charge": CHARGE,
                    "multiplicity": mult,
                    "method": METHOD,
                    "dispersion": DISPERSION,
                    "basis": BASIS,
                    "nprocs": NPROCS,
                    "maxcore_MB_per_process": MAXCORE_MB,
                    "run_type": "Opt Freq",
                    **stats,
                }
                write_metadata(calc_dir / "metadata.json", metadata)
                manifest.append({"path": str(calc_dir), **metadata})

    # Final refinement folder template. analyze_results.py will fill it automatically for the best candidate.
    final_dir = ROOT / "final_refinement_template"
    final_dir.mkdir(parents=True, exist_ok=True)
    (final_dir / "input.inp").write_text(final_refinement_template_text(), encoding="utf-8")
    (Path("final_refinement_template.inp")).write_text(final_refinement_template_text(), encoding="utf-8")

    write_metadata(ROOT / "manifest.json", {"calculations": manifest})
    write_readme(Path("."))

    print(f"Created {len(manifest)} ORCA input folders under {ROOT}")
    print(f"Created example XYZ files under {examples_dir}")
    print("Next steps:")
    print("  chmod +x run_all.sh")
    print("  tmux new -s b12_orca")
    print("  ./run_all.sh")
    print("  python3 analyze_results.py")


if __name__ == "__main__":
    main()
