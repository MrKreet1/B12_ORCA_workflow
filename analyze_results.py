#!/usr/bin/env python3
"""
analyze_results.py

Parser and summarizer for the B12-only ORCA 6.1 workflow.

It creates:
- calculations/B12/summary.csv
- calculations/B12/best_structure_B12.txt
- optimized.xyz in each completed calculation folder, when coordinates are available
- frequencies.txt in each completed calculation folder, when frequencies are available
- normal_modes_raw.txt in each completed calculation folder, when normal mode text is available
- calculations/B12/final_refinement/input.inp and best_start.xyz, when a true minimum is found

No energies, frequencies, coordinates, spectra, or normal modes are invented.
Values are extracted only from existing ORCA output / Hessian files.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = Path("calculations") / "B12"
CLUSTER = "B12"
N_ATOMS = 12
EXPECTED_VIBRATIONAL_MODES = 3 * N_ATOMS - 6
HARTREE_TO_EV = 27.211386245988
IMAGINARY_THRESHOLD_CM = -1.0
ZERO_MODE_ABS_THRESHOLD_CM = 1.0

SUMMARY_HEADERS = [
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

Atom = Tuple[str, float, float, float]


@dataclass
class ModeData:
    mode_number: int
    frequency_cm: Optional[float] = None
    ir_intensity: Optional[float] = None
    reduced_mass: Optional[float] = None
    force_constant: Optional[float] = None
    displacements: Optional[List[Tuple[float, float, float]]] = None


@dataclass
class CalcResult:
    path: Path
    output_path: Path
    structure_name: str = ""
    initial_distance_A: Optional[float] = None
    charge: Optional[int] = None
    multiplicity: Optional[int] = None
    method: str = ""
    basis: str = ""
    total_energy_Eh: Optional[float] = None
    relative_energy_eV: Optional[float] = None
    terminated_normally: bool = False
    scf_converged: bool = False
    has_imaginary_frequencies: bool = False
    number_of_imaginary_frequencies: int = 0
    min_real_frequency_cm: Optional[float] = None
    number_of_vibrational_modes: int = 0
    optimized_xyz_path: str = ""
    frequencies: List[float] = field(default_factory=list)
    vibrational_frequencies: List[float] = field(default_factory=list)
    modes: List[ModeData] = field(default_factory=list)
    ir_by_mode: Dict[int, float] = field(default_factory=dict)
    reduced_masses_by_mode: Dict[int, float] = field(default_factory=dict)
    force_constants_by_mode: Dict[int, float] = field(default_factory=dict)
    normal_modes_raw: str = ""
    disk_usage_bytes: int = 0
    parse_notes: List[str] = field(default_factory=list)

    def is_true_minimum(self) -> bool:
        return (
            self.terminated_normally
            and self.scf_converged
            and self.total_energy_Eh is not None
            and not self.has_imaginary_frequencies
            and self.number_of_vibrational_modes in (0, EXPECTED_VIBRATIONAL_MODES)
        )


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def safe_float(text: str) -> Optional[float]:
    try:
        return float(text.replace("D", "E").replace("d", "E"))
    except Exception:
        return None


def parse_path_metadata(calc_dir: Path) -> Dict[str, object]:
    """Parse structure/distance/multiplicity from calculations/B12/<structure>/d_x/mult_y."""
    parts = calc_dir.parts
    metadata: Dict[str, object] = {}
    try:
        idx = parts.index("B12")
        metadata["structure_name"] = parts[idx + 1]
        metadata["initial_distance_A"] = float(parts[idx + 2].replace("d_", ""))
        metadata["multiplicity"] = int(parts[idx + 3].replace("mult_", ""))
    except Exception:
        pass

    meta_path = calc_dir / "metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(read_text(meta_path))
            metadata.update(meta)
        except Exception:
            pass
    return metadata


def parse_input_metadata(input_path: Path) -> Dict[str, object]:
    meta: Dict[str, object] = {"charge": None, "multiplicity": None, "method": "", "basis": ""}
    if not input_path.exists():
        return meta
    text = read_text(input_path)
    simple_lines = [line.strip() for line in text.splitlines() if line.strip().startswith("!")]
    if simple_lines:
        tokens: List[str] = []
        for line in simple_lines:
            tokens.extend(line[1:].split())
        # Keep method and basis simple and explicit for summary.
        upper = [t.upper() for t in tokens]
        for candidate in ["PBE0", "B3LYP", "R2SCAN-3C", "BP86", "TPSS"]:
            if candidate.upper() in upper:
                meta["method"] = candidate
                break
        for token in tokens:
            if token.lower().startswith("def2-"):
                meta["basis"] = token
                break
    m = re.search(r"^\s*\*\s*xyzfile\s+(-?\d+)\s+(\d+)\s+\S+", text, re.IGNORECASE | re.MULTILINE)
    if not m:
        m = re.search(r"^\s*\*\s*xyz\s+(-?\d+)\s+(\d+)\s*$", text, re.IGNORECASE | re.MULTILINE)
    if m:
        meta["charge"] = int(m.group(1))
        meta["multiplicity"] = int(m.group(2))
    return meta


def parse_final_energy(text: str) -> Optional[float]:
    energies = re.findall(r"FINAL\s+SINGLE\s+POINT\s+ENERGY\s+(-?\d+\.\d+(?:[EeDd][+-]?\d+)?)", text)
    if not energies:
        return None
    return safe_float(energies[-1])


def parse_scf_converged(text: str) -> bool:
    if re.search(r"SCF\s+NOT\s+CONVERGED", text, re.IGNORECASE):
        return False
    patterns = [
        r"SCF\s+CONVERGED\s+AFTER",
        r"SCF\s+CONVERGENCE\s+HAS\s+BEEN\s+ACHIEVED",
        r"TOTAL\s+SCF\s+TIME",
    ]
    return any(re.search(p, text, re.IGNORECASE) for p in patterns)


def parse_cartesian_blocks(text: str) -> List[List[Atom]]:
    """Extract all CARTESIAN COORDINATES (ANGSTROEM) blocks from ORCA output."""
    lines = text.splitlines()
    blocks: List[List[Atom]] = []
    i = 0
    header_re = re.compile(r"CARTESIAN\s+COORDINATES\s+\(ANGSTROEM\)", re.IGNORECASE)
    atom_re = re.compile(
        r"^\s*([A-Z][a-z]?)\s+(-?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)\s+"
        r"(-?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)\s+"
        r"(-?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)"
    )
    while i < len(lines):
        if header_re.search(lines[i]):
            i += 1
            block: List[Atom] = []
            while i < len(lines):
                line = lines[i]
                if not line.strip() or set(line.strip()) <= {"-"}:
                    i += 1
                    # Skip leading blank/dash lines before coordinates.
                    if not block:
                        continue
                    break
                m = atom_re.match(line)
                if m:
                    vals = [safe_float(m.group(k)) for k in (2, 3, 4)]
                    if all(v is not None for v in vals):
                        block.append((m.group(1), vals[0], vals[1], vals[2]))  # type: ignore[arg-type]
                elif block:
                    break
                i += 1
            if block:
                blocks.append(block)
        else:
            i += 1
    return blocks


def parse_orca_xyz_file(calc_dir: Path) -> Optional[List[Atom]]:
    """Read likely ORCA-generated final XYZ files, if present."""
    candidates = [calc_dir / "input.xyz"] + sorted(calc_dir.glob("*.xyz"))
    for path in candidates:
        if not path.exists() or path.name == "start.xyz" or path.name == "optimized.xyz":
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            if len(lines) < N_ATOMS + 2:
                continue
            n = int(lines[0].strip())
            if n != N_ATOMS:
                continue
            atoms: List[Atom] = []
            for line in lines[2 : 2 + n]:
                parts = line.split()
                if len(parts) >= 4 and parts[0] == "B":
                    vals = [safe_float(parts[1]), safe_float(parts[2]), safe_float(parts[3])]
                    if all(v is not None for v in vals):
                        atoms.append((parts[0], vals[0], vals[1], vals[2]))  # type: ignore[arg-type]
            if len(atoms) == N_ATOMS:
                return atoms
        except Exception:
            continue
    return None


def write_xyz(path: Path, atoms: Sequence[Atom], comment: str) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment.strip() + "\n")
        for el, x, y, z in atoms:
            f.write(f"{el:<2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")


def extract_optimized_xyz(calc_dir: Path, output_text: str) -> Tuple[str, List[str]]:
    notes: List[str] = []
    atoms = parse_orca_xyz_file(calc_dir)
    if atoms is None:
        blocks = parse_cartesian_blocks(output_text)
        b12_blocks = [b for b in blocks if len(b) == N_ATOMS and all(a[0] == "B" for a in b)]
        if b12_blocks:
            atoms = b12_blocks[-1]
        elif blocks:
            notes.append(f"Cartesian coordinate blocks found, but no complete {N_ATOMS}-atom B-only block.")
    if atoms is None:
        return "", notes
    out = calc_dir / "optimized.xyz"
    write_xyz(out, atoms, f"Optimized coordinates extracted from {calc_dir / 'output.out'}")
    return str(out), notes


def parse_frequencies_from_output(text: str) -> List[float]:
    lines = text.splitlines()
    freqs: List[float] = []
    in_section = False
    for line in lines:
        if "VIBRATIONAL FREQUENCIES" in line.upper():
            in_section = True
            continue
        if in_section:
            upper = line.upper()
            if "NORMAL MODES" in upper or "IR SPECTRUM" in upper or "THERMOCHEMISTRY" in upper:
                break
            m = re.match(r"\s*\d+\s*:\s*(-?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?)\s*cm", line, re.IGNORECASE)
            if m:
                val = safe_float(m.group(1))
                if val is not None:
                    freqs.append(val)
    return freqs


def parse_simple_table_section(text: str, title_patterns: Sequence[str]) -> Dict[int, List[float]]:
    """Extract numeric rows from a named text table. Heuristic but robust."""
    lines = text.splitlines()
    table: Dict[int, List[float]] = {}
    in_section = False
    title_re = re.compile("|".join(title_patterns), re.IGNORECASE)
    next_section_re = re.compile(
        r"^(\s*-{5,}\s*$|\s*[A-Z][A-Z0-9 /_()\-]{8,}\s*$)",
        re.IGNORECASE,
    )
    for line in lines:
        if title_re.search(line):
            in_section = True
            continue
        if not in_section:
            continue
        if next_section_re.match(line) and table:
            break
        nums = re.findall(r"-?\d+(?:\.\d*)?(?:[EeDd][+-]?\d+)?", line)
        if len(nums) >= 2:
            # First number in ORCA tables is often the mode index.
            idx_float = safe_float(nums[0])
            if idx_float is None:
                continue
            idx = int(round(idx_float))
            vals: List[float] = []
            for n in nums[1:]:
                v = safe_float(n)
                if v is not None:
                    vals.append(v)
            if vals:
                table[idx] = vals
    return table


def extract_hess_section(hess_text: str, section_name: str) -> str:
    pattern = re.compile(rf"^\${re.escape(section_name)}\s*$", re.IGNORECASE | re.MULTILINE)
    m = pattern.search(hess_text)
    if not m:
        return ""
    start = m.end()
    next_m = re.search(r"^\$[A-Za-z_]+\s*$", hess_text[start:], re.MULTILINE)
    end = start + next_m.start() if next_m else len(hess_text)
    return hess_text[start:end].strip("\n")


def parse_hess_frequencies(hess_text: str) -> List[float]:
    section = extract_hess_section(hess_text, "vibrational_frequencies")
    if not section:
        return []
    lines = [ln.strip() for ln in section.splitlines() if ln.strip()]
    if not lines:
        return []
    freqs: List[float] = []
    for line in lines[1:]:  # first line usually count
        parts = line.split()
        if len(parts) >= 2:
            v = safe_float(parts[-1])
            if v is not None:
                freqs.append(v)
    return freqs


def parse_hess_ir(hess_text: str) -> Dict[int, float]:
    section = extract_hess_section(hess_text, "ir_spectrum")
    result: Dict[int, float] = {}
    if not section:
        return result
    lines = [ln.strip() for ln in section.splitlines() if ln.strip()]
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 3:
            idx_val = safe_float(parts[0])
            # ORCA .hess IR columns may include frequency, intensity and transition dipole components.
            # Store the second numeric column after index as a conservative intensity-like value.
            intensity_val = safe_float(parts[2]) if len(parts) >= 3 else None
            if idx_val is not None and intensity_val is not None:
                result[int(round(idx_val))] = intensity_val
    return result


def parse_hess_raw_normal_modes(calc_dir: Path) -> str:
    hess_files = sorted(calc_dir.glob("*.hess"))
    for hp in hess_files:
        text = read_text(hp)
        section = extract_hess_section(text, "normal_modes")
        if section:
            return f"# Source: {hp}\n$normal_modes\n{section}\n"
    return ""


def parse_frequency_related_data(calc_dir: Path, output_text: str) -> Tuple[List[float], Dict[int, float], Dict[int, float], Dict[int, float], str, List[str]]:
    notes: List[str] = []
    freqs = parse_frequencies_from_output(output_text)

    hess_files = sorted(calc_dir.glob("*.hess"))
    hess_text = ""
    if hess_files:
        try:
            hess_text = read_text(hess_files[0])
        except Exception as exc:
            notes.append(f"Could not read Hessian file {hess_files[0]}: {exc}")

    if not freqs and hess_text:
        freqs = parse_hess_frequencies(hess_text)
        if freqs:
            notes.append("Frequencies were read from .hess file.")

    ir_by_mode: Dict[int, float] = {}
    reduced_by_mode: Dict[int, float] = {}
    force_by_mode: Dict[int, float] = {}

    # Output IR table heuristic.
    ir_table = parse_simple_table_section(output_text, [r"IR\s+SPECTRUM"])
    for idx, vals in ir_table.items():
        # In many ORCA outputs, an IR table row has frequency and intensity among the first columns.
        if len(vals) >= 2:
            ir_by_mode[idx] = vals[1]
        elif vals:
            ir_by_mode[idx] = vals[-1]

    # Thermochemistry/frequency block may print reduced masses and force constants in some ORCA modes.
    # These regexes intentionally accept several common headings.
    red_table = parse_simple_table_section(output_text, [r"REDUCED\s+MASSES", r"RED\.\s*MASSES"])
    for idx, vals in red_table.items():
        if vals:
            reduced_by_mode[idx] = vals[-1]

    force_table = parse_simple_table_section(output_text, [r"FORCE\s+CONSTANTS", r"FORCE\s+CONST"])
    for idx, vals in force_table.items():
        if vals:
            force_by_mode[idx] = vals[-1]

    if hess_text:
        hess_ir = parse_hess_ir(hess_text)
        # Do not overwrite output-derived values unless absent.
        for k, v in hess_ir.items():
            ir_by_mode.setdefault(k, v)

    normal_modes_raw = ""
    # Prefer explicitly printed normal mode section from output.
    m = re.search(r"(NORMAL\s+MODES.*?)(?:\n\s*[A-Z][A-Z0-9 /_()\-]{8,}\s*\n|\Z)", output_text, re.IGNORECASE | re.DOTALL)
    if m:
        normal_modes_raw = m.group(1).strip()
    if not normal_modes_raw:
        normal_modes_raw = parse_hess_raw_normal_modes(calc_dir)

    return freqs, ir_by_mode, reduced_by_mode, force_by_mode, normal_modes_raw, notes


def select_vibrational_frequencies(all_freqs: Sequence[float]) -> List[float]:
    """
    ORCA frequency lists may include 3N entries with 6 translational/rotational near-zero modes.
    For nonlinear B12, report 30 vibrational modes. If only 30 are printed, keep them.
    """
    freqs = list(all_freqs)
    if len(freqs) >= 3 * N_ATOMS:
        return freqs[-EXPECTED_VIBRATIONAL_MODES:]
    if len(freqs) > EXPECTED_VIBRATIONAL_MODES:
        # Drop the smallest-by-absolute-value modes until 30 remain.
        indexed = list(enumerate(freqs))
        indexed_sorted = sorted(indexed, key=lambda p: abs(p[1]))
        drop = {idx for idx, _ in indexed_sorted[: len(freqs) - EXPECTED_VIBRATIONAL_MODES]}
        return [v for idx, v in indexed if idx not in drop]
    return freqs


def write_frequencies_file(calc_dir: Path, freqs: Sequence[float], vib_freqs: Sequence[float], result: CalcResult) -> None:
    out = calc_dir / "frequencies.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write(f"Source output: {calc_dir / 'output.out'}\n")
        f.write(f"All frequencies parsed: {len(freqs)}\n")
        f.write(f"Vibrational frequencies used for B12 check: {len(vib_freqs)}\n")
        f.write(f"Expected nonlinear B12 vibrational modes: {EXPECTED_VIBRATIONAL_MODES}\n\n")
        f.write("# mode_number,frequency_cm-1,imaginary,ir_intensity,reduced_mass,force_constant\n")
        for i, freq in enumerate(vib_freqs, start=1):
            mode_index_guess = i + max(0, len(freqs) - len(vib_freqs)) - 1
            ir = result.ir_by_mode.get(mode_index_guess, result.ir_by_mode.get(i, None))
            rm = result.reduced_masses_by_mode.get(mode_index_guess, result.reduced_masses_by_mode.get(i, None))
            fc = result.force_constants_by_mode.get(mode_index_guess, result.force_constants_by_mode.get(i, None))
            f.write(
                f"{i},{freq:.8f},{str(freq < IMAGINARY_THRESHOLD_CM).lower()},"
                f"{'' if ir is None else f'{ir:.8f}'},"
                f"{'' if rm is None else f'{rm:.8f}'},"
                f"{'' if fc is None else f'{fc:.8f}'}\n"
            )


def directory_size_bytes(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            fp = Path(root) / name
            try:
                total += fp.stat().st_size
            except OSError:
                pass
    return total


def parse_calculation(calc_dir: Path) -> CalcResult:
    output_path = calc_dir / "output.out"
    result = CalcResult(path=calc_dir, output_path=output_path)

    path_meta = parse_path_metadata(calc_dir)
    input_meta = parse_input_metadata(calc_dir / "input.inp")

    result.structure_name = str(path_meta.get("structure_name", ""))
    result.initial_distance_A = path_meta.get("initial_distance_A") if isinstance(path_meta.get("initial_distance_A"), float) else None
    result.charge = input_meta.get("charge") if isinstance(input_meta.get("charge"), int) else None
    result.multiplicity = input_meta.get("multiplicity") if isinstance(input_meta.get("multiplicity"), int) else None
    if result.multiplicity is None and isinstance(path_meta.get("multiplicity"), int):
        result.multiplicity = int(path_meta["multiplicity"])
    result.method = str(input_meta.get("method") or path_meta.get("method") or "")
    result.basis = str(input_meta.get("basis") or path_meta.get("basis") or "")
    result.disk_usage_bytes = directory_size_bytes(calc_dir)

    if not output_path.exists():
        result.parse_notes.append("output.out not found")
        return result

    text = read_text(output_path)
    result.terminated_normally = "ORCA TERMINATED NORMALLY" in text
    result.scf_converged = parse_scf_converged(text)
    result.total_energy_Eh = parse_final_energy(text)

    opt_path, coord_notes = extract_optimized_xyz(calc_dir, text)
    result.optimized_xyz_path = opt_path
    result.parse_notes.extend(coord_notes)

    freqs, ir, red, force, normal_modes_raw, freq_notes = parse_frequency_related_data(calc_dir, text)
    result.frequencies = freqs
    result.vibrational_frequencies = select_vibrational_frequencies(freqs)
    result.ir_by_mode = ir
    result.reduced_masses_by_mode = red
    result.force_constants_by_mode = force
    result.normal_modes_raw = normal_modes_raw
    result.parse_notes.extend(freq_notes)

    result.number_of_vibrational_modes = len(result.vibrational_frequencies)
    imag = [v for v in result.vibrational_frequencies if v < IMAGINARY_THRESHOLD_CM]
    result.number_of_imaginary_frequencies = len(imag)
    result.has_imaginary_frequencies = len(imag) > 0
    real_freqs = [v for v in result.vibrational_frequencies if v > ZERO_MODE_ABS_THRESHOLD_CM]
    result.min_real_frequency_cm = min(real_freqs) if real_freqs else None

    for i, freq in enumerate(result.vibrational_frequencies, start=1):
        mode_index_guess = i + max(0, len(result.frequencies) - len(result.vibrational_frequencies)) - 1
        result.modes.append(
            ModeData(
                mode_number=i,
                frequency_cm=freq,
                ir_intensity=result.ir_by_mode.get(mode_index_guess, result.ir_by_mode.get(i)),
                reduced_mass=result.reduced_masses_by_mode.get(mode_index_guess, result.reduced_masses_by_mode.get(i)),
                force_constant=result.force_constants_by_mode.get(mode_index_guess, result.force_constants_by_mode.get(i)),
            )
        )

    if result.frequencies:
        write_frequencies_file(calc_dir, result.frequencies, result.vibrational_frequencies, result)
    if result.normal_modes_raw:
        (calc_dir / "normal_modes_raw.txt").write_text(result.normal_modes_raw + "\n", encoding="utf-8")

    if result.number_of_vibrational_modes not in (0, EXPECTED_VIBRATIONAL_MODES):
        result.parse_notes.append(
            f"Expected {EXPECTED_VIBRATIONAL_MODES} vibrational modes for nonlinear B12; parsed {result.number_of_vibrational_modes}."
        )

    return result


def format_optional_float(value: Optional[float], ndigits: int = 10) -> str:
    if value is None:
        return ""
    return f"{value:.{ndigits}f}"


def result_to_row(result: CalcResult) -> Dict[str, object]:
    return {
        "cluster": CLUSTER,
        "structure_name": result.structure_name,
        "initial_distance_A": "" if result.initial_distance_A is None else f"{result.initial_distance_A:.1f}",
        "charge": "" if result.charge is None else result.charge,
        "multiplicity": "" if result.multiplicity is None else result.multiplicity,
        "method": result.method,
        "basis": result.basis,
        "total_energy_Eh": format_optional_float(result.total_energy_Eh, 12),
        "relative_energy_eV": format_optional_float(result.relative_energy_eV, 8),
        "terminated_normally": result.terminated_normally,
        "scf_converged": result.scf_converged,
        "has_imaginary_frequencies": result.has_imaginary_frequencies,
        "number_of_imaginary_frequencies": result.number_of_imaginary_frequencies,
        "min_real_frequency_cm-1": format_optional_float(result.min_real_frequency_cm, 4),
        "number_of_vibrational_modes": result.number_of_vibrational_modes,
        "optimized_xyz_path": result.optimized_xyz_path,
        "output_path": str(result.output_path),
    }


def write_summary(results: Sequence[CalcResult], root: Path) -> None:
    path = root / "summary.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_HEADERS)
        writer.writeheader()
        for r in results:
            writer.writerow(result_to_row(r))


def load_xyz_text(path: str) -> str:
    if not path:
        return ""
    p = Path(path)
    if not p.exists():
        return ""
    return p.read_text(encoding="utf-8", errors="replace")


def create_final_refinement_input(best: CalcResult, root: Path) -> None:
    if not best.optimized_xyz_path:
        return
    final_dir = root / "final_refinement"
    final_dir.mkdir(parents=True, exist_ok=True)
    dst_xyz = final_dir / "best_start.xyz"
    shutil.copyfile(best.optimized_xyz_path, dst_xyz)
    mult = best.multiplicity if best.multiplicity is not None else 1
    input_text = f"""# Final refinement for the best preliminary B12 true minimum.
# Source calculation: {best.path}
# Preliminary energy: {best.total_energy_Eh:.12f} Eh
# Preliminary relative energy: {best.relative_energy_eV:.8f} eV
# Run with: orca input.inp > output.out

! PBE0 D3BJ def2-TZVP Opt Freq TightSCF Grid5 NoAutoStart XYZFile

%pal
nprocs 8
end

%maxcore 2500

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

* xyzfile 0 {mult} best_start.xyz
"""
    (final_dir / "input.inp").write_text(input_text, encoding="utf-8")
    metadata = {
        "purpose": "final PBE0-D3BJ/def2-TZVP Opt Freq refinement",
        "source_calculation": str(best.path),
        "source_optimized_xyz": best.optimized_xyz_path,
        "charge": best.charge,
        "multiplicity": best.multiplicity,
        "preliminary_total_energy_Eh": best.total_energy_Eh,
        "preliminary_relative_energy_eV": best.relative_energy_eV,
    }
    (final_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_best_report(results: Sequence[CalcResult], root: Path) -> None:
    best_candidates = [r for r in results if r.is_true_minimum()]
    best = min(best_candidates, key=lambda r: r.total_energy_Eh) if best_candidates else None

    out = root / "best_structure_B12.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write("B12 best structure report\n")
        f.write("=========================\n\n")
        f.write("Selection rule: lowest total energy among true minima only.\n")
        f.write("True minimum: ORCA terminated normally, SCF converged, no imaginary frequencies.\n")
        f.write(f"Expected nonlinear B12 vibrational modes: {EXPECTED_VIBRATIONAL_MODES}\n\n")

        if best is None:
            f.write("No true minimum was found among parsed calculations.\n")
            f.write("Possible reasons: calculations are not finished, SCF did not converge, output parsing failed, or all completed structures have imaginary frequencies.\n\n")
            f.write("Completed calculations with energies:\n")
            for r in sorted([x for x in results if x.total_energy_Eh is not None], key=lambda x: x.total_energy_Eh):
                f.write(
                    f"- {r.path}: E={r.total_energy_Eh:.12f} Eh, "
                    f"imag={r.number_of_imaginary_frequencies}, modes={r.number_of_vibrational_modes}, "
                    f"terminated={r.terminated_normally}, scf={r.scf_converged}\n"
                )
            return

        create_final_refinement_input(best, root)

        f.write("Best preliminary B12 true minimum\n")
        f.write("---------------------------------\n")
        f.write(f"Structure name: {best.structure_name}\n")
        f.write(f"Calculation path: {best.path}\n")
        f.write(f"Initial distance d: {best.initial_distance_A} A\n")
        f.write(f"Charge: {best.charge}\n")
        f.write(f"Multiplicity: {best.multiplicity}\n")
        f.write(f"Method: {best.method}\n")
        f.write(f"Basis: {best.basis}\n")
        f.write(f"Total energy: {best.total_energy_Eh:.12f} Eh\n")
        f.write(f"Relative energy: {best.relative_energy_eV:.8f} eV\n")
        f.write(f"ORCA terminated normally: {best.terminated_normally}\n")
        f.write(f"SCF converged: {best.scf_converged}\n")
        f.write(f"Imaginary frequencies present: {best.has_imaginary_frequencies}\n")
        f.write(f"Number of imaginary frequencies: {best.number_of_imaginary_frequencies}\n")
        f.write(f"Number of vibrational modes parsed: {best.number_of_vibrational_modes}\n")
        f.write(f"Optimized XYZ path: {best.optimized_xyz_path}\n")
        f.write(f"Output path: {best.output_path}\n")
        f.write(f"Disk used by calculation folder: {best.disk_usage_bytes / (1024 ** 2):.2f} MB\n")
        f.write("Final refinement input created at: calculations/B12/final_refinement/input.inp\n\n")

        f.write("Optimized coordinates XYZ\n")
        f.write("-------------------------\n")
        xyz_text = load_xyz_text(best.optimized_xyz_path)
        f.write(xyz_text if xyz_text else "Optimized coordinates were not extracted.\n")
        if xyz_text and not xyz_text.endswith("\n"):
            f.write("\n")
        f.write("\n")

        f.write("Vibrational modes\n")
        f.write("-----------------\n")
        f.write("mode, frequency_cm-1, IR_intensity_if_available, reduced_mass_if_available, force_constant_if_available\n")
        for mode in best.modes:
            f.write(
                f"{mode.mode_number:3d} "
                f"{'' if mode.frequency_cm is None else f'{mode.frequency_cm:14.6f}'} "
                f"IR={'' if mode.ir_intensity is None else f'{mode.ir_intensity:.8f}'} "
                f"red_mass={'' if mode.reduced_mass is None else f'{mode.reduced_mass:.8f}'} "
                f"force_const={'' if mode.force_constant is None else f'{mode.force_constant:.8f}'}\n"
            )
        f.write("\n")

        if best.normal_modes_raw:
            f.write("Normal-mode displacement data\n")
            f.write("-----------------------------\n")
            f.write("Raw normal-mode data were extracted to normal_modes_raw.txt in the calculation folder.\n")
            f.write(f"Path: {best.path / 'normal_modes_raw.txt'}\n\n")
        else:
            f.write("Normal-mode displacement data were not found in output/.hess files.\n\n")

        if best.parse_notes:
            f.write("Parser notes\n")
            f.write("------------\n")
            for note in best.parse_notes:
                f.write(f"- {note}\n")
            f.write("\n")

        f.write("Conclusion\n")
        f.write("----------\n")
        f.write("This structure is the best preliminary true minimum among the parsed PBE0-D3BJ/def2-SVP Opt Freq calculations. ")
        f.write("It must still be refined with the generated PBE0-D3BJ/def2-TZVP Opt Freq input before final reporting.\n")


def write_disk_report(results: Sequence[CalcResult], root: Path) -> None:
    total = directory_size_bytes(root)
    out = root / "disk_usage_report.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write("Disk usage report for B12 calculations\n")
        f.write("======================================\n")
        f.write(f"Total under {root}: {total / (1024 ** 3):.3f} GB\n")
        f.write("VPS disk reference: 400 GB SSD\n\n")
        for r in sorted(results, key=lambda x: x.disk_usage_bytes, reverse=True):
            f.write(f"{r.disk_usage_bytes / (1024 ** 2):10.2f} MB  {r.path}\n")


def main() -> None:
    if not ROOT.exists():
        raise SystemExit("calculations/B12 not found. Run python3 generate_inputs.py first.")

    input_files = sorted(p for p in ROOT.rglob("input.inp") if "final_refinement_template" not in p.parts and "final_refinement" not in p.parts)
    if not input_files:
        raise SystemExit("No input.inp files found under calculations/B12. Run python3 generate_inputs.py first.")

    results: List[CalcResult] = []
    for input_path in input_files:
        results.append(parse_calculation(input_path.parent))

    # Relative energies: zero is the lowest true minimum if available; otherwise the lowest normally terminated converged energy.
    true_minima = [r for r in results if r.is_true_minimum()]
    if true_minima:
        reference_energy = min(r.total_energy_Eh for r in true_minima if r.total_energy_Eh is not None)
    else:
        finished = [r for r in results if r.terminated_normally and r.scf_converged and r.total_energy_Eh is not None]
        reference_energy = min((r.total_energy_Eh for r in finished), default=None)
    if reference_energy is not None:
        for r in results:
            if r.total_energy_Eh is not None:
                r.relative_energy_eV = (r.total_energy_Eh - reference_energy) * HARTREE_TO_EV

    write_summary(results, ROOT)
    write_best_report(results, ROOT)
    write_disk_report(results, ROOT)

    print(f"Parsed {len(results)} calculation folders.")
    print(f"Wrote {ROOT / 'summary.csv'}")
    print(f"Wrote {ROOT / 'best_structure_B12.txt'}")
    print(f"Wrote {ROOT / 'disk_usage_report.txt'}")
    if true_minima:
        best = min(true_minima, key=lambda r: r.total_energy_Eh)
        print(f"Best preliminary true minimum: {best.path}")
        print(f"Energy: {best.total_energy_Eh:.12f} Eh")
        print(f"Final refinement input: {ROOT / 'final_refinement' / 'input.inp'}")
    else:
        print("No true minimum found yet. Finish ORCA jobs or inspect failed/imaginary-frequency cases.")


if __name__ == "__main__":
    main()
