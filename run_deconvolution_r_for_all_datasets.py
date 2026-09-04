#!/usr/bin/env python3
"""Run deconvolution_methods.R (CIBERSORTx + marker-based proxy scoring) as
a subprocess for every dataset that has a deduplicated expression matrix,
writing each dataset's output into a deconvolution_results/ subfolder next
to its expression_analysis_results/.

R equivalent of run_deconvolution_for_all_datasets.py, used because it's
often convenient to keep the R and Python deconvolution pipelines callable
side by side. Each dataset is deconvolved against a marker-based binary
signature matrix: PBMC datasets use markers/pbmc_cell_markers.csv, Whole
Blood datasets use markers/blood_cell_markers.csv, mirroring the
marker-file convention already used elsewhere in this project.
"""

from __future__ import annotations

import subprocess
import traceback
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SKIP_DIR_NAMES = {".venv", "node_modules", ".git"}

R_SCRIPT_PATH = PROJECT_ROOT / "deconvolution_methods.R"
PBMC_MARKERS_CSV = PROJECT_ROOT / "markers" / "pbmc_cell_markers.csv"
BLOOD_MARKERS_CSV = PROJECT_ROOT / "markers" / "blood_cell_markers.csv"


def find_dataset_expression_files(project_root: Path) -> list[Path]:
    """Return each deduplicated_expression_dataset.csv under an
    expression_analysis_results/ directory, skipping vendored directories."""
    expression_files = []
    for expr_file in sorted(
        project_root.rglob("expression_analysis_results/deduplicated_expression_dataset.csv")
    ):
        if SKIP_DIR_NAMES & set(expr_file.parts):
            continue
        expression_files.append(expr_file)
    return expression_files


def choose_markers_csv_for_dataset(dataset_dir: Path, project_root: Path) -> Path:
    dataset_label = str(dataset_dir.relative_to(project_root))
    if dataset_label.startswith("PBMC"):
        return PBMC_MARKERS_CSV
    return BLOOD_MARKERS_CSV


def run_deconvolution_r(
    bulk_file: str | Path,
    markers_csv: str | Path,
    output_dir: str | Path,
    r_script_path: str | Path = R_SCRIPT_PATH,
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        str(r_script_path),
        f"--bulk_file={bulk_file}",
        f"--markers_csv={markers_csv}",
        f"--output_dir={output_dir}",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"deconvolution_methods.R failed (exit code {result.returncode}):\n{result.stderr}"
        )

    expected = output_dir / "cibersortx_cell_type_proportions.csv"
    if not expected.exists():
        raise FileNotFoundError(
            f"Expected R script to produce {expected}, but it wasn't found. "
            f"R stderr:\n{result.stderr}"
        )


def run_deconvolution_r_for_all_datasets(
    project_root: str | Path = PROJECT_ROOT,
) -> dict[str, str]:
    project_root = Path(project_root)
    expression_files = find_dataset_expression_files(project_root)
    if not expression_files:
        print("No expression_analysis_results/deduplicated_expression_dataset.csv files were found.")
        return {}

    outcomes: dict[str, str] = {}
    for index, expression_file in enumerate(expression_files):
        dataset_dir = expression_file.parent.parent
        dataset_label = str(dataset_dir.relative_to(project_root))
        output_dir = dataset_dir / "deconvolution_results"
        markers_csv = choose_markers_csv_for_dataset(dataset_dir, project_root)

        print(f"[{index + 1}/{len(expression_files)}] Running R deconvolution for {dataset_label}")
        try:
            run_deconvolution_r(
                bulk_file=expression_file,
                markers_csv=markers_csv,
                output_dir=output_dir,
            )
            outcomes[dataset_label] = "success"
        except Exception:
            outcomes[dataset_label] = "failed"
            print(f"R deconvolution failed for {dataset_label}:")
            traceback.print_exc()

    succeeded = [label for label, outcome in outcomes.items() if outcome == "success"]
    failed = [label for label, outcome in outcomes.items() if outcome == "failed"]
    print(f"\nR deconvolution complete: {len(succeeded)} succeeded, {len(failed)} failed.")
    if failed:
        print("Failed datasets:")
        for label in failed:
            print(f"- {label}")

    return outcomes


if __name__ == "__main__":
    run_deconvolution_r_for_all_datasets()
