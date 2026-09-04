"""Thin wrapper that runs qc_normalize.R as a subprocess and returns the
path to the normalized/batch-corrected expression CSV it produces.

Drop-in replacement for qc_normalize.run_qc_normalization(), used because
the ComBat step was too memory-hungry / unstable in the Python 'combat'
package. R's sva::ComBat is the reference implementation and is generally
more memory-efficient.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

# Adjust this if qc_normalize.R lives somewhere else relative to your project.
R_SCRIPT_PATH = Path(__file__).resolve().parent / "qc_normalize.R"


def run_qc_normalization(
    expression_file: str | Path,
    metadata_file: str | Path,
    output_dir: str | Path,
    metadata_sample_column: str,
    phenotype_column: str,
    batch_column: str,
    r_script_path: str | Path = R_SCRIPT_PATH,
) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        str(r_script_path),
        f"--expression_file={expression_file}",
        f"--metadata_file={metadata_file}",
        f"--output_dir={output_dir}",
        f"--metadata_sample_column={metadata_sample_column}",
        f"--phenotype_column={phenotype_column}",
        f"--batch_column={batch_column}",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"qc_normalize.R failed (exit code {result.returncode}):\n{result.stderr}"
        )

    output_path = output_dir / "qc_normalized_expression.csv"
    if not output_path.exists():
        raise FileNotFoundError(
            f"Expected R script to produce {output_path}, but it wasn't found. "
            f"R stderr:\n{result.stderr}"
        )
    return output_path
