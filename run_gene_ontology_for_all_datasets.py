#!/usr/bin/env python3
"""Run GO/KEGG enrichment for every dataset that has differential-expression
results, writing each dataset's output into a gene_ontology/ subfolder next
to its expression_analysis_results/.

This is intentionally kept separate from each dataset's main.py pipeline: it
calls the public Enrichr API, which is rate-limited, so a batch run needs a
delay between datasets and should keep going (and report) if one dataset
fails rather than aborting the rest.
"""

from __future__ import annotations

import time
import traceback
from pathlib import Path

from analyze_gene_ontology import run_gene_ontology_analysis

PROJECT_ROOT = Path(__file__).resolve().parent
SKIP_DIR_NAMES = {".venv", "node_modules", ".git"}
DELAY_BETWEEN_DATASETS_SECONDS = 5.0


def find_dataset_result_dirs(project_root: Path) -> list[Path]:
    """Return each expression_analysis_results/ directory that has both
    top_up_genes.csv and top_down_genes.csv, skipping vendored directories."""
    result_dirs = []
    for up_genes_path in sorted(project_root.rglob("expression_analysis_results/top_up_genes.csv")):
        if SKIP_DIR_NAMES & set(up_genes_path.parts):
            continue
        results_dir = up_genes_path.parent
        if (results_dir / "top_down_genes.csv").exists():
            result_dirs.append(results_dir)
    return result_dirs


def run_gene_ontology_for_all_datasets(
    project_root: str | Path = PROJECT_ROOT,
    delay_between_datasets_seconds: float = DELAY_BETWEEN_DATASETS_SECONDS,
) -> dict[str, str]:
    project_root = Path(project_root)
    result_dirs = find_dataset_result_dirs(project_root)
    if not result_dirs:
        print("No expression_analysis_results directories with up/down gene lists were found.")
        return {}

    outcomes: dict[str, str] = {}
    for index, results_dir in enumerate(result_dirs):
        dataset_dir = results_dir.parent
        dataset_label = str(dataset_dir.relative_to(project_root))
        output_dir = dataset_dir / "gene_ontology"

        background_file = results_dir / "deduplicated_expression_dataset.csv"

        print(f"[{index + 1}/{len(result_dirs)}] Running gene ontology analysis for {dataset_label}")
        try:
            run_gene_ontology_analysis(
                up_genes_file=results_dir / "top_up_genes.csv",
                down_genes_file=results_dir / "top_down_genes.csv",
                background_file=background_file if background_file.exists() else None,
                background_gene_column="Gene Symbol",
                output_dir=output_dir,
            )
            outcomes[dataset_label] = "success"
        except Exception:
            outcomes[dataset_label] = "failed"
            print(f"Gene ontology analysis failed for {dataset_label}:")
            traceback.print_exc()

        if index < len(result_dirs) - 1:
            time.sleep(delay_between_datasets_seconds)

    succeeded = [label for label, outcome in outcomes.items() if outcome == "success"]
    failed = [label for label, outcome in outcomes.items() if outcome == "failed"]
    print(f"\nGene ontology analysis complete: {len(succeeded)} succeeded, {len(failed)} failed.")
    if failed:
        print("Failed datasets:")
        for label in failed:
            print(f"- {label}")

    return outcomes


if __name__ == "__main__":
    run_gene_ontology_for_all_datasets()
