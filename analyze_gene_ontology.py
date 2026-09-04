#!/usr/bin/env python3
"""Run GO and KEGG over-representation analysis on a dataset's up/down gene lists.

This script is designed to be reusable across different GEO datasets: point it
at the top_up_genes.csv / top_down_genes.csv produced by
analyze_depression_expression.py (or any CSV with a gene symbol column) and it
queries the Enrichr API for GO (Biological Process, Molecular Function,
Cellular Component) and KEGG pathway enrichment, separately for the
up-regulated and down-regulated gene sets.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

DEFAULT_GO_LIBRARIES = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
]
DEFAULT_KEGG_LIBRARY = "KEGG_2021_Human"

EMPTY_RESULTS_COLUMNS = [
    "Gene_set",
    "Term",
    "P-value",
    "Adjusted P-value",
    "Odds Ratio",
    "Combined Score",
    "Genes",
]


def load_gene_list(
    path: str | Path,
    gene_column: str | None = None,
) -> list[str]:
    """Read a CSV and return a cleaned, de-duplicated list of gene symbols."""
    df = pd.read_csv(path, dtype=str)
    if df.empty:
        return []

    if gene_column is None:
        candidates = [
            "gene_symbol",
            "Gene Symbol",
            "Gene_Symbol",
            "SYMBOL",
            "Symbol",
            "top_up_genes",
            "top_down_genes",
        ]
        gene_column = next((col for col in candidates if col in df.columns), df.columns[0])

    if gene_column not in df.columns:
        raise ValueError(f"Gene column not found: {gene_column}")

    genes = df[gene_column].dropna().astype(str).str.strip()
    genes = genes[genes != ""]
    return genes.drop_duplicates().tolist()


def run_ontology_enrichment(
    genes: list[str],
    gene_sets: list[str] | str,
    organism: str,
    background: list[str] | None,
) -> pd.DataFrame:
    """Query Enrichr for one or more gene-set libraries and return the raw
    results table (empty with the expected columns if there are no genes)."""
    if not genes:
        return pd.DataFrame(columns=EMPTY_RESULTS_COLUMNS)

    enrichment = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism=organism,
        background=background,
        outdir=None,
        no_plot=True,
    )
    return enrichment.results


def save_enrichment_barplot(
    results: pd.DataFrame,
    title: str,
    output_path: Path,
    top_n: int,
) -> None:
    if results.empty:
        plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(0.5, 0.5, "No enrichment results available.", ha="center", va="center", fontsize=12)
        plt.savefig(output_path, dpi=220)
        plt.close()
        return

    top_terms = results.sort_values("Adjusted P-value").head(top_n).copy()
    top_terms["neg_log10_fdr"] = -np.log10(
        np.clip(top_terms["Adjusted P-value"].to_numpy(dtype=float), 1e-300, None)
    )
    top_terms = top_terms.sort_values("neg_log10_fdr")

    plt.figure(figsize=(10, max(4, math.ceil(len(top_terms) * 0.4))))
    sns.barplot(
        data=top_terms,
        x="neg_log10_fdr",
        y="Term",
        color="#4C72B0",
    )
    plt.axvline(-np.log10(0.05), color="grey", linestyle="--", linewidth=1)
    plt.title(title)
    plt.xlabel("-log10(Adjusted P-value)")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def write_ontology_report(
    results_by_key: dict[str, pd.DataFrame],
    output_dir: Path,
    fdr_threshold: float,
    top_n: int,
) -> None:
    lines = ["Gene ontology (GO + KEGG) enrichment analysis", ""]
    for key, results in results_by_key.items():
        significant = results[results["Adjusted P-value"] < fdr_threshold]
        lines.append(f"{key}: {len(significant)} terms with FDR < {fdr_threshold}")
        if not significant.empty:
            display_columns = [
                col
                for col in ("Term", "Overlap", "Adjusted P-value", "Combined Score")
                if col in significant.columns
            ]
            top_table = significant.sort_values("Adjusted P-value").head(top_n)[display_columns]
            lines.append(top_table.to_string(index=False))
        lines.append("")
    (output_dir / "gene_ontology_report.txt").write_text("\n".join(lines), encoding="utf-8")


def run_gene_ontology_analysis(
    up_genes_file: str | Path,
    down_genes_file: str | Path,
    output_dir: str | Path,
    background_file: str | Path | None = None,
    gene_column: str | None = None,
    background_gene_column: str | None = None,
    organism: str = "human",
    go_libraries: list[str] | None = None,
    kegg_library: str = DEFAULT_KEGG_LIBRARY,
    fdr_threshold: float = 0.05,
    top_n_terms: int = 15,
) -> dict[str, pd.DataFrame]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    go_libraries = go_libraries or DEFAULT_GO_LIBRARIES

    up_genes = load_gene_list(up_genes_file, gene_column)
    down_genes = load_gene_list(down_genes_file, gene_column)

    background = None
    if background_file is not None:
        background = load_gene_list(background_file, background_gene_column) or None

    gene_sets = go_libraries + [kegg_library]
    results_by_key: dict[str, pd.DataFrame] = {}

    for direction, genes in (("up", up_genes), ("down", down_genes)):
        results = run_ontology_enrichment(genes, gene_sets, organism, background)
        for gene_set_name, group in (
            list(results.groupby("Gene_set")) if not results.empty else [(name, pd.DataFrame(columns=EMPTY_RESULTS_COLUMNS)) for name in gene_sets]
        ):
            key = f"{gene_set_name}_{direction}"
            results_by_key[key] = group.reset_index(drop=True)
            group.to_csv(output_dir / f"{key}_enrichment.csv", index=False)
            save_enrichment_barplot(
                group,
                title=f"{gene_set_name} enrichment ({direction}-regulated genes)",
                output_path=output_dir / f"{key}_enrichment.png",
                top_n=top_n_terms,
            )

    write_ontology_report(results_by_key, output_dir, fdr_threshold, top_n_terms)
    print(f"Saved gene ontology analysis outputs to {output_dir}")
    return results_by_key
