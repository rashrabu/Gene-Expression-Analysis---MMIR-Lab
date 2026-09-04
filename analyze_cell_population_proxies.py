"""Estimate blood cell-population proxy scores from GEO-style expression datasets.

This script is designed to be reusable across different GEO datasets by:
1) flexibly identifying gene and sample columns,
2) inferring phenotype groups from either sample-name prefixes or metadata, and
3) loading blood-cell marker genes from an external CSV in multiple formats.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests




def _split_keywords(text: str) -> list[str]:
    return [token.strip().lower() for token in text.split(",") if token.strip()]


def _maybe_split_keywords(text: str | None) -> list[str]:
    if not text:
        return []
    return _split_keywords(text)


def _tokenize_markers(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    tokens = re.split(r"[;,|]", text)
    return [token.strip() for token in tokens if token.strip()]


def _tokenize_gene_identifier(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    # GEO annotations often store aliases as "GENE1 /// GENE2".
    tokens = re.split(r"\s*///\s*|\s*[;,|]\s*", text)
    return [token.strip() for token in tokens if token.strip()]


def infer_geo_platform_from_dataset_path(input_path: Path) -> str | None:
    current = input_path.resolve().parent
    for parent in [current, *current.parents]:
        for raw_file in parent.glob("original_data/*_series_matrix.txt"):
            try:
                with raw_file.open("r", encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if line.startswith("!Series_platform_id"):
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) >= 2:
                                platform = parts[1].strip().strip('"')
                                return platform or None
                        if line.startswith("!series_matrix_table_begin"):
                            break
            except OSError:
                continue
    return None


def choose_gene_identifier_column(df: pd.DataFrame) -> str:
    candidates = [
        "Gene Symbol",
        "Gene_Symbol",
        "gene_symbol",
        "SYMBOL",
        "Symbol",
        "ID_REF",
        "ID",
        "Probe Set ID",
    ]
    for candidate in candidates:
        if candidate in df.columns:
            if candidate in ("ID_REF", "ID", "Probe Set ID"):
                continue
            series = df[candidate].fillna("").astype(str).str.strip()
            if (series != "").sum() > 0:
                return candidate

    for fallback in ("ID_REF", "ID", "Probe Set ID"):
        if fallback in df.columns:
            return fallback

    raise ValueError(
        "Could not identify a gene identifier column. "
        "Expected Gene Symbol/Gene_Symbol/gene_symbol/ID_REF/ID."
    )


def find_sample_columns_by_prefix(
    columns: pd.Index,
    control_prefix: str | None,
    case_prefix: str | None,
) -> tuple[list[str], list[str]]:
    if not control_prefix or not case_prefix:
        raise ValueError("Group prefixes were not provided")
    control_cols = [col for col in columns if str(col).startswith(control_prefix)]
    case_cols = [col for col in columns if str(col).startswith(case_prefix)]
    if not control_cols or not case_cols:
        raise ValueError("Could not identify both control and case sample columns by prefix")
    return control_cols, case_cols


def find_sample_columns_from_metadata(
    columns: pd.Index,
    metadata: pd.DataFrame,
    sample_column: str,
    phenotype_column: str | None,
    control_keywords: list[str],
    case_keywords: list[str],
) -> tuple[list[str], list[str]]:
    if sample_column not in metadata.columns:
        raise ValueError(f"Metadata sample column not found: {sample_column}")

    column_set = {str(col) for col in columns}
    aligned = metadata[metadata[sample_column].astype(str).isin(column_set)].copy()
    if aligned.empty:
        raise ValueError("No metadata samples match expression column names")

    aligned["_sample_name"] = aligned[sample_column].astype(str).str.lower()
    if phenotype_column and phenotype_column in aligned.columns:
        aligned["_phenotype"] = aligned[phenotype_column].fillna("").astype(str).str.lower()
    else:
        aligned["_phenotype"] = ""

    control_mask = np.zeros(len(aligned), dtype=bool)
    case_mask = np.zeros(len(aligned), dtype=bool)
    for key in control_keywords:
        control_mask |= aligned["_phenotype"].str.contains(key, na=False)
        control_mask |= aligned["_sample_name"].str.contains(key, na=False)
    for key in case_keywords:
        case_mask |= aligned["_phenotype"].str.contains(key, na=False)
        case_mask |= aligned["_sample_name"].str.contains(key, na=False)

    control_cols = aligned.loc[control_mask, sample_column].astype(str).tolist()
    case_cols = aligned.loc[case_mask, sample_column].astype(str).tolist()

    if not control_cols or not case_cols:
        if phenotype_column and phenotype_column in aligned.columns:
            values = aligned[phenotype_column].fillna("").astype(str).str.strip()
            unique_values = [value for value in sorted(values.unique()) if value]
            if len(unique_values) == 2:
                control_cols = aligned.loc[
                    values == unique_values[0], sample_column
                ].astype(str).tolist()
                case_cols = aligned.loc[
                    values == unique_values[1], sample_column
                ].astype(str).tolist()

    if not control_cols or not case_cols:
        raise ValueError(
            "Could not infer both control and case groups from metadata and sample names"
        )
    return control_cols, case_cols


def infer_metadata_sample_column(
    metadata: pd.DataFrame,
    expression_columns: pd.Index,
) -> str:
    best_col = None
    best_overlap = 0
    expr_set = {str(col) for col in expression_columns}

    for col in metadata.columns:
        values = metadata[col].dropna().astype(str)
        if values.empty:
            continue
        overlap = int(values.isin(expr_set).sum())
        if overlap > best_overlap:
            best_overlap = overlap
            best_col = col

    if best_col is None or best_overlap == 0:
        raise ValueError(
            "Could not infer metadata sample column. Provide metadata_sample_column explicitly."
        )
    return best_col


def normalize_marker_dictionary(marker_dict: dict[str, list[str]]) -> dict[str, list[str]]:
    normalized: dict[str, list[str]] = {}
    for cell_type, genes in marker_dict.items():
        unique = []
        seen = set()
        for gene in genes:
            g = str(gene).strip()
            if not g:
                continue
            key = g.upper()
            if key in seen:
                continue
            seen.add(key)
            unique.append(g)
        if unique:
            normalized[str(cell_type).strip()] = unique
    return normalized


def determine_blood_cell_markers_from_csv(
    marker_csv: Path,
    celltype_column: str = "auto",
    gene_column: str = "auto",
) -> dict[str, list[str]]:
    """Determine blood-cell markers from a CSV in long or wide format.

    Supported formats:
    1) Long format: one row per marker with columns like cell_type, gene.
    2) Wide row format: first column cell type, remaining columns marker genes.
    3) Wide column format: each column is a cell type, rows are marker genes.
    """
    marker_df = pd.read_csv(marker_csv, low_memory=False)
    if marker_df.empty:
        raise ValueError(f"Marker CSV is empty: {marker_csv}")

    columns = list(marker_df.columns)
    lower_map = {col.lower().strip(): col for col in columns}

    if celltype_column != "auto" and celltype_column not in marker_df.columns:
        raise ValueError(f"Marker cell-type column not found: {celltype_column}")
    if gene_column != "auto" and gene_column not in marker_df.columns:
        raise ValueError(f"Marker gene column not found: {gene_column}")

    if celltype_column == "auto":
        for candidate in ("cell_type", "celltype", "cell", "population", "immune_cell"):
            if candidate in lower_map:
                celltype_column = lower_map[candidate]
                break
    if gene_column == "auto":
        for candidate in ("gene", "gene_symbol", "symbol", "marker", "marker_gene"):
            if candidate in lower_map:
                gene_column = lower_map[candidate]
                break

    marker_dict: dict[str, list[str]] = {}

    # Long format.
    if celltype_column != "auto" and gene_column != "auto":
        grouped = marker_df.groupby(celltype_column, dropna=True)
        for cell_type, group in grouped:
            genes = []
            for value in group[gene_column]:
                genes.extend(_tokenize_markers(value))
            marker_dict[str(cell_type)] = genes
        return normalize_marker_dictionary(marker_dict)

    # Wide row format: first column is cell type, the rest are marker genes.
    first_col = marker_df.columns[0]
    if marker_df.shape[1] >= 2:
        non_null_first = marker_df[first_col].notna().sum()
        if non_null_first >= max(2, int(0.5 * len(marker_df))):
            temp_dict: dict[str, list[str]] = {}
            for row in marker_df.itertuples(index=False):
                cell_type = str(row[0]).strip()
                if not cell_type or cell_type.lower() in ("nan", "none"):
                    continue
                genes = []
                for value in row[1:]:
                    genes.extend(_tokenize_markers(value))
                if genes:
                    temp_dict[cell_type] = genes
            normalized = normalize_marker_dictionary(temp_dict)
            if normalized:
                return normalized

    # Wide column format: each column is a cell type, rows are marker genes.
    temp_dict = {}
    for col in marker_df.columns:
        genes = []
        for value in marker_df[col].tolist():
            genes.extend(_tokenize_markers(value))
        if genes:
            temp_dict[str(col)] = genes

    normalized = normalize_marker_dictionary(temp_dict)
    if not normalized:
        raise ValueError(
            "Could not determine marker format from CSV. Provide explicit "
            "--marker-celltype-column and --marker-gene-column for long format."
        )
    return normalized


def build_expression_matrices(
    df: pd.DataFrame,
    gene_col: str,
    control_cols: list[str],
    case_cols: list[str],
) -> pd.DataFrame:
    sample_cols = control_cols + case_cols
    matrix = df[[gene_col] + sample_cols].copy()
    matrix = matrix.set_index(gene_col)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]
    return matrix


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    z = matrix.sub(matrix.mean(axis=1), axis=0)
    z = z.div(matrix.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0.0)


def compute_pca(
    matrix: pd.DataFrame,
    control_cols: list[str],
    metadata: pd.DataFrame | None,
    metadata_sample_column: str,
    control_label: str,
    case_label: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    transposed = matrix.T
    transposed = transposed.apply(lambda col: col.fillna(col.mean()), axis=0).fillna(0.0)
    scaled = StandardScaler().fit_transform(transposed)
    pca_model = PCA(n_components=2, random_state=0)
    components = pca_model.fit_transform(scaled)
    pca_df = pd.DataFrame(
        {
            "sample": transposed.index,
            "PC1": components[:, 0],
            "PC2": components[:, 1],
            "phenotype": [
                control_label if sample in control_cols else case_label
                for sample in transposed.index
            ],
        }
    )

    if metadata is not None and metadata_sample_column in metadata.columns:
        pca_df = pca_df.merge(
            metadata,
            how="left",
            left_on="sample",
            right_on=metadata_sample_column,
        )

    return pca_df, pca_model.explained_variance_ratio_


def compute_cell_population_scores(
    matrix: pd.DataFrame,
    blood_cell_markers: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_records = []
    marker_records = []

    # Case-insensitive matching to support mixed gene-symbol casing.
    upper_to_gene: dict[str, str] = {}
    for gene in matrix.index.astype(str):
        upper_to_gene.setdefault(gene.upper(), gene)
        for token in _tokenize_gene_identifier(gene):
            upper_to_gene.setdefault(token.upper(), gene)

    for cell_type, markers in blood_cell_markers.items():
        used_markers = []
        for marker in markers:
            key = str(marker).upper()
            if key in upper_to_gene:
                used_markers.append(upper_to_gene[key])

        used_markers = list(dict.fromkeys(used_markers))
        if not used_markers:
            continue

        z_markers = zscore_rows(matrix.loc[used_markers])
        cell_scores = z_markers.mean(axis=0)
        for sample, score in cell_scores.items():
            score_records.append(
                {
                    "sample": sample,
                    "cell_type": cell_type,
                    "proxy_score": score,
                }
            )

        marker_records.append(
            {
                "cell_type": cell_type,
                "markers_requested": len(markers),
                "markers_used": len(used_markers),
                "coverage_proportion": len(used_markers) / len(markers) if markers else 0.0,
                "marker_genes_used": "; ".join(used_markers),
            }
        )
        

    return pd.DataFrame(score_records), pd.DataFrame(marker_records)


def summarize_associations(score_df: pd.DataFrame, pca_df: pd.DataFrame) -> pd.DataFrame:
    if score_df.empty:
        return pd.DataFrame(
            columns=[
                "cell_type",
                "pc1_spearman_rho",
                "pc1_p_value",
                "pc2_spearman_rho",
                "pc2_p_value",
                "control_mean_proxy",
                "case_mean_proxy",
                "mean_difference",
                "case_control_t_statistic",
                "case_control_p_value",
            ]
        )

    merged = score_df.merge(pca_df, how="left", on="sample")
    rows = []
    for cell_type, group in merged.groupby("cell_type"):
        pc1_rho, pc1_p = stats.spearmanr(group["proxy_score"], group["PC1"])
        pc2_rho, pc2_p = stats.spearmanr(group["proxy_score"], group["PC2"])
        case_scores = group.loc[group["phenotype"] == "Case", "proxy_score"]
        control_scores = group.loc[group["phenotype"] == "Control", "proxy_score"]
        t_stat, p_value = stats.ttest_ind(
            case_scores,
            control_scores,
            equal_var=False,
            nan_policy="omit",
        )
        rows.append(
            {
                "cell_type": cell_type,
                "pc1_spearman_rho": pc1_rho,
                "pc1_p_value": pc1_p,
                "pc2_spearman_rho": pc2_rho,
                "pc2_p_value": pc2_p,
                "control_mean_proxy": control_scores.mean(),
                "case_mean_proxy": case_scores.mean(),
                "mean_difference": case_scores.mean() - control_scores.mean(),
                "case_control_t_statistic": t_stat,
                "case_control_p_value": p_value,
            }
        )

    result = pd.DataFrame(rows)
    result["case_control_fdr"] = multipletests(
        result["case_control_p_value"].fillna(1.0), method="fdr_bh"
    )[1]
    return result.sort_values(
        by=["pc1_spearman_rho"],
        key=lambda col: col.abs(),
        ascending=False,
    )


def save_proxy_boxplot(score_df: pd.DataFrame, pca_df: pd.DataFrame, output_dir: Path) -> None:
    if score_df.empty:
        plt.figure(figsize=(12, 6))
        plt.axis("off")
        plt.text(
            0.5,
            0.5,
            "No cell-population proxy scores available",
            ha="center",
            va="center",
            fontsize=13,
        )
        plt.title("Blood Cell-Population Proxy Scores by Group")
        plt.tight_layout()
        plt.savefig(output_dir / "cell_population_proxy_scores_by_group.png", dpi=220)
        plt.close()
        return

    plot_df = score_df.merge(pca_df[["sample", "phenotype"]], how="left", on="sample")
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        data=plot_df,
        x="cell_type",
        y="proxy_score",
        hue="phenotype",
        palette="Set2",
    )
    plt.xticks(rotation=30, ha="right")
    plt.xlabel("")
    plt.ylabel("Average marker-gene z-score")
    plt.title("Blood Cell-Population Proxy Scores by Group")
    plt.tight_layout()
    plt.savefig(output_dir / "cell_population_proxy_scores_by_group.png", dpi=220)
    plt.close()


def select_top_cell_types(
    associations: pd.DataFrame,
    top_n: int,
) -> tuple[list[str], list[str]]:
    """Rank cell types by case-vs-control mean_difference in proxy score and
    split into the top_n most increased ("up") and top_n most decreased
    ("down") in cases, mirroring the up/down gene selection used for the
    differential-expression heatmap."""
    if associations.empty:
        return [], []

    ranked = associations.copy()
    ranked["abs_mean_difference"] = ranked["mean_difference"].abs()

    up_cell_types = (
        ranked[ranked["mean_difference"] > 0]
        .sort_values("abs_mean_difference", ascending=False)
        .head(top_n)["cell_type"]
        .tolist()
    )
    down_cell_types = (
        ranked[ranked["mean_difference"] < 0]
        .sort_values("abs_mean_difference", ascending=False)
        .head(top_n)["cell_type"]
        .tolist()
    )
    return up_cell_types, down_cell_types


def save_proxy_heatmap(
    score_df: pd.DataFrame,
    pca_df: pd.DataFrame,
    associations: pd.DataFrame,
    output_dir: Path,
    top_n: int = 10,
) -> tuple[list[str], list[str]]:
    if score_df.empty:
        return [], []

    up_cell_types, down_cell_types = select_top_cell_types(associations, top_n)
    selected_cell_types = [
        cell_type
        for cell_type in dict.fromkeys(up_cell_types + down_cell_types)
        if cell_type in score_df["cell_type"].unique()
    ]
    if not selected_cell_types:
        selected_cell_types = score_df["cell_type"].unique().tolist()

    merge_cols = ["sample", "phenotype"]
    if "gender" in pca_df.columns:
        merge_cols.append("gender")

    plot_df = score_df[score_df["cell_type"].isin(selected_cell_types)].merge(
        pca_df[merge_cols], how="left", on="sample"
    )

    order_cols = ["phenotype", "sample"]
    if "gender" in plot_df.columns:
        order_cols = ["phenotype", "gender", "sample"]

    order_df = plot_df[["sample"] + [col for col in order_cols if col != "sample"]].drop_duplicates()
    order_df = order_df.sort_values(order_cols)

    wide_proxy_scores = plot_df.pivot(index="cell_type", columns="sample", values="proxy_score")
    heatmap_df = zscore_rows(wide_proxy_scores).reindex(
        index=selected_cell_types, columns=order_df["sample"]
    )

    plt.figure(figsize=(14, max(6, math.ceil(len(selected_cell_types) * 0.4))))
    sns.heatmap(
        heatmap_df,
        cmap="vlag",
        center=0,
        xticklabels=False,
        yticklabels=True,
        cbar_kws={"label": "Cell-type-wise z-score of proxy score"},
    )
    plt.title("Top Up and Down Blood Cell-Population Proxy Scores Across Samples")
    plt.xlabel("Samples ordered by group, optional sex, and sample name")
    plt.ylabel("Cell population")
    plt.tight_layout()
    plt.savefig(output_dir / "cell_population_proxy_heatmap.png", dpi=220)
    plt.close()

    return up_cell_types, down_cell_types


def save_pca_colored_by_top_proxies(
    score_df: pd.DataFrame,
    pca_df: pd.DataFrame,
    explained: np.ndarray,
    associations: pd.DataFrame,
    output_dir: Path,
) -> None:
    if score_df.empty or associations.empty:
        return

    top_cell_types = associations.head(4)["cell_type"].tolist()
    merged = score_df.merge(pca_df, how="left", on="sample")
    subset = merged[merged["cell_type"].isin(top_cell_types)].copy()
    if subset.empty:
        return

    subset["panel_title"] = subset["cell_type"]
    grid = sns.FacetGrid(
        subset,
        col="panel_title",
        col_wrap=2,
        height=4,
        sharex=True,
        sharey=True,
    )
    grid.map_dataframe(
        sns.scatterplot,
        x="PC1",
        y="PC2",
        hue="proxy_score",
        palette="viridis",
        s=70,
    )
    grid.set_titles("{col_name}")
    grid.set_axis_labels(
        f"PC1 ({explained[0] * 100:.1f}% variance)",
        f"PC2 ({explained[1] * 100:.1f}% variance)",
    )
    grid.add_legend()
    grid.figure.subplots_adjust(top=0.88)
    grid.figure.suptitle("PCA Colored by Top Cell-Population Proxy Scores")
    grid.savefig(output_dir / "pca_colored_by_top_cell_proxies.png", dpi=220)
    plt.close(grid.figure)


def report_pc1_correlation_by_cell_type(associations: pd.DataFrame) -> pd.DataFrame:
    """Return PC1 vs. proxy-score Spearman correlations, one row per cell type,
    ranked by correlation strength (strongest association first)."""
    if associations.empty:
        return pd.DataFrame(columns=["cell_type", "pc1_spearman_rho", "pc1_p_value"])

    return (
        associations[["cell_type", "pc1_spearman_rho", "pc1_p_value"]]
        .sort_values("pc1_spearman_rho", key=lambda col: col.abs(), ascending=False)
        .reset_index(drop=True)
    )


def write_report(
    marker_df: pd.DataFrame,
    associations: pd.DataFrame,
    output_dir: Path,
    marker_source: str,
) -> None:
    pc1_correlations = report_pc1_correlation_by_cell_type(associations)
    lines = [
        "Cell-population proxy analysis",
        "",
        f"Marker source: {marker_source}",
        "",
        "Method:",
        "Marker-based proxy scores were computed for each configured cell type",
        "from the bulk expression matrix and tested against PCA axes.",
        "",
        "Marker coverage:",
        marker_df.to_string(index=False) if not marker_df.empty else "No markers matched.",
        "",
        "PC1 correlation by cell type (Spearman, ranked by |rho|):",
        pc1_correlations.to_string(index=False) if not pc1_correlations.empty else "No associations available.",
        "",
        "Top PCA associations:",
        associations.head(8).to_string(index=False) if not associations.empty else "No associations available.",
    ]
    (output_dir / "cell_population_proxy_report.txt").write_text("\n".join(lines), encoding="utf-8")


def run_cell_proxy_analysis(
    input_file: str | Path,
    metadata_file: str | Path | None = None,
    output_dir: str | Path | None = None,
    markers_csv: str | Path | None = None,
    marker_celltype_column: str | None = None,
    marker_gene_column: str | None = None,
    control_prefix: str | None = None,
    case_prefix: str | None = None,
    metadata_sample_column: str | None = None,
    metadata_phenotype_column: str | None = None,
    control_label_keywords: str | None = None,
    case_label_keywords: str | None = None,
    control_label: str = "Group_A",
    case_label: str = "Group_B",
) -> None:
    input_path = Path(input_file)
    output_dir = Path(output_dir) if output_dir else input_path.parent / "cell_population_proxies"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    metadata = None
    if metadata_file:
        metadata_path = Path(metadata_file)
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
        metadata = pd.read_csv(metadata_path)

    if not markers_csv:
        raise ValueError(
            "markers_csv is required for generic GEO processing. "
            "Provide a marker definition CSV (long or wide format)."
        )
    marker_path = Path(markers_csv)
    if not marker_path.exists():
        raise FileNotFoundError(f"Marker CSV not found: {marker_path}")

    blood_cell_markers = determine_blood_cell_markers_from_csv(
        marker_path,
        celltype_column=marker_celltype_column or "auto",
        gene_column=marker_gene_column or "auto",
    )
    marker_source = str(marker_path)

    control_keywords = _maybe_split_keywords(control_label_keywords)
    case_keywords = _maybe_split_keywords(case_label_keywords)

    sns.set_theme(style="whitegrid", context="talk")

    df = pd.read_csv(input_path, low_memory=False)
    gene_col = choose_gene_identifier_column(df)

    if metadata is not None and not metadata_sample_column:
        metadata_sample_column = infer_metadata_sample_column(metadata, df.columns)

    try:
        control_cols, case_cols = find_sample_columns_by_prefix(
            df.columns,
            control_prefix,
            case_prefix,
        )
    except ValueError:
        if metadata is None:
            raise ValueError(
                "Could not infer two phenotype groups from prefixes and no metadata was provided"
            )
        if not metadata_sample_column:
            raise ValueError(
                "metadata_sample_column is required when metadata-based grouping is used"
            )
        control_cols, case_cols = find_sample_columns_from_metadata(
            df.columns,
            metadata,
            sample_column=metadata_sample_column,
            phenotype_column=metadata_phenotype_column,
            control_keywords=control_keywords,
            case_keywords=case_keywords,
        )

    matrix = build_expression_matrices(df, gene_col, control_cols, case_cols)
    pca_df, explained = compute_pca(
        matrix,
        control_cols,
        metadata,
        metadata_sample_column=metadata_sample_column,
        control_label=control_label,
        case_label=case_label,
    )

    score_df, marker_df = compute_cell_population_scores(matrix, blood_cell_markers)
    total_markers_used = int(marker_df["markers_used"].sum()) if not marker_df.empty else 0
    if total_markers_used == 0:
        preview = (
            pd.Series(df[gene_col])
            .dropna()
            .astype(str)
            .str.strip()
        )
        preview = preview[preview != ""]
        preview_values = preview.drop_duplicates().head(5).tolist()
        platform = infer_geo_platform_from_dataset_path(input_path)
        platform_hint = (
            f"Detected platform: {platform}. Add a matching annotation CSV in a folder named "
            f"{platform}/ (for example, {platform}/{platform}_annotation.csv), then rerun preprocessing. "
            if platform
            else "Provide the correct platform annotation CSV to preprocessing (annotation_path), then rerun. "
        )
        raise ValueError(
            "No configured marker genes matched the expression gene identifiers. "
            f"Detected gene column '{gene_col}' with examples: {preview_values}. "
            "This usually means the matrix is probe-level (or unannotated) and not mapped "
            f"to gene symbols. {platform_hint}"
        )
    associations = summarize_associations(score_df, pca_df)

    if score_df.empty:
        wide_scores = pd.DataFrame(columns=["sample"])
    else:
        wide_scores = score_df.pivot(index="sample", columns="cell_type", values="proxy_score")
        wide_scores = wide_scores.reset_index()

    if metadata is not None and metadata_sample_column in metadata.columns:
        wide_scores = wide_scores.merge(
            metadata,
            how="left",
            left_on="sample",
            right_on=metadata_sample_column,
        )

    wide_scores.to_csv(output_dir / "cell_population_proxy_scores.csv", index=False)
    marker_df.to_csv(output_dir / "cell_population_marker_coverage.csv", index=False)
    associations.to_csv(output_dir / "cell_population_pca_associations.csv", index=False)

    pc1_correlations = report_pc1_correlation_by_cell_type(associations)
    pc1_correlations.to_csv(output_dir / "pc1_correlation_by_cell_type.csv", index=False)

    save_proxy_boxplot(score_df, pca_df, output_dir)
    up_cell_types, down_cell_types = save_proxy_heatmap(score_df, pca_df, associations, output_dir)
    pd.DataFrame({"top_up_cell_types": up_cell_types}).to_csv(
        output_dir / "top_up_cell_types.csv", index=False
    )
    pd.DataFrame({"top_down_cell_types": down_cell_types}).to_csv(
        output_dir / "top_down_cell_types.csv", index=False
    )
    save_pca_colored_by_top_proxies(score_df, pca_df, explained, associations, output_dir)
    write_report(marker_df, associations, output_dir, marker_source=marker_source)

    print(f"Saved cell-population proxy outputs to {output_dir}")
    print("PC1 correlation by cell type (Spearman, ranked by |rho|):")
    print(pc1_correlations.to_string(index=False) if not pc1_correlations.empty else "No associations available.")