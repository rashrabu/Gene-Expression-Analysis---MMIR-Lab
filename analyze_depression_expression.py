#!/usr/bin/env python3
"""Prepare and explore the GSE98793 depression expression dataset.

This script performs three tasks on the annotated, renamed, unlogged matrix:
1. Resolves duplicate probes by keeping the probe with the highest mean signal
   for each gene symbol.
2. Runs exploratory analyses comparing depressed (case_*) vs non-depressed
   (control_*) samples.
3. Saves figures as PNG files and tabular summaries as CSV files.
"""

from __future__ import annotations

import math
from pathlib import Path
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


def _split_keywords(text: str) -> list[str]:
    return [token.strip().lower() for token in text.split(",") if token.strip()]


def _maybe_split_keywords(text: str | None) -> list[str]:
    if not text:
        return []
    return _split_keywords(text)

# this function will help to figure out which column means what. It is done instead of hard coding the columns
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
        "Could not identify a gene identifier column. Expected Gene Symbol/Gene_Symbol/"
        "gene_symbol/ID_REF/ID."
    )


def find_sample_columns(
    columns: pd.Index,
    control_prefix: str | tuple[str, ...],
    case_prefix: str | tuple[str, ...],
) -> tuple[list[str], list[str]]:
    control_cols = [col for col in columns if str(col).startswith(control_prefix)]
    case_cols = [col for col in columns if str(col).startswith(case_prefix)]
    if not control_cols or not case_cols:
        raise ValueError(
            "Could not find both control and case sample columns using prefixes"
        )
    return control_cols, case_cols

# backup to extract from metadata. If we are unable to determine the difference based on prefixes.
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
        control_mask |= aligned["_phenotype"].str.contains(re.escape(key), na=False)
        control_mask |= aligned["_sample_name"].str.contains(re.escape(key), na=False)
    for key in case_keywords:
        case_mask |= aligned["_phenotype"].str.contains(re.escape(key), na=False)
        case_mask |= aligned["_sample_name"].str.contains(re.escape(key), na=False)

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


def choose_gene_symbol_column(columns: pd.Index) -> str:
    for candidate in ("Gene Symbol", "Gene_Symbol", "gene_symbol", "SYMBOL", "Symbol"):
        if candidate in columns:
            return candidate
    raise ValueError(
        "Could not find a gene symbol column. Expected one of: "
        "Gene Symbol, Gene_Symbol, gene_symbol, SYMBOL, Symbol"
    )


def clean_gene_symbols(series: pd.Series) -> pd.Series:
    cleaned = series.fillna("").astype(str).str.strip()
    cleaned = cleaned.replace({"nan": "", "---": ""})
    return cleaned


def deduplicate_probes(
    df: pd.DataFrame, gene_col: str, sample_cols: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    working = df.copy()
    working[gene_col] = clean_gene_symbols(working[gene_col])
    working = working[working[gene_col] != ""].copy()
    numeric_samples = working[sample_cols].apply(pd.to_numeric, errors="coerce")
    working = pd.concat(
        [working.drop(columns=sample_cols), numeric_samples],
        axis=1,
    )
    working = working.assign(
        probe_mean_signal=numeric_samples.mean(axis=1, skipna=True)
    )
    working = working.sort_values(
        by=[gene_col, "probe_mean_signal"], ascending=[True, False]
    )
    duplicates = (
        working.groupby(gene_col, as_index=False)
        .agg(
            probes_for_gene=(gene_col, "size"),
            retained_probe_mean_signal=("probe_mean_signal", "first"),
        )
        .sort_values("probes_for_gene", ascending=False)
    )
    deduplicated = working.drop_duplicates(subset=gene_col, keep="first").copy()
    deduplicated = deduplicated.sort_values(gene_col).reset_index(drop=True)
    return deduplicated, duplicates


def build_expression_matrices(
    df: pd.DataFrame, gene_col: str, control_cols: list[str], case_cols: list[str]
) -> pd.DataFrame:
    matrix = df[[gene_col] + control_cols + case_cols].copy()
    matrix = matrix.set_index(gene_col)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.loc[~matrix.index.duplicated(keep="first")]
    return matrix


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    z = matrix.sub(matrix.mean(axis=1), axis=0)
    z = z.div(matrix.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0.0)


def compute_differential_expression(
    matrix: pd.DataFrame, control_cols: list[str], case_cols: list[str]
) -> pd.DataFrame:
    control = matrix[control_cols].to_numpy(dtype=float)
    case = matrix[case_cols].to_numpy(dtype=float)

    control_mean = np.nanmean(control, axis=1)
    case_mean = np.nanmean(case, axis=1)
    mean_difference = case_mean - control_mean
    log2_fold_change = np.log2(case_mean + 1.0) - np.log2(control_mean + 1.0)

    t_stat, p_value = stats.ttest_ind(
        case,
        control,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )

    diff = pd.DataFrame(
        {
            "gene_symbol": matrix.index,
            "control_mean": control_mean,
            "case_mean": case_mean,
            "mean_difference": mean_difference,
            "log2_fold_change": log2_fold_change,
            "t_statistic": t_stat,
            "p_value": p_value,
        }
    )
    diff["abs_mean_difference"] = diff["mean_difference"].abs()
    diff["fdr"] = multipletests(diff["p_value"].fillna(1.0), method="fdr_bh")[1]
    diff["neg_log10_p_value"] = -np.log10(
        np.clip(diff["p_value"].fillna(1.0).to_numpy(), 1e-300, None)
    )
    diff = diff.sort_values(
        by=["p_value", "abs_mean_difference"], ascending=[True, False]
    ).reset_index(drop=True)
    return diff


def compute_group_scores(
    matrix: pd.DataFrame,
    up_genes: list[str],
    down_genes: list[str],
    control_cols: list[str],
    case_label: str,
    control_label: str,
) -> pd.DataFrame:
    selected_genes = [
        gene for gene in dict.fromkeys(up_genes + down_genes) if gene in matrix.index
    ]
    z_selected = zscore_rows(matrix.loc[selected_genes]) if selected_genes else matrix.loc[[]]
    up_index = z_selected.index.intersection(up_genes)
    down_index = z_selected.index.intersection(down_genes)

    sample_groups = []
    control_set = set(control_cols)
    for sample in matrix.columns:
        sample_groups.append(
            {
                "sample": sample,
                "phenotype": control_label if sample in control_set else case_label,
                "up_group_score": z_selected.loc[up_index, sample].mean(),
                "down_group_score": z_selected.loc[down_index, sample].mean(),
            }
        )

    scores = pd.DataFrame(sample_groups)
    scores["contrast_score"] = scores["up_group_score"] - scores["down_group_score"]
    return scores


def _describe_regression_term(
    term: str,
    control_label: str,
    case_label: str,
    covariate_map: dict[str, str],
) -> str:
    if term == "Intercept":
        return f"Intercept ({control_label} baseline)"
    if term.startswith("C(group"):
        return f"{case_label} vs {control_label} (group)"
    for safe_name, original_name in covariate_map.items():
        if term == safe_name:
            return original_name
        if term.startswith(f"C({safe_name})"):
            level = term.split("[T.")[1].rstrip("]") if "[T." in term else ""
            return f"{original_name}: {level}" if level else original_name
    return term


def compute_regression_contrast_score(
    scores: pd.DataFrame,
    metadata: pd.DataFrame | None,
    metadata_sample_column: str | None,
    control_label: str,
    case_label: str,
    covariate_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Fit contrast_score ~ group (+ covariates) via OLS and report the
    case-vs-control regression coefficient adjusted for available covariates
    (e.g. age, gender), alongside the covariate effects themselves."""
    regression_df = scores[["sample", "phenotype", "contrast_score"]].copy()
    covariate_map: dict[str, str] = {}

    if metadata is not None and metadata_sample_column and metadata_sample_column in metadata.columns:
        candidate_covariates = covariate_columns or ["age", "gender"]
        available_covariates = [col for col in candidate_covariates if col in metadata.columns]
        if available_covariates:
            covariate_slice = metadata[[metadata_sample_column] + available_covariates].copy()
            regression_df = regression_df.merge(
                covariate_slice,
                how="left",
                left_on="sample",
                right_on=metadata_sample_column,
            )
            for index, column in enumerate(available_covariates):
                safe_name = f"covariate_{index}"
                regression_df = regression_df.rename(columns={column: safe_name})
                covariate_map[safe_name] = column

    regression_df = regression_df.dropna(
        subset=["contrast_score", "phenotype", *covariate_map.keys()]
    )
    empty_result = pd.DataFrame(
        columns=[
            "term", "coefficient", "std_error", "t_value", "p_value",
            "conf_int_low", "conf_int_high", "n_samples", "r_squared",
            "covariates_included",
        ]
    )
    if {control_label, case_label} - set(regression_df["phenotype"].unique()):
        return empty_result

    regression_df["group"] = pd.Categorical(
        regression_df["phenotype"], categories=[control_label, case_label]
    )

    formula_terms = [f"C(group, Treatment(reference='{control_label}'))"]
    for safe_name in covariate_map:
        if pd.api.types.is_numeric_dtype(regression_df[safe_name]):
            formula_terms.append(safe_name)
        else:
            formula_terms.append(f"C({safe_name})")
    formula = "contrast_score ~ " + " + ".join(formula_terms)

    model = smf.ols(formula, data=regression_df).fit()
    conf_int = model.conf_int()

    rows = []
    for term in model.params.index:
        rows.append(
            {
                "term": _describe_regression_term(term, control_label, case_label, covariate_map),
                "coefficient": model.params[term],
                "std_error": model.bse[term],
                "t_value": model.tvalues[term],
                "p_value": model.pvalues[term],
                "conf_int_low": conf_int.loc[term, 0],
                "conf_int_high": conf_int.loc[term, 1],
            }
        )

    result = pd.DataFrame(rows)
    result["n_samples"] = len(regression_df)
    result["r_squared"] = model.rsquared
    result["covariates_included"] = ", ".join(covariate_map.values()) if covariate_map else "none"
    return result


def save_distribution_plot(
    matrix: pd.DataFrame,
    control_cols: list[str],
    case_cols: list[str],
    output_dir: Path,
    control_label: str,
    case_label: str,
) -> None:
    long_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "sample": control_cols,
                    "phenotype": control_label,
                    "median_expression": matrix[control_cols].median(axis=0).values,
                }
            ),
            pd.DataFrame(
                {
                    "sample": case_cols,
                    "phenotype": case_label,
                    "median_expression": matrix[case_cols].median(axis=0).values,
                }
            ),
        ],
        ignore_index=True,
    )

    plt.figure(figsize=(8, 5))
    sns.boxplot(
        data=long_df,
        x="phenotype",
        y="median_expression",
        hue="phenotype",
        palette="Set2",
        legend=False,
    )
    sns.stripplot(
        data=long_df,
        x="phenotype",
        y="median_expression",
        color="black",
        alpha=0.6,
        size=4,
    )
    plt.title("Sample-Level Median Expression by Phenotype")
    plt.xlabel("")
    plt.ylabel("Median expression across deduplicated genes")
    plt.tight_layout()
    plt.savefig(output_dir / "sample_median_expression.png", dpi=200)
    plt.close()


def save_pca_plot(
    matrix: pd.DataFrame,
    control_cols: list[str],
    case_cols: list[str],
    output_dir: Path,
    metadata: pd.DataFrame | None = None,
    metadata_sample_column: str = "analysis_sample",
    control_label: str = "Control",
    case_label: str = "Case",
) -> None:
    transposed = matrix.T
    transposed = transposed.apply(lambda col: col.fillna(col.mean()), axis=0).fillna(0.0)
    scaled = StandardScaler().fit_transform(transposed)
    pca_model = PCA(n_components=2, random_state=0)
    components = pca_model.fit_transform(scaled)
    explained = pca_model.explained_variance_ratio_

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
        metadata_cols = [metadata_sample_column]
        if "gender" in metadata.columns:
            metadata_cols.append("gender")
        metadata_slice = metadata[metadata_cols].copy()
        if "gender" in metadata_slice.columns:
            metadata_slice["gender"] = metadata_slice["gender"].fillna("Unknown")
        pca_df = pca_df.merge(
            metadata_slice,
            how="left",
            left_on="sample",
            right_on=metadata_sample_column,
        )
        if "gender" in pca_df.columns:
            pca_df["gender"] = pca_df["gender"].fillna("Unknown")

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=pca_df,
        x="PC1",
        y="PC2",
        hue="phenotype",
        style="phenotype",
        s=80,
        palette="Set1",
    )
    plt.title("PCA of Deduplicated Expression Profiles")
    plt.xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
    plt.ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
    plt.tight_layout()
    plt.savefig(output_dir / "pca_samples.png", dpi=200)
    plt.close()

    if metadata is not None and "gender" in pca_df.columns:
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            data=pca_df,
            x="PC1",
            y="PC2",
            hue="gender",
            style="phenotype",
            s=80,
            palette="Set2",
        )
        plt.title("PCA of Deduplicated Expression Profiles Colored by Sex")
        plt.xlabel(f"PC1 ({explained[0] * 100:.1f}% variance)")
        plt.ylabel(f"PC2 ({explained[1] * 100:.1f}% variance)")
        plt.tight_layout()
        plt.savefig(output_dir / "pca_samples_by_sex.png", dpi=200)
        plt.close()


def save_volcano_plot(diff: pd.DataFrame, output_dir: Path, top_n: int) -> None:
    plt.figure(figsize=(9, 6))
    sns.scatterplot(
        data=diff,
        x="log2_fold_change",
        y="neg_log10_p_value",
        alpha=0.65,
        edgecolor=None,
    )
    highlight = diff.nsmallest(top_n, "p_value")
    for row in highlight.itertuples(index=False):
        plt.text(
            row.log2_fold_change,
            row.neg_log10_p_value,
            row.gene_symbol,
            fontsize=7,
            alpha=0.8,
        )
    plt.axvline(0.0, color="grey", linestyle="--", linewidth=1)
    plt.title("Differential Signal: Depressed vs Control")
    plt.xlabel("log2(case mean + 1) - log2(control mean + 1)")
    plt.ylabel("-log10 Welch t-test p-value")
    plt.tight_layout()
    plt.savefig(output_dir / "volcano_differential_expression.png", dpi=220)
    plt.close()


def save_heatmap(
    matrix: pd.DataFrame,
    diff: pd.DataFrame,
    control_cols: list[str],
    case_cols: list[str],
    output_dir: Path,
    top_n: int,
    fdr_threshold: float = 0.05, # default value for FDR threshold
) -> tuple[list[str], list[str]]:
    significant = diff[diff["fdr"] < fdr_threshold].copy()
    significant["abs_log2_fold_change"] = significant["log2_fold_change"].abs()

    up_genes = (
        significant[significant["log2_fold_change"] > 0]
        .sort_values("abs_log2_fold_change", ascending=False)
        .head(top_n)["gene_symbol"]
        .tolist()
    )
    down_genes = (
        significant[significant["log2_fold_change"] < 0]
        .sort_values("abs_log2_fold_change", ascending=False)
        .head(top_n)["gene_symbol"]
        .tolist()
    )

    selected_genes = [gene for gene in up_genes + down_genes if gene in matrix.index]

    if not selected_genes:
        plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(
            0.5, 0.5,
            f"No genes passed FDR < {fdr_threshold}; heatmap skipped.",
            ha="center", va="center", fontsize=12,
        )
        plt.savefig(output_dir / "heatmap_top_differential_genes.png", dpi=220)
        plt.close()
        return up_genes, down_genes

    ordered_cols = control_cols + case_cols
    heatmap_df = zscore_rows(matrix.loc[selected_genes, ordered_cols])

    plt.figure(figsize=(12, max(8, math.ceil(len(selected_genes) * 0.18))))
    sns.heatmap(
        heatmap_df,
        cmap="vlag",
        center=0,
        xticklabels=False,
        yticklabels=True,
        cbar_kws={"label": "Gene-wise z-score"},
    )
    plt.title(
        f"Top {top_n} Up and Top {top_n} Down Differential Genes by Mean Difference"
    )
    plt.xlabel("Samples (controls first, then depressed cases)")
    plt.ylabel("Genes")
    plt.tight_layout()
    plt.savefig(output_dir / "heatmap_top_differential_genes.png", dpi=220)
    plt.close()

    return up_genes, down_genes


def save_gene_group_plot(scores: pd.DataFrame, output_dir: Path) -> None:
    melted = scores.melt(
        id_vars=["sample", "phenotype"],
        value_vars=["up_group_score", "down_group_score", "contrast_score"],
        var_name="group_metric",
        value_name="score",
    )

    if melted["score"].notna().sum() == 0:
        plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(
            0.5, 0.5,
            "No differentially expressed genes available to compute group scores.",
            ha="center", va="center", fontsize=12,
        )
        plt.savefig(output_dir / "gene_group_scores.png", dpi=220)
        plt.close()
        return
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(
        data=melted,
        x="group_metric",
        y="score",
        hue="phenotype",
        palette="Set2",
    )
    plt.title("Gene-Group Scores Derived From Differentially Expressed Genes")
    plt.xlabel("")
    plt.ylabel("Average gene-wise z-score")
    plt.tight_layout()
    plt.savefig(output_dir / "gene_group_scores.png", dpi=220)
    plt.close()


def save_regression_contrast_plot(
    regression_contrast: pd.DataFrame, output_dir: Path
) -> None:
    if regression_contrast.empty:
        plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(
            0.5, 0.5,
            "Regression could not be fit.",
            ha="center", va="center", fontsize=12,
        )
        plt.savefig(output_dir / "regression_contrast_score.png", dpi=220)
        plt.close()
        return

    plot_df = regression_contrast.iloc[::-1].reset_index(drop=True)
    significant = plot_df["p_value"] < 0.05
    sig_color, nonsig_color = "#3B6FA0", "#9CA3AF"

    fig, ax = plt.subplots(figsize=(10, max(4, 0.6 * len(plot_df) + 2)))
    for i, (_, row) in enumerate(plot_df.iterrows()):
        color = sig_color if significant.iloc[i] else nonsig_color
        ax.errorbar(
            row["coefficient"], i,
            xerr=[[row["coefficient"] - row["conf_int_low"]],
                  [row["conf_int_high"] - row["coefficient"]]],
            fmt="o", color=color, ecolor=color, elinewidth=1.5, capsize=4,
            markersize=8, markeredgecolor="white", markeredgewidth=0.8, zorder=3,
        )

    ax.axvline(0, color="#444444", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["term"])
    ax.set_ylim(-0.5, len(plot_df) - 0.5)

    covariates_note = (
        plot_df["covariates_included"].iloc[0]
        if "covariates_included" in plot_df.columns
        else "none"
    )
    xlabel = "Coefficient (95% CI)"
    if covariates_note and covariates_note != "none":
        xlabel += f" — model adjusted for {covariates_note}"
    ax.set_xlabel(xlabel)
    ax.set_title("Regression Contrast Score Coefficients")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=sig_color,
               markersize=9, label="p < 0.05"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=nonsig_color,
               markersize=9, label="p ≥ 0.05"),
    ]
    ax.legend(handles=legend_elements, loc="best", frameon=True)

    plt.tight_layout()
    plt.savefig(output_dir / "regression_contrast_score.png", dpi=220)
    plt.close()


def save_individual_gene_panels(
    matrix: pd.DataFrame,
    diff: pd.DataFrame,
    output_dir: Path,
    panel_count: int = 6,
    control_cols: list[str] | None = None,
    control_label: str = "Control",
    case_label: str = "Case",
) -> None:
    selected = diff.nsmallest(panel_count, "p_value")["gene_symbol"].tolist()
    records = []
    control_set = set(control_cols or [])
    for gene in selected:
        if gene not in matrix.index:
            continue
        for sample, value in matrix.loc[gene].items():
            records.append(
                {
                    "gene_symbol": gene,
                    "sample": sample,
                    "expression": value,
                    "phenotype": control_label if sample in control_set else case_label,
                }
            )

    gene_df = pd.DataFrame(records)
    if gene_df.empty:
        return

    g = sns.catplot(
        data=gene_df,
        x="phenotype",
        y="expression",
        col="gene_symbol",
        kind="box",
        hue="phenotype",
        col_wrap=3,
        sharey=False,
        height=3.5,
        aspect=1.05,
        palette="Set2",
        legend=False,
    )
    g.map_dataframe(
        sns.stripplot,
        x="phenotype",
        y="expression",
        color="black",
        alpha=0.55,
        size=2.8,
    )
    g.set_titles("{col_name}")
    g.set_xlabels("")
    g.set_ylabels("Expression")
    g.figure.subplots_adjust(top=0.88)
    g.figure.suptitle("Top Individual Genes Distinguishing Depressed vs Control")
    g.savefig(output_dir / "individual_gene_panels.png", dpi=220)
    plt.close(g.figure)


def save_summary_tables(
    deduplicated: pd.DataFrame,
    duplicates: pd.DataFrame,
    diff: pd.DataFrame,
    scores: pd.DataFrame,
    up_genes: list[str],
    down_genes: list[str],
    regression_contrast: pd.DataFrame,
    output_dir: Path,
) -> None:
    deduplicated.to_csv(output_dir / "deduplicated_expression_dataset.csv", index=False)
    duplicates.to_csv(output_dir / "duplicate_probe_summary.csv", index=False)
    diff.to_csv(output_dir / "differential_expression_summary.csv", index=False)
    scores.to_csv(output_dir / "gene_group_scores.csv", index=False)
    pd.DataFrame({"top_up_genes": up_genes}).to_csv(
        output_dir / "top_up_genes.csv", index=False
    )
    pd.DataFrame({"top_down_genes": down_genes}).to_csv(
        output_dir / "top_down_genes.csv", index=False
    )
    regression_contrast.to_csv(output_dir / "regression_contrast_score.csv", index=False)


def write_report(
    deduplicated: pd.DataFrame,
    duplicates: pd.DataFrame,
    diff: pd.DataFrame,
    regression_contrast: pd.DataFrame,
    control_cols: list[str],
    case_cols: list[str],
    output_dir: Path,
    top_n: int,
    control_label: str,
    case_label: str,
) -> None:
    duplicate_gene_count = int((duplicates["probes_for_gene"] > 1).sum())
    top_table = diff.head(top_n).loc[
        :, [
            "gene_symbol",
            "control_mean",
            "case_mean",
            "mean_difference",
            "log2_fold_change",
            "p_value",
        ]
    ]

    lines = [
        "Depression expression exploratory analysis",
        "",
        f"Samples: {len(control_cols)} {control_label}, {len(case_cols)} {case_label}",
        f"Genes after deduplication: {len(deduplicated)}",
        f"Genes with duplicate probes resolved: {duplicate_gene_count}",
        "",
        f"Top {top_n} genes ranked by Welch t-test p-value:",
        top_table.to_string(index=False),
        "",
        f"Regression contrast score (contrast_score ~ group"
        + (
            f" + {regression_contrast['covariates_included'].iloc[0]}"
            if not regression_contrast.empty
            and regression_contrast["covariates_included"].iloc[0] != "none"
            else ""
        )
        + "):",
        regression_contrast.drop(columns=["covariates_included"]).to_string(index=False)
        if not regression_contrast.empty
        else "Regression could not be fit.",
        "",
        "Generated figures:",
        "- sample_median_expression.png",
        "- pca_samples.png",
        "- pca_samples_by_sex.png",
        "- volcano_differential_expression.png",
        "- heatmap_top_differential_genes.png",
        "- gene_group_scores.png",
        "- individual_gene_panels.png",
        "- regression_contrast_score.png",
    ]
    (output_dir / "analysis_report.txt").write_text("\n".join(lines))


def run_expression_analysis(
    input_file: str | Path,
    output_dir: str | Path | None = None,
    top_genes: int = 25,
    metadata_file: str | Path | None = None,
    control_prefix: str | tuple[str, ...] | None = None,
    case_prefix: str | tuple[str, ...] | None = None,
    metadata_sample_column: str | None = None,
    metadata_phenotype_column: str | None = None,
    control_label_keywords: str | None = None,
    case_label_keywords: str | None = None,
    control_label: str = "Group_A",
    case_label: str = "Group_B",
    covariate_columns: list[str] | None = None,
) -> None:
    input_path = Path(input_file)
    output_dir = Path(output_dir) if output_dir else input_path.parent / "expression_analysis_results"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    metadata = None
    if metadata_file:
        metadata_path = Path(metadata_file)
        if metadata_path.exists():
            metadata = pd.read_csv(metadata_path)
        else:
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    sns.set_theme(style="whitegrid", context="talk")

    df = pd.read_csv(input_path)
    if metadata is not None and not metadata_sample_column:
        metadata_sample_column = infer_metadata_sample_column(metadata, df.columns)

    gene_col = choose_gene_identifier_column(df)
    try:
        if not control_prefix or not case_prefix:
            raise ValueError("Group prefixes were not provided")
        control_cols, case_cols = find_sample_columns(df.columns, control_prefix, case_prefix)
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
            control_keywords=_maybe_split_keywords(control_label_keywords),
            case_keywords=_maybe_split_keywords(case_label_keywords),
        )

    if "Gene Symbol" not in df.columns and gene_col in df.columns:
        df = df.copy()
        df["Gene Symbol"] = df[gene_col].astype(str)
        gene_col = "Gene Symbol"

    sample_cols = control_cols + case_cols

    deduplicated, duplicates = deduplicate_probes(df, gene_col, sample_cols)
    matrix = build_expression_matrices(
        deduplicated, gene_col, control_cols, case_cols
    )
    diff = compute_differential_expression(matrix, control_cols, case_cols)
    up_genes, down_genes = save_heatmap(
        matrix, diff, control_cols, case_cols, output_dir, top_genes
    )
    scores = compute_group_scores(
        matrix,
        up_genes,
        down_genes,
        control_cols,
        case_label,
        control_label,
    )
    regression_contrast = compute_regression_contrast_score(
        scores,
        metadata,
        metadata_sample_column,
        control_label,
        case_label,
        covariate_columns,
    )
    print("Regression contrast score (contrast_score ~ group + covariates):")
    print(
        regression_contrast.drop(columns=["covariates_included"]).to_string(index=False)
        if not regression_contrast.empty
        else "Regression could not be fit."
    )

    save_distribution_plot(
        matrix,
        control_cols,
        case_cols,
        output_dir,
        control_label,
        case_label,
    )
    save_pca_plot(
        matrix,
        control_cols,
        case_cols,
        output_dir,
        metadata,
        metadata_sample_column=metadata_sample_column,
        control_label=control_label,
        case_label=case_label,
    )
    save_volcano_plot(diff, output_dir, min(15, top_genes))
    save_gene_group_plot(scores, output_dir)
    save_regression_contrast_plot(regression_contrast, output_dir)
    save_individual_gene_panels(
        matrix,
        diff,
        output_dir,
        control_cols=control_cols,
        control_label=control_label,
        case_label=case_label,
    )
    save_summary_tables(
        deduplicated,
        duplicates,
        diff,
        scores,
        up_genes,
        down_genes,
        regression_contrast,
        output_dir,
    )
    write_report(
        deduplicated,
        duplicates,
        diff,
        regression_contrast,
        control_cols,
        case_cols,
        output_dir,
        min(15, top_genes),
        control_label,
        case_label,
    )

    print(f"Saved analysis outputs to {output_dir}")


def main() -> None:
    run_expression_analysis(
        input_file="Analysis/GSE98793_annotated_renamed_unlogged.csv",
        output_dir="Analysis/results",
        top_genes=25,
        metadata_file="Analysis/GSE98793_sample_metadata.csv",
        control_prefix="control_",
        case_prefix="case_",
        metadata_sample_column="analysis_sample",
        metadata_phenotype_column="depression_status",
        control_label_keywords="control,healthy",
        case_label_keywords="case,depressed,mdd,disease",
        control_label="Control",
        case_label="Case",
    )


if __name__ == "__main__":
    main()