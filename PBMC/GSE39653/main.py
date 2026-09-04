from pathlib import Path
import sys

def _find_project_root(start: Path) -> Path:
    current = start.resolve()
    for parent in [current, *current.parents]:
        if (parent / "preprocessing.py").exists():
            return parent
    raise RuntimeError("Could not locate project root containing preprocessing.py")


PROJECT_ROOT = _find_project_root(Path(__file__).parent)
TEST_DIR = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from preprocessing import run_preprocessor
from analyze_depression_expression import run_expression_analysis
from analyze_cell_population_proxies import run_cell_proxy_analysis
from run_qc_normalize_r import run_qc_normalization

# Study includes MDD, bipolar disorder (BD1/BD2), and healthy control (HC)
# subjects. Only MDD vs HC is scored as Case/Control here; BD samples are
# preprocessed and kept in the metadata but left unlabeled (phenotype blank)
# so they are excluded from the two-group differential-expression/cell-proxy
# steps rather than being folded into either group.
mapping_text = """
GSM973762	MDD-376
GSM973763	HC-104
GSM973764	HC-261
GSM973765	MDD-465
GSM973766	HC-513
GSM973767	BD2-231
GSM973768	HC-365
GSM973769	HC-3
GSM973770	MDD-456
GSM973771	HC-284
GSM973772	BD2-282
GSM973773	MDD-425
GSM973774	MDD-397
GSM973775	BD1-351
GSM973776	MDD-329
GSM973777	MDD-395
GSM973778	MDD-305
GSM973779	HC-45
GSM973780	HC-209
GSM973781	HC-592
GSM973782	BD1-588
GSM973783	HC-418
GSM973784	MDD-341
GSM973785	HC-84
GSM973786	HC-296
GSM973787	MDD-328
GSM973788	HC-172
GSM973789	HC-216
GSM973790	BD1-362
GSM973791	MDD-347
GSM973792	HC-321
GSM973793	HC-285
GSM973794	MDD-343
GSM973795	MDD-369
GSM973796	MDD-286
GSM973797	HC-350
GSM973798	MDD-398
GSM973799	HC-253
GSM973800	BD2-244
GSM973801	HC-262
GSM973802	MDD-580
GSM973803	MDD-407
GSM973804	MDD-573
GSM973805	MDD-402
GSM973806	BD2-338
GSM973807	HC-115
GSM973808	HC-204
GSM973809	BD1-568
GSM973810	HC-569
GSM973811	HC-279
GSM973812	MDD-566
GSM973813	HC-88
GSM973814	MDD-433
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE39653_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data" / "GPL10558_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="disease",
    phenotype_value_map={
        "healthy control": "Control",
        "major depressive disorder/ bipolar disorder": "Case",
    },
)

qc_dir = TEST_DIR / "qc_normalized"
qc_normalized_csv = run_qc_normalization(
    expression_file=TEST_DIR / "processed_data" / "annotated_renamed_unlogged.csv",
    metadata_file=metadata_csv,
    output_dir=qc_dir,
    metadata_sample_column="analysis_sample",
    phenotype_column="depression_status",
    batch_column="batch",
)

results_dir = TEST_DIR / "expression_analysis_results"
run_expression_analysis(
    input_file=qc_normalized_csv,
    output_dir=results_dir,
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_label_keywords="control",
    case_label_keywords="case",
    control_label="Control",
    case_label="Case",
    control_prefix="HC",
    case_prefix=("MDD", "BD"),
)

run_cell_proxy_analysis(
    input_file=results_dir / "deduplicated_expression_dataset.csv",
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    output_dir=TEST_DIR / "cell_population_proxies",
    markers_csv=PROJECT_ROOT / "markers" / "pbmc_cell_markers.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_label_keywords="control",
    case_label_keywords="case",
    control_label="Control",
    case_label="Case",
    control_prefix="HC",
    case_prefix=("MDD", "BD"),
)

run_cell_proxy_analysis(
    input_file=results_dir / "deduplicated_expression_dataset.csv",
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    output_dir=TEST_DIR / "cell_population_proxies_neutrophil",
    markers_csv=PROJECT_ROOT / "markers" / "blood_cell_markers.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_label_keywords="control",
    case_label_keywords="case",
    control_label="Control",
    case_label="Case",
    control_prefix="HC",
    case_prefix=("MDD", "BD"),
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
