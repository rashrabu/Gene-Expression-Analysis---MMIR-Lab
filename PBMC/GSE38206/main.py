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

mapping_text = """
GSM936718	PBMC_P1-0w_rep1
GSM936719	PBMC_P1-8w_rep1
GSM936720	PBMC_P2-0w_rep1
GSM936721	PBMC_P2-8w_rep1
GSM936722	PBMC_P3-0w_rep1
GSM936723	PBMC_P3-8w_rep1
GSM936724	PBMC_P4-0w_rep1
GSM936725	PBMC_P4-8w_rep1
GSM936726	PBMC_P5-0w_rep1
GSM936727	PBMC_P5-8w_rep1
GSM936728	PBMC_P6-0w_rep1
GSM936729	PBMC_P6-8w_rep1
GSM936730	PBMC_P7-0w_rep1
GSM936731	PBMC_P7-8w_rep1
GSM936732	PBMC_P8-0w_rep1
GSM936733	PBMC_P8-8w_rep1
GSM936734	PBMC_P9-0w_rep1
GSM936735	PBMC_P9-8w_rep1
GSM936736	PBMC_C1-0w_rep1
GSM936737	PBMC_C1-8w_rep1
GSM936738	PBMC_C2-0w_rep1
GSM936739	PBMC_C2-8w_rep1
GSM936740	PBMC_C3-0w_rep1
GSM936741	PBMC_C3-8w_rep1
GSM936742	PBMC_C4-0w_rep1
GSM936743	PBMC_C4-8w_rep1
GSM936744	PBMC_C5-0w_rep1
GSM936745	PBMC_C5-8w_rep1
GSM936746	PBMC_C6-0w_rep1
GSM936747	PBMC_C6-8w_rep1
GSM936748	PBMC_C7-0w_rep1
GSM936749	PBMC_C7-8w_rep1
GSM936750	PBMC_C8-0w_rep1
GSM936751	PBMC_C8-8w_rep1
GSM936752	PBMC_C9-0w_rep1
GSM936753	PBMC_C9-8w_rep1
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE38206_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data" / "GPL13607_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="status",
    phenotype_value_map={
        "control": "Control",
        "major depressive episode": "Case",
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
    control_prefix="PBMC_C",
    case_prefix="PBMC_P",
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
    control_prefix="PBMC_C",
    case_prefix="PBMC_P",
)

# to see with neutrophil in PBMC
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
    control_prefix="PBMC_C",
    case_prefix="PBMC_P",
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
