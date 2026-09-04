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
GSM7473078	PBMC RNA, control 1
GSM7473079	PBMC RNA, control 2
GSM7473080	PBMC RNA, control 3
GSM7473081	PBMC RNA, control 4
GSM7473082	PBMC RNA, control 5
GSM7473083	PBMC RNA, control 6
GSM7473084	PBMC RNA, control 7
GSM7473085	PBMC RNA, control 8
GSM7473086	PBMC RNA, control 9
GSM7473087	PBMC RNA, control 10
GSM7473088	PBMC RNA, control 11
GSM7473089	PBMC RNA, control 12
GSM7473090	PBMC RNA, ecf 1, 1st visit
GSM7473091	PBMC RNA, ecf 2, 1st visit
GSM7473092	PBMC RNA, ecf 3, 1st visit
GSM7473093	PBMC RNA, ecf 4, 1st visit
GSM7473094	PBMC RNA, ecf 1, 2nd visit
GSM7473095	PBMC RNA, ecf 2, 2nd visit
GSM7473096	PBMC RNA, ecf 3, 2nd visit
GSM7473097	PBMC RNA, ecf 4, 2nd visit
GSM7473098	PBMC RNA, ecf 5, 2nd visit
GSM7473099	PBMC RNA, kat 1, 1st visit
GSM7473100	PBMC RNA, kat 2, 1st visit
GSM7473101	PBMC RNA, kat 3, 1st visit
GSM7473102	PBMC RNA, kat 4, 1st visit
GSM7473103	PBMC RNA, kat 1, 2nd visit
GSM7473104	PBMC RNA, kat 2, 2nd visit
GSM7473105	PBMC RNA, kat 3, 2nd visit
GSM7473106	PBMC RNA, kat 4, 2nd visit
GSM7473107	PBMC RNA, sri 1, 1st visit
GSM7473108	PBMC RNA, sri 2, 1st visit
GSM7473109	PBMC RNA, sri 1, 2nd visit
GSM7473110	PBMC RNA, sri 2, 2nd visit
GSM7473111	PBMC RNA, sri 3, 2nd visit
GSM7473112	PBMC RNA, sri 4, 2nd visit
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE234669_series_matrix.txt",
    TEST_DIR / "processed_data",
    supplementary_tar_path=TEST_DIR / "original_data" / "GSE234669_RAW.tar",
    supplementary_value_format="star_gene_counts",
    supplementary_annotation_path=PROJECT_ROOT / "Ensembl" / "ensembl_gene_id_to_symbol.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="treatment",
    phenotype_value_map={
        "None": "Control",
        "electroconvulsive therapy": "Case",
        "ketamine": "Case",
        "ssri": "Case",
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
    control_prefix="PBMC RNA, control",
    case_prefix=("PBMC RNA, ecf", "PBMC RNA, kat", "PBMC RNA, sri"),
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
    control_prefix="PBMC RNA, control",
    case_prefix=("PBMC RNA, ecf", "PBMC RNA, kat", "PBMC RNA, sri"),
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
    control_prefix="PBMC RNA, control",
    case_prefix=("PBMC RNA, ecf", "PBMC RNA, kat", "PBMC RNA, sri"),
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
