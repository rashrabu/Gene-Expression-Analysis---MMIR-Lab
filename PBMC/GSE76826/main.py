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
GSM2038455	HC_1
GSM2038456	HC_2
GSM2038457	HC_3
GSM2038458	HC_4
GSM2038459	HC_5
GSM2038460	HC_6
GSM2038461	HC_7
GSM2038462	HC_8
GSM2038463	HC_9
GSM2038464	HC_10
GSM2038465	HC_11
GSM2038466	HC_12
GSM2038467	LOD-dep_1
GSM2038468	LOD-dep_2
GSM2038469	LOD-dep_3
GSM2038470	LOD-dep_4
GSM2038471	LOD-dep_5
GSM2038472	LOD-dep_6
GSM2038473	LOD-dep_7
GSM2038474	LOD-dep_8
GSM2038475	LOD-dep_9
GSM2038476	LOD-dep_10
GSM2038477	LOD-rem_1
GSM2038478	LOD-rem_2
GSM2038479	LOD-rem_3
GSM2038480	LOD-rem_4
GSM2038481	LOD-rem_5
GSM2038482	LOD-rem_6
GSM2038483	LOD-rem_7
GSM2038484	LOD-rem_8
GSM2038485	LOD-rem_9
GSM2038486	LOD-rem_10
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE76826_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data" / "GPL17077_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="diagnosis",
    phenotype_value_map={
        "Healthy": "Control",
        "Major depressive disorder": "Case",
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
    case_prefix="LOD",
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
    case_prefix="LOD",
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
    case_prefix="LOD",
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
