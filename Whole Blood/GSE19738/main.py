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
GSM492850	whole blood-control-basal-CON01
GSM492851	whole blood-control-basal-CON02
GSM492852	whole blood-control-basal-CON03
GSM492853	whole blood-control-basal-CON04
GSM492854	whole blood-control-basal-CON05
GSM492855	whole blood-control-basal-CON06
GSM492856	whole blood-control-basal-CON07
GSM492857	whole blood-control-basal-CON08
GSM492858	whole blood-control-basal-CON09
GSM492859	whole blood-control-basal-CON10
GSM492860	whole blood-control-basal-CON11
GSM492861	whole blood-control-basal-CON12
GSM492862	whole blood-control-basal-CON13
GSM492863	whole blood-control-basal-CON14
GSM492864	whole blood-control-basal-CON15
GSM492865	whole blood-control-basal-CON16
GSM492866	whole blood-control-basal-CON17
GSM492867	whole blood-control-basal-CON18
GSM492868	whole blood-control-basal-CON19
GSM492869	whole blood-control-basal-CON20
GSM492870	whole blood-control-basal-CON21
GSM492871	whole blood-control-basal-CON22
GSM492872	whole blood-control-basal-CON23
GSM492873	whole blood-control-basal-CON24
GSM492874	whole blood-control-basal-CON25
GSM492875	whole blood-control-basal-CON26
GSM492876	whole blood-control-basal-CON27
GSM492877	whole blood-control-basal-CON28
GSM492878	whole blood-control-basal-CON29
GSM492879	whole blood-control-basal-CON30
GSM492880	whole blood-control-basal-CON31
GSM492881	whole blood-control-basal-CON32
GSM492882	whole blood-control-basal-CON36
GSM492883	whole blood-control-basal-CON37
GSM492884	whole blood-MDD-basal-MDD01
GSM492885	whole blood-MDD-basal-MDD02
GSM492886	whole blood-MDD-basal-MDD03
GSM492887	whole blood-MDD-basal-MDD04
GSM492888	whole blood-MDD-basal-MDD05
GSM492889	whole blood-MDD-basal-MDD06
GSM492890	whole blood-MDD-basal-MDD07
GSM492891	whole blood-MDD-basal-MDD08
GSM492892	whole blood-MDD-basal-MDD09
GSM492893	whole blood-MDD-basal-MDD10
GSM492894	whole blood-MDD-basal-MDD11
GSM492895	whole blood-MDD-basal-MDD12
GSM492896	whole blood-MDD-basal-MDD13
GSM492897	whole blood-MDD-basal-MDD14
GSM492898	whole blood-MDD-basal-MDD15
GSM492899	whole blood-MDD-basal-MDD16
GSM492900	whole blood-MDD-basal-MDD17
GSM492901	whole blood-MDD-basal-MDD18
GSM492902	whole blood-MDD-basal-MDD19
GSM492903	whole blood-MDD-basal-MDD20
GSM492904	whole blood-MDD-basal-MDD21
GSM492905	whole blood-MDD-basal-MDD22
GSM492906	whole blood-MDD-basal-MDD23
GSM492907	whole blood-MDD-basal-MDD24
GSM492908	whole blood-MDD-basal-MDD25
GSM492909	whole blood-MDD-basal-MDD26
GSM492910	whole blood-MDD-basal-MDD27
GSM492911	whole blood-MDD-basal-MDD28
GSM492912	whole blood-MDD-basal-MDD29
GSM492913	whole blood-MDD-basal-MDD30
GSM492914	whole blood-MDD-basal-MDD31
GSM492915	whole blood-MDD-basal-MDD32
GSM492916	whole blood-MDD-basal-MDD35
GSM492917	whole blood-control-LPS-± 6 h-CON01
GSM492918	whole blood-control-LPS-± 6 h-CON02
GSM492919	whole blood-control-LPS-± 6 h-CON03
GSM492920	whole blood-control-LPS-± 6 h-CON04
GSM492921	whole blood-control-LPS-± 6 h-CON05
GSM492922	whole blood-control-LPS-± 6 h-CON06
GSM492923	whole blood-control-LPS-± 6 h-CON07
GSM492924	whole blood-control-LPS-± 6 h-CON08
GSM492925	whole blood-control-LPS-± 6 h-CON09
GSM492926	whole blood-control-LPS-± 6 h-CON10
GSM492927	whole blood-control-LPS-± 6 h-CON11
GSM492928	whole blood-control-LPS-± 6 h-CON12
GSM492929	whole blood-control-LPS-± 6 h-CON13
GSM492930	whole blood-control-LPS-± 6 h-CON14
GSM492931	whole blood-control-LPS-± 6 h-CON15
GSM492932	whole blood-control-LPS-± 6 h-CON16
GSM492933	whole blood-control-LPS-± 6 h-CON17
GSM492934	whole blood-control-LPS-± 6 h-CON18
GSM492935	whole blood-control-LPS-± 6 h-CON19
GSM492936	whole blood-control-LPS-± 6 h-CON20
GSM492938	whole blood-control-LPS-± 6 h-CON22
GSM492939	whole blood-control-LPS-± 6 h-CON23
GSM492940	whole blood-control-LPS-± 6 h-CON24
GSM492941	whole blood-control-LPS-± 6 h-CON25
GSM492942	whole blood-control-LPS-± 6 h-CON26
GSM492943	whole blood-control-LPS-± 6 h-CON27
GSM492944	whole blood-control-LPS-± 6 h-CON28
GSM492945	whole blood-control-LPS-± 6 h-CON29
GSM492947	whole blood-control-LPS-± 6 h-CON31
GSM492948	whole blood-control-LPS-± 6 h-CON32
GSM492949	whole blood-control-LPS-± 6 h-CON36
GSM492950	whole blood-control-LPS-± 6 h-CON37
GSM492951	whole blood-MDD-LPS-± 6 h-MDD01
GSM492952	whole blood-MDD-LPS-± 6 h-MDD02
GSM492953	whole blood-MDD-LPS-± 6 h-MDD03
GSM492954	whole blood-MDD-LPS-± 6 h-MDD04
GSM492955	whole blood-MDD-LPS-± 6 h-MDD05
GSM492956	whole blood-MDD-LPS-± 6 h-MDD06
GSM492957	whole blood-MDD-LPS-± 6 h-MDD07
GSM492958	whole blood-MDD-LPS-± 6 h-MDD08
GSM492959	whole blood-MDD-LPS-± 6 h-MDD09
GSM492960	whole blood-MDD-LPS-± 6 h-MDD10
GSM492961	whole blood-MDD-LPS-± 6 h-MDD11
GSM492962	whole blood-MDD-LPS-± 6 h-MDD12
GSM492963	whole blood-MDD-LPS-± 6 h-MDD13
GSM492964	whole blood-MDD-LPS-± 6 h-MDD14
GSM492965	whole blood-MDD-LPS-± 6 h-MDD15
GSM492966	whole blood-MDD-LPS-± 6 h-MDD16
GSM492967	whole blood-MDD-LPS-± 6 h-MDD17
GSM492968	whole blood-MDD-LPS-± 6 h-MDD18
GSM492969	whole blood-MDD-LPS-± 6 h-MDD19
GSM492970	whole blood-MDD-LPS-± 6 h-MDD20
GSM492971	whole blood-MDD-LPS-± 6 h-MDD21
GSM492972	whole blood-MDD-LPS-± 6 h-MDD22
GSM492973	whole blood-MDD-LPS-± 6 h-MDD23
GSM492974	whole blood-MDD-LPS-± 6 h-MDD24
GSM492975	whole blood-MDD-LPS-± 6 h-MDD25
GSM492976	whole blood-MDD-LPS-± 6 h-MDD26
GSM492977	whole blood-MDD-LPS-± 6 h-MDD27
GSM492978	whole blood-MDD-LPS-± 6 h-MDD28
GSM492979	whole blood-MDD-LPS-± 6 h-MDD29
GSM492980	whole blood-MDD-LPS-± 6 h-MDD30
GSM492981	whole blood-MDD-LPS-± 6 h-MDD31
GSM492982	whole blood-MDD-LPS-± 6 h-MDD32
GSM492983	whole blood-MDD-LPS-± 6 h-MDD35
"""



expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE19738_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data/GPL6848_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    control_prefix="whole blood-control-",
    case_prefix="whole blood-MDD-",
    control_label="Control",
    case_label="Case",
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
    control_prefix="whole blood-control-",
    case_prefix="whole blood-MDD-",
    control_label="Control",
    case_label="Case",
)

run_cell_proxy_analysis(
    input_file=results_dir / "deduplicated_expression_dataset.csv",
    metadata_file=TEST_DIR / "processed_data" / "sample_metadata.csv",
    output_dir=TEST_DIR / "cell_population_proxies",
    markers_csv=PROJECT_ROOT / "markers" / "blood_cell_markers.csv",
    metadata_sample_column="analysis_sample",
    metadata_phenotype_column="depression_status",
    control_prefix="whole blood-control-",
    case_prefix="whole blood-MDD-",
    control_label_keywords="control,healthy",
    case_label_keywords="case,depressed,mdd,disease",
    control_label="Control",
    case_label="Case",
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")