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
GSM2612096	whole_blood_control_1
GSM2612097	whole_blood_control_2
GSM2612098	whole_blood_control_3
GSM2612099	whole_blood_control_4
GSM2612100	whole_blood_control_5
GSM2612101	whole_blood_control_6
GSM2612102	whole_blood_control_7
GSM2612103	whole_blood_control_8
GSM2612104	whole_blood_control_9
GSM2612105	whole_blood_control_10
GSM2612106	whole_blood_control_11
GSM2612107	whole_blood_control_12
GSM2612108	whole_blood_control_13
GSM2612109	whole_blood_control_14
GSM2612110	whole_blood_control_15
GSM2612111	whole_blood_control_16
GSM2612112	whole_blood_control_17
GSM2612113	whole_blood_control_18
GSM2612114	whole_blood_control_19
GSM2612115	whole_blood_control_20
GSM2612116	whole_blood_control_21
GSM2612117	whole_blood_control_22
GSM2612118	whole_blood_control_23
GSM2612119	whole_blood_control_24
GSM2612120	whole_blood_control_25
GSM2612121	whole_blood_control_26
GSM2612122	whole_blood_control_27
GSM2612123	whole_blood_control_28
GSM2612124	whole_blood_control_29
GSM2612125	whole_blood_control_30
GSM2612126	whole_blood_control_31
GSM2612127	whole_blood_control_32
GSM2612128	whole_blood_case_1
GSM2612129	whole_blood_case_2
GSM2612130	whole_blood_case_3
GSM2612131	whole_blood_case_4
GSM2612132	whole_blood_case_5
GSM2612133	whole_blood_case_6
GSM2612134	whole_blood_case_7
GSM2612135	whole_blood_case_8
GSM2612136	whole_blood_case_9
GSM2612137	whole_blood_case_10
GSM2612138	whole_blood_case_11
GSM2612139	whole_blood_case_12
GSM2612140	whole_blood_case_13
GSM2612141	whole_blood_case_14
GSM2612142	whole_blood_case_15
GSM2612143	whole_blood_case_16
GSM2612144	whole_blood_case_17
GSM2612145	whole_blood_case_18
GSM2612146	whole_blood_case_19
GSM2612147	whole_blood_case_20
GSM2612148	whole_blood_case_21
GSM2612149	whole_blood_case_22
GSM2612150	whole_blood_case_23
GSM2612151	whole_blood_case_24
GSM2612152	whole_blood_case_25
GSM2612153	whole_blood_case_26
GSM2612154	whole_blood_case_27
GSM2612155	whole_blood_case_28
GSM2612156	whole_blood_case_29
GSM2612157	whole_blood_case_30
GSM2612158	whole_blood_case_31
GSM2612159	whole_blood_case_32
GSM2612160	whole_blood_case_33
GSM2612161	whole_blood_case_34
GSM2612162	whole_blood_case_35
GSM2612163	whole_blood_case_36
GSM2612164	whole_blood_case_37
GSM2612165	whole_blood_case_38
GSM2612166	whole_blood_case_39
GSM2612167	whole_blood_case_40
GSM2612168	whole_blood_case_41
GSM2612169	whole_blood_case_42
GSM2612170	whole_blood_case_43
GSM2612171	whole_blood_case_44
GSM2612172	whole_blood_case_45
GSM2612173	whole_blood_case_46
GSM2612174	whole_blood_case_47
GSM2612175	whole_blood_case_48
GSM2612176	whole_blood_case_49
GSM2612177	whole_blood_case_50
GSM2612178	whole_blood_case_51
GSM2612179	whole_blood_case_52
GSM2612180	whole_blood_case_53
GSM2612181	whole_blood_case_54
GSM2612182	whole_blood_case_55
GSM2612183	whole_blood_case_56
GSM2612184	whole_blood_case_57
GSM2612185	whole_blood_case_58
GSM2612186	whole_blood_case_59
GSM2612187	whole_blood_case_60
GSM2612188	whole_blood_case_61
GSM2612189	whole_blood_case_62
GSM2612190	whole_blood_case_63
GSM2612191	whole_blood_case_64
GSM2612192	whole_blood_control_33
GSM2612193	whole_blood_case_65
GSM2612194	whole_blood_case_66
GSM2612195	whole_blood_case_67
GSM2612196	whole_blood_case_68
GSM2612197	whole_blood_control_34
GSM2612198	whole_blood_case_69
GSM2612199	whole_blood_control_35
GSM2612200	whole_blood_case_70
GSM2612201	whole_blood_case_71
GSM2612202	whole_blood_case_72
GSM2612203	whole_blood_control_36
GSM2612204	whole_blood_case_73
GSM2612205	whole_blood_case_74
GSM2612206	whole_blood_control_37
GSM2612207	whole_blood_case_75
GSM2612208	whole_blood_case_76
GSM2612209	whole_blood_case_77
GSM2612210	whole_blood_case_78
GSM2612211	whole_blood_case_79
GSM2612212	whole_blood_case_80
GSM2612213	whole_blood_case_81
GSM2612214	whole_blood_case_82
GSM2612215	whole_blood_case_83
GSM2612216	whole_blood_control_38
GSM2612217	whole_blood_control_39
GSM2612218	whole_blood_control_40
GSM2612219	whole_blood_case_84
GSM2612220	whole_blood_case_85
GSM2612221	whole_blood_case_86
GSM2612222	whole_blood_case_87
GSM2612223	whole_blood_case_88
GSM2612224	whole_blood_case_89
GSM2612225	whole_blood_case_90
GSM2612226	whole_blood_case_91
GSM2612227	whole_blood_control_41
GSM2612228	whole_blood_case_92
GSM2612229	whole_blood_case_93
GSM2612230	whole_blood_control_42
GSM2612231	whole_blood_case_94
GSM2612232	whole_blood_case_95
GSM2612233	whole_blood_control_43
GSM2612234	whole_blood_control_44
GSM2612235	whole_blood_case_96
GSM2612236	whole_blood_control_45
GSM2612237	whole_blood_control_46
GSM2612238	whole_blood_control_47
GSM2612239	whole_blood_control_48
GSM2612240	whole_blood_case_97
GSM2612241	whole_blood_control_49
GSM2612242	whole_blood_case_98
GSM2612243	whole_blood_case_99
GSM2612244	whole_blood_case_100
GSM2612245	whole_blood_control_50
GSM2612246	whole_blood_case_101
GSM2612247	whole_blood_case_102
GSM2612248	whole_blood_control_51
GSM2612249	whole_blood_case_103
GSM2612250	whole_blood_control_52
GSM2612251	whole_blood_control_53
GSM2612252	whole_blood_control_54
GSM2612253	whole_blood_case_104
GSM2612254	whole_blood_case_105
GSM2612255	whole_blood_case_106
GSM2612256	whole_blood_case_107
GSM2612257	whole_blood_case_108
GSM2612258	whole_blood_case_109
GSM2612259	whole_blood_control_55
GSM2612260	whole_blood_case_110
GSM2612261	whole_blood_case_111
GSM2612262	whole_blood_case_112
GSM2612263	whole_blood_case_113
GSM2612264	whole_blood_case_114
GSM2612265	whole_blood_case_115
GSM2612266	whole_blood_case_116
GSM2612267	whole_blood_case_117
GSM2612268	whole_blood_control_56
GSM2612269	whole_blood_case_118
GSM2612270	whole_blood_case_119
GSM2612271	whole_blood_case_120
GSM2612272	whole_blood_control_57
GSM2612273	whole_blood_case_121
GSM2612274	whole_blood_case_122
GSM2612275	whole_blood_control_58
GSM2612276	whole_blood_case_123
GSM2612277	whole_blood_case_124
GSM2612278	whole_blood_control_59
GSM2612279	whole_blood_control_60
GSM2612280	whole_blood_control_61
GSM2612281	whole_blood_case_125
GSM2612282	whole_blood_case_126
GSM2612283	whole_blood_case_127
GSM2612284	whole_blood_control_62
GSM2612285	whole_blood_case_128
GSM2612286	whole_blood_control_63
GSM2612287	whole_blood_control_64
"""

expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE98793_series_matrix.txt",
    TEST_DIR / "processed_data",
    annotation_path=TEST_DIR / "original_data/GPL570_annotation.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    control_prefix="whole_blood_control_",
    case_prefix="whole_blood_case_",
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
    control_prefix="whole_blood_control_",
    case_prefix="whole_blood_case_",
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
    control_prefix="whole_blood_control_",
    case_prefix="whole_blood_case_",
    control_label_keywords="control,healthy",
    case_label_keywords="case,depressed,mdd,disease",
    control_label="Control",
    case_label="Case",
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")