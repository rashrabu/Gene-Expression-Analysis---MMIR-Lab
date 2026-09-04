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
GSM8785021	total RNA-seq from sample MUC23510_S328
GSM8785022	total RNA-seq from sample MUC23422_S240
GSM8785023	total RNA-seq from sample MUC23353_S171
GSM8785024	total RNA-seq from sample MUC23268_S86
GSM8785025	total RNA-seq from sample MUC23419_S237
GSM8785026	total RNA-seq from sample MUC23247_S65
GSM8785027	total RNA-seq from sample MUC23471_S289
GSM8785028	total RNA-seq from sample MUC23252_S70
GSM8785029	total RNA-seq from sample MUC23432_S250
GSM8785030	total RNA-seq from sample MUC23259_S77
GSM8785031	total RNA-seq from sample MUC23544_S362
GSM8785032	total RNA-seq from sample MUC23251_S69
GSM8785033	total RNA-seq from sample MUC23307_S125
GSM8785034	total RNA-seq from sample MUC23481_S299
GSM8785035	total RNA-seq from sample MUC23536_S354
GSM8785036	total RNA-seq from sample MUC23225_S43
GSM8785037	total RNA-seq from sample MUC23201_S19
GSM8785038	total RNA-seq from sample MUC23489_S307
GSM8785039	total RNA-seq from sample MUC23211_S29
GSM8785040	total RNA-seq from sample MUC23406_S224
GSM8785041	total RNA-seq from sample MUC23498_S316
GSM8785042	total RNA-seq from sample MUC23313_S131
GSM8785043	total RNA-seq from sample MUC23506_S324
GSM8785044	total RNA-seq from sample MUC23452_S270
GSM8785045	total RNA-seq from sample MUC23237_S55
GSM8785046	total RNA-seq from sample MUC23202_S20
GSM8785047	total RNA-seq from sample MUC23449_S267
GSM8785048	total RNA-seq from sample MUC23390_S208
GSM8785049	total RNA-seq from sample MUC23330_S148
GSM8785050	total RNA-seq from sample MUC23547_S365
GSM8785051	total RNA-seq from sample MUC23187_S5
GSM8785052	total RNA-seq from sample MUC23192_S10
GSM8785053	total RNA-seq from sample MUC23340_S158
GSM8785054	total RNA-seq from sample MUC23335_S153
GSM8785055	total RNA-seq from sample MUC23497_S315
GSM8785056	total RNA-seq from sample MUC23426_S244
GSM8785057	total RNA-seq from sample MUC23311_S129
GSM8785058	total RNA-seq from sample MUC23545_S363
GSM8785059	total RNA-seq from sample MUC23283_S101
GSM8785060	total RNA-seq from sample MUC23370_S188
GSM8785061	total RNA-seq from sample MUC23222_S40
GSM8785062	total RNA-seq from sample MUC23372_S190
GSM8785063	total RNA-seq from sample MUC23559_S377
GSM8785064	total RNA-seq from sample MUC23468_S286
GSM8785065	total RNA-seq from sample MUC23408_S226
GSM8785066	total RNA-seq from sample MUC23229_S47
GSM8785067	total RNA-seq from sample MUC23513_S331
GSM8785068	total RNA-seq from sample MUC23381_S199
GSM8785069	total RNA-seq from sample MUC23206_S24
GSM8785070	total RNA-seq from sample MUC23280_S98
GSM8785071	total RNA-seq from sample MUC23336_S154
GSM8785072	total RNA-seq from sample MUC23236_S54
GSM8785073	total RNA-seq from sample MUC23217_S35
GSM8785074	total RNA-seq from sample MUC23246_S64
GSM8785075	total RNA-seq from sample MUC23453_S271
GSM8785076	total RNA-seq from sample MUC23509_S327
GSM8785077	total RNA-seq from sample MUC23443_S261
GSM8785078	total RNA-seq from sample MUC23546_S364
GSM8785079	total RNA-seq from sample MUC23474_S292
GSM8785080	total RNA-seq from sample MUC23548_S366
GSM8785081	total RNA-seq from sample MUC23279_S97
GSM8785082	total RNA-seq from sample MUC23485_S303
GSM8785083	total RNA-seq from sample MUC23238_S56
GSM8785084	total RNA-seq from sample MUC23413_S231
GSM8785085	total RNA-seq from sample MUC23344_S162
GSM8785086	total RNA-seq from sample MUC23433_S251
GSM8785087	total RNA-seq from sample MUC23529_S347
GSM8785088	total RNA-seq from sample MUC23429_S247
GSM8785089	total RNA-seq from sample MUC23436_S254
GSM8785090	total RNA-seq from sample MUC23368_S186
GSM8785091	total RNA-seq from sample MUC23418_S236
GSM8785092	total RNA-seq from sample MUC23508_S326
GSM8785093	total RNA-seq from sample MUC23188_S6
GSM8785094	total RNA-seq from sample MUC23387_S205
GSM8785095	total RNA-seq from sample MUC23405_S223
GSM8785096	total RNA-seq from sample MUC23440_S258
GSM8785097	total RNA-seq from sample MUC23235_S53
GSM8785098	total RNA-seq from sample MUC23321_S139
GSM8785099	total RNA-seq from sample MUC23329_S147
GSM8785100	total RNA-seq from sample MUC23400_S218
GSM8785101	total RNA-seq from sample MUC23257_S75
GSM8785102	total RNA-seq from sample MUC23382_S200
GSM8785103	total RNA-seq from sample MUC23376_S194
GSM8785104	total RNA-seq from sample MUC23505_S323
GSM8785105	total RNA-seq from sample MUC23480_S298
GSM8785106	total RNA-seq from sample MUC23241_S59
GSM8785107	total RNA-seq from sample MUC23343_S161
GSM8785108	total RNA-seq from sample MUC23521_S339
GSM8785109	total RNA-seq from sample MUC23560_S378
GSM8785110	total RNA-seq from sample MUC23261_S79
GSM8785111	total RNA-seq from sample MUC23286_S104
GSM8785112	total RNA-seq from sample MUC23534_S352
GSM8785113	total RNA-seq from sample MUC23384_S202
GSM8785114	total RNA-seq from sample MUC23324_S142
GSM8785115	total RNA-seq from sample MUC23447_S265
GSM8785116	total RNA-seq from sample MUC23231_S49
GSM8785117	total RNA-seq from sample MUC23345_S163
GSM8785118	total RNA-seq from sample MUC23386_S204
GSM8785119	total RNA-seq from sample MUC23379_S197
GSM8785120	total RNA-seq from sample MUC23278_S96
GSM8785121	total RNA-seq from sample MUC23563_S381
GSM8785122	total RNA-seq from sample MUC23269_S87
GSM8785123	total RNA-seq from sample MUC23403_S221
GSM8785124	total RNA-seq from sample MUC23535_S353
GSM8785125	total RNA-seq from sample MUC23439_S257
GSM8785126	total RNA-seq from sample MUC23541_S359
GSM8785127	total RNA-seq from sample MUC23328_S146
GSM8785128	total RNA-seq from sample MUC23460_S278
GSM8785129	total RNA-seq from sample MUC23437_S255
GSM8785130	total RNA-seq from sample MUC23319_S137
GSM8785131	total RNA-seq from sample MUC23256_S74
GSM8785132	total RNA-seq from sample MUC23250_S68
GSM8785133	total RNA-seq from sample MUC23215_S33
GSM8785134	total RNA-seq from sample MUC23219_S37
GSM8785135	total RNA-seq from sample MUC23504_S322
GSM8785136	total RNA-seq from sample MUC23520_S338
GSM8785137	total RNA-seq from sample MUC23484_S302
GSM8785138	total RNA-seq from sample MUC23446_S264
GSM8785139	total RNA-seq from sample MUC23378_S196
GSM8785140	total RNA-seq from sample MUC23184_S2
GSM8785141	total RNA-seq from sample MUC23438_S256
GSM8785142	total RNA-seq from sample MUC23431_S249
GSM8785143	total RNA-seq from sample MUC23310_S128
GSM8785144	total RNA-seq from sample MUC23226_S44
GSM8785145	total RNA-seq from sample MUC23190_S8
GSM8785146	total RNA-seq from sample MUC23483_S301
GSM8785147	total RNA-seq from sample MUC23216_S34
GSM8785148	total RNA-seq from sample MUC23320_S138
GSM8785149	total RNA-seq from sample MUC23500_S318
GSM8785150	total RNA-seq from sample MUC23374_S192
GSM8785151	total RNA-seq from sample MUC23553_S371
GSM8785152	total RNA-seq from sample MUC23199_S17
GSM8785153	total RNA-seq from sample MUC23200_S18
GSM8785154	total RNA-seq from sample MUC23467_S285
GSM8785155	total RNA-seq from sample MUC23416_S234
GSM8785156	total RNA-seq from sample MUC23424_S242
GSM8785157	total RNA-seq from sample MUC23463_S281
GSM8785158	total RNA-seq from sample MUC23248_S66
GSM8785159	total RNA-seq from sample MUC23341_S159
GSM8785160	total RNA-seq from sample MUC23315_S133
GSM8785161	total RNA-seq from sample MUC23493_S311
GSM8785162	total RNA-seq from sample MUC23294_S112
GSM8785163	total RNA-seq from sample MUC23465_S283
GSM8785164	total RNA-seq from sample MUC23212_S30
GSM8785165	total RNA-seq from sample MUC23410_S228
GSM8785166	total RNA-seq from sample MUC23281_S99
GSM8785167	total RNA-seq from sample MUC23186_S4
GSM8785168	total RNA-seq from sample MUC23260_S78
GSM8785169	total RNA-seq from sample MUC23472_S290
GSM8785170	total RNA-seq from sample MUC23531_S349
GSM8785171	total RNA-seq from sample MUC23220_S38
GSM8785172	total RNA-seq from sample MUC23218_S36
GSM8785173	total RNA-seq from sample MUC23197_S15
GSM8785174	total RNA-seq from sample MUC23317_S135
GSM8785175	total RNA-seq from sample MUC23264_S82
GSM8785176	total RNA-seq from sample MUC23363_S181
GSM8785177	total RNA-seq from sample MUC23293_S111
GSM8785178	total RNA-seq from sample MUC23331_S149
GSM8785179	total RNA-seq from sample MUC23503_S321
GSM8785180	total RNA-seq from sample MUC23409_S227
GSM8785181	total RNA-seq from sample MUC23242_S60
GSM8785182	total RNA-seq from sample MUC23492_S310
GSM8785183	total RNA-seq from sample MUC23398_S216
GSM8785184	total RNA-seq from sample MUC23297_S115
GSM8785185	total RNA-seq from sample MUC23457_S275
GSM8785186	total RNA-seq from sample MUC23558_S376
GSM8785187	total RNA-seq from sample MUC23477_S295
GSM8785188	total RNA-seq from sample MUC23308_S126
GSM8785189	total RNA-seq from sample MUC23425_S243
GSM8785190	total RNA-seq from sample MUC23475_S293
GSM8785191	total RNA-seq from sample MUC23309_S127
GSM8785192	total RNA-seq from sample MUC23334_S152
GSM8785193	total RNA-seq from sample MUC23517_S335
GSM8785194	total RNA-seq from sample MUC23528_S346
GSM8785195	total RNA-seq from sample MUC23232_S50
GSM8785196	total RNA-seq from sample MUC23407_S225
GSM8785197	total RNA-seq from sample MUC23448_S266
GSM8785198	total RNA-seq from sample MUC23333_S151
GSM8785199	total RNA-seq from sample MUC23423_S241
GSM8785200	total RNA-seq from sample MUC23348_S166
GSM8785201	total RNA-seq from sample MUC23228_S46
GSM8785202	total RNA-seq from sample MUC23496_S314
GSM8785203	total RNA-seq from sample MUC23527_S345
GSM8785204	total RNA-seq from sample MUC23359_S177
GSM8785205	total RNA-seq from sample MUC23209_S27
GSM8785206	total RNA-seq from sample MUC23537_S355
GSM8785207	total RNA-seq from sample MUC23249_S67
GSM8785208	total RNA-seq from sample MUC23267_S85
GSM8785209	total RNA-seq from sample MUC23193_S11
GSM8785210	total RNA-seq from sample MUC23499_S317
GSM8785211	total RNA-seq from sample MUC23522_S340
GSM8785212	total RNA-seq from sample MUC23542_S360
GSM8785213	total RNA-seq from sample MUC23287_S105
GSM8785214	total RNA-seq from sample MUC23358_S176
GSM8785215	total RNA-seq from sample MUC23476_S294
GSM8785216	total RNA-seq from sample MUC23263_S81
GSM8785217	total RNA-seq from sample MUC23430_S248
GSM8785218	total RNA-seq from sample MUC23532_S350
GSM8785219	total RNA-seq from sample MUC23411_S229
GSM8785220	total RNA-seq from sample MUC23213_S31
GSM8785221	total RNA-seq from sample MUC23332_S150
GSM8785222	total RNA-seq from sample MUC23208_S26
GSM8785223	total RNA-seq from sample MUC23183_S1
GSM8785224	total RNA-seq from sample MUC23354_S172
GSM8785225	total RNA-seq from sample MUC23551_S369
GSM8785226	total RNA-seq from sample MUC23487_S305
GSM8785227	total RNA-seq from sample MUC23325_S143
GSM8785228	total RNA-seq from sample MUC23414_S232
GSM8785229	total RNA-seq from sample MUC23360_S178
GSM8785230	total RNA-seq from sample MUC23501_S319
GSM8785231	total RNA-seq from sample MUC23198_S16
GSM8785232	total RNA-seq from sample MUC23478_S296
GSM8785233	total RNA-seq from sample MUC23339_S157
GSM8785234	total RNA-seq from sample MUC23210_S28
GSM8785235	total RNA-seq from sample MUC23445_S263
GSM8785236	total RNA-seq from sample MUC23361_S179
GSM8785237	total RNA-seq from sample MUC23262_S80
GSM8785238	total RNA-seq from sample MUC23486_S304
GSM8785239	total RNA-seq from sample MUC23479_S297
GSM8785240	total RNA-seq from sample MUC23322_S140
GSM8785241	total RNA-seq from sample MUC23357_S175
GSM8785242	total RNA-seq from sample MUC23265_S83
GSM8785243	total RNA-seq from sample MUC23346_S164
GSM8785244	total RNA-seq from sample MUC23458_S276
GSM8785245	total RNA-seq from sample MUC23355_S173
GSM8785246	total RNA-seq from sample MUC23254_S72
GSM8785247	total RNA-seq from sample MUC23323_S141
GSM8785248	total RNA-seq from sample MUC23191_S9
GSM8785249	total RNA-seq from sample MUC23230_S48
GSM8785250	total RNA-seq from sample MUC23412_S230
GSM8785251	total RNA-seq from sample MUC23524_S342
GSM8785252	total RNA-seq from sample MUC23326_S144
GSM8785253	total RNA-seq from sample MUC23543_S361
GSM8785254	total RNA-seq from sample MUC23462_S280
GSM8785255	total RNA-seq from sample MUC23196_S14
GSM8785256	total RNA-seq from sample MUC23494_S312
GSM8785257	total RNA-seq from sample MUC23312_S130
GSM8785258	total RNA-seq from sample MUC23455_S273
GSM8785259	total RNA-seq from sample MUC23538_S356
GSM8785260	total RNA-seq from sample MUC23224_S42
GSM8785261	total RNA-seq from sample MUC23373_S191
GSM8785262	total RNA-seq from sample MUC23282_S100
GSM8785263	total RNA-seq from sample MUC23417_S235
GSM8785264	total RNA-seq from sample MUC23397_S215
GSM8785265	total RNA-seq from sample MUC23470_S288
GSM8785266	total RNA-seq from sample MUC23388_S206
GSM8785267	total RNA-seq from sample MUC23338_S156
GSM8785268	total RNA-seq from sample MUC23285_S103
GSM8785269	total RNA-seq from sample MUC23377_S195
GSM8785270	total RNA-seq from sample MUC23253_S71
GSM8785271	total RNA-seq from sample MUC23451_S269
GSM8785272	total RNA-seq from sample MUC23415_S233
GSM8785273	total RNA-seq from sample MUC23302_S120
GSM8785274	total RNA-seq from sample MUC23561_S379
GSM8785275	total RNA-seq from sample MUC23223_S41
GSM8785276	total RNA-seq from sample MUC23349_S167
GSM8785277	total RNA-seq from sample MUC23350_S168
GSM8785278	total RNA-seq from sample MUC23482_S300
GSM8785279	total RNA-seq from sample MUC23495_S313
GSM8785280	total RNA-seq from sample MUC23461_S279
GSM8785281	total RNA-seq from sample MUC23533_S351
GSM8785282	total RNA-seq from sample MUC23362_S180
GSM8785283	total RNA-seq from sample MUC23364_S182
GSM8785284	total RNA-seq from sample MUC23556_S374
GSM8785285	total RNA-seq from sample MUC23185_S3
GSM8785286	total RNA-seq from sample MUC23507_S325
GSM8785287	total RNA-seq from sample MUC23512_S330
GSM8785288	total RNA-seq from sample MUC23530_S348
GSM8785289	total RNA-seq from sample MUC23434_S252
GSM8785290	total RNA-seq from sample MUC23292_S110
GSM8785291	total RNA-seq from sample MUC23304_S122
GSM8785292	total RNA-seq from sample MUC23557_S375
GSM8785293	total RNA-seq from sample MUC23288_S106
GSM8785294	total RNA-seq from sample MUC23464_S282
GSM8785295	total RNA-seq from sample MUC23375_S193
GSM8785296	total RNA-seq from sample MUC23435_S253
GSM8785297	total RNA-seq from sample MUC23207_S25
GSM8785298	total RNA-seq from sample MUC23421_S239
GSM8785299	total RNA-seq from sample MUC23233_S51
GSM8785300	total RNA-seq from sample MUC23420_S238
GSM8785301	total RNA-seq from sample MUC23456_S274
GSM8785302	total RNA-seq from sample MUC23244_S62
GSM8785303	total RNA-seq from sample MUC23515_S333
GSM8785304	total RNA-seq from sample MUC23347_S165
GSM8785305	total RNA-seq from sample MUC23502_S320
GSM8785306	total RNA-seq from sample MUC23399_S217
GSM8785307	total RNA-seq from sample MUC23490_S308
GSM8785308	total RNA-seq from sample MUC23337_S155
GSM8785309	total RNA-seq from sample MUC23396_S214
GSM8785310	total RNA-seq from sample MUC23245_S63
GSM8785311	total RNA-seq from sample MUC23564_S382
GSM8785312	total RNA-seq from sample MUC23314_S132
GSM8785313	total RNA-seq from sample MUC23550_S368
GSM8785314	total RNA-seq from sample MUC23391_S209
GSM8785315	total RNA-seq from sample MUC23365_S183
GSM8785316	total RNA-seq from sample MUC23380_S198
GSM8785317	total RNA-seq from sample MUC23195_S13
GSM8785318	total RNA-seq from sample MUC23392_S210
GSM8785319	total RNA-seq from sample MUC23356_S174
"""


expression_csv, metadata_csv = run_preprocessor(
    mapping_text,
    TEST_DIR / "original_data/GSE289146_series_matrix.txt",
    TEST_DIR / "processed_data",
    supplementary_tar_path=TEST_DIR / "original_data" / "GSE289146_RAW.tar",
    supplementary_value_format="rsem_gene_results",
    supplementary_value_column="TPM",
    supplementary_annotation_path=PROJECT_ROOT / "Ensembl" / "ensembl_gene_id_to_symbol.csv",
    mapping_strip_prefix=None,
    sample_title_prefix=None,
    phenotype_column_name="depression_status",
    phenotype_characteristic_name="status",
    phenotype_value_map={
        "unaffected": "Control",
        "affected": "Case",
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
)

print(f"Expression output: {expression_csv}")
print(f"Metadata output:   {metadata_csv}")
print(f"Analysis outputs:  {TEST_DIR / 'cell_population_proxies'}")
