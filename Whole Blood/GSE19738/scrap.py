# import pandas as pd

# soft_file = "Whole Blood/GSE19738/GPL6848_family.soft"
# output_file = "Whole Blood/GSE19738/GPL6848_annotation.tsv"

# # Read SOFT file
# with open(soft_file, "r", encoding="utf-8") as f:
#     lines = f.readlines()

# # Locate annotation table
# start = lines.index("!platform_table_begin\n") + 1
# end = lines.index("!platform_table_end\n")

# # Extract table
# annotation_lines = lines[start:end]

# # Save temporary file
# with open("GPL6848_annotation.txt", "w", encoding="utf-8") as f:
#     f.writelines(annotation_lines)

# # Read and save as TSV
# df = pd.read_csv(
#     "GPL6848_annotation.txt",
#     sep="\t",
#     dtype=str
# )

# df.to_csv(
#     output_file,
#     sep="\t",
#     index=False
# )

# print(f"Created: {output_file}")
# print(df.head())
# print(df.columns.tolist())

import pandas as pd

# Read the TSV file
df = pd.read_csv("Whole Blood/GSE19738/GPL6848_annotation.tsv", sep="\t")

# Save as a CSV file
df.to_csv("Whole Blood/GSE19738/GPL6848_annotation.csv", index=False)
