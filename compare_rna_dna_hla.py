# compare_rna_dna_hla.py

# Load libraries
import pandas as pd

# -------------------------------
# INPUT FILES
# -------------------------------

sample_sheet = "/workspace/lab/houlahanlab/projects/pancancer-geb-egfr/data/alignment/rna/gdc_sample_sheet.2026-07-24_female.tsv"

rna_hla = "/workspace/lab/houlahanlab/oneillks/projects/hla_pipeline/hla/merged/GBM.genotypes.tsv"

dna_hla = "/workspace/lab/houlahanlab/resources/external_datasets/tcga_pancancer_atlas/Shukla_Wu_Getz_Polysolver_HLA_Types_2015.tsv"

# -------------------------------
# LOAD SAMPLE SHEET
# -------------------------------

samples = pd.read_csv(sample_sheet, sep="\t")
#store it as sample, read this and seperate it by tabs

# Keep only GBM
# == will respond T/F
samples = samples[samples["Project.ID"] == "TCGA-GBM"]

# Extract UUID from File.Name
#Take the File.Name column, remove everything after the first period, and save the result as a new column called subject
samples["subject"] = samples["File.Name"].str.split(".").str[0]

# Keep only columns needed for mapping
samples = samples[
    [
        "subject",
        "Case.ID",
        "Sample.ID"
    ]
]

# -------------------------------
# LOAD RNA HLA
# -------------------------------

rna = pd.read_csv(rna_hla, sep="\t")

# Attach the TCGA files to the RNA data
# Merge RNA HLA with sample sheet
rna = rna.merge(
    samples,
    on="subject",
    how="left"
)

# -------------------------------
# LOAD DNA HLA
# -------------------------------

#white space delimited 
dna = pd.read_csv(
    dna_hla,
    sep=r"\s+",
    header=None
)

# Adding column names
dna.columns = [
    "DNA_ID",
    "DNA_A1",
    "DNA_A2",
    "DNA_B1",
    "DNA_B2",
    "DNA_C1",
    "DNA_C2"
]

# Convert:
# GBM-TCGA-12-3648
# ->
# TCGA-12-3648

# Where GBM is replace it with blank space
dna["Case.ID"] = dna["DNA_ID"].str.replace(
    "GBM-",
    "",
    regex=False
)

# -------------------------------
# MERGE RNA + DNA
# -------------------------------

# Keep only ones that are in both
merged = rna.merge(
    dna,
    on="Case.ID",
    how="inner"
)

# -------------------------------
# SAVE RESULT
# -------------------------------

outfile = "GBM_RNA_DNA_HLA_mapping.tsv"

merged.to_csv(
    outfile,
    sep="\t",
    index=False
)

print(f"Merged samples: {len(merged)}")
print(f"Output written to: {outfile}")
