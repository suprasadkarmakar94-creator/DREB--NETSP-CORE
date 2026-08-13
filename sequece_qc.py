from Bio import SeqIO
from pathlib import Path
import csv

# ============================================================
# DREB-NETSP CORE
# Protein Sequence Quality Control
# ============================================================

print("=" * 60)
print("DREB-NETSP CORE - PROTEIN SEQUENCE QUALITY CONTROL")
print("=" * 60)

# ------------------------------------------------------------
# 1. Locate the UniProt FASTA file
# ------------------------------------------------------------

fasta_file = Path(r"C:\Users\supra\Downloads\DREB_proteins.fasta")

if not fasta_file.exists():
    print("\nERROR: FASTA file not found.")
    print(f"Expected location: {fasta_file}")
    raise SystemExit(1)

# ------------------------------------------------------------
# 2. Read FASTA sequences
# ------------------------------------------------------------

records = list(SeqIO.parse(fasta_file, "fasta"))

print(f"\nInput file: {fasta_file.name}")
print(f"Number of sequences: {len(records)}")

# ------------------------------------------------------------
# 3. Define accepted amino-acid characters
# ------------------------------------------------------------

valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")

# ------------------------------------------------------------
# 4. Perform sequence quality control
# ------------------------------------------------------------

qc_results = []
sequence_dictionary = {}

for record in records:

    sequence = str(record.seq).upper()
    sequence_id = record.id
    length = len(sequence)

    # Check for invalid amino-acid characters
    invalid_characters = sorted(
        set(sequence) - valid_amino_acids
    )

    # Check whether sequence is empty
    empty_sequence = length == 0

    # Check for duplicate sequences
    if sequence in sequence_dictionary:
        duplicate_of = sequence_dictionary[sequence]
    else:
        sequence_dictionary[sequence] = sequence_id
        duplicate_of = ""

    # Determine overall QC status
    if empty_sequence:
        status = "FAIL - Empty sequence"

    elif invalid_characters:
        status = "FAIL - Invalid characters"

    elif duplicate_of:
        status = f"WARNING - Duplicate of {duplicate_of}"

    else:
        status = "PASS"

    qc_results.append({
        "Protein_ID": sequence_id,
        "Protein_Length": length,
        "Invalid_Characters": (
            ",".join(invalid_characters)
            if invalid_characters else "None"
        ),
        "Duplicate_Of": (
            duplicate_of
            if duplicate_of else "None"
        ),
        "QC_Status": status
    })

# ------------------------------------------------------------
# 5. Display individual results
# ------------------------------------------------------------

print("\nSEQUENCE QC RESULTS")
print("-" * 60)

for result in qc_results:

    print(
        f"{result['Protein_ID']:12s} | "
        f"{result['Protein_Length']:4d} aa | "
        f"{result['QC_Status']}"
    )

# ------------------------------------------------------------
# 6. Calculate summary statistics
# ------------------------------------------------------------

total_sequences = len(qc_results)

passed_sequences = sum(
    result["QC_Status"] == "PASS"
    for result in qc_results
)

invalid_sequences = sum(
    result["Invalid_Characters"] != "None"
    for result in qc_results
)

duplicate_sequences = sum(
    result["Duplicate_Of"] != "None"
    for result in qc_results
)

empty_sequences = sum(
    result["Protein_Length"] == 0
    for result in qc_results
)

lengths = [
    result["Protein_Length"]
    for result in qc_results
    if result["Protein_Length"] > 0
]

# ------------------------------------------------------------
# 7. Display summary
# ------------------------------------------------------------

print("\nSUMMARY")
print("-" * 60)

print(f"Total sequences        : {total_sequences}")
print(f"Sequences passing QC   : {passed_sequences}")
print(f"Invalid sequences      : {invalid_sequences}")
print(f"Duplicate sequences    : {duplicate_sequences}")
print(f"Empty sequences        : {empty_sequences}")

if lengths:

    print(f"Shortest sequence      : {min(lengths)} aa")
    print(f"Longest sequence       : {max(lengths)} aa")
    print(f"Average length         : {sum(lengths) / len(lengths):.1f} aa")

# ------------------------------------------------------------
# 8. Save QC report
# ------------------------------------------------------------

output_file = Path(r"C:\Users\supra\Downloads\sequence_qc_report.csv")

with open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as csvfile:

    fieldnames = [
        "Protein_ID",
        "Protein_Length",
        "Invalid_Characters",
        "Duplicate_Of",
        "QC_Status"
    ]

    writer = csv.DictWriter(
        csvfile,
        fieldnames=fieldnames
    )

    writer.writeheader()
    writer.writerows(qc_results)

print("\nQC report saved to:")
print(output_file)

print("\n" + "=" * 60)
print("SEQUENCE QC COMPLETE")
print("=" * 60)
