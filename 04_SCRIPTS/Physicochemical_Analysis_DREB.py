"""
DREB-NETSP CORE
Physicochemical Analysis of Rice DREB Proteins

Input:
    DREB_MSA.aln-fasta

Output:
    DREB_physicochemical_properties.csv

Notes:
    - The input may be a gapped multiple-sequence alignment.
    - Alignment gaps ('-' and '.') are removed before physicochemical
      properties are calculated.
    - Calculations use Biopython ProteinAnalysis.
    - Aliphatic index is calculated explicitly because ProteinAnalysis
      does not provide an aliphatic_index() method.
"""

import csv
import math
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ---------------------------------------------------------------------
# FILES
# ---------------------------------------------------------------------

INPUT_FASTA = Path(__file__).with_name("DREB_MSA.aln-fasta")
OUTPUT_CSV = Path(__file__).with_name("DREB_physicochemical_properties.csv")


# ---------------------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------------------

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------

def clean_protein_sequence(raw_sequence):
    """
    Remove alignment gaps and whitespace and convert to uppercase.
    Returns a clean amino-acid sequence.
    """
    sequence = str(raw_sequence).upper().replace("-", "").replace(".", "")
    sequence = "".join(sequence.split())

    invalid = sorted(set(sequence) - VALID_AA)

    if invalid:
        raise ValueError(
            f"Invalid amino-acid characters found: {', '.join(invalid)}"
        )

    if not sequence:
        raise ValueError("Sequence is empty after removing alignment gaps.")

    return sequence


def calculate_aliphatic_index(analysis, sequence_length):
    """
    Aliphatic index:

        X_A + 2.9*X_V + 3.9*(X_I + X_L)

    where X values are mole percentages.

    Because ProteinAnalysis.count_amino_acids() returns residue counts,
    convert the counts to percentages first.
    """
    counts = analysis.count_amino_acids()

    x_a = (counts["A"] / sequence_length) * 100
    x_v = (counts["V"] / sequence_length) * 100
    x_i = (counts["I"] / sequence_length) * 100
    x_l = (counts["L"] / sequence_length) * 100

    return x_a + (2.9 * x_v) + (3.9 * (x_i + x_l))


def safe_round(value, digits=3):
    """Return a rounded numeric value."""
    if value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return round(float(value), digits)


# ---------------------------------------------------------------------
# MAIN ANALYSIS
# ---------------------------------------------------------------------

def main():
    print("=" * 64)
    print("DREB-NETSP CORE - PHYSICOCHEMICAL ANALYSIS")
    print("=" * 64)
    print()

    if not INPUT_FASTA.exists():
        print("ERROR: Input FASTA file not found!")
        print(f"Expected file: {INPUT_FASTA}")
        print()
        print("Place DREB_MSA.aln-fasta in the same folder as this script.")
        sys.exit(1)

    records = list(SeqIO.parse(str(INPUT_FASTA), "fasta"))

    if not records:
        print("ERROR: No FASTA sequences were found.")
        sys.exit(1)

    print("Sequences loaded successfully!")
    print(f"Number of sequences: {len(records)}")
    print()

    results = []

    for record in records:
        try:
            clean_sequence = clean_protein_sequence(record.seq)
            analysis = ProteinAnalysis(clean_sequence)

            length_aa = len(clean_sequence)

            # Core physicochemical properties
            molecular_weight = analysis.molecular_weight()
            theoretical_pI = analysis.isoelectric_point()
            gravy = analysis.gravy()
            instability_index = analysis.instability_index()
            aromaticity = analysis.aromaticity()
            aliphatic_index = calculate_aliphatic_index(
                analysis, length_aa
            )

            # Composition-derived properties
            counts = analysis.count_amino_acids()

            acidic_residues = counts["D"] + counts["E"]
            basic_residues = counts["K"] + counts["R"]

            acidic_percentage = (acidic_residues / length_aa) * 100
            basic_percentage = (basic_residues / length_aa) * 100

            results.append({
                "Protein_ID": record.id,
                "Description": record.description,
                "Length_aa": length_aa,
                "Molecular_Weight_Da": safe_round(molecular_weight, 2),
                "Theoretical_pI": safe_round(theoretical_pI, 3),
                "GRAVY": safe_round(gravy, 3),
                "Instability_Index": safe_round(instability_index, 3),
                "Aliphatic_Index": safe_round(aliphatic_index, 3),
                "Aromaticity": safe_round(aromaticity, 4),
                "Acidic_Residues_DE": acidic_residues,
                "Basic_Residues_KR": basic_residues,
                "Acidic_Percentage": safe_round(acidic_percentage, 2),
                "Basic_Percentage": safe_round(basic_percentage, 2),
            })

        except Exception as exc:
            print(f"WARNING: Could not analyse {record.id}")
            print(f"Reason: {exc}")
            print()

    if not results:
        print("ERROR: No sequences could be analysed.")
        sys.exit(1)

    # -----------------------------------------------------------------
    # SAVE CSV
    # -----------------------------------------------------------------

    fieldnames = list(results[0].keys())

    with OUTPUT_CSV.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    # -----------------------------------------------------------------
    # DISPLAY SUMMARY
    # -----------------------------------------------------------------

    print("PHYSICOCHEMICAL SUMMARY")
    print("-" * 64)

    for result in results:
        print(f"\n{result['Protein_ID']}")
        print(f"  Length              : {result['Length_aa']} aa")
        print(f"  Molecular weight    : {result['Molecular_Weight_Da']} Da")
        print(f"  Theoretical pI      : {result['Theoretical_pI']}")
        print(f"  GRAVY               : {result['GRAVY']}")
        print(f"  Instability index   : {result['Instability_Index']}")
        print(f"  Aliphatic index     : {result['Aliphatic_Index']}")
        print(f"  Aromaticity         : {result['Aromaticity']}")
        print(f"  Acidic residues DE  : {result['Acidic_Residues_DE']}")
        print(f"  Basic residues KR   : {result['Basic_Residues_KR']}")

    print()
    print("=" * 64)
    print("RESULTS SAVED")
    print("=" * 64)
    print(f"CSV file: {OUTPUT_CSV}")
    print()
    print("PHYSICOCHEMICAL ANALYSIS COMPLETE")
    print("=" * 64)


if __name__ == "__main__":
    main()
