from Bio import AlignIO
from pathlib import Path
import csv

# ============================================================
# DREB-NETSP CORE
# Conserved Motif Analysis
# ============================================================

print("=" * 60)
print("DREB-NETSP CORE - CONSERVED MOTIF ANALYSIS")
print("=" * 60)

# ------------------------------------------------------------
# 1. Locate MSA file
# ------------------------------------------------------------

msa_file = Path(
    r"C:\Users\supra\Downloads\DREB_MSA.aln-fasta"
)

if not msa_file.exists():
    print("\nERROR: MSA file not found.")
    print(f"Expected location: {msa_file}")
    raise SystemExit

# ------------------------------------------------------------
# 2. Load alignment
# ------------------------------------------------------------

alignment = AlignIO.read(msa_file, "fasta")

print("\nAlignment loaded successfully!")
print(f"Number of sequences : {len(alignment)}")
print(f"Alignment length    : {alignment.get_alignment_length()}")

# ------------------------------------------------------------
# 3. Identify conserved positions
# ------------------------------------------------------------

conservation_threshold = 0.80

conserved_positions = []

for column_index in range(alignment.get_alignment_length()):

    column = alignment[:, column_index]

    # Ignore gaps
    residues = [aa for aa in column if aa != "-"]

    if len(residues) == 0:
        continue

    counts = {}

    for aa in residues:
        counts[aa] = counts.get(aa, 0) + 1

    most_common = max(counts, key=counts.get)
    frequency = counts[most_common] / len(residues)

    if frequency >= conservation_threshold:

        conserved_positions.append({
            "alignment_position": column_index + 1,
            "residue": most_common,
            "conservation": frequency
        })

# ------------------------------------------------------------
# 4. Print conserved positions
# ------------------------------------------------------------

print("\nCONSERVATION SUMMARY")
print("-" * 60)

print(
    f"Conserved positions (>= "
    f"{conservation_threshold * 100:.0f}%): "
    f"{len(conserved_positions)}"
)

for item in conserved_positions[:30]:

    print(
        f"Position {item['alignment_position']:>3} | "
        f"Residue {item['residue']} | "
        f"Conservation "
        f"{item['conservation'] * 100:.1f}%"
    )

# ------------------------------------------------------------
# 5. Save results
# ------------------------------------------------------------

output_file = Path(
    r"C:\Users\supra\Downloads\DREB_conserved_positions.csv"
)

with open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "alignment_position",
            "residue",
            "conservation"
        ]
    )

    writer.writeheader()
    writer.writerows(conserved_positions)

print("\nRESULTS SAVED")
print("-" * 60)
print(output_file)

print("\n" + "=" * 60)
print("CONSERVATION ANALYSIS COMPLETE")
print("=" * 60)


# ------------------------------------------------------------
# 6. Identify contiguous conserved regions
# ------------------------------------------------------------

minimum_region_length = 4

conserved_regions = []

current_region = []

for item in conserved_positions:

    if not current_region:
        current_region.append(item)

    elif item["alignment_position"] == current_region[-1]["alignment_position"] + 1:
        current_region.append(item)

    else:

        if len(current_region) >= minimum_region_length:
            conserved_regions.append(current_region)

        current_region = [item]

# Check final region
if len(current_region) >= minimum_region_length:
    conserved_regions.append(current_region)


# ------------------------------------------------------------
# 7. Print conserved regions
# ------------------------------------------------------------

print("\nCONSERVED REGIONS")
print("-" * 60)

for i, region in enumerate(conserved_regions, start=1):

    start = region[0]["alignment_position"]
    end = region[-1]["alignment_position"]

    motif = "".join(
        item["residue"] for item in region
    )

    average_conservation = sum(
        item["conservation"] for item in region
    ) / len(region)

    print(
        f"Region {i}: "
        f"{start}-{end} | "
        f"Length: {len(region)} | "
        f"Motif: {motif} | "
        f"Mean conservation: "
        f"{average_conservation * 100:.1f}%"
    )


# ------------------------------------------------------------
# 8. Save conserved regions
# ------------------------------------------------------------

regions_file = Path(
    r"C:\Users\supra\Downloads\DREB_conserved_regions.csv"
)

with open(
    regions_file,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(f)

    writer.writerow([
        "Region",
        "Start",
        "End",
        "Length",
        "Consensus_Motif",
        "Mean_Conservation"
    ])

    for i, region in enumerate(conserved_regions, start=1):

        start = region[0]["alignment_position"]
        end = region[-1]["alignment_position"]

        motif = "".join(
            item["residue"] for item in region
        )

        average_conservation = sum(
            item["conservation"] for item in region
        ) / len(region)

        writer.writerow([
            i,
            start,
            end,
            len(region),
            motif,
            round(average_conservation, 3)
        ])

print("\nCONSERVED REGIONS SAVED")
print("-" * 60)
print(regions_file)