from Bio import Phylo
from pathlib import Path
import csv

# ============================================================
# DREB-NETSP CORE
# Phylogenetic Tree Analysis
# ============================================================

print("=" * 65)
print("DREB-NETSP CORE - PHYLOGENETIC ANALYSIS")
print("=" * 65)


# ------------------------------------------------------------
# 1. Locate input phylogenetic tree
# ------------------------------------------------------------

tree_file = Path(
    r"C:\Users\supra\Downloads\DREB_Phylogenetic_TREE.phylotree"
)

if not tree_file.exists():
    print("\nERROR: Phylogenetic tree file not found.")
    print(f"Expected location: {tree_file}")
    raise SystemExit


# ------------------------------------------------------------
# 2. Read the phylogenetic tree
# ------------------------------------------------------------

try:
    tree = Phylo.read(tree_file, "newick")
except Exception as error:
    print("\nERROR: Could not read the phylogenetic tree.")
    print(f"Details: {error}")
    raise SystemExit


print("\nPhylogenetic tree loaded successfully!")


# ------------------------------------------------------------
# 3. Basic tree information
# ------------------------------------------------------------

terminals = tree.get_terminals()
non_terminals = tree.get_nonterminals()

print("\nTREE INFORMATION")
print("-" * 65)

print(f"Number of sequences : {len(terminals)}")
print(f"Internal nodes      : {len(non_terminals)}")


# ------------------------------------------------------------
# 4. List sequences and branch lengths
# ------------------------------------------------------------

print("\nSEQUENCE / BRANCH INFORMATION")
print("-" * 65)

tree_records = []

for terminal in terminals:

    name = terminal.name
    branch_length = terminal.branch_length

    if branch_length is None:
        branch_length = 0.0

    print(
        f"{name} | branch length: {branch_length:.5f}"
    )

    tree_records.append({
        "sequence": name,
        "branch_length": branch_length
    })


# ------------------------------------------------------------
# 5. Calculate basic tree statistics
# ------------------------------------------------------------

branch_lengths = [
    record["branch_length"]
    for record in tree_records
]

if branch_lengths:

    mean_branch_length = (
        sum(branch_lengths) / len(branch_lengths)
    )

    longest_branch = max(branch_lengths)
    shortest_branch = min(branch_lengths)

else:

    mean_branch_length = 0.0
    longest_branch = 0.0
    shortest_branch = 0.0


print("\nTREE STATISTICS")
print("-" * 65)

print(
    f"Mean terminal branch length : "
    f"{mean_branch_length:.5f}"
)

print(
    f"Shortest terminal branch    : "
    f"{shortest_branch:.5f}"
)

print(
    f"Longest terminal branch     : "
    f"{longest_branch:.5f}"
)


# ------------------------------------------------------------
# 6. Identify the most closely separated terminal pair
# ------------------------------------------------------------

print("\nPAIRWISE TREE DISTANCES")
print("-" * 65)

distance_records = []

for i in range(len(terminals)):

    for j in range(i + 1, len(terminals)):

        seq1 = terminals[i]
        seq2 = terminals[j]

        distance = tree.distance(seq1, seq2)

        distance_records.append({
            "sequence_1": seq1.name,
            "sequence_2": seq2.name,
            "tree_distance": distance
        })

        print(
            f"{seq1.name} vs {seq2.name} "
            f"| distance: {distance:.5f}"
        )


if distance_records:

    closest_pair = min(
        distance_records,
        key=lambda x: x["tree_distance"]
    )

    print("\nClosest pair in the tree:")
    print(
        f"{closest_pair['sequence_1']} "
        f"vs "
        f"{closest_pair['sequence_2']}"
    )

    print(
        f"Tree distance: "
        f"{closest_pair['tree_distance']:.5f}"
    )


# ------------------------------------------------------------
# 7. Save sequence / branch information
# ------------------------------------------------------------

output_file = Path(
    r"C:\Users\supra\Downloads\DREB_phylogenetic_summary.csv"
)

with open(
    output_file,
    "w",
    newline="",
    encoding="utf-8"
) as csvfile:

    writer = csv.DictWriter(
        csvfile,
        fieldnames=[
            "sequence",
            "branch_length"
        ]
    )

    writer.writeheader()
    writer.writerows(tree_records)


# ------------------------------------------------------------
# 8. Save pairwise tree distances
# ------------------------------------------------------------

distance_file = Path(
    r"C:\Users\supra\Downloads\DREB_tree_distances.csv"
)

with open(
    distance_file,
    "w",
    newline="",
    encoding="utf-8"
) as csvfile:

    writer = csv.DictWriter(
        csvfile,
        fieldnames=[
            "sequence_1",
            "sequence_2",
            "tree_distance"
        ]
    )

    writer.writeheader()
    writer.writerows(distance_records)


# ------------------------------------------------------------
# 9. Export a standard Newick copy
# ------------------------------------------------------------

newick_file = Path(
    r"C:\Users\supra\Downloads\DREB_phylogenetic_tree.newick"
)

Phylo.write(
    tree,
    newick_file,
    "newick"
)


# ------------------------------------------------------------
# 10. Final report
# ------------------------------------------------------------

print("\nRESULTS SAVED")
print("-" * 65)

print(f"Phylogenetic summary : {output_file}")
print(f"Tree distances       : {distance_file}")
print(f"Newick tree          : {newick_file}")

print("\n" + "=" * 65)
print("PHYLOGENETIC ANALYSIS COMPLETE")
print("=" * 65)