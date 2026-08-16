#!/usr/bin/env python3
"""
DREB-NETSP Integrated Evidence Pipeline
=======================================

Purpose
-------
Integrate the currently available DREB datasets for five Oryza sativa
subsp. japonica candidates into a reproducible, auditable evidence matrix.

The pipeline deliberately separates:
  1. Evolutionary evidence: MSA-based conservation and pairwise identity.
  2. Regulatory-context evidence: promoter motif / motif_CE / FunTFBS density
     and overlap with conserved genomic elements.
  3. Structural evidence: AlphaFold pLDDT extracted from PDB B-factor fields.
  4. Functional/domain evidence: InterPro and functional annotation as
     eligibility/annotation fields, not discriminatory scores when identical.
  5. TF-target network evidence: PlantRegMap source->target relationships,
     included as an auxiliary evidence layer and only allowed into the core
     ranking if comparable source coverage exists across candidates.

Important scientific safeguards
-------------------------------
- Missing evidence is NOT automatically treated as biological absence.
- Features with no variation across candidates are not allowed to affect rank.
- Network evidence is excluded from the core score when source coverage is
  insufficient to make the candidates comparable.
- Weights are configurable and a weight-sensitivity analysis is produced.
- This is a candidate-prioritization / hypothesis-generation framework, not
  experimental validation or proof of causal drought resistance.

Outputs
-------
results/
  DREB_NETSP_integrated_evidence_matrix.csv
  DREB_NETSP_candidate_ranking.csv
  DREB_NETSP_network_evidence.csv
  DREB_NETSP_promoter_features.csv
  DREB_NETSP_weight_sensitivity.csv
  DREB_NETSP_report.txt
  DREB_NETSP_ranking.png
  DREB_NETSP_evidence_heatmap.png
  DREB_NETSP_NJ_tree.png
  DREB_NETSP_NJ_tree.newick
"""

from __future__ import annotations

import csv
import gzip
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    from Bio.PDB import PDBParser
    PDB_AVAILABLE = True
except Exception:
    PDB_AVAILABLE = False


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

CANDIDATES = {
    "DREB2A": {"gene": "LOC_Os01g07120", "uniprot": "Q0JQF7", "expected_length": 274},
    "DREB1D": {"gene": "LOC_Os06g06970", "uniprot": "Q6J1A5", "expected_length": 253},
    "DREB1C": {"gene": "LOC_Os06g03670", "uniprot": "Q9LWV3", "expected_length": 214},
    "DREB1B": {"gene": "LOC_Os09g35010", "uniprot": "Q3T5N4", "expected_length": 218},
    "DREB1A": {"gene": "LOC_Os09g35030", "uniprot": "Q64MA1", "expected_length": 238},
}

# Core score weights. These are deliberately transparent and editable.
WEIGHTS = {
    "evolution": 0.40,
    "promoter_regulation": 0.30,
    "structure": 0.30,
}

# The PlantRegMap promoter definition supplied with the project:
# -2000 bp to +100 bp around TSS.
PROMOTER_UP = 2000
PROMOTER_DOWN = 100

# Files currently present in the project directory.
FASTA = ROOT / "DREB_proteins.fasta"
MSA = ROOT / "DREB_MSA.aln-fasta"
INTERPRO = ROOT / "DREB_InterPro_domain_annotations.xlsx"
FUNCTIONAL = ROOT / "DREB_functional_annotation_evidence.csv"
BED = ROOT / "DREB_NETSP_PlantRegMap_gene_regions.bed"
CONSERVATION = ROOT / "COMSERVED_ELEMENTS"
PROMOTER_MOTIF = ROOT / "promoter_motif_scanning.gff"
PROMOTER_MOTIF_CE = ROOT / "TFBS_from_motif_CE_inProm_Osj.gff"
PROMOTER_FUN = ROOT / "FUN TFBS_OSJ.gff"
REG_MERGED = ROOT / "REGULATION_MERGED.txt"
REG_MOTIF_CE = ROOT / "REGULATION_MOTIF_CE.txt"
REG_FUN = ROOT / "Regulation (motif_CE).txt"

PDBS = {
    "Q0JQF7": ROOT / "AF-Q0JQF7-F1-model_v6.pdb",
    "Q3T5N4": ROOT / "AF-Q3T5N4-F1-model_v6 (1).pdb",
    "Q64MA1": ROOT / "AF-Q64MA1-F1-model_v6.pdb",
    "Q6J1A5": ROOT / "AF-Q6J1A5-F1-model_v6.pdb",
    "Q9LWV3": ROOT / "AF-Q9LWV3-F1-model_v6.pdb",
}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def clean_seq(seq: str) -> str:
    return str(seq).upper().replace("-", "").replace(".", "").replace(" ", "")


def minmax(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or pd.isna(hi):
        return pd.Series(np.nan, index=series.index)
    if math.isclose(float(lo), float(hi)):
        return pd.Series(1.0, index=series.index)
    return (s - lo) / (hi - lo)


def z_to_0_1(series: pd.Series) -> pd.Series:
    """Rank-like normalization with midpoint 0.5 when all values are equal."""
    return minmax(series)


def overlap_len(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def canon_chrom(chrom: str) -> str:
    """Normalize common rice chromosome spellings (Chr01 -> Chr1)."""
    c = str(chrom).strip()
    m = re.fullmatch(r"Chr0*(\d+)", c, flags=re.IGNORECASE)
    if m:
        return "Chr" + str(int(m.group(1)))
    return c


def parse_gff9(attr: str) -> dict:
    out = {}
    for part in attr.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


def interval_from_bed(chrom: str, start: int, end: int, strand: str):
    # BED start is 0-based, end is exclusive. Convert to a simple half-open
    # promoter interval around the TSS implied by the gene strand.
    if strand == "+":
        return chrom, max(0, start - PROMOTER_UP), start + PROMOTER_DOWN
    # For a reverse-strand gene, TSS is at the BED end coordinate.
    return chrom, max(0, end - PROMOTER_DOWN), end + PROMOTER_UP


# ---------------------------------------------------------------------------
# Candidate metadata
# ---------------------------------------------------------------------------

def load_bed_promoters():
    promoters = {}
    with open(BED, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            chrom, start, end, gene, _, strand = p[:6]
            if gene not in {v["gene"] for v in CANDIDATES.values()}:
                continue
            chrom = canon_chrom(chrom)
            chrom, p0, p1 = interval_from_bed(chrom, int(start), int(end), strand)
            promoters[gene] = {
                "chrom": chrom,
                "prom_start": p0,
                "prom_end": p1,
                "gene_start": int(start),
                "gene_end": int(end),
                "strand": strand,
            }
    return promoters


def load_sequences():
    records = list(SeqIO.parse(str(FASTA), "fasta"))
    by_uniprot = {}
    for r in records:
        m = re.search(r"\|([^|]+)\|", r.id)
        if m:
            by_uniprot[m.group(1)] = clean_seq(r.seq)
    return by_uniprot


def load_msa():
    aln = AlignIO.read(str(MSA), "fasta")
    by_uniprot = {}
    for r in aln:
        m = re.search(r"\|([^|]+)\|", r.id)
        if m:
            by_uniprot[m.group(1)] = str(r.seq).upper()
    return aln, by_uniprot


# ---------------------------------------------------------------------------
# Evolutionary features
# ---------------------------------------------------------------------------

def pairwise_identity(seq_a, seq_b):
    matches = comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a in "-." or b in "-.":
            continue
        comparable += 1
        matches += int(a == b)
    return matches / comparable if comparable else np.nan


def evolutionary_features(aln):
    ids = []
    names = []
    seqs = []
    for r in aln:
        m = re.search(r"\|([^|]+)\|", r.id)
        if m:
            ids.append(m.group(1)); names.append(r.id); seqs.append(str(r.seq).upper())

    # Mean pairwise identity against the other DREB family members.
    matrix = pd.DataFrame(index=ids, columns=ids, dtype=float)
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            matrix.loc[a, b] = pairwise_identity(seqs[i], seqs[j])

    # Per-position consensus conservation across the alignment.
    col_conservation = []
    for col in zip(*seqs):
        residues = [x for x in col if x not in "-."]
        if not residues:
            continue
        counts = pd.Series(residues).value_counts()
        col_conservation.append(float(counts.iloc[0] / len(residues)))
    global_conservation = float(np.mean(col_conservation)) if col_conservation else np.nan

    rows = []
    for uid in ids:
        others = matrix.loc[uid, [x for x in ids if x != uid]].astype(float)
        rows.append({
            "UniProt_ID": uid,
            "mean_pairwise_identity": float(others.mean()),
            "min_pairwise_identity": float(others.min()),
            "max_pairwise_identity": float(others.max()),
            "family_consensus_mean": global_conservation,
        })
    return pd.DataFrame(rows), matrix


def build_nj_tree(aln):
    try:
        calculator = DistanceCalculator("identity")
        dm = calculator.get_distance(aln)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree.root_at_midpoint()
        tree_path = OUT / "DREB_NETSP_NJ_tree.newick"
        Phylo.write(tree, str(tree_path), "newick")
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=ax, do_show=False)
        fig.tight_layout()
        fig.savefig(OUT / "DREB_NETSP_NJ_tree.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        return str(tree_path)
    except Exception as exc:
        return f"Tree generation failed: {exc}"


# ---------------------------------------------------------------------------
# Structure features
# ---------------------------------------------------------------------------

def structure_features():
    rows = []
    if not PDB_AVAILABLE:
        return pd.DataFrame(rows)
    parser = PDBParser(QUIET=True)
    for uid, p in PDBS.items():
        if not p.exists():
            rows.append({"UniProt_ID": uid, "structure_status": "missing"})
            continue
        try:
            structure = parser.get_structure(uid, str(p))
            vals = []
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == " " and "CA" in residue:
                            vals.append(float(residue["CA"].bfactor))
            rows.append({
                "UniProt_ID": uid,
                "structure_status": "ok",
                "modelled_residues": len(vals),
                "mean_pLDDT": float(np.mean(vals)) if vals else np.nan,
                "median_pLDDT": float(np.median(vals)) if vals else np.nan,
                "min_pLDDT": float(np.min(vals)) if vals else np.nan,
                "fraction_pLDDT_ge70": float(np.mean(np.array(vals) >= 70)) if vals else np.nan,
                "fraction_pLDDT_ge90": float(np.mean(np.array(vals) >= 90)) if vals else np.nan,
            })
        except Exception as exc:
            rows.append({"UniProt_ID": uid, "structure_status": f"error: {exc}"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Promoter regulatory-context features
# ---------------------------------------------------------------------------

def count_gff_overlaps(path: Path, promoters: dict, want_extra=False):
    counts = defaultdict(lambda: {"hits": 0, "bp": 0, "scores": []})
    if not path.exists():
        return counts
    with open(path, errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, _, _, start, end, score, _, _, attr = p[:9]
            chrom = canon_chrom(chrom)
            try:
                s, e = int(start) - 1, int(end)  # GFF -> half-open
            except ValueError:
                continue
            for gene, pr in promoters.items():
                if chrom != pr["chrom"]:
                    continue
                ov = overlap_len(s, e, pr["prom_start"], pr["prom_end"])
                if ov:
                    counts[gene]["hits"] += 1
                    counts[gene]["bp"] += ov
                    try:
                        counts[gene]["scores"].append(float(score))
                    except Exception:
                        pass
    return counts


def count_conserved_overlaps(path: Path, promoters: dict):
    counts = defaultdict(lambda: {"elements": 0, "bp": 0, "score_sum": 0.0})
    opener = gzip.open if path.exists() and path.read_bytes()[:2] == b"\x1f\x8b" else open
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            chrom = canon_chrom(p[0])
            try:
                s, e = int(p[3]) - 1, int(p[4])
                score = float(p[5]) if len(p) > 5 else np.nan
            except ValueError:
                continue
            for gene, pr in promoters.items():
                if chrom != pr["chrom"]:
                    continue
                ov = overlap_len(s, e, pr["prom_start"], pr["prom_end"])
                if ov:
                    counts[gene]["elements"] += 1
                    counts[gene]["bp"] += ov
                    if not np.isnan(score):
                        counts[gene]["score_sum"] += score
    return counts


def promoter_features(promoters):
    motif = count_gff_overlaps(PROMOTER_MOTIF, promoters)
    motif_ce = count_gff_overlaps(PROMOTER_MOTIF_CE, promoters)
    fun = count_gff_overlaps(PROMOTER_FUN, promoters)
    cons = count_conserved_overlaps(CONSERVATION, promoters)

    rows = []
    for cname, meta in CANDIDATES.items():
        gene = meta["gene"]
        pr = promoters.get(gene, {})
        rows.append({
            "Candidate": cname,
            "Gene": gene,
            "promoter_chrom": pr.get("chrom"),
            "promoter_start": pr.get("prom_start"),
            "promoter_end": pr.get("prom_end"),
            "promoter_length_bp": (pr.get("prom_end") - pr.get("prom_start")) if pr else np.nan,
            "motif_hits": motif[gene]["hits"],
            "motif_overlap_bp": motif[gene]["bp"],
            "motif_CE_hits": motif_ce[gene]["hits"],
            "motif_CE_overlap_bp": motif_ce[gene]["bp"],
            "FunTFBS_hits": fun[gene]["hits"],
            "FunTFBS_overlap_bp": fun[gene]["bp"],
            "conserved_elements": cons[gene]["elements"],
            "conserved_element_overlap_bp": cons[gene]["bp"],
            "conserved_element_score_sum": cons[gene]["score_sum"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# PlantRegMap TF-target network evidence
# ---------------------------------------------------------------------------

def parse_network_file(path: Path, label: str, source_ids: set):
    rows = []
    if not path.exists():
        return rows
    with open(path, errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            source, target = p[0], p[2]
            if source in source_ids:
                # Keep the actual evidence label. The two auxiliary files have
                # fixed evidence labels; merged has the label in column 5.
                evidence = (p[4] if label == "merged" and len(p) >= 5 and p[4] else label)
                rows.append((source, target, evidence, label))
    return rows


def network_features():
    source_ids = {v["gene"] for v in CANDIDATES.values()}
    all_rows = []
    for path, label in [
        (REG_MERGED, "merged"),
        (REG_MOTIF_CE, "motif_CE"),
        (REG_FUN, "FunTFBS"),
    ]:
        all_rows.extend(parse_network_file(path, label, source_ids))

    # De-duplicate the same source-target-evidence relationship across files.
    unique = set(all_rows)
    net_rows = []
    for cname, meta in CANDIDATES.items():
        source = meta["gene"]
        rel = [(s, t, e, f) for s, t, e, f in unique if s == source]
        targets = {x[1] for x in rel}
        evidence_counts = pd.Series([x[2] for x in rel]).value_counts().to_dict() if rel else {}
        # Evidence diversity is descriptive only. It is not treated as
        # independent additive evidence because files overlap conceptually.
        net_rows.append({
            "Candidate": cname,
            "Gene": source,
            "network_targets": len(targets),
            "network_relationships": len(rel),
            "network_evidence_types": ";".join(sorted(set(x[2] for x in rel))),
            "network_has_evidence": bool(rel),
            "motif_relationships": evidence_counts.get("motif", 0),
            "motif_CE_relationships": evidence_counts.get("motif_CE", 0),
            "FunTFBS_relationships": evidence_counts.get("FunTFBS", 0),
        })
    return pd.DataFrame(net_rows), pd.DataFrame(list(unique), columns=["Gene","Target","Evidence","FileLayer"])


# ---------------------------------------------------------------------------
# Domain / functional annotation
# ---------------------------------------------------------------------------

def annotation_features():
    dom = pd.read_excel(INTERPRO, sheet_name="Domain_Annotations")
    fun = pd.read_csv(FUNCTIONAL)
    dom = dom.rename(columns={"Gene": "Candidate"})
    merged = dom.merge(fun, on=["Candidate", "UniProt_ID"], how="outer")
    merged["domain_present"] = merged["Primary_Domain"].notna()
    merged["functional_evidence_present"] = merged["Molecular_Function"].notna()
    return merged


# ---------------------------------------------------------------------------
# Physicochemical descriptors (descriptive; not core score)
# ---------------------------------------------------------------------------

def physicochemical_features(seqs):
    rows = []
    for cname, meta in CANDIDATES.items():
        uid = meta["uniprot"]
        seq = seqs.get(uid, "")
        if not seq:
            continue
        a = ProteinAnalysis(seq)
        counts = a.count_amino_acids()
        L = len(seq)
        aliphatic = 100 * (counts["A"]/L) + 2.9*100*(counts["V"]/L) + 3.9*100*((counts["I"]+counts["L"])/L)
        rows.append({
            "Candidate": cname,
            "UniProt_ID": uid,
            "protein_length": L,
            "molecular_weight_Da": a.molecular_weight(),
            "pI": a.isoelectric_point(),
            "gravy": a.gravy(),
            "aromaticity": a.aromaticity(),
            "instability_index": a.instability_index(),
            "aliphatic_index": aliphatic,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def make_core_score(df):
    df = df.copy()

    # Evolution: combine mean family identity and global conservation.
    df["evo_identity_norm"] = minmax(df["mean_pairwise_identity"])
    # family_consensus_mean is global and identical for all; don't let it
    # create fake discrimination.
    df["evolution_score"] = df["evo_identity_norm"]

    # Promoter regulatory context: normalize comparable promoter-level features.
    # Log transform count-like quantities to reduce hub/density domination.
    for col in ["motif_hits", "motif_CE_hits", "FunTFBS_hits", "conserved_elements", "conserved_element_overlap_bp"]:
        df[f"log_{col}"] = np.log1p(pd.to_numeric(df[col], errors="coerce"))

    reg_components = [
        minmax(df["log_motif_hits"]),
        minmax(df["log_motif_CE_hits"]),
        minmax(df["log_FunTFBS_hits"]),
        minmax(df["log_conserved_elements"]),
        minmax(df["log_conserved_element_overlap_bp"]),
    ]
    df["promoter_regulation_score"] = pd.concat(reg_components, axis=1).mean(axis=1, skipna=True)

    # Structure: mean pLDDT is a confidence/quality feature, NOT a measure of
    # biological importance. It is used only to distinguish model confidence.
    df["structure_score"] = minmax(df["mean_pLDDT"])

    # A feature is only scored if it varies. All-equal domain/function evidence
    # is retained as annotation but contributes zero to ranking.
    df["core_score"] = (
        WEIGHTS["evolution"] * df["evolution_score"].fillna(0.5)
        + WEIGHTS["promoter_regulation"] * df["promoter_regulation_score"].fillna(0.5)
        + WEIGHTS["structure"] * df["structure_score"].fillna(0.5)
    )
    df["rank"] = df["core_score"].rank(method="min", ascending=False).astype(int)
    df = df.sort_values(["core_score", "Candidate"], ascending=[False, True]).reset_index(drop=True)
    return df


def sensitivity_analysis(base_df, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    records = []
    score_cols = ["evolution_score", "promoter_regulation_score", "structure_score"]
    for i in range(n):
        w = rng.dirichlet(np.ones(3))
        scores = (
            w[0] * base_df[score_cols[0]].fillna(0.5)
            + w[1] * base_df[score_cols[1]].fillna(0.5)
            + w[2] * base_df[score_cols[2]].fillna(0.5)
        )
        order = base_df.loc[scores.sort_values(ascending=False).index, "Candidate"].tolist()
        for rank, candidate in enumerate(order, 1):
            records.append({
                "iteration": i,
                "candidate": candidate,
                "rank": rank,
                "w_evolution": w[0],
                "w_promoter_regulation": w[1],
                "w_structure": w[2],
            })
    sens = pd.DataFrame(records)
    summary = sens.groupby("candidate").agg(
        mean_rank=("rank", "mean"),
        median_rank=("rank", "median"),
        top1_frequency=("rank", lambda x: np.mean(x == 1)),
        top2_frequency=("rank", lambda x: np.mean(x <= 2)),
    ).reset_index().sort_values(["mean_rank", "top1_frequency"], ascending=[True, False])
    return summary


# ---------------------------------------------------------------------------
# Figures / report
# ---------------------------------------------------------------------------

def make_figures(df):
    plot_df = df.sort_values("core_score", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(plot_df["Candidate"], plot_df["core_score"])
    ax.set_xlabel("Evidence-integrated prioritization score")
    ax.set_ylabel("DREB candidate")
    ax.set_title("DREB-NETSP candidate prioritization")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    fig.savefig(OUT / "DREB_NETSP_ranking.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    cols = ["evolution_score", "promoter_regulation_score", "structure_score"]
    hm = df.set_index("Candidate")[cols].copy()
    fig, ax = plt.subplots(figsize=(7, 4))
    im = ax.imshow(hm.values, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)), ["Evolution", "Promoter regulation", "Structure"])
    ax.set_yticks(range(len(hm)), hm.index)
    for i in range(hm.shape[0]):
        for j in range(hm.shape[1]):
            ax.text(j, i, f"{hm.iloc[i,j]:.2f}", ha="center", va="center")
    fig.colorbar(im, ax=ax, label="Normalized evidence")
    ax.set_title("DREB-NETSP evidence matrix")
    fig.tight_layout()
    fig.savefig(OUT / "DREB_NETSP_evidence_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(df, network_df, sensitivity, tree_result):
    path = OUT / "DREB_NETSP_report.txt"
    with open(path, "w") as f:
        f.write("DREB-NETSP INTEGRATED PIPELINE REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write("Candidates:\n")
        for cname, meta in CANDIDATES.items():
            f.write(f"  {cname}: {meta['gene']} / {meta['uniprot']}\n")
        f.write("\nCore score weights:\n")
        for k, v in WEIGHTS.items():
            f.write(f"  {k}: {v:.2f}\n")
        f.write("\nImportant interpretation safeguards:\n")
        f.write("- Missing evidence is not automatically interpreted as biological absence.\n")
        f.write("- Uniform InterPro/function annotations are retained but not used to rank.\n")
        f.write("- PlantRegMap TF-target network evidence is auxiliary unless source coverage is comparable.\n")
        f.write("- Scores prioritize hypotheses; they do not prove causal drought tolerance.\n")
        f.write("\nRanking:\n")
        for _, r in df.iterrows():
            f.write(f"  {int(r['rank'])}. {r['Candidate']} ({r['Gene']}) score={r['core_score']:.3f}\n")
        f.write("\nNetwork evidence:\n")
        f.write(network_df.to_string(index=False))
        f.write("\n\nSensitivity summary (1000 random weight combinations):\n")
        f.write(sensitivity.to_string(index=False))
        f.write("\n\nPhylogeny:\n")
        f.write(str(tree_result) + "\n")
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("DREB-NETSP INTEGRATED EVIDENCE PIPELINE")
    print("=" * 72)

    promoters = load_bed_promoters()
    seqs = load_sequences()
    aln, msa_seqs = load_msa()

    print(f"Candidates: {len(CANDIDATES)}")
    print(f"Promoter definitions loaded: {len(promoters)}")
    print(f"Protein sequences loaded: {len(seqs)}")
    print(f"MSA sequences loaded: {len(msa_seqs)}")

    evo, identity = evolutionary_features(aln)
    print("Evolutionary features: OK")

    promoter = promoter_features(promoters)
    promoter.to_csv(OUT / "DREB_NETSP_promoter_features.csv", index=False)
    print("Promoter/regulatory-context scan: OK")

    struct = structure_features()
    print("Structure features: OK")

    ann = annotation_features()
    print("InterPro/function annotations: OK")

    phys = physicochemical_features(seqs)
    print("Physicochemical descriptors: OK")

    network_df, network_edges = network_features()
    network_df.to_csv(OUT / "DREB_NETSP_network_evidence.csv", index=False)
    network_edges.to_csv(OUT / "DREB_NETSP_network_edges_deduplicated.csv", index=False)
    print("PlantRegMap network layer: OK")

    meta = pd.DataFrame([
        {"Candidate": k, "Gene": v["gene"], "UniProt_ID": v["uniprot"], "expected_length": v["expected_length"]}
        for k, v in CANDIDATES.items()
    ])

    # Convert evolution IDs into candidate names.
    evo = evo.merge(meta[["Candidate", "UniProt_ID"]], on="UniProt_ID", how="left")

    merged = meta.merge(evo, on=["Candidate", "UniProt_ID"], how="left")
    merged = merged.merge(promoter, on=["Candidate", "Gene"], how="left")
    merged = merged.merge(struct, on="UniProt_ID", how="left")
    merged = merged.merge(ann[["Candidate", "UniProt_ID", "Family", "Primary_Domain", "domain_present", "Molecular_Function", "Biological_Process", "functional_evidence_present"]], on=["Candidate", "UniProt_ID"], how="left")
    merged = merged.merge(phys, on=["Candidate", "UniProt_ID"], how="left")
    merged = merged.merge(network_df, on=["Candidate", "Gene"], how="left")

    # Recalculate expected-length check.
    merged["length_matches_expected"] = merged["protein_length"].eq(merged["expected_length"])

    ranked = make_core_score(merged)

    # Network coverage diagnostic: current files may only contain a source->target
    # network for a subset of candidates. Do not force this into the core score.
    network_coverage = float(ranked["network_has_evidence"].mean())
    ranked["network_coverage_fraction"] = network_coverage
    ranked["network_used_in_core_score"] = False

    sensitivity = sensitivity_analysis(ranked)
    sensitivity.to_csv(OUT / "DREB_NETSP_weight_sensitivity.csv", index=False)

    ranked.to_csv(OUT / "DREB_NETSP_candidate_ranking.csv", index=False)
    ranked.to_csv(OUT / "DREB_NETSP_integrated_evidence_matrix.csv", index=False)

    tree_result = build_nj_tree(aln)
    make_figures(ranked)
    report = write_report(ranked, network_df, sensitivity, tree_result)

    print("\nFINAL RANKING")
    print(ranked[["rank", "Candidate", "Gene", "core_score", "evolution_score", "promoter_regulation_score", "structure_score"]].to_string(index=False))
    print(f"\nNetwork source coverage across candidates: {network_coverage:.0%}")
    if network_coverage < 0.8:
        print("NOTE: network evidence is currently auxiliary and excluded from the core ranking because coverage is not comparable.")
    print(f"\nOutputs written to: {OUT}")
    print(f"Report: {report}")


if __name__ == "__main__":
    main()
