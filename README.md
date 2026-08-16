# DREB-NETSP Core 🧬🌱

## An Evolutionary and Computational Framework for Prioritizing DREB Candidates in Rice

> **DREB-NETSP** is a modular computational framework developed to integrate
> evolutionary, sequence, protein, functional, and regulatory evidence for
> systematic prioritization of DREB candidates associated with drought-response
> biology in rice.

---

## 🌱 Overview

Drought is one of the major environmental stresses affecting crop productivity,
and the **Dehydration-Responsive Element-Binding (DREB)** family represents an
important group of transcription factors involved in plant stress responses.

DREB-NETSP was developed as a computational approach to investigate DREB
candidates in rice by bringing together multiple evidence layers rather than
relying on a single sequence or annotation-based criterion.

The framework currently integrates:

- sequence quality assessment
- multiple sequence alignment
- pairwise sequence identity
- evolutionary and phylogenetic analysis
- conserved-position and conserved-region analysis
- protein physicochemical characterization
- domain characterization
- functional annotation
- regulatory evidence
- integrated evidence scoring
- network-based evidence representation
- candidate prioritization
- weight-sensitivity analysis

The goal is not simply to identify DREB proteins, but to establish a
**structured evidence-integration framework** through which candidates can be
ranked for further biological investigation.

---

## 🎯 Research Question

**Can heterogeneous computational evidence—including sequence conservation,
evolutionary relationships, conserved regions, protein characteristics,
functional annotation, and regulatory evidence—be systematically integrated
to prioritize DREB candidates for further investigation of drought-response
biology in rice?**

---

## 🧬 Why DREB?

DREB transcription factors are important components of plant stress-response
networks.

Their evolutionary conservation, sequence characteristics, conserved
regions, domain architecture, functional annotation, and regulatory context
provide complementary information that can potentially distinguish candidates
with stronger evidence of biological relevance.

DREB-NETSP therefore approaches candidate prioritization as an
**evidence-integration problem**, rather than depending on one individual
feature.

---

# 🔬 DREB-NETSP Workflow

The current framework follows a modular computational workflow:

```text
Literature & Study Rationale
             │
             ▼
      DREB Candidate Curation
             │
             ▼
    Sequence Quality Assessment
             │
             ▼
     Multiple Sequence Alignment
             │
             ▼
 Sequence Identity & Evolutionary Analysis
             │
             ▼
  Conserved Positions / Regions / Motifs
             │
             ▼
    Protein Characterization
             │
             ▼
   Domain & Functional Annotation
             │
             ▼
      Regulatory Evidence
             │
             ▼
     Evidence Integration
             │
             ▼
       Network Evidence
             │
             ▼
    Candidate Prioritization
             │
             ▼
    Weight Sensitivity Analysis
             │
             ▼
     Prioritized DREB Candidates
```

The framework is designed to remain modular so that additional evidence layers
can be incorporated as the project develops.

---

# 📌 Current Objectives

The current implementation of DREB-NETSP focuses on:

1. Curating a rice DREB candidate dataset.
2. Performing sequence quality assessment.
3. Performing multiple sequence alignment.
4. Characterizing sequence identity and evolutionary relationships.
5. Identifying conserved positions and regions.
6. Characterizing DREB protein physicochemical properties.
7. Investigating domain architecture.
8. Integrating functional annotation evidence.
9. Incorporating regulatory evidence.
10. Constructing an integrated evidence matrix.
11. Developing a candidate prioritization approach.
12. Evaluating ranking sensitivity to evidence weighting.
13. Generating reproducible computational outputs and visualizations.

---

# 🏆 Candidate Prioritization

One of the central components of DREB-NETSP is the integration of multiple
evidence layers into a candidate prioritization framework.

Rather than treating one computational result as definitive, the framework
combines complementary evidence to generate an **evidence-integrated
prioritization score**.

The current candidate-prioritization module contains:

- integrated evidence matrix
- evidence heatmap
- candidate ranking
- network evidence
- deduplicated network edges
- promoter-related features
- weight-sensitivity analysis
- Neighbor-Joining tree representation
- candidate ranking report

All current candidate-prioritization outputs are organized under:

```text
03_RESULTS/
└── 06_CANDIDATE_PRIORITISATION/
```

### Current candidate-prioritization outputs

| File | Description |
|---|---|
| `DREB_NETSP_NJ_tree.newick` | Neighbor-Joining tree data |
| `DREB_NETSP_NJ_tree.png` | Visualization of the NJ tree |
| `DREB_NETSP_candidate_ranking.csv` | Candidate ranking table |
| `DREB_NETSP_evidence_heatmap.png` | Integrated evidence heatmap |
| `DREB_NETSP_integrated_evidence_matrix.csv` | Integrated evidence matrix |
| `DREB_NETSP_network_edges_deduplicated.csv` | Deduplicated network edge table |
| `DREB_NETSP_network_evidence.csv` | Network-associated evidence |
| `DREB_NETSP_promoter_features.csv` | Promoter-related features |
| `DREB_NETSP_ranking.png` | Candidate prioritization visualization |
| `DREB_NETSP_report.txt` | Candidate prioritization report |
| `DREB_NETSP_weight_sensitivity.csv` | Sensitivity analysis of evidence weighting |

---

# 📊 Evidence Integration

The current implementation explores multiple dimensions of evidence, including:

### Evolutionary Evidence

- sequence similarity
- pairwise identity
- multiple sequence alignment
- phylogenetic relationships
- conserved positions
- conserved regions

### Protein Evidence

- physicochemical properties
- domain architecture
- sequence characteristics

### Functional Evidence

- functional annotations
- biological process information
- molecular function information

### Regulatory Evidence

- regulatory mapping
- promoter-associated features
- plant regulatory resources

### Integrated Evidence

These evidence layers are brought together into a common computational
framework for candidate prioritization.

The current evidence matrix provides a transparent representation of how
different evidence dimensions contribute to the prioritization process.

---

# 📈 Robustness and Weight Sensitivity

Candidate ranking can depend on how individual evidence layers are weighted.

DREB-NETSP therefore includes a **weight-sensitivity analysis** to examine
whether candidate prioritization remains relatively stable when evidence
weights are varied.

This is important because a prioritization framework should not depend
entirely on one arbitrary weighting configuration.

The current sensitivity-analysis output is available at:

```text
03_RESULTS/06_CANDIDATE_PRIORITISATION/
└── DREB_NETSP_weight_sensitivity.csv
```

This component will be further developed as additional datasets and evidence
layers are incorporated.

---

# 📁 Repository Structure

```text
DREB--NETSP-CORE/
│
├── 01_LITERATURE REVIEW/
│   └── Research gap, study rationale and literature resources
│
├── 02_INDEX/
│   └── 02-DATA/
│       └── 02_DATA/
│           ├── curated DREB datasets
│           ├── candidate registry
│           ├── protein sequences
│           ├── domain annotations
│           └── data index
│
├── 03_RESULTS/
│   │
│   ├── 01_SEQUENCE_EVOLUTION/
│   │   ├── MSA
│   │   ├── pairwise identity
│   │   ├── phylogenetic analysis
│   │   ├── conserved positions
│   │   ├── conserved regions
│   │   └── evolutionary summaries
│   │
│   ├── 02_FUNCTIONAL_ANNOTATION/
│   │   └── functional annotation evidence
│   │
│   ├── 03_PROTEIN_CHARACTERIZATION/
│   │   └── protein physicochemical characterization
│   │
│   ├── 04_REGULATORY_EVIDENCE/
│   │   └── regulatory evidence and PlantRegMap mapping
│   │
│   ├── 05_QUALITY_CONTROL/
│   │   └── sequence quality-control outputs
│   │
│   └── 06_CANDIDATE_PRIORITISATION/
│       ├── integrated evidence matrix
│       ├── candidate ranking
│       ├── evidence heatmap
│       ├── network evidence
│       ├── promoter features
│       ├── NJ tree
│       ├── sensitivity analysis
│       └── prioritization report
│
├── 04_SCRIPTS/
│   ├── MSA_analysis.py
│   ├── Phylogenetic_Analysis (1).py
│   ├── Physicochemical_Analysis_DREB.py
│   ├── conserved_motif_analysis (1).py
│   └── sequence_qc.py
│
├── 05_FIGURES/
│   ├── MSA figures
│   ├── phylogenetic figures
│   ├── conserved-region figures
│   ├── functional evidence figures
│   ├── domain architecture figures
│   └── structural analysis figures
│
├── 06_POSTER/
│
├── 07_DOCUMENTATION/
│
└── README.md
```

---

# 🛠️ Tools & Resources

## Databases and Resources

The current project makes use of publicly available biological resources,
including:

- PlantTFDB
- NCBI
- Ensembl Plants
- InterPro
- Pfam
- PlantRegMap

## Computational Tools

The framework is primarily implemented using:

- Python
- Biopython
- MAFFT
- sequence-analysis utilities
- publicly available plant regulatory resources

## Version Control

- Git
- GitHub

---

# 🧪 Computational Implementation

The framework is implemented as a collection of modular Python scripts.

The scripts are separated from the generated results to make the project
easier to inspect, reproduce, modify, and extend.

Current analysis scripts include:

```text
04_SCRIPTS/
├── MSA_analysis.py
├── Phylogenetic_Analysis (1).py
├── Physicochemical_Analysis_DREB.py
├── conserved_motif_analysis (1).py
└── sequence_qc.py
```

The modular structure allows additional analyses to be incorporated without
rebuilding the complete framework from the beginning.

---

# ♻️ Reproducibility

DREB-NETSP is organized around a separation between:

```text
Input / Indexed Data
        ↓
Computational Scripts
        ↓
Intermediate Results
        ↓
Integrated Evidence
        ↓
Figures & Visualizations
        ↓
Candidate Prioritization
```

The repository therefore preserves both the computational workflow and the
corresponding outputs generated during the current development stage.

As the framework develops, additional reproducibility components such as
environment specifications, dependency management, and automated workflows
may be incorporated.

---

# 🚧 Current Project Status

## Conference-Stage Computational Framework

The current repository represents an **initial implementation and
proof-of-concept of the DREB-NETSP computational framework**.

The present version demonstrates the integration of multiple computational
evidence layers and provides an initial candidate-prioritization workflow
using a curated dataset.

### Important scientific note

The current implementation is **not presented as a final biological
validation study or publication-complete analysis**.

The present goal is to establish and demonstrate the computational framework,
its modular architecture, evidence integration strategy, and candidate
prioritization concept.

Following the conference stage, the framework is intended to be expanded using
broader datasets, additional evidence layers, and independent validation.

---

# ⚠️ Current Limitations

The current implementation has several limitations that will be addressed in
future development:

- The initial dataset is limited in scope.
- Computational evidence is dependent on the quality and availability of
  public annotations and databases.
- Candidate prioritization is computational and does not constitute
  experimental validation.
- Evidence integration involves weighting assumptions that require further
  evaluation.
- Additional datasets are required to test the generalizability of the
  prioritization framework.
- Independent biological evidence will be required to strengthen candidate
  prioritization.

These limitations are considered part of the current development stage rather
than hidden from the interpretation of the framework.

---

# 🚀 Future Development

The next stage of DREB-NETSP will focus on expanding both the dataset and the
evidence architecture.

Planned development includes:

### Dataset Expansion

- expansion of the DREB candidate dataset
- inclusion of broader genomic resources
- additional rice lineages/accessions where appropriate

### Evidence Expansion

- additional regulatory evidence
- expanded functional evidence
- transcriptomic/expression evidence where appropriate
- additional sequence and evolutionary features

### Framework Development

- improved evidence normalization
- evaluation of alternative weighting strategies
- expanded robustness analysis
- independent candidate validation
- improved network-based integration

### Biological Interpretation

The ultimate objective is to determine whether candidates prioritized by the
computational framework are consistently supported by independent biological
evidence and therefore represent stronger candidates for future experimental
investigation.

---

# 📚 Research Philosophy

DREB-NETSP is being developed with a simple principle:

> **A biological candidate should not be prioritized because of one feature
> alone; stronger prioritization should emerge when independent evidence
> converges.**

The framework is therefore designed around evidence integration rather than
single-metric prediction.

At the current stage, the framework is a computational foundation. Its value
will ultimately depend on how well it performs when challenged with broader
datasets and independent biological evidence.

---

# 🌱 Why This Project Exists

DREB-NETSP started as an attempt to understand how bioinformatics could move
beyond individual analyses and become a connected computational workflow.

Instead of performing sequence alignment, phylogenetics, protein
characterization, annotation, and regulatory analysis as isolated exercises,
this project attempts to connect these analyses into a single framework.

The current implementation is the first version of that idea.

It is expected to evolve.

---

# 👨‍💻 Author

**Surya Prasad Karmakar**

B.Tech Biotechnology  
Lovely Professional University, India

Interested in:

- Computational Biology
- Bioinformatics
- Evolutionary Biology
- Plant Molecular Biology
- Genomics
- AI/ML applications in Biology
- Research and scientific discovery

---

# 🧬 DREB-NETSP

**DREB-NETSP Core**  
*An evolutionary and computational framework for prioritizing DREB
candidates in rice.*

Built as an ongoing computational research project.

> *This repository documents the development of the framework — including its
> current capabilities, limitations, and future direction.*

---

**Surya Prasad Karmakar**  
*Building, learning, questioning, and iterating — one biological problem at a time.* 🌱🧬💻
