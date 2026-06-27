# BioAlign Extensions: Clustering Coefficient & Degree Matching

Extensions to [BioAlign](https://github.com/original/BioAlign) that modify the Stage-2 topology step. Everything else — Stage-1 biological similarity, remote homology, secondary structure — is identical to the original.

---

## What Was Changed

BioAlign's original topology step scores candidate pairs by counting mutually-aligned neighbours. Two alternative strategies are provided here:

| File | Method | What it does |
|------|--------|--------------|
| `ClusterCoefficientBioAlign.py` | Neighbourhood Expansion | Aligns residual nodes by expanding outward from already-aligned clusters using local clustering patterns |
| `degree_matching_technique.py` | Degree Centrality Matching | Aligns residual nodes by matching proteins with the closest degree (number of interactions) |

---

## Setup

```bash
pip install numpy
```

Place your data files in the same directory structure as original BioAlign:

```
./3D-structure-similarity/   ← .tmscore files
./global-sequence-similarity/ ← .blast files
./local-sequence-similarity/  ← .local files
./networks/                   ← .interaction files
./homologs/                   ← per-species homolog .txt files
./sec/                        ← per-species secondary structure .txt files
./new_alignments/             ← output alignments written here (create this folder)
```

---

## Run

Both scripts take the same three arguments: `specie1 specie2 choice`

`choice` controls which stages run:
- `b` — Stage-1 + Remote Homology + Secondary Structure (no topology)
- `t` — Stage-1 + Topology only
- `bt` — Stage-1 + Remote Homology + Secondary Structure + Topology

```bash
# Clustering Coefficient
python ClusterCoefficientBioAlign.py mouse human bt
python ClusterCoefficientBioAlign.py yeast human t

# Degree Matching
python degree_matching_technique.py mouse human bt
python degree_matching_technique.py yeast human t
```

Output alignment is saved to `./new_alignments/{specie1}-{specie2}-{choice}.alignment`

---

## How Each Method Works

### Clustering Coefficient (CC)
Extends the alignment using neighbourhood expansion. For each unaligned node pair (u, v), it checks whether their already-aligned neighbours are mutually aligned — same as the original — but applies this as a looser threshold test rather than a strict ranking. It additionally computes the local clustering coefficient of every node (fraction of its neighbours that are connected to each other) to characterise local network density. The result is that CC expands the alignment further into peripheral regions of the network that the original topology step would skip.

**Effect:** Coverage increases to ~95–100% on t and bt variants. AFS drops slightly (~0.003–0.006) because more peripheral, less-annotated proteins are included.

### Degree Matching (DM)
Replaces topology entirely with degree-based greedy matching. For each unaligned protein in network 1 (sorted by degree, highest first), it finds the unaligned protein in network 2 whose degree is closest. Hub proteins are matched first. No neighbourhood topology is required.

**Effect:** AFS improves by ~0.006–0.013 on MF because hub proteins carry richer GO annotations. Coverage drops slightly (~1.5–2.5%) because some peripheral pairs that topology-propagation would catch are missed.

---

## Results Summary (bt variant, averaged over 5 species pairs)

| Method | MF AFS | BP AFS | MF Coverage | BP Coverage |
|--------|--------|--------|-------------|-------------|
| BioAlign (original) | 0.2516 | 0.1479 | ~75.8% | ~76.5% |
| Clustering Coefficient | 0.2477 | 0.1470 | ~99.2% | ~99.2% |
| Degree Matching | **0.2579** | **0.1508** | ~74.2% | ~74.2% |

Use **CC** when coverage matters. Use **DM** when alignment quality (AFS) matters.

---

## Author

Asiya Jahangir
