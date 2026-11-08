#!/usr/bin/env python3
"""
Build a list of generic "housekeeping" GO term IDs to filter from GO
enrichment results. Hand-picked anchors, read against go_ancestor_map.tsv
(build_go_ancestor_map.py):
  SUBTREE_SEEDS   seed and all its descendants excluded (bounded machinery
                  subtrees: translation, proteasome, DNA replication/repair,
                  splicing, central carbon metabolism, protein folding)
  TERM_ONLY_SEEDS seed term only; descendants kept (ontology roots and
                  near-root MF/CC terms whose children are informative)

Usage: python3 build_housekeeping_go_exclusion_list.py [ancestor_map] [out_path]
  ancestor_map: path to go_ancestor_map.tsv (default: the genome dir copy)
  out_path: default BASE_DIR/genome/go_terms_excluded_housekeeping.txt
"""
import sys
from collections import defaultdict

BASE_DIR = "/private/groups/meyerlab/esteiger"

DEFAULT_ANCESTOR_MAP = f"{BASE_DIR}/genome/go_ancestor_map.tsv"
DEFAULT_OUT = f"{BASE_DIR}/genome/go_terms_excluded_housekeeping.txt"

# (GO ID, name, short rationale); seed and every descendant excluded
SUBTREE_SEEDS = [
    # Core translation / ribosome machinery
    ("GO:0006412", "translation", "core protein synthesis machinery"),
    ("GO:0042254", "ribosome biogenesis", "core ribosome assembly machinery"),
    ("GO:0006364", "rRNA processing", "core ribosome assembly machinery"),
    ("GO:0005840", "ribosome", "constitutive translation machinery (CC)"),
    ("GO:0003735", "structural constituent of ribosome", "constitutive translation machinery (MF)"),
    # Core protein degradation machinery
    ("GO:0043161", "proteasome-mediated ubiquitin-dependent protein catabolic process",
     "bulk protein turnover machinery"),
    ("GO:0000502", "proteasome complex", "bulk protein turnover machinery (CC)"),
    # Core DNA replication / repair
    ("GO:0006260", "DNA replication", "core genome maintenance machinery"),
    ("GO:0006281", "DNA repair", "core genome maintenance machinery"),
    # Core RNA splicing / processing
    ("GO:0008380", "RNA splicing", "core spliceosome machinery"),
    ("GO:0006397", "mRNA processing", "core mRNA maturation machinery"),
    # Core central carbon / energy metabolism
    ("GO:0006096", "glycolytic process", "core central carbon metabolism"),
    ("GO:0006099", "tricarboxylic acid cycle", "core central carbon metabolism"),
    ("GO:0006119", "oxidative phosphorylation", "core cellular energy production"),
    # Core protein folding
    ("GO:0006457", "protein folding", "constitutive chaperone machinery"),
]

# (GO ID, name, short rationale); seed term only, descendants kept
TERM_ONLY_SEEDS = [
    ("GO:0005515", "protein binding", "too generic to be informative alone"),
    ("GO:0003824", "catalytic activity", "too generic to be informative alone"),
    ("GO:0005488", "binding", "too generic to be informative alone"),
    ("GO:0016787", "hydrolase activity", "too generic to be informative alone"),
    ("GO:0016740", "transferase activity", "too generic to be informative alone"),
    ("GO:0003676", "nucleic acid binding", "too generic to be informative alone"),
    ("GO:0097159", "organic cyclic compound binding", "too generic to be informative alone"),
    # Ontology roots and near-root CC terms (near-universal gene coverage)
    ("GO:0008150", "biological_process", "ontology root — no differential meaning"),
    ("GO:0003674", "molecular_function", "ontology root — no differential meaning"),
    ("GO:0005575", "cellular_component", "ontology root — no differential meaning"),
    ("GO:0110165", "cellular anatomical entity", "near-root CC term, near-universal coverage"),
    ("GO:0043226", "organelle", "near-root CC term, near-universal coverage"),
    ("GO:0043227", "membrane-bounded organelle", "near-root CC term, near-universal coverage"),
    ("GO:0043229", "intracellular organelle", "near-root CC term, near-universal coverage"),
    ("GO:0043231", "intracellular membrane-bounded organelle", "near-root CC term, near-universal coverage"),
    ("GO:0032991", "protein-containing complex", "near-root CC term, near-universal coverage"),
]


def main():
    ancestor_map = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ANCESTOR_MAP
    out_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUT

    # descendants[ancestor] = set of terms that have `ancestor` in their ancestor list
    descendants = defaultdict(set)
    with open(ancestor_map) as f:
        next(f)  # header
        for line in f:
            term, ancestor = line.rstrip("\n").split("\t")
            descendants[ancestor].add(term)

    excluded = set()
    for sid, name, rationale in SUBTREE_SEEDS:
        n_desc = len(descendants.get(sid, set()))
        print(f"{sid} ({name}): seed + {n_desc} descendant terms — {rationale}")
        excluded.add(sid)
        excluded |= descendants.get(sid, set())
    for sid, name, rationale in TERM_ONLY_SEEDS:
        print(f"{sid} ({name}): seed only (descendants kept) — {rationale}")
        excluded.add(sid)

    with open(out_path, "w") as out:
        for go_id in sorted(excluded):
            out.write(f"{go_id}\thousekeeping\n")

    n_seeds = len(SUBTREE_SEEDS) + len(TERM_ONLY_SEEDS)
    print(f"\nWrote {len(excluded)} excluded GO term IDs (from {n_seeds} seed anchors) to {out_path}")


if __name__ == "__main__":
    main()
