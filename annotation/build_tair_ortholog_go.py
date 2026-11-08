#!/usr/bin/env python3
"""
Transfer Arabidopsis (TAIR) GO annotation onto Deinandra genes via 1:1
orthologs from OrthoFinder (run_orthofinder.sh); one input to
merge_go_sources.sh. 1:1 = single-copy orthogroups in Orthogroups.tsv. The
Orthologues/*.tsv pairwise files are not used because one Deinandra gene can
appear in several "1:1" rows within large paralog families.

Inputs:
  - OrthoFinder's Orthogroups/Orthogroups.tsv
  - the reference GFF, to resolve Deinandra transcript IDs to gene IDs
  - tair_geneid_go.tsv (dump_tair_go.R); AT proteome is keyed by AGI gene ID

Usage: python3 build_tair_ortholog_go.py [output.tsv]
"""
import sys
import os
import re
import glob

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
ORTHOFINDER_RESULTS = f"{BASE_DIR}/genome/annotation/orthofinder_results"
TAIR_GO = f"{BASE_DIR}/genome/annotation/tair_geneid_go.tsv"
DEFAULT_OUT = f"{BASE_DIR}/genome/annotated_sequences_orthofinder_tair_geneid_go.tsv"


def find_orthogroups_file():
    hits = glob.glob(os.path.join(ORTHOFINDER_RESULTS, "Results_*", "Orthogroups", "Orthogroups.tsv"))
    if not hits:
        sys.exit(f"Could not find Orthogroups.tsv under {ORTHOFINDER_RESULTS}")
    return hits[0]


def load_tx2gene():
    tx2gene = {}
    with open(GFF) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "mRNA":
                continue
            parent = re.search(r"Parent=([^;]+)", fields[8])
            mid = re.search(r"ID=([^;]+)", fields[8])
            if parent and mid:
                tx2gene[mid.group(1)] = parent.group(1)
    return tx2gene


def load_tair_go():
    gene2go = {}
    with open(TAIR_GO) as f:
        next(f)  # header
        for line in f:
            gene_id, go_id = line.rstrip("\n").split("\t")
            gene2go.setdefault(gene_id, set()).add(go_id)
    return gene2go


def main():
    out_file = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT

    orthogroups_file = find_orthogroups_file()
    print(f"Using OrthoFinder orthogroups file: {orthogroups_file}")

    tx2gene = load_tx2gene()
    print(f"Loaded {len(tx2gene)} transcript -> gene mappings from {GFF}")

    gene2go = load_tair_go()
    print(f"Loaded GO terms for {len(gene2go)} AT genes from {TAIR_GO}")

    n_orthogroups = 0
    n_single_copy = 0
    n_deinandra_tx_unmapped = 0
    n_at_no_go = 0
    pairs = set()
    deinandra_genes_with_ortholog = set()
    deinandra_genes_with_go = set()

    with open(orthogroups_file) as f:
        header = next(f).rstrip("\n").split("\t")
        deinandra_idx = header.index("Deinandra")
        at_idx = header.index("Arabidopsis_thaliana")

        for line in f:
            n_orthogroups += 1
            cols = line.rstrip("\n").split("\t")
            deinandra_field = cols[deinandra_idx].strip()
            at_field = cols[at_idx].strip()
            if not deinandra_field or not at_field:
                continue
            # Single-copy orthogroup: exactly one gene per species (no comma).
            if "," in deinandra_field or "," in at_field:
                continue
            n_single_copy += 1

            deinandra_tx, at_gene = deinandra_field, at_field

            gene_id = tx2gene.get(deinandra_tx)
            if gene_id is None:
                n_deinandra_tx_unmapped += 1
                continue

            deinandra_genes_with_ortholog.add(gene_id)

            go_terms = gene2go.get(at_gene)
            if not go_terms:
                n_at_no_go += 1
                continue

            deinandra_genes_with_go.add(gene_id)
            for go_id in go_terms:
                pairs.add((gene_id, go_id))

    pairs = sorted(pairs)
    with open(out_file, "w") as out:
        out.write("gene_id\tgo_term\n")
        for gene_id, go_id in pairs:
            out.write(f"{gene_id}\t{go_id}\n")

    print(f"Orthogroups: {n_orthogroups}, single-copy (1:1) Deinandra<->AT orthogroups: {n_single_copy} "
          f"({n_deinandra_tx_unmapped} Deinandra tx not found in GFF)")
    print(f"1:1-ortholog Deinandra genes: {len(deinandra_genes_with_ortholog)}, "
          f"of which {len(deinandra_genes_with_go)} had >=1 AT GO term "
          f"({n_at_no_go} AT orthologs had no GO annotation)")
    print(f"Wrote {len(pairs)} unique gene-GO pairs to {out_file}")


if __name__ == "__main__":
    main()
