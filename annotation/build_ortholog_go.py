#!/usr/bin/env python3
"""
Transfer one species' GO annotation onto Deinandra genes via 1:1 orthologs
from OrthoFinder (one input to merge_go_sources.sh). A row of Orthogroups.tsv
counts as 1:1 if the Deinandra and <species> columns each hold exactly one
gene; the other species' columns are ignored.

Inputs:
  - OrthoFinder's Orthogroups/Orthogroups.tsv
  - the reference GFF, to resolve Deinandra transcript IDs to gene IDs
  - <species>_geneid_go.tsv: species gene_id -> GO term (from
    dump_biomart_go.sh or dump_tair_go.R). The species proteome must be keyed
    by gene id, as written by prep_tair_proteome.sh / prep_ensembl_plant_proteome.sh.

Usage: python3 build_ortholog_go.py <species_column> <species_go_file> <output.tsv>
  species_column: the species' column name in Orthogroups.tsv, e.g. Helianthus_annuus
  species_go_file: path to that species' gene_id/go_term TSV
  output.tsv: where to write the Deinandra gene_id/go_term result

Example:
  python3 build_ortholog_go.py Helianthus_annuus \\
      BASE_DIR/genome/annotation/Helianthus_annuus_geneid_go.tsv \\
      BASE_DIR/genome/annotated_sequences_orthofinder_hannuus_geneid_go.tsv
"""
import sys
import os
import re
import glob

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
ORTHOFINDER_RESULTS = f"{BASE_DIR}/genome/annotation/orthofinder_results"


def find_orthogroups_file():
    hits = glob.glob(os.path.join(ORTHOFINDER_RESULTS, "Results_*", "Orthogroups", "Orthogroups.tsv"))
    if not hits:
        sys.exit(f"Could not find Orthogroups.tsv under {ORTHOFINDER_RESULTS}")
    # use the most recent run if there are several
    hits.sort(key=os.path.getmtime, reverse=True)
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


def load_species_go(path):
    gene2go = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            gene_id, go_id = line.rstrip("\n").split("\t")
            gene2go.setdefault(gene_id, set()).add(go_id)
    return gene2go


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: python3 build_ortholog_go.py <species_column> <species_go_file> <output.tsv>")
    species_column, species_go_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

    orthogroups_file = find_orthogroups_file()
    print(f"Using OrthoFinder orthogroups file: {orthogroups_file}")
    print(f"Species column: {species_column}")

    tx2gene = load_tx2gene()
    print(f"Loaded {len(tx2gene)} transcript -> gene mappings from {GFF}")

    gene2go = load_species_go(species_go_file)
    print(f"Loaded GO terms for {len(gene2go)} {species_column} genes from {species_go_file}")

    n_orthogroups = 0
    n_single_copy = 0
    n_deinandra_tx_unmapped = 0
    n_species_no_go = 0
    pairs = set()
    deinandra_genes_with_ortholog = set()
    deinandra_genes_with_go = set()

    with open(orthogroups_file) as f:
        header = next(f).rstrip("\n").split("\t")
        deinandra_idx = header.index("Deinandra")
        species_idx = header.index(species_column)

        for line in f:
            n_orthogroups += 1
            cols = line.rstrip("\n").split("\t")
            deinandra_field = cols[deinandra_idx].strip()
            species_field = cols[species_idx].strip()
            if not deinandra_field or not species_field:
                continue
            # 1:1 for this species pair only; other columns ignored
            if "," in deinandra_field or "," in species_field:
                continue
            n_single_copy += 1

            deinandra_tx, species_gene = deinandra_field, species_field

            gene_id = tx2gene.get(deinandra_tx)
            if gene_id is None:
                n_deinandra_tx_unmapped += 1
                continue

            deinandra_genes_with_ortholog.add(gene_id)

            go_terms = gene2go.get(species_gene)
            if not go_terms:
                n_species_no_go += 1
                continue

            deinandra_genes_with_go.add(gene_id)
            for go_id in go_terms:
                pairs.add((gene_id, go_id))

    pairs = sorted(pairs)
    with open(out_file, "w") as out:
        out.write("gene_id\tgo_term\n")
        for gene_id, go_id in pairs:
            out.write(f"{gene_id}\t{go_id}\n")

    print(f"Orthogroups: {n_orthogroups}, single-copy (1:1) Deinandra<->{species_column} pairs: {n_single_copy} "
          f"({n_deinandra_tx_unmapped} Deinandra tx not found in GFF)")
    print(f"1:1-ortholog Deinandra genes: {len(deinandra_genes_with_ortholog)}, "
          f"of which {len(deinandra_genes_with_go)} had >=1 {species_column} GO term "
          f"({n_species_no_go} {species_column} orthologs had no GO annotation)")
    print(f"Wrote {len(pairs)} unique gene-GO pairs to {out_file}")


if __name__ == "__main__":
    main()
