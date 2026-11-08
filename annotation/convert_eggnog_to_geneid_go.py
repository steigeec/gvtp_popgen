#!/usr/bin/env python3
"""
Reshape eggNOG-mapper .emapper.annotations into two-column gene_id / go_term
format. emapper's query column is the transcript id, so transcripts are
mapped to genes from the reference GFF's mRNA lines (ID=, Parent=).

Also writes eggnog_annotation_coverage.jpg next to the output TSV: genes with
no eggNOG hit, a hit but no GO term, or a GO term.

Usage: python3 convert_eggnog_to_geneid_go.py <prefix>.emapper.annotations [output.tsv]
"""
import sys
import os
import re

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: convert_eggnog_to_geneid_go.py <prefix>.emapper.annotations [output.tsv]")
    ann_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) > 2 else \
        f"{BASE_DIR}/genome/annotated_sequences_eggnog_geneid_go.tsv"

    # --- transcript_id -> gene_id, from the GFF's own mRNA feature lines ---
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
    print(f"Loaded {len(tx2gene)} transcript -> gene mappings from {GFF}")

    # --- parse emapper.annotations: find header row, then query/GOs columns ---
    header_cols = None
    query_idx = go_idx = None
    n_rows = 0
    n_unmapped_query = 0
    pairs = []  # (gene_id, go_id)
    genes_with_hit = set()   # any emapper.annotations row at all (ortholog found)
    genes_with_go = set()    # that row's GOs field was non-empty

    with open(ann_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                # header row, e.g. "#query  seed_ortholog  evalue ... GOs ..."
                header_cols = line.lstrip("#").split("\t")
                query_idx = header_cols.index("query")
                go_idx = header_cols.index("GOs")
                continue
            if header_cols is None:
                continue  # shouldn't happen, but don't crash on unexpected format
            cols = line.split("\t")
            if len(cols) <= max(query_idx, go_idx):
                continue
            n_rows += 1
            query = cols[query_idx]
            go_field = cols[go_idx]

            gene_id = tx2gene.get(query)
            if gene_id is None:
                # not in GFF: fall back to the query id with suffix stripped
                gene_id = query.split(".")[0]
                n_unmapped_query += 1

            genes_with_hit.add(gene_id)
            if go_field not in ("-", ""):
                genes_with_go.add(gene_id)
                for go_id in go_field.split(","):
                    go_id = go_id.strip()
                    if go_id:
                        pairs.append((gene_id, go_id))

    pairs = sorted(set(pairs))
    with open(out_file, "w") as out:
        out.write("gene_id\tgo_term\n")
        for gene_id, go_id in pairs:
            out.write(f"{gene_id}\t{go_id}\n")

    print(f"Parsed {n_rows} queries with an eggNOG hit, {len(genes_with_go)} genes with at least one GO term "
          f"({n_unmapped_query} queries not found in GFF, used dot-strip fallback)")
    print(f"Wrote {len(pairs)} unique gene-GO pairs across {len(genes_with_go)} genes to {out_file}")

    # --- annotation coverage plot ---
    all_genes = set(tx2gene.values())
    n_total = len(all_genes)
    n_has_go = len(genes_with_go & all_genes)
    n_hit_no_go = len((genes_with_hit & all_genes) - genes_with_go)
    n_no_hit = n_total - len(genes_with_hit & all_genes)

    plot_coverage(n_total, n_no_hit, n_hit_no_go, n_has_go,
                  os.path.join(os.path.dirname(os.path.abspath(out_file)) or ".",
                               "eggnog_annotation_coverage.jpg"))


def plot_coverage(n_total, n_no_hit, n_hit_no_go, n_has_go, plot_file):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    categories = ["No eggNOG hit", "Hit, no GO term", "GO term assigned"]
    counts = [n_no_hit, n_hit_no_go, n_has_go]
    pcts = [100 * c / n_total for c in counts]
    colors = ["#b3b3b3", "#f0a04b", "#4c9a5b"]

    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(categories, counts, color=colors)
    for bar, count, pct in zip(bars, counts, pcts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{count}\n({pct:.1f}%)", ha="center", va="bottom")
    ax.set_ylabel("Number of genes")
    ax.set_title(f"eggNOG-mapper annotation coverage (n={n_total} genes)")
    ax.set_ylim(0, max(counts) * 1.2)
    fig.tight_layout()
    fig.savefig(plot_file, format="jpg", dpi=150)
    plt.close(fig)
    print(f"Wrote coverage plot to {plot_file}")
    print(f"  No hit: {n_no_hit} ({100*n_no_hit/n_total:.1f}%)  |  "
          f"Hit, no GO: {n_hit_no_go} ({100*n_hit_no_go/n_total:.1f}%)  |  "
          f"GO assigned: {n_has_go} ({100*n_has_go/n_total:.1f}%)")


if __name__ == "__main__":
    main()
