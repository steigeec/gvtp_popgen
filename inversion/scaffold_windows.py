#!/usr/bin/env python3
"""
Sliding-window gene density and repeat content for any list of scaffolds, with
the same inputs and conventions as sections 2-3 of inversion_analysis.py
(features counted only if entirely inside a window; RepeatMasker hits flagged
"*" excluded; genome-wide reference = all other scaffolds).

Writes to OUT_DIR, per scaffold:
  gene_density_windows_<scaf>.tsv   start, end, n_genes, density (genes/Mb), mid
  repeat_windows_<scaf>.tsv         start, end, mid, <repeat class %>..., total
  genome_wide_refs_<scaf>.tsv       gene_density, total_repeat, repeat_<class> (rest of genome)

Usage: scaffold_windows.py DEINVI_00003 DEINVI_00008 ...
"""

import sys

import pandas as pd

BASE_DIR = "/private/groups/meyerlab/esteiger"

GFF      = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"
FASTA    = f"{BASE_DIR}/genome/DEINVI_genome.fasta.masked"
RM_OUT   = f"{BASE_DIR}/inversions/gene_density/DEINVI_genome.fasta.LTRfix.out"
OUT_DIR  = f"{BASE_DIR}/inversions/gene_density"
WIN_SIZE = 500_000

CLASS_ORDER = ["LTR/Gypsy", "LTR/Copia", "LTR/Other", "LINE", "SINE",
               "DNA", "Satellite", "Simple_repeat", "Low_complexity", "Unknown", "Other"]


def broad_class(cls):
    """Same lumping as inversion_analysis.py."""
    c = cls.split("/")[0]
    if c == "LTR":
        sub = cls.split("/")[1] if "/" in cls else "Unknown"
        return "LTR/Gypsy" if "Gypsy" in sub else "LTR/Copia" if "Copia" in sub else "LTR/Other"
    if c in ("Simple_repeat", "Low_complexity", "Satellite", "LINE", "SINE", "DNA", "Unknown"):
        return c
    if c in ("RC", "PLE"):
        return "DNA"
    return "Other"


def main(scaffolds):
    print("Parsing genome sizes...", flush=True)
    sizes, cur, n = {}, None, 0
    with open(FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    sizes[cur] = n
                cur, n = line[1:].split()[0], 0
            else:
                n += len(line.rstrip())
    sizes[cur] = n

    print("Parsing genes...", flush=True)
    genes = []
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) >= 9 and p[2] == "gene":
                genes.append((p[0], int(p[3]), int(p[4])))
    genes = pd.DataFrame(genes, columns=["chrom", "start", "end"])

    print("Parsing RepeatMasker .out...", flush=True)
    rec = []
    with open(RM_OUT) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip()
            if i < 3 or not line or line.endswith("*"):
                continue
            p = line.split()
            if len(p) >= 11:
                rec.append((p[4], int(p[5]), int(p[6]), p[10]))
    rm = pd.DataFrame(rec, columns=["chrom", "start", "end", "class"])
    rm["length"] = rm["end"] - rm["start"] + 1
    rm["broad"] = rm["class"].map(broad_class)

    for scaf in scaffolds:
        L = sizes[scaf]
        g = genes[genes.chrom == scaf]
        r = rm[rm.chrom == scaf]
        gd_rows, rep_rows = [], []
        for ws in range(0, L, WIN_SIZE):
            we = min(ws + WIN_SIZE, L)
            n_g = ((g.start >= ws) & (g.end <= we)).sum()
            gd_rows.append({"start": ws, "end": we, "n_genes": n_g,
                            "density": n_g / ((we - ws) / 1e6), "mid": (ws + we) / 2})
            sub = r[(r.start >= ws) & (r.end <= we)]
            row = {"start": ws, "end": we, "mid": (ws + we) / 2}
            for cls in CLASS_ORDER:
                row[cls] = sub.loc[sub.broad == cls, "length"].sum() / (we - ws) * 100
            row["total"] = sub["length"].sum() / (we - ws) * 100
            rep_rows.append(row)
        pd.DataFrame(gd_rows).to_csv(f"{OUT_DIR}/gene_density_windows_{scaf}.tsv", sep="\t", index=False)
        pd.DataFrame(rep_rows).to_csv(f"{OUT_DIR}/repeat_windows_{scaf}.tsv", sep="\t", index=False)

        other_bp = sum(v for k, v in sizes.items() if k != scaf)
        other_rm = rm[rm.chrom != scaf]
        frac = other_rm.groupby("broad")["length"].sum().reindex(CLASS_ORDER, fill_value=0) / other_bp * 100
        refs = {"gene_density": (genes.chrom != scaf).sum() / (other_bp / 1e6),
                "total_repeat": frac.sum()}
        refs.update({f"repeat_{c}": frac[c] for c in CLASS_ORDER})
        pd.DataFrame({"measure": list(refs), "genome_wide": list(refs.values())}).to_csv(
            f"{OUT_DIR}/genome_wide_refs_{scaf}.tsv", sep="\t", index=False)
        print(f"  {scaf}: {len(gd_rows)} windows written", flush=True)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
