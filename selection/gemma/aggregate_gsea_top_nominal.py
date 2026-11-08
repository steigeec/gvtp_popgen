#!/usr/bin/env python3
"""
Combines cluster_gsea_go_results.R output (<var>/GSEA_GO_top_nominal_<var>.tsv)
across environmental variables into one table, adding a Variable column and
the GO term definition (from go_term_definitions.tsv).

Usage: python3 aggregate_gsea_top_nominal.py <var1> [<var2> ...]
Writes BASE_DIR/selection/gemma/GSEA_GO_top_nominal_all_variables.tsv
"""
import csv
import sys

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA_DIR = f"{BASE_DIR}/selection/gemma"
DEFINITIONS_FILE = f"{BASE_DIR}/genome/go_term_definitions.tsv"
OUT_FILE = f"{GEMMA_DIR}/GSEA_GO_top_nominal_all_variables.tsv"


def load_definitions():
    defs = {}
    with open(DEFINITIONS_FILE) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            defs[row["id"]] = row["definition"]
    return defs


def main():
    variables = sys.argv[1:]
    if not variables:
        sys.exit("Usage: aggregate_gsea_top_nominal.py <var1> [<var2> ...]")

    definitions = load_definitions()
    print(f"Loaded {len(definitions)} GO term definitions", file=sys.stderr)

    all_rows = []
    fieldnames = None
    for var in variables:
        path = f"{GEMMA_DIR}/{var}/GSEA_GO_top_nominal_{var}.tsv"
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if fieldnames is None:
                fieldnames = reader.fieldnames
            for row in reader:
                row["Variable"] = var
                row["Definition"] = definitions.get(row["pathway"], "NA")
                all_rows.append(row)
        print(f"{var}: {sum(1 for r in all_rows if r['Variable'] == var)} rows", file=sys.stderr)

    all_rows.sort(key=lambda r: float(r["pval"]))

    out_fieldnames = ["Variable"] + fieldnames + ["Definition"]
    with open(OUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Wrote {len(all_rows)} rows across {len(variables)} variables to {OUT_FILE}", file=sys.stderr)


if __name__ == "__main__":
    main()
