#!/usr/bin/env python3
"""
Merges gene_summary_<var>.tsv across variables into the master
gene_summary_all.tsv, tagging each row with its source Variable.

Column sources:
  - Gene_Name, Gene_Coordinates, Orthologs, BLAST_Gene_Names, Acronym,
    Functional_Description, Mutation_Types, Impact_Levels: from
    gene_summary_<var>.tsv, joined on Gene_ID (Gene_Name is a display label only).
  - SNP_Count, Min_P_Value: from gene_snp_stats_<var>.tsv
    (compute_gene_snp_stats.py), joined on Gene_ID.
  - Significant_GO_Terms: "GO:ID(term); GO:ID(term)" for any FDR<0.05
    category that lists this gene, from
    GO_enrichment_results_<var>_unfiltered_<GROUP>_goseq.tsv if present.
  - Putative_Adaptive_Role, Citation: reused from the current
    gene_summary_all.tsv or its pre-*.tsv backups for the same
    (Variable, Gene_Coordinates) pair if the Acronym is unchanged;
    otherwise "NA", to be filled from the literature.

After writing, two consistency checks flag rows whose role text names a gene
identity not supported by the row's current Orthologs/Acronym/Functional_Description.

Usage: python3 build_gene_summary_all.py <var1> [<var2> ...]
Writes gene_summary_all.tsv (backs up the existing one first).
"""
import csv
import glob
import os
import re
import shutil
import sys
from datetime import datetime

BASE_DIR = "/private/groups/meyerlab/esteiger"

GEMMA_DIR = f"{BASE_DIR}/selection/gemma"
MASTER = os.path.join(GEMMA_DIR, "gene_summary_all.tsv")
OUT_COLUMNS = [
    "Variable", "Gene_Name", "Gene_Coordinates", "Orthologs", "BLAST_Gene_Names",
    "Acronym", "Functional_Description", "Significant_GO_Terms", "SNP_Count",
    "Min_P_Value", "Mutation_Types", "Impact_Levels", "Putative_Adaptive_Role",
    "Citation",
]


def load_gene_summary(var):
    path = os.path.join(GEMMA_DIR, var, f"gene_summary_{var}.tsv")
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows[row["Gene_ID"]] = row
    return rows


def load_snp_stats(var):
    path = os.path.join(GEMMA_DIR, var, f"gene_snp_stats_{var}.tsv")
    stats = {}
    if not os.path.exists(path):
        return stats
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            stats[row["gene_id"]] = (row["snp_count"], row["min_p"])
    return stats


def load_significant_go(var):
    """gene_id -> ["GO:ID(term)", ...] from FDR<0.05 goseq unfiltered results."""
    sig = {}
    pattern = os.path.join(GEMMA_DIR, var, f"GO_enrichment_results_{var}_unfiltered_*_goseq.tsv")
    for path in glob.glob(pattern):
        if path.endswith("_reduced.tsv"):
            continue
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                label = f"GO:{row['category'].split(':')[-1]}({row['term']})" if not row["category"].startswith("GO:") else f"{row['category']}({row['term']})"
                for gid in row.get("Genes", "").split(","):
                    gid = gid.strip()
                    if gid:
                        sig.setdefault(gid, []).append(label)
    return sig


def load_prior_adaptive_role():
    """(Variable, Gene_Coordinates) -> (role, citation, acronym_then), first
    non-empty hit across the current master sheet and every pre-*.tsv
    backup. acronym_then is the Acronym stored with that row; the caller checks
    it against the current Acronym before reusing role/citation."""
    files = [MASTER] + sorted(glob.glob(os.path.join(GEMMA_DIR, "gene_summary_all.pre-*.tsv")))
    found = {}
    for path in files:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header is None:
                continue
            try:
                var_i = header.index("Variable")
                coord_i = header.index("Gene_Coordinates")
                role_i = header.index("Putative_Adaptive_Role")
                cit_i = header.index("Citation")
                acr_i = header.index("Acronym")
            except ValueError:
                continue
            for row in reader:
                if len(row) <= max(var_i, coord_i, role_i, cit_i, acr_i):
                    continue
                key = (row[var_i], row[coord_i])
                role = row[role_i]
                if key in found or role in ("NA", "", None):
                    continue
                found[key] = (role, row[cit_i], row[acr_i])
    return found


IDENTITY_OPENER_RE = re.compile(r'^([A-Z][A-Za-z0-9]{1,9})\s+(?:is|encodes?|acts?\s+as|functions?\s+as|belongs?\s+to)\b')
# Sentence openers that match IDENTITY_OPENER_RE but are not gene symbols
NON_IDENTITY_OPENERS = {"This", "The", "It", "A", "An", "These", "That"}


def verify_role_identity_consistency(out_rows):
    """Flags output rows whose role text opens by asserting a specific identity
    ("GIS3 is a transcription factor...") not found in the row's current
    Orthologs/Acronym/Functional_Description. Heuristic; prints a report only."""
    flagged = []
    for row in out_rows:
        role = row["Putative_Adaptive_Role"]
        if role == "NA":
            continue
        m = IDENTITY_OPENER_RE.match(role)
        if not m or m.group(1) in NON_IDENTITY_OPENERS:
            continue
        claimed = m.group(1)
        current_text = " ".join([row["Orthologs"], row["Acronym"], row["Functional_Description"]])
        if claimed not in current_text:
            flagged.append((row["Variable"], row["Gene_Name"], claimed, row["Acronym"], role[:70]))
    if flagged:
        print(f"\nWARNING: {len(flagged)} row(s) have a Putative_Adaptive_Role that opens by "
              f"asserting a specific identity not found in that row's current Orthologs/Acronym/"
              f"Functional_Description -- possible stale reuse, review before trusting:",
              file=sys.stderr)
        for f in flagged:
            print(f"  {f}", file=sys.stderr)
    else:
        print("\nIdentity consistency check: no suspicious rows found.", file=sys.stderr)


ROLE_TOKEN_RE = re.compile(
    r'\b([A-Z][a-z]+[A-Z][A-Za-z0-9]*|[A-Z]{2,10}[0-9]*|[A-Z][A-Za-z]*[0-9]+[A-Za-z]*)\b'
)
# Generic abbreviations that look like gene symbols, ignored by verify_role_acronym_tokens
ROLE_TOKEN_GENERIC = {
    "BLAST", "UNIPROT", "LOC", "GO", "ER", "PI", "PQ", "GPI", "GDI", "RAB",
    "COPII", "DNA", "RNA", "ATP", "ADP", "GTP", "ABA", "UV", "ROS", "CO2",
    "PH", "NADPH", "NADH", "FAD", "IAA", "GA", "JA", "SA", "QTL", "SNP",
    "LD", "GWAS", "PPR", "RING", "CCCH", "AGC", "PIN", "VIII", "ACA",
    "ABCC", "MRP", "F-BOX", "SCF", "E3", "GLCN",
}


def verify_role_acronym_tokens(out_rows):
    """For rows with Acronym == "NA", flags role text whose first 20 words name
    a gene-symbol-like token (e.g. "RxLR", "STY13") not found in the row's
    current Orthologs/Functional_Description. Heuristic; prints a report only."""
    flagged = []
    for row in out_rows:
        role = row["Putative_Adaptive_Role"]
        if role == "NA" or row["Acronym"] != "NA":
            continue
        scope = " ".join(role.split()[:20])
        tokens = {m.group(1) for m in ROLE_TOKEN_RE.finditer(scope)}
        tokens = {t for t in tokens if t.upper() not in ROLE_TOKEN_GENERIC}
        if not tokens:
            continue
        current_text = " ".join([row["Orthologs"], row["Acronym"], row["Functional_Description"]])
        if not any(t in current_text for t in tokens):
            flagged.append((row["Variable"], row["Gene_Name"], sorted(tokens), role[:70]))
    if flagged:
        print(f"\nWARNING: {len(flagged)} row(s) with Acronym=NA have a Putative_Adaptive_Role "
              f"naming a specific gene family/motif token not found anywhere in that row's current "
              f"Orthologs/Functional_Description -- possible stale reuse via the both-NA reuse "
              f"branch, review before trusting:", file=sys.stderr)
        for f in flagged:
            print(f"  {f}", file=sys.stderr)
    else:
        print("Acronym-token consistency check: no suspicious rows found.", file=sys.stderr)


def main():
    variables = sys.argv[1:]
    if not variables:
        sys.exit("Usage: build_gene_summary_all.py <var1> [<var2> ...]")

    prior_roles = load_prior_adaptive_role()
    print(f"Loaded {len(prior_roles)} prior (Variable, Coordinates) -> "
          f"Adaptive_Role/Citation pairs from master + backups", file=sys.stderr)

    out_rows = []
    n_reused = 0
    n_placeholder = 0
    n_invalidated = 0
    for var in variables:
        gene_rows = load_gene_summary(var)
        snp_stats = load_snp_stats(var)
        sig_go = load_significant_go(var)

        for gid, row in sorted(gene_rows.items()):
            snp_count, min_p = snp_stats.get(gid, ("NA", "NA"))
            go_terms = "; ".join(sig_go.get(gid, [])) or "NA"

            key = (var, row["Coordinates"])
            prior = prior_roles.get(key)
            # Reuse only if the stored Acronym matches the current one (both NA counts as a match)
            if prior is not None:
                role, citation, acronym_then = prior
                current_acronym = row["Acronym"] if row["Acronym"] != "NA" else None
                acronym_then = acronym_then if acronym_then not in ("NA", "", None) else None
                if acronym_then == current_acronym:
                    n_reused += 1
                else:
                    role, citation = "NA", "NA"
                    n_invalidated += 1
                    n_placeholder += 1
            else:
                role, citation = "NA", "NA"
                n_placeholder += 1

            out_rows.append({
                "Variable": var,
                "Gene_Name": row["Gene_Name"],
                "Gene_Coordinates": row["Coordinates"],
                "Orthologs": row["Orthologs"],
                "BLAST_Gene_Names": row["BLAST_Gene_Names"],
                "Acronym": row["Acronym"],
                "Functional_Description": row["Function_Description"],
                "Significant_GO_Terms": go_terms,
                "SNP_Count": snp_count,
                "Min_P_Value": min_p,
                "Mutation_Types": row["Mutation_Types"],
                "Impact_Levels": row["Impact_Levels"],
                "Putative_Adaptive_Role": role,
                "Citation": citation,
            })

    if os.path.exists(MASTER):
        stamp = datetime.now().strftime("%Y%m%d-%H%M")
        backup = os.path.join(GEMMA_DIR, f"gene_summary_all.pre-merge-{stamp}.tsv")
        shutil.copy2(MASTER, backup)
        print(f"Backed up existing master sheet to {backup}", file=sys.stderr)

    with open(MASTER, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows to {MASTER} "
          f"({n_reused} with reused Adaptive_Role/Citation, "
          f"{n_placeholder} pending literature lookup, "
          f"{n_invalidated} of those invalidated by an Acronym mismatch "
          f"vs. their historical entry)", file=sys.stderr)

    verify_role_identity_consistency(out_rows)
    verify_role_acronym_tokens(out_rows)


if __name__ == "__main__":
    main()
