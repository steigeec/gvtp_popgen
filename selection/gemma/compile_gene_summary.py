#!/usr/bin/env python3
"""
Compile per-gene summary from GEMMA outlier analysis outputs.
Usage: python3 compile_gene_summary.py <var_name>
"""
import sys
import os
import re
import csv
import time
from collections import defaultdict
import urllib.request
import urllib.parse

sys.path.insert(0, "/private/groups/meyerlab/esteiger/gvtp_popgen/selection/gemma")
from curate_intergenic_by_ld import GENE_TOKEN_RE, load_gff_resolver, resolve_to_gene

var_name = sys.argv[1]
GFF = "/private/groups/meyerlab/esteiger/genome/remove_duplicates_rm_internalstops_longest_isoform.gff"

# Shared token-extraction + resolution logic (same as curate_intergenic_by_ld.py
# and compute_gene_snp_stats.py use to build the gene list itself), so ANN
# parsing below resolves compound "geneA-geneB" intergenic annotations to
# BOTH real genes independently, instead of a fragile substring scan.
_gene_ids_gff, _tx2gene_shared = load_gff_resolver()

# --- 1. Read gene list ---
gene_list_file = os.path.join(var_name, f"gene_list_{var_name}_unfiltered.txt")
gene_ids = []
with open(gene_list_file) as f:
    for line in f:
        gid = line.strip().rstrip('.')
        if gid:
            gene_ids.append(gid)
gene_id_set = set(gene_ids)

# Genes also supported by a 40000_3 (multi-SNP peak) outlier, if that scope
# was run this pass — absent file just means "not run yet", not "no peaks".
peak_40000_3_file = os.path.join(var_name, f"gene_list_{var_name}_40000_3.txt")
genes_in_40000_3_peak = set()
if os.path.exists(peak_40000_3_file):
    with open(peak_40000_3_file) as f:
        for line in f:
            gid = line.strip().rstrip('.')
            if gid:
                genes_in_40000_3_peak.add(gid)

# --- 2. Parse ANN INFO file ---
gene_ann_name = {}
gene_mutations = defaultdict(set)
gene_impacts = defaultdict(set)

ann_file = os.path.join(var_name, f"outliers_{var_name}.ann.INFO")
if os.path.exists(ann_file) and os.path.getsize(ann_file) > 0:
    with open(ann_file) as f:
        f.readline()
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 5 or cols[4] == '.':
                continue
            for entry in cols[4].split(','):
                parts = entry.split('|')
                if len(parts) < 5:
                    continue
                annotation = parts[1]
                impact = parts[2]
                gene_id_ann = parts[4]
                # gene_id_ann is "geneA-geneB" for intergenic_region (both
                # flanking genes named at once) and a single id otherwise.
                # Resolve and credit EVERY token independently -- both
                # flanking genes are legitimate, separately-LD-verified
                # candidates (see curate_intergenic_by_ld.py), not one
                # combined entity, so neither should be skipped.
                for tok in GENE_TOKEN_RE.findall(gene_id_ann):
                    gid = resolve_to_gene(tok, _gene_ids_gff, _tx2gene_shared)
                    if gid is None or gid not in gene_id_set:
                        continue
                    gene_mutations[gid].add(annotation)
                    gene_impacts[gid].add(impact)
                    if gid not in gene_ann_name:
                        gene_ann_name[gid] = tok

# --- 3. Build gene_id <-> coordinates mapping from GFF (moved ahead of the
# orthologs read below, since matching orthologs by gene now needs tx2gene) ---
# gene_list_{var_name}.txt sometimes carries entries that aren't a clean gene
# id: stray transcript ids (e.g. "g14399" whose real gene is "g4525" — see
# resolve_gene_ids.R for the same issue upstream) and SnpEff's compound
# "geneA-geneB" names for intergenic SNPs. resolve_coord() below handles both
# so the Coordinates column isn't NA more often than it has to be.
gene_id_to_coord_all = {}
tx2gene = {}
gene_id_pattern = re.compile(r'gene_id=([^;]+)')
mrna_id_pattern = re.compile(r'ID=([^;]+)')
mrna_parent_pattern = re.compile(r'Parent=([^;]+)')
with open(GFF) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue
        # tx2gene from the mRNA lines' own ID=/Parent= attributes, NOT the
        # gene_id=/transcript_id= attributes (which can be stale — see
        # lookup_by_gene_id.py header for the same underlying issue).
        if fields[2] == 'mRNA':
            tm = mrna_id_pattern.search(fields[8])
            pm = mrna_parent_pattern.search(fields[8])
            if tm and pm and tm.group(1) not in tx2gene:
                tx2gene[tm.group(1)] = pm.group(1)
        if fields[2] == 'gene':
            gm = gene_id_pattern.search(fields[8])
            if gm:
                gid = gm.group(1)
                if gid not in gene_id_to_coord_all:
                    gene_id_to_coord_all[gid] = f"{fields[0]}:{fields[3]}-{fields[4]}"

# --- 4. Read orthologs ---
# orthologs_{var}.txt (build_gene_orthology.py) is already keyed directly by
# the real resolved Gene_ID -- no transcript/prefix resolution needed here,
# unlike the old lookup_by_gene_id.py-based file. It also carries a direct
# Acronym/Function_Description for the OrthoFinder:Arabidopsis_thaliana and
# eggNOG tiers (curated symbol + description, no UniProt query needed for
# those genes -- see section 6). BLAST_Gene_Names comes from the same file
# now too (the old separate blast_{var}.filtered.txt was stale/orphaned --
# nothing in the active pipeline regenerated it).
gene_orthologs = {}
gene_ortholog_source = {}
gene_acronym_direct = {}
gene_function_direct = {}
gene_blast_names = {}

ortho_file = os.path.join(var_name, f"orthologs_{var_name}.txt")
if os.path.exists(ortho_file):
    with open(ortho_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gid = row["Gene_ID"].rstrip('.')
            if gid not in gene_id_set:
                continue
            if row["Orthologs"] != "NA":
                gene_orthologs[gid] = row["Orthologs"]
                gene_ortholog_source[gid] = row["Ortholog_Source"]
            if row["Acronym"] != "NA":
                gene_acronym_direct[gid] = row["Acronym"]
            if row["Function_Description"] != "NA":
                gene_function_direct[gid] = row["Function_Description"]
            if row["BLAST_Gene_Names"] != "NA":
                gene_blast_names[gid] = row["BLAST_Gene_Names"]

def resolve_coord(name):
    name = name.rstrip(".")
    if name in gene_id_to_coord_all:
        return gene_id_to_coord_all[name]
    for suffix in (".t1", ".1"):
        cand = tx2gene.get(name + suffix)
        if cand and cand in gene_id_to_coord_all:
            return gene_id_to_coord_all[cand]
    for tx, gid in tx2gene.items():
        if tx.startswith(name + ".") and gid in gene_id_to_coord_all:
            return gene_id_to_coord_all[gid]
    if "-" in name:
        found = [resolve_coord(p) for p in name.split("-")]
        found = [x for x in found if x]
        if found:
            return "; ".join(dict.fromkeys(found))
    return None

gene_id_to_coord = {gid: (resolve_coord(gid) or "NA") for gid in gene_ids}

# --- 5. Query UniProt for plant functional descriptions ---
# Only needed as a fallback for genes without a direct Acronym/Function_Description
# already supplied by build_gene_orthology.py (see section 4) -- i.e. genes
# whose only hit was the strict-evalue sunflower BLAST tier, or with no
# ortholog call at all.
# Prefer Asteraceae, then flowering plants, then all green plants
TAXONOMY_TIERS = [
    ("Asteraceae", 4210),
    ("Magnoliophyta", 3398),
    ("Viridiplantae", 33090),
]

def query_uniprot_taxon(search_term, taxon_id, reviewed_only=False):
    """Query UniProt for a gene function + primary gene symbol within a given taxonomy.
    reviewed_only restricts to curated Swiss-Prot entries, which are far more
    likely to have a gene_primary symbol populated than automatic/TrEMBL ones."""
    try:
        q = f'{search_term} AND taxonomy_id:{taxon_id}'
        if reviewed_only:
            q += ' AND reviewed:true'
        query = urllib.parse.quote(q)
        url = (f'https://rest.uniprot.org/uniprotkb/search?query={query}'
               f'&fields=protein_name,gene_primary,cc_function,organism_name&format=tsv&size=1')
        req = urllib.request.Request(url, headers={'User-Agent': 'Python/gemma_pipeline'})
        with urllib.request.urlopen(req, timeout=15) as resp:
            lines = resp.read().decode().strip().split('\n')
            if len(lines) > 1:
                cols = lines[1].split('\t')
                protein_name = cols[0].strip() if cols[0] else ""
                acronym = cols[1].strip() if len(cols) > 1 and cols[1] else ""
                # NCBI LOC###### placeholders aren't real gene symbols — just
                # another anonymous ID — so don't surface them as an "acronym".
                if re.match(r'^LOC\d+$', acronym, re.IGNORECASE):
                    acronym = ""
                function = ""
                organism = cols[3].strip() if len(cols) > 3 and cols[3] else ""
                if len(cols) > 2 and cols[2]:
                    function = re.sub(r'FUNCTION:\s*', '', cols[2])
                    function = re.sub(r'\{[^}]*\}', '', function).strip().rstrip('.')
                if function:
                    desc = f"{protein_name}: {function}"
                elif protein_name:
                    desc = protein_name
                else:
                    return None
                if organism:
                    desc += f" [{organism}]"
                return {"desc": desc, "acronym": acronym}
    except Exception:
        pass
    return None

def query_uniprot(search_term):
    """Query UniProt across plant taxonomic tiers, preferring Asteraceae.
    Two full passes across ALL tiers: reviewed-only first (best shot at a
    real gene symbol), then unrestricted (for at least a description) only
    if no tier had a reviewed hit. A narrow tier's unreviewed/uncurated hit
    (e.g. an unnamed Asteraceae relative) must not shadow a broader tier's
    well-curated one (e.g. Arabidopsis) that actually has a gene symbol —
    which is what tier-by-tier reviewed-then-unreviewed fallback used to do."""
    for taxon_name, taxon_id in TAXONOMY_TIERS:
        result = query_uniprot_taxon(search_term, taxon_id, reviewed_only=True)
        if result:
            return result
        time.sleep(0.15)
    for taxon_name, taxon_id in TAXONOMY_TIERS:
        result = query_uniprot_taxon(search_term, taxon_id, reviewed_only=False)
        if result:
            return result
        time.sleep(0.15)
    return None

def clean_for_search(text):
    """Strip NCBI-specific 'isoform X#' suffixes before querying UniProt.
    These are RefSeq auto-annotation artifacts absent from curated entries;
    left in, they skew UniProt's ranking toward whichever poorly-annotated
    species' record happens to contain that literal phrase, instead of the
    best-curated match (e.g. "protein GRAVITROPIC IN THE LIGHT 1 isoform X2"
    matches an uncharacterized squash LOC gene instead of Arabidopsis GIL1).
    Only affects the search term — display columns keep the original text."""
    return re.sub(r'\s+isoform\s+\S+$', '', text, flags=re.IGNORECASE).strip()

# Collect unique search terms and batch lookup -- skip genes that already
# have a direct Acronym or Function_Description from build_gene_orthology.py
# (OrthoFinder:Arabidopsis_thaliana / eggNOG tiers): querying UniProt for
# those would just add a slower, less-curated second opinion.
function_cache = {}
search_terms_per_gene = {}

for gid in gene_ids:
    if gid in gene_acronym_direct or gid in gene_function_direct:
        continue
    terms = []
    ortho_text = gene_orthologs.get(gid)
    if ortho_text:
        cleaned = re.sub(r'^[A-Z]{2}_\d+\.\d+\s*', '', ortho_text)
        cleaned = re.sub(r'\s*\[.*?\]\s*$', '', cleaned)
        cleaned = re.sub(r',?\s*partial$', '', cleaned)
        cleaned = re.sub(r'^PREDICTED:\s*', '', cleaned)
        if cleaned and 'uncharacterized' not in cleaned.lower() and 'hypothetical' not in cleaned.lower():
            terms.append(cleaned)
    if not terms:
        bn = gene_blast_names.get(gid)
        if bn:
            cleaned = re.sub(r'^[A-Z]{2}_\d+\.\d+\s*', '', bn)
            if cleaned and 'uncharacterized' not in cleaned.lower() and 'hypothetical' not in cleaned.lower():
                terms.append(cleaned)
    if terms:
        search_terms_per_gene[gid] = clean_for_search(terms[0])

n_direct = len(set(gene_acronym_direct) | set(gene_function_direct))
unique_terms = set(search_terms_per_gene.values())
print(f"{n_direct} gene(s) already have a direct OrthoFinder/eggNOG description "
      f"(no UniProt query needed). Querying UniProt for {len(unique_terms)} unique "
      f"gene names for the rest...")

for i, term in enumerate(sorted(unique_terms)):
    if term not in function_cache:
        result = query_uniprot(term)
        function_cache[term] = result if result else None
        if (i + 1) % 10 == 0:
            print(f"  ... queried {i + 1}/{len(unique_terms)}")
        time.sleep(0.3)

# --- 6. Write summary ---
out_file = os.path.join(var_name, f"gene_summary_{var_name}.tsv")
with open(out_file, 'w') as out:
    # Gene_ID is the real resolved key (same id used in gene_list_*_unfiltered.txt
    # and gene_snp_stats_*.tsv) -- Gene_Name is a display label only, so
    # downstream joins must use Gene_ID, never Gene_Name. Ortholog_Source
    # records which tier (BLAST:Helianthus_annuus / OrthoFinder:Arabidopsis_
    # thaliana / eggNOG / none) produced the Orthologs value, for transparency.
    out.write("Gene_ID\tCoordinates\tGene_Name\tIn_40000_3_Peak\tOrthologs\tOrtholog_Source\t"
               "BLAST_Gene_Names\tAcronym\tFunction_Description\tMutation_Types\tImpact_Levels\n")
    for gid in sorted(gene_ids):
        coord = gene_id_to_coord.get(gid, "NA")
        name = gene_ann_name.get(gid, gid)
        in_peak = "Yes" if gid in genes_in_40000_3_peak else "No"

        ortho = gene_orthologs.get(gid, "NA")
        ortho_source = gene_ortholog_source.get(gid, "none")
        blast = gene_blast_names.get(gid, "NA")

        if gid in gene_acronym_direct or gid in gene_function_direct:
            acronym = gene_acronym_direct.get(gid, "NA")
            func = gene_function_direct.get(gid, "NA")
        else:
            search_term = search_terms_per_gene.get(gid)
            result = function_cache.get(search_term) if search_term else None
            func = result["desc"] if result else "NA"
            acronym = result["acronym"] if result and result.get("acronym") else "NA"

        muts = "; ".join(sorted(gene_mutations.get(gid, set()))) if gene_mutations.get(gid) else "NA"
        imps = "; ".join(sorted(gene_impacts.get(gid, set()))) if gene_impacts.get(gid) else "NA"

        out.write(f"{gid}\t{coord}\t{name}\t{in_peak}\t{ortho}\t{ortho_source}\t"
                  f"{blast}\t{acronym}\t{func}\t{muts}\t{imps}\n")

print(f"Gene summary: {len(gene_ids)} genes written to {out_file}")
