#!/usr/bin/env python3
"""
Collect the CENH3 coding sequences used by cenh3_dnds.py into
OUT_DIR/cenh3_cds_raw.fa (one sequence per entry, stop codon removed).

Sources (all stored in CDS_DIR; see cenh3_dnds.py's header for how each CENH3
was identified -- best match to the Deinandra CENH3 domain over canonical H3,
and the CENH3-specific C-terminus):
  Deinandra         genome/remove_duplicates_rm_internalstops_longest_isoform.cds,
                    IDmodified-mRNA-22870 (merged Argonaute + CENH3 model; the
                    CENH3 part is cut out in cenh3_dnds.py)
  Hannuus_Chr14/15  Ensembl Plants HanXRQr2.0-SUNRISE CDS: HanXRQr2_Chr14g0649211,
                    HanXRQr2_Chr15g0703621
  Hanomalus_Chr14/15 NCBI MFS8009295.1, MFS8018859.1 (CDS via efetch fasta_cds_na)
  Hdebilis_Chr14    NCBI KAM0005100.1 (H. debilis subsp. tardiflorus, GCA_052426805.1;
                    KAM0038730.1 is an identical CDS on another scaffold -- an
                    assembly duplicate -- and is not used; no Chr15-type paralogue
                    was found among the 83,904 annotated H. debilis proteins)
  Layia             L. glandulosa alternate-haplotype assembly GCA_046128875.1
                    (the primary assembly lacks CENH3), gene model by miniprot
                    (Deinandra CENH3 protein as query): lalt000001,
                    JBJXVY010000114.1:879,022-882,030 (-), 9 exons. A second copy,
                    lalt000002 (JBJXVY010000140.1:1,500,933-1,504,009, +), differs
                    by 3 nt (1 aa) and is not used (recent duplicate or assembly
                    duplication).
  Lsativa_chr2/chr5 NCBI KAJ0221726.1, KAJ0205668.1 (chr5 is fragmentary and is
                    dropped in cenh3_dnds.py; chr2 is used for the histone fold only)
  Athaliana         Ensembl Plants TAIR10 CDS, AT1G01370.1
Tomato is not included (CDS not retrievable).

Usage: cenh3_collect_cds.py
"""

import gzip

from Bio.Seq import Seq

BASE_DIR = "/private/groups/meyerlab/esteiger"

OUT_DIR = f"{BASE_DIR}/inversions/centromere/cenh3_dnds"
CDS_DIR = f"{OUT_DIR}/cds"
DEIN_CDS = f"{BASE_DIR}/genome/remove_duplicates_rm_internalstops_longest_isoform.cds"


def read(path, keep=lambda h: True):
    op = gzip.open if path.endswith(".gz") else open
    out, cur = {}, None
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                cur = l[1:].strip() if keep(l[1:].strip()) else None
                if cur:
                    out[cur] = []
            elif cur:
                out[cur].append(l.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def only(d):
    assert len(d) == 1, list(d)
    return next(iter(d.values()))


def main():
    seqs = {"Deinandra": only(read(DEIN_CDS, lambda h: h.split()[0] == "IDmodified-mRNA-22870"))}
    ha = read(f"{CDS_DIR}/Helianthus_annuus.HanXRQr2.0-SUNRISE.cds.all.fa.gz",
              lambda h: "HanXRQr2_Chr14g0649211" in h or "HanXRQr2_Chr15g0703621" in h)
    for k, v in ha.items():
        seqs["Hannuus_" + ("Chr14" if "Chr14" in k else "Chr15")] = v
    for acc, name in [("MFS8009295.1", "Hanomalus_Chr14"), ("MFS8018859.1", "Hanomalus_Chr15"),
                      ("KAM0005100.1", "Hdebilis_Chr14"),
                      ("KAJ0221726.1", "Lsativa_chr2"), ("KAJ0205668.1", "Lsativa_chr5")]:
        seqs[name] = only(read(f"{CDS_DIR}/{acc}.fa"))
    seqs["Layia"] = read(f"{CDS_DIR}/Layia_glandulosa_alt_CENH3_miniprot.fa")["lalt000001"]
    seqs["Athaliana"] = only(read(f"{CDS_DIR}/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz",
                                  lambda h: h.startswith("AT1G01370.1 ")))
    with open(f"{OUT_DIR}/cenh3_cds_raw.fa", "w") as out:
        for k, v in seqs.items():
            v = v[:len(v) - len(v) % 3]
            if str(Seq(v[-3:]).translate()) == "*":
                v = v[:-3]
            p = str(Seq(v).translate())
            assert "*" not in p, (k, p)
            print(f"{k:16s} {len(v):4d} nt  {p}")
            out.write(f">{k}\n{v}\n")


if __name__ == "__main__":
    main()
