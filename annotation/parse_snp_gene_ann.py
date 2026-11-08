#!/usr/bin/env python3
"""
Read "CHROM\tPOS\tANN" lines on stdin and write "scaffold\tposition\tgene_id"
to stdout, one row per SnpEff ANN entry other than intergenic_region. A SNP
can map to more than one gene (e.g. upstream of one, downstream of another).
Output is not deduplicated; build_snp_gene_map.sh runs `sort -u` on it.

Usage: bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' in.vcf | python3 parse_snp_gene_ann.py
"""
import sys

EXCLUDED_ANNOTATION = "intergenic_region"


def main():
    out = sys.stdout
    for line in sys.stdin:
        line = line.rstrip("\n")
        try:
            chrom, pos, ann = line.split("\t", 2)
        except ValueError:
            continue
        if ann == "." or not ann:
            continue
        for entry in ann.split(","):
            fields = entry.split("|")
            if len(fields) < 5:
                continue
            annotation, gene_id = fields[1], fields[4]
            if annotation == EXCLUDED_ANNOTATION or not gene_id:
                continue
            out.write(f"{chrom}\t{pos}\t{gene_id}\n")


if __name__ == "__main__":
    main()
