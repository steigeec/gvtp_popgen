#!/bin/bash

#SBATCH --job-name=bp_reads
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=4gb
#SBATCH --time=4:00:00
#SBATCH --output=/private/groups/meyerlab/esteiger/inversions/breakpoints/logs/bp_reads_%A_%a.log

# Read-level inversion-breakpoint evidence on DEINVI_00006 for ONE sample
# (array task N = line N of samples.txt; written by breakpoint_samples.py's
# output via: cut -f3 breakpoint_samples.tsv | tail -n +2 | sort -u).
#
# Output: reads/<sample>.tsv, one line per read (MAPQ >= 20, primary only):
#   type  pos  strand  other_pos  other_strand  qname
#   PAIR  -- paired read whose mate maps on DEINVI_00006 > MIN_DIST away on the
#            SAME strand (inward-facing pairs become same-strand across an
#            inversion breakpoint). Both mates are written; a pair counts
#            only if both appear (i.e. both have MAPQ >= 20).
#   SPLIT -- read with a supplementary/secondary part (SA tag, part MAPQ >= 20)
#            on DEINVI_00006 > MIN_DIST away on the OPPOSITE strand.
# plus reads/<sample>.count: number of primary MAPQ >= 20 reads on the scaffold
# (for depth normalization).
#
# bwa was run with -P (independent mate placement, no proper-pair flags) and
# -M (split parts flagged secondary rather than supplementary), so pairs are
# judged from strand + position only, and split parts are taken from SA tags.

set -euo pipefail

OUT=/private/groups/meyerlab/esteiger/inversions/breakpoints
BAM_DIR=/private/groups/meyerlab/esteiger/mapping/merged
CHROM=DEINVI_00006
MIN_DIST=20000

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${OUT}/samples.txt)
BAM=${BAM_DIR}/${SAMPLE}
mkdir -p ${OUT}/reads

samtools view -c -q 20 -F 0x904 ${BAM} ${CHROM} > ${OUT}/reads/${SAMPLE}.count

samtools view -q 20 -F 0x904 ${BAM} ${CHROM} | awk -v OFS='\t' -v chrom=${CHROM} -v d=${MIN_DIST} '
{
    rev = (and($2, 16) > 0) ? "-" : "+"
    # same-strand distant pair on the same scaffold
    if (and($2, 1) && !and($2, 8) && $7 == "=") {
        mrev = (and($2, 32) > 0) ? "-" : "+"
        dist = $8 - $4; if (dist < 0) dist = -dist
        if (dist > d && rev == mrev) print "PAIR", $4, rev, $8, mrev, $1
    }
    # split read: SA:Z:rname,pos,strand,CIGAR,mapQ,NM;...
    for (i = 12; i <= NF; i++) {
        if (substr($i, 1, 5) != "SA:Z:") continue
        n = split(substr($i, 6), parts, ";")
        for (j = 1; j <= n; j++) {
            if (parts[j] == "") continue
            split(parts[j], f, ",")
            dist = f[2] - $4; if (dist < 0) dist = -dist
            if (f[1] == chrom && f[5] >= 20 && dist > d && f[3] != rev)
                print "SPLIT", $4, rev, f[2], f[3], $1
        }
    }
}' > ${OUT}/reads/${SAMPLE}.tsv

echo "$(date): ${SAMPLE}: $(cut -f1 ${OUT}/reads/${SAMPLE}.tsv | sort | uniq -c | tr '\n' ' ') of $(cat ${OUT}/reads/${SAMPLE}.count) reads"
