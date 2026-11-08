#!/bin/bash
# run_PCAngsd.sh
# Write per-subset sample filter files (1 = include, in bamlist_4x.txt order) and
# submit PCAngsd on the full LD-pruned dataset and each subset.
# Usage: bash run_PCAngsd.sh
# Conda env: PCAngsd
#   conda create -n PCAngsd -c conda-forge python numpy scipy cython \
#     c-compiler cxx-compiler scikit-allel plotly python-kaleido
#   pip install -e <clone of https://github.com/Rosemeis/pcangsd>   # PCAngsd 1.36.0
export BASE_DIR=/private/groups/meyerlab/esteiger

SCRIPT_DIR="${BASE_DIR}/gvtp_popgen/structure"
cd ${BASE_DIR}/structure/PCAngsd

Rscript - <<'EOF'
pop <- read.table(file.path(Sys.getenv("BASE_DIR"), "field_data.txt"), header=T, sep="\t")
pop$file <- paste0(pop$SampleID, "_rmdup.bam")
pop <- unique(pop[,c("SampleID", "Cov", "index", "file", "Location", "Year")])

bam4x <- read.table(file.path(Sys.getenv("BASE_DIR"), "mapping/merged/bamlist_4x.txt"), header=F)
bam4x$bam_index <- seq(1, nrow(bam4x))

merged <- merge(bam4x, pop, by.x="V1", by.y="file")
merged <- merged[order(merged$bam_index),]

year1_bams <- read.table(file.path(Sys.getenv("BASE_DIR"), "mapping/merged/bamlist_year1.txt"), header=F)$V1
strauss_bams <- read.table(file.path(Sys.getenv("BASE_DIR"), "mapping/merged/bamlist_strauss.txt"), header=F)$V1
cdfw_bams <- read.table(file.path(Sys.getenv("BASE_DIR"), "mapping/merged/bamlist_CDFW.txt"), header=F)$V1
coastal_bams <- read.table(file.path(Sys.getenv("BASE_DIR"), "mapping/merged/bamlist_coastal.txt"), header=F)$V1

# > 4x for Year 1
merged$year1 <- ifelse(merged$V1 %in% year1_bams & merged$Cov > 4, 1, 0)
write.table(merged$year1, "year1_4x_subset.txt", row.names=F, col.names=F, quote=F)
cat("year1 >4x:", sum(merged$year1), "\n")

# > 5x for Strauss
merged$strauss <- ifelse(merged$V1 %in% strauss_bams & merged$Cov > 5, 1, 0)
write.table(merged$strauss, "strauss_4x_subset.txt", row.names=F, col.names=F, quote=F)
cat("strauss >5x:", sum(merged$strauss), "\n")

# > 2x for CDFW
merged$cdfw <- ifelse(merged$V1 %in% cdfw_bams & merged$Cov > 2, 1, 0)
write.table(merged$cdfw, "CDFW_4x_subset.txt", row.names=F, col.names=F, quote=F)
cat("CDFW >2x:", sum(merged$cdfw), "\n")

# > 3x for coastal (Pt Conception, Dangermond, Lions Head, Pt Arguello)
merged$coastal <- ifelse(merged$V1 %in% coastal_bams & merged$Cov > 3, 1, 0)
write.table(merged$coastal, "coastal_4x_subset.txt", row.names=F, col.names=F, quote=F)
cat("coastal >3x:", sum(merged$coastal), "\n")
EOF

BEAGLE="postParalog_4x_LD"

# PCAngsd (no filter) on full 4x LD-pruned dataset
sbatch ${SCRIPT_DIR}/PCAngsd_nofilter.sh ${BEAGLE}

# PCAngsd with locality filters
sbatch ${SCRIPT_DIR}/PCAngsd_filter.sh ${BEAGLE} year1_4x
sbatch ${SCRIPT_DIR}/PCAngsd_filter.sh ${BEAGLE} strauss_4x
sbatch ${SCRIPT_DIR}/PCAngsd_filter.sh ${BEAGLE} CDFW_4x
sbatch ${SCRIPT_DIR}/PCAngsd_filter.sh ${BEAGLE} coastal_4x
