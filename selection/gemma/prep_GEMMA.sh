#!/bin/bash
# Prepare GEMMA inputs for >4x individuals: PCs and environmental variables (PC_env_creation.R),
# PLINK bed/bim/fam with renumbered scaffolds and SNP IDs, and scaffold_snp_map.txt.

export BASE_DIR=/private/groups/meyerlab/esteiger

GITHUB="${BASE_DIR}/gvtp_popgen"

# ---------------------------------------------------------------------------
echo ""
echo "=== Create PCs and env variables for >4x individuals ==="

# Requires PCAngsd results for >4x samples.
Rscript ${GITHUB}/selection/gemma/PC_env_creation.R

# ---------------------------------------------------------------------------
echo ""
echo "=== PLINK files of maf-filtered sites for >4x individuals ==="
echo "=== for selection work and other things ==="

ANGSD_DIR="${BASE_DIR}/variant_calling/angsd"
GEMMA_DIR="${BASE_DIR}/selection/gemma"
mkdir -p "$GEMMA_DIR"

VCF="${ANGSD_DIR}/postParalog_MAF_4x_GLimp.vcf"
PLINK="${BASE_DIR}/software/plink"
CLEAN_BASE="${GEMMA_DIR}/$(basename ${VCF%.vcf})_clean"

# Strip "_rmdup.bam" from sample names (underscores cause problems in PLINK)
bcftools query -l ${VCF} | sed 's/_rmdup\.bam//g' > "${GEMMA_DIR}/new_sample_names.txt"
bcftools reheader --samples "${GEMMA_DIR}/new_sample_names.txt" ${VCF} -o "${CLEAN_BASE}.vcf"
rm "${GEMMA_DIR}/new_sample_names.txt"

${PLINK} --vcf "${CLEAN_BASE}.vcf" --make-bed \
         --double-id \
         --remove <(printf '%s\n' 4548 4549 4550 4551 4552 4553 4554 4555 4556 4557 4559 4561 4562 | awk '{print $1, $1}') \
         --mac 1 --allow-extra-chr --allow-no-sex --out "${CLEAN_BASE}"
rm "${CLEAN_BASE}.vcf"

# Map GEMMA ids back to scaffolds; must run before the .bim is renamed
# (columns: GEMMA chr, SNP id, scaffold, position)
awk 'BEGIN{OFS="\t"} !($1 in id){id[$1]=++n} {print id[$1], "snp"NR, $1, $4}' \
    "${CLEAN_BASE}.bim" > "${GEMMA_DIR}/scaffold_snp_map.txt"

# Re-name the bim file so that snps are snp1, snp2, etc.
awk 'BEGIN { snp=1 } { $2 = "snp"snp++; print }' OFS='\t' "${CLEAN_BASE}.bim" \
| awk '{
    if (!($1 in seen)) {
        seen[$1] = ++count
    }
    $1 = seen[$1]
    print
}' OFS='\t' > "${GEMMA_DIR}/tmp_renamed.bim"
mv "${GEMMA_DIR}/tmp_renamed.bim" "${CLEAN_BASE}.bim"

# Fix family file
export CLEAN_BASE
Rscript - <<'EOF'
fam <- read.table(paste0(Sys.getenv("CLEAN_BASE"), ".fam"), header=F)
fam$index <- seq(1, nrow(fam))
fam <- fam[,-2]
fam$V1 <- sub("_rmdup\\.bam$", "", fam$V1)

pop <- read.table(file.path(Sys.getenv("BASE_DIR"), "field_data.txt"), header=T, sep="\t")
pop <- pop[,c("SampleID", "Clusters")]

fam <- merge(pop, fam, by.x="SampleID", by.y="V1")
fam <- fam[order(fam$index),]
fam <- fam[,c("Clusters", "SampleID", "V3", "V4", "V5", "V6")]
write.table(fam, paste0(Sys.getenv("CLEAN_BASE"), ".fam"), row.names=F, col.names=F, quote=F, sep=' ')
EOF
