#!/bin/bash

GITHUB="/private/groups/meyerlab/esteiger/gvtp_popgen"

# ---------------------------------------------------------------------------
echo ""
echo "=== Create PCs and env variables for >4x individuals ==="

# PCANGSD should have already been run on samples >4x... 
Rscript ${GITHUB}/selection/gemma/PC_env_creation.R

# ---------------------------------------------------------------------------
echo ""
echo "=== PLINK files of maf-filtered sites for >4x individuals ==="
echo "=== for selection work and other things ==="

ANGSD_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd"
GEMMA_DIR="/private/groups/meyerlab/esteiger/selection/gemma"
mkdir -p "$GEMMA_DIR"

VCF="${ANGSD_DIR}/postParalog_MAF_4x_GLimp.vcf"
PLINK="/private/groups/meyerlab/esteiger/software/plink"
CLEAN_BASE="${GEMMA_DIR}/$(basename ${VCF%.vcf})_clean"

# Make sure the VCF doesn't contain underscores, which causes issues with PLINK
# Create sample names file
bcftools query -l ${VCF} | sed 's/_rmdup\.bam//g' > "${GEMMA_DIR}/new_sample_names.txt"
# Reheader with the file
bcftools reheader --samples "${GEMMA_DIR}/new_sample_names.txt" ${VCF} -o "${CLEAN_BASE}.vcf"
# Clean up
rm "${GEMMA_DIR}/new_sample_names.txt"

# Now use PLINK on the cleaned VCF
${PLINK} --vcf "${CLEAN_BASE}.vcf" --make-bed \
         --double-id \
         --remove <(printf '%s\n' 4548 4549 4550 4551 4552 4553 4554 4555 4556 4557 4559 4561 4562 | awk '{print $1, $1}') \
         --mac 1 --allow-extra-chr --allow-no-sex --out "${CLEAN_BASE}"
rm "${CLEAN_BASE}.vcf"

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

pop <- read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop <- pop[,c("SampleID", "Clusters")]

fam <- merge(pop, fam, by.x="SampleID", by.y="V1")
fam <- fam[order(fam$index),]
fam <- fam[,c("Clusters", "SampleID", "V3", "V4", "V5", "V6")]
write.table(fam, paste0(Sys.getenv("CLEAN_BASE"), ".fam"), row.names=F, col.names=F, quote=F, sep=' ')
EOF


# Create mapping of old to new names (new_chr, new_snp, orig_chr, orig_pos)
# Must be done before renaming the .bim
awk 'BEGIN { snp=1 } {
    orig_chr = $1
    orig_snp = $4
    new_snp = "snp"snp++
    print orig_chr"\t"orig_snp"\t"orig_chr"\t"new_snp
}' OFS='\t' "${CLEAN_BASE}.bim" > "${GEMMA_DIR}/mapping_step1.txt"
awk '{
    orig_chr = $1
    orig_snp = $2
    renamed_chr = $3
    renamed_snp = $4
    if (!(renamed_chr in seen)) {
        seen[renamed_chr] = ++count
    }
    new_chr = seen[renamed_chr]
    print new_chr"\t"renamed_snp"\t"orig_chr"\t"orig_snp
}' OFS='\t' "${GEMMA_DIR}/mapping_step1.txt" > "${GEMMA_DIR}/scaffold_snp_map.txt"
rm "${GEMMA_DIR}/mapping_step1.txt"
