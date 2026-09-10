#!/bin/bash
#SBATCH --job-name=merge_filter_snps
#SBATCH --partition=long
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=esteiger@ucsc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=64gb
#SBATCH --time=168:00:00
#SBATCH --output=merge_filter_snps_%j.log

# Post-processing after variant_call_noMAF.sh completes.
# 1. Merge per-chunk outputs into noMAF_preParalog.* (all SNPs)
# 2. Build noMAF_preParalog.sites (sorted, indexed) for the doSaf/TajD pipeline
# 3. Remove per-chunk files
# NOTE: GEA subsetting (MAF filter + beagle/BCF filter) happens in create_file_subsets.sh
#       AFTER paralog filtering with ngsParalog, since GEA requires paralog-clean sites.

SNPS_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd"
CHUNKS_DIR="/private/groups/meyerlab/esteiger/variant_calling/angsd/chromo_30_chunks"
ANGSD="/private/groups/meyerlab/esteiger/software/angsd/angsd"

export TMPDIR="/private/groups/meyerlab/esteiger/tmp/"
ulimit -n 100000

cd "$SNPS_DIR"

# ---------------------------------------------------------------------------
echo "=== Pre-check: verifying all chunks completed successfully ==="
missing=0
for chunk in "$CHUNKS_DIR"/input_split_*; do
  name=$(basename "$chunk")
  arg="${SNPS_DIR}/${name}.arg"
  if [ ! -f "$arg" ]; then
    echo "  MISSING .arg: ${name}"
    missing=$((missing + 1))
  elif ! tail -1 "$arg" | grep -q "\[ALL done\] walltime used ="; then
    echo "  INCOMPLETE: ${name}"
    missing=$((missing + 1))
  fi
done
if [ "$missing" -gt 0 ]; then
  echo "ERROR: $missing chunk(s) did not complete successfully. Aborting merge."
  exit 1
fi
echo "  All 23 chunks completed successfully. Proceeding."

# ---------------------------------------------------------------------------
echo "=== Step 1: Merging per-chunk outputs ==="

# mafs.gz and geno.gz: cols 1-2 are chromo and position (tab-separated)
for type in mafs geno; do
  if [ -f "noMAF_preParalog.${type}.gz" ]; then
    echo "  ${type}.gz already exists, skipping."
  else
    echo "  ${type}.gz..."
    header_file=$(ls input_split_*.${type}.gz | sort | head -1)
    zcat "$header_file" | head -1 > tmp_header.txt
    find . -name "input_split_*.${type}.gz" | sort | xargs zcat \
      | grep -v "^chromo" \
      | sort -k1,1 -k2,2n | uniq > tmp_data.txt
    cat tmp_header.txt tmp_data.txt | gzip > "noMAF_preParalog.${type}.gz"
    rm tmp_header.txt tmp_data.txt
    echo "    -> $(zcat noMAF_preParalog.${type}.gz | wc -l) lines (including header)"
  fi
done

# beagle.gz: col 1 is chromo_pos; split on last underscore for sorting
if [ -f "noMAF_preParalog.beagle.gz" ]; then
  echo "  beagle.gz already exists, skipping."
else
  echo "  beagle.gz..."
  header_file=$(ls input_split_*.beagle.gz | sort | head -1)
  zcat "$header_file" | head -1 > tmp_header.txt
  find . -name "input_split_*.beagle.gz" | sort | xargs zcat \
    | grep -v "^marker" \
    | awk '{
        n = split($1, parts, "_")
        left = parts[1]; for (i=2; i<n; i++) left = left "_" parts[i]
        print left "\t" parts[n] "\t" $0
      }' \
    | sort -k1,1 -k2,2n \
    | cut -f3- \
    | uniq > tmp_data.txt
  cat tmp_header.txt tmp_data.txt | gzip > "noMAF_preParalog.beagle.gz"
  rm tmp_header.txt tmp_data.txt
  echo "    -> $(zcat noMAF_preParalog.beagle.gz | wc -l) lines (including header)"
fi

# BCF: index each chunk, concat, sort, deduplicate
if [ -f "noMAF_preParalog.bcf" ]; then
  echo "  BCF already exists, skipping."
else
  echo "  BCF..."
  for f in input_split_*.bcf; do bcftools index "$f"; done
  bcftools concat -a $(ls input_split_*.bcf) -o tmp_concat.bcf
  bcftools sort tmp_concat.bcf -O b -o noMAF_preParalog.bcf -T "$TMPDIR"
  bcftools index noMAF_preParalog.bcf
  rm tmp_concat.bcf input_split_*.bcf.csi
  echo "    -> $(bcftools view -H noMAF_preParalog.bcf | wc -l) sites in noMAF_preParalog.bcf"
fi

# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Building noMAF_preParalog.sites (all SNPs) ==="
if [ -f "noMAF_preParalog.sites" ]; then
  echo "  noMAF_preParalog.sites already exists, skipping."
else
  zcat noMAF_preParalog.mafs.gz \
    | grep -v "^chromo" | cut -f1-2 | sort -k1,1 -k2,2n | uniq > noMAF_preParalog.sites
  echo "  Total SNPs: $(wc -l < noMAF_preParalog.sites)"
  $ANGSD sites index noMAF_preParalog.sites
fi

# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: Building preParalog_0.05maf.sites (MAF > 0.05 sites that still have paralog problem sites) ==="
MAF_SITES="/private/groups/meyerlab/esteiger/variant_calling/ngsParalog/preParalog_0.05maf.sites"
if [ -f "$MAF_SITES" ]; then
  echo "  preParalog_0.05maf.sites already exists, skipping."
else
  zcat noMAF_preParalog.mafs.gz \
    | awk 'NR>1 && $6 > 0.05 {print $1 "\t" $2}' \
    > "$MAF_SITES"
  echo "  Sites with MAF > 0.05: $(wc -l < "$MAF_SITES")"
fi

# ---------------------------------------------------------------------------
echo ""
echo "=== Step 4: Removing per-chunk files ==="
rm input_split_*.mafs.gz input_split_*.beagle.gz \
   input_split_*.bcf input_split_*.geno.gz input_split_*.hwe.gz 2>/dev/null
rm input_split_*.log input_split_*.arg input_split_*.snpStat.gz 2>/dev/null
echo "  Done."

echo ""
echo "=== Complete ==="
echo "  ${SNPS_DIR}/noMAF_preParalog.sites (+ .mafs.gz, .beagle.gz, .bcf)"
