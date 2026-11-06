#!/bin/bash

#SBATCH --job-name=windowed_coverage
#SBATCH --partition=long    
#SBATCH --nodes=1    	
#SBATCH --ntasks=1                  
#SBATCH --mem=1gb                     # 1gb per task with winPCA 
#SBATCH --time=168:00:00              # long max is 14 days * 24 hrs = 168 hrs
#SBATCH --output=serial_test_%j.log   # Standard output and error log

cd /private/groups/meyerlab/esteiger/inversions

VCF="/private/groups/meyerlab/esteiger/variant_calling/angsd/paralog_filtered/input_split_13_chunk_postNGSparalog.bcf"

# Extract sample names from the VCF
SAMPLES=$(bcftools query -l ${VCF} | tr '\n' '\t')

# Process each chromosome separately
bcftools query -f '%CHROM\t%POS\t[%DP\t]\n' ${VCF} | awk -v OFS="\t" -v WIN=500000000 -v STEP=10000 -v SAMPLES="$SAMPLES" '
BEGIN {
    split(SAMPLES, sample_names, "\t");
}

{
    chrom = $1;
    pos = $2;

    # If on a new chromosome, initialize tracking
    if (chrom != prev_chrom) {
        if (NR > 1) close(output_file);  # Close previous file
        start = 0;
        output_file = chrom ".depth.tsv";
        print "pos\t" SAMPLES > output_file;  # Write header
        prev_chrom = chrom;
        delete depth_values;
    }

    # Store depths for this position
    for (i = 3; i <= NF; i++) {
        depth_values[pos, i] = $i;
    }

    # Move the window and calculate depth averages
    while (pos >= start + STEP) {
        midpoint = start + (WIN / 2);
        sum_depths = ""; count = 0;

        # Sum depths over the full 500Mb window
        for (p = start; p < start + WIN; p++) {
            if ((p, 3) in depth_values) {
                count++;
                for (i = 3; i <= NF; i++) {
                    sum[i] += depth_values[p, i];
                }
            }
        }

        # Compute averages and write to file
        for (i = 3; i <= NF; i++) {
            avg_depth = (count > 0) ? sum[i] / count : "NA";
            sum_depths = sum_depths (i > 3 ? OFS : "") avg_depth;
            sum[i] = 0;  # Reset sum for next window
        }

        print midpoint, sum_depths >> output_file;

        start += STEP;  # Slide the window
    }
}
'