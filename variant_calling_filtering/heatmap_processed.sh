#!/bin/bash
# Bin an ngsLD output into 100 kb x 100 kb windows (second-highest r2 per bin
# pair) and plot an LD heatmap.
# Usage: heatmap_processed.sh <file.ld>  (run by launch_heatmap_plot.sh)
FILE="$1"
BIN_SIZE=100000
HALF_BIN=50000
OUTPUT="${FILE%.ld}_${BIN_SIZE}window.txt"
PLOT="${FILE%.ld}_${BIN_SIZE}window.jpg"

if [ -f "$OUTPUT" ]; then
    echo "Skipping binning - $OUTPUT already exists"
else
    echo "Binning and filtering..."
    awk -v bin=$BIN_SIZE -v half=$HALF_BIN '
BEGIN {OFS="\t"}
NR==1 {next}
$1 != $2 && $3 != "inf" && $3 != "nan" && $3 ~ /^[0-9]/ {
    bin1 = int($1/bin)*bin + half;
    bin2 = int($2/bin)*bin + half;
    key = bin1 "_" bin2;
    r2 = $3 + 0;
    if (key in max1) {
        if (r2 > max1[key]) {
            max2[key] = max1[key];
            max1[key] = r2;
        } else if (r2 > max2[key] || !(key in max2)) {
            max2[key] = r2;
        }
    } else {
        max1[key] = r2;
        max2[key] = r2;
    }
}
END {
    print "site1_bin\tsite2_bin\tr2";
    for (key in max2) {
        split(key, bins, "_");
        print bins[1] "\t" bins[2] "\t" max2[key];
    }
}' "$FILE" > "$OUTPUT"
fi

echo "Creating plot in R..."
Rscript - "$OUTPUT" "$PLOT" <<'EOF'


args <- commandArgs(trailingOnly=TRUE)
library(ggplot2)

data <- read.table(args[1], header=TRUE, sep="\t")

jpeg(args[2], width=1200, height=1200, res=150)

# Calculate breaks for major and minor grids
x_range <- range(data$site1_bin)
y_range <- range(data$site2_bin)
total_range <- max(x_range[2] - x_range[1], y_range[2] - y_range[1])

# Major gridlines (fewer, bolder)
if (total_range > 50e6) {
  major_spacing <- 10e6
  minor_spacing <- 2e6
  label_format <- function(x) paste0(round(x/1e6, 1), " Mb")
} else if (total_range > 10e6) {
  major_spacing <- 5e6
  minor_spacing <- 1e6
  label_format <- function(x) paste0(round(x/1e6, 1), " Mb")
} else if (total_range > 1e6) {
  major_spacing <- 1e6
  minor_spacing <- 200e3
  label_format <- function(x) paste0(round(x/1e6, 1), " Mb")
} else {
  major_spacing <- 500e3
  minor_spacing <- 100e3
  label_format <- function(x) paste0(round(x/1e3), " kb")
}

# Major breaks (labeled)
x_major <- seq(floor(x_range[1]/major_spacing) * major_spacing,
               ceiling(x_range[2]/major_spacing) * major_spacing,
               by = major_spacing)
y_major <- seq(floor(y_range[1]/major_spacing) * major_spacing,
               ceiling(y_range[2]/major_spacing) * major_spacing,
               by = major_spacing)

# Minor breaks (unlabeled)
x_minor <- seq(floor(x_range[1]/minor_spacing) * minor_spacing,
               ceiling(x_range[2]/minor_spacing) * minor_spacing,
               by = minor_spacing)
y_minor <- seq(floor(y_range[1]/minor_spacing) * minor_spacing,
               ceiling(y_range[2]/minor_spacing) * minor_spacing,
               by = minor_spacing)

ggplot(data, aes(x=site1_bin, y=site2_bin, fill=r2)) +
  geom_tile() +
  # Add custom gridlines as geom_vline/hline for more control
  geom_vline(xintercept = x_minor, color="white", size=0.2, alpha=0.5) +
  geom_hline(yintercept = y_minor, color="white", size=0.2, alpha=0.5) +
  geom_vline(xintercept = x_major, color="white", size=0.6) +
  geom_hline(yintercept = y_major, color="white", size=0.6) +
  scale_fill_gradientn(
    colors=c("navy","blue","cyan","yellow","orange","red"),
    name=expression(LD~(r^2)),
    limits=c(0,1)
  ) +
  scale_x_continuous(
    breaks = x_major,
    labels = label_format,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = y_major,
    labels = label_format,
    expand = c(0, 0)
  ) +
  labs(x="Position", y="Position") +
  coord_fixed() +  # Equal aspect ratio for easier reading
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove default grid since we added custom
    axis.text.x = element_text(angle=45, hjust=1, size=11, face="bold", color="black"),
    axis.text.y = element_text(size=11, face="bold", color="black"),
    axis.title = element_text(size=13, face="bold"),
    panel.background = element_rect(fill="gray98", color=NA),
    panel.border = element_rect(fill=NA, color="black", size=1),
    plot.background = element_rect(fill="white"),
    legend.position = "right",
    legend.key.height = unit(2.5, "cm"),
    legend.title = element_text(size=12, face="bold")
  )

dev.off()
EOF

echo "Done! Created $PLOT"

exit 0
