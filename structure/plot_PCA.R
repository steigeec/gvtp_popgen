#!/usr/bin/env Rscript

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(ggplot2)
library(sf)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_PCA.R <path_to.log>")
}
log_file <- args[1]

gvtp_theme <- theme_bw() +
  theme(
    text = element_text(size = 50),
    axis.text = element_text(angle = 45, hjust=1),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
  )

popcolors <- c("CDFW" = "#e69f00", "CDFW_MR" = "#e69f00", "Strauss"="#009e73",
               "Point_Conception"="#0072b2", "VAFB_Lion_Head"="#f0e442",
               "VAFB_Arguello"="#d55e00", "Dangermond"="#56b4e9")

clustercolors <- c("71"="#EE8866", "72"="#DDDDDD", "73"="#AAAA00", "74"="#AAAA00",
                    "75"="#AAAA00", "76"="#99DDFF", "77"="#EEDD88", "78"="#99DDFF",
                    "79"="#77AADD", "80"="#77AADD", "81"="#77AADD", "82"="#77AADD",
                    "83"="#77AADD", "84"="#77AADD", "85"="#77AADD", "86"="#77AADD",
                    "87"="#77AADD", "88"="#77AADD")

# --- Parse PCAngsd eigenvector count from log ---
parse_n_eigenvectors <- function(log_lines) {
  ev_line <- grep("eigenvectors", log_lines, value=TRUE)
  if (length(ev_line) == 0) return(NULL)
  m <- regmatches(ev_line, regexpr("[0-9]+(?= eigenvector)", ev_line, perl=TRUE))
  if (length(m) > 0) return(as.integer(m)) else return(NULL)
}

# --- Parse log file ---
log_lines <- readLines(log_file)
out_dir <- dirname(log_file)

beagle_path <- trimws(sub(".*--beagle\\s+", "", grep("--beagle", log_lines, value=TRUE)))
out_name <- trimws(sub(".*--out\\s+", "", grep("--out", log_lines, value=TRUE)))
filter_lines <- grep("--filter", log_lines, value=TRUE)
has_filter <- length(filter_lines) > 0

cov_file <- file.path(out_dir, paste0(out_name, ".cov"))

# Derive subset name from output name (strip _filter suffix and angsd_paralog_LD_ prefix)
subset_name <- sub("_filter$", "", out_name)
subset_name <- sub("^angsd_paralog_LD_", "", subset_name)
out_name <- sub("^angsd_paralog_LD_", "", out_name)

# Derive bamlist: try matching subset name, fall back to beagle name, then 4x
beagle_base <- sub("\\.beagle\\.gz$", "", basename(beagle_path))
bamlist_dir <- "/private/groups/meyerlab/esteiger/mapping/merged"
beagle_suffix <- sub("^(angsd_paralog_LD_|postParalog_)", "", beagle_base)
beagle_suffix <- sub("_miss[0-9]+$", "", beagle_suffix)
beagle_suffix_nold <- sub("_LD$", "", beagle_suffix)
bamlist_candidates <- c(
  file.path(bamlist_dir, paste0("bamlist_", beagle_suffix, ".txt")),
  file.path(bamlist_dir, paste0("bamlist_", beagle_suffix_nold, ".txt")),
  file.path(bamlist_dir, paste0("bamlist_", sub("s$", "", beagle_suffix), ".txt")),
  file.path(bamlist_dir, paste0("bamlist_", subset_name, ".txt")),
  file.path(bamlist_dir, paste0("bamlist_", sub("_4x$", "", subset_name), ".txt"))
)
bamlist_path <- bamlist_candidates[file.exists(bamlist_candidates)][1]
if (is.na(bamlist_path)) stop("Could not find matching bamlist for: ", beagle_base)

cat("Subset:", subset_name, "\n")
cat("Cov file:", cov_file, "\n")
cat("Filtered:", has_filter, "\n")
cat("Bamlist:", bamlist_path, "\n")

# Load population metadata
pop <- read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop$Location <- as.factor(pop$Location)
pop$Year <- as.factor(pop$Year)
pop$Cov <- as.numeric(pop$Cov)

# Load study area shapefile for maps
CA <- st_read("/private/groups/meyerlab/esteiger/env/study_area/GVTP_clip.shp", quiet=TRUE)

# Read bamlist to get sample order
sample_order <- read.table(bamlist_path, header=F, sep="\t")
sample_order$bam_index <- seq(1, nrow(sample_order))
sample_order$V1 <- gsub("_rmdup\\.bam$", "", sample_order$V1)

# Merge pop with sample order
my_pop <- merge(pop, sample_order, by.y="V1", by.x="SampleID", all=FALSE)
my_pop <- my_pop[order(my_pop$bam_index),]

# Apply coverage filter if present
if (has_filter) {
  filter_path <- trimws(sub(".*--filter\\s+", "", filter_lines))
  filter_file <- file.path(out_dir, filter_path)
  filter_vec <- read.table(filter_file, header=FALSE)$V1
  my_pop$keep <- filter_vec[my_pop$bam_index]
  my_pop <- my_pop[my_pop$keep == 1,]
}

cat("Samples:", nrow(my_pop), "\n")

# Eigen decomposition
C <- as.matrix(read.table(cov_file))
e <- eigen(C)

# Use PCAngsd's eigenvector count, minimum 2 for PC1 vs PC2
pcangsd_k <- parse_n_eigenvectors(log_lines)
if (is.null(pcangsd_k)) pcangsd_k <- min(4, ncol(C))
n_pcs <- max(pcangsd_k, 2)
cat("PCAngsd eigenvectors:", pcangsd_k, "- plotting", n_pcs, "PCs\n")

# Scree plot with elbow marked
scree <- data.frame(PC = seq_along(e$values), Eigenvalue = e$values)
scree_sub <- scree[1:min(20, nrow(scree)),]
jpeg(file.path(out_dir, paste0(out_name, '_scree.jpg')), height=1000, width=1000)
print(ggplot(scree_sub, aes(x = PC, y = Eigenvalue)) +
  geom_point(size = 3) + geom_line() +
  geom_vline(xintercept = n_pcs, linetype = "dashed", color = "red") +
  annotate("text", x = n_pcs + 0.5, y = max(scree_sub$Eigenvalue) * 0.9,
           label = paste("k =", n_pcs), hjust = 0, color = "red", size = 5) +
  labs(x = "Principal Component", y = "Eigenvalue", title = paste(out_name, "Scree Plot")) +
  theme_bw())
dev.off()

# Variance explained
var_exp <- round((e$values / sum(e$values))[1:n_pcs], 4) * 100

# Extract PCs and merge with metadata
myPCs <- as.data.frame(e$vectors[, 1:n_pcs])
names(myPCs) <- paste0("PC", 1:n_pcs)
data2 <- cbind(my_pop, myPCs)
data2$Cov <- as.numeric(data2$Cov)

# --- Determine color scheme ---
base_name <- sub("_4x$", "", subset_name)
is_year1 <- (base_name == "year1")

bad_coords <- c("4548", "4549", "4550", "4551", "4552", "4553",
                "4554", "4555", "4556", "4557", "4559", "4561", "4562")
if (!is_year1) {
  data2 <- data2[!data2$SampleID %in% bad_coords,]
  cat("After excluding bad_coords:", nrow(data2), "samples for plotting\n")
}
is_strauss <- (base_name == "strauss")

bangwong <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
              "#0072B2", "#D55E00", "#CC79A7", "#000000",
              "#332288", "#88CCEE")

ohana_colors <- c(
  "1"="#E69F00", "2"="#56B4E9", "10"="#009E73",
  "25"="#F0E442", "36"="#DDDDDD",
  "39"="#0072B2", "47"="#D55E00", "51"="#CC79A7",
  "54"="#000000", "57"="#332288", "66"="#88CCEE",
  "488"="#E69F00", "504"="#56B4E9"
)

if (is_year1) {
  data2$color_var <- as.factor(data2$Clusters)
  color_scale <- scale_fill_manual(values=popcolors, guide="none")
} else {
  data2$color_var <- as.factor(data2$OHANA_order)
  ohana_vals <- as.character(sort(unique(data2$OHANA_order)))
  if (is_strauss) {
    color_map <- clustercolors[ohana_vals[ohana_vals %in% names(clustercolors)]]
    unmapped <- ohana_vals[!ohana_vals %in% names(clustercolors)]
    if (length(unmapped) > 0) {
      avail <- setdiff(bangwong, unname(color_map))
      extra_cols <- setNames(rep_len(avail, length(unmapped)), unmapped)
      color_map <- c(color_map, extra_cols)
    }
  } else {
    color_map <- ohana_colors[ohana_vals]
    unmapped <- ohana_vals[is.na(color_map)]
    if (length(unmapped) > 0) {
      avail <- setdiff(bangwong, unname(color_map[!is.na(color_map)]))
      extra <- setNames(rep_len(avail, length(unmapped)), unmapped)
      color_map[unmapped] <- extra
    }
    names(color_map) <- ohana_vals
  }
  color_scale <- scale_fill_manual(values=color_map, breaks=ohana_vals)
}

# --- Helper to make a PCA plot ---
make_pca_plot <- function(pc_x, pc_y, ve_x, ve_y) {
  ggplot(data2, aes(x=.data[[pc_x]], y=.data[[pc_y]], fill=color_var)) +
    geom_point(colour="black", size=9, shape=21) +
    color_scale +
    labs(x=paste0(pc_x, " (", ve_x, "% variance)"),
         y=paste0(pc_y, " (", ve_y, "% variance)"),
         fill=ifelse(is_year1, "Locality", "OHANA")) +
    gvtp_theme + coord_cartesian(expand = TRUE)
}

# --- PCA colored by locality/OHANA ---
pc_pairs <- list()
for (i in 1:(n_pcs - 1)) {
  pc_x <- paste0("PC", i)
  pc_y <- paste0("PC", i + 1)
  p <- make_pca_plot(pc_x, pc_y, var_exp[i], var_exp[i+1])
  pc_pairs[[length(pc_pairs) + 1]] <- list(plot=p, name=paste0(i, i+1))

  jpeg(file.path(out_dir, paste0(out_name, '_PC', i, i+1, '.jpg')), height=1000, width=1000)
  print(p)
  dev.off()
}

if (length(pc_pairs) >= 2) {
  jpeg(file.path(out_dir, paste0(out_name, '_PC_combined.jpg')), height=1000, width=2000)
  print(pc_pairs[[1]]$plot | pc_pairs[[length(pc_pairs)]]$plot)
  dev.off()
}

# --- PCA colored by coverage ---
make_cov_plot <- function(pc_x, pc_y, ve_x, ve_y) {
  ggplot(data2, aes(x=.data[[pc_x]], y=.data[[pc_y]], fill=Cov)) +
    geom_point(colour="black", size=2, shape=21) +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=mean(data2$Cov),
                         guide=guide_colorbar(title="Coverage")) +
    labs(x=paste0(pc_x, " (", ve_x, "% variance)"),
         y=paste0(pc_y, " (", ve_y, "% variance)")) +
    theme_bw(base_size=12)
}

for (i in 1:(n_pcs - 1)) {
  pc_x <- paste0("PC", i)
  pc_y <- paste0("PC", i + 1)
  jpeg(file.path(out_dir, paste0(out_name, '_cov_PC', i, i+1, '.jpg')), units="in", height=5, width=6, res=300)
  print(make_cov_plot(pc_x, pc_y, var_exp[i], var_exp[i+1]))
  dev.off()
}

# --- Map plots with PC gradients ---
x_range <- diff(range(data2$Long, na.rm=TRUE))
y_range <- diff(range(data2$Lat, na.rm=TRUE))
buf <- max(x_range, y_range) * 0.1 + 0.005
jit_w <- x_range * 0.02 + 0.001
jit_h <- y_range * 0.02 + 0.001
xlims <- range(data2$Long, na.rm=TRUE) + c(-buf, buf)
ylims <- range(data2$Lat, na.rm=TRUE) + c(-buf, buf)

for (pc in paste0("PC", 1:n_pcs)) {
  jpeg(file.path(out_dir, paste0(out_name, '_map_', pc, '.jpg')), units="in", height=5, width=6, res=300)
  print(ggplot() +
    geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", size = 0.5) +
    geom_jitter(data=data2, aes(x=Long, y=Lat, fill=.data[[pc]]), colour="black", size=1.5, shape=21, width=jit_w, height=jit_h) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         guide = guide_colorbar(title=pc)) +
    coord_sf(xlim=xlims, ylim=ylims) +
    theme_bw(base_size=12))
  dev.off()
}

# --- Isolation by distance: geographic distance vs PC distance ---
bad_coords <- c("4548", "4549", "4550", "4551", "4552", "4553",
                "4554", "4555", "4556", "4557", "4559", "4561", "4562")
ibd_data <- data2[!data2$SampleID %in% bad_coords,]

if (nrow(ibd_data) >= 3) {
  deg2rad <- function(d) d * pi / 180
  haversine <- function(lat1, lon1, lat2, lon2) {
    R <- 6371
    dlat <- deg2rad(lat2 - lat1)
    dlon <- deg2rad(lon2 - lon1)
    a <- sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2)^2
    R * 2 * atan2(sqrt(a), sqrt(1-a))
  }

  pc_cols <- paste0("PC", 1:n_pcs)
  n <- nrow(ibd_data)
  pairs <- expand.grid(i = 1:n, j = 1:n)
  pairs <- pairs[pairs$i < pairs$j,]

  geo_dist <- mapply(function(i, j) {
    haversine(ibd_data$Lat[i], ibd_data$Long[i], ibd_data$Lat[j], ibd_data$Long[j])
  }, pairs$i, pairs$j)

  pc_dist <- mapply(function(i, j) {
    sqrt(sum((ibd_data[i, pc_cols] - ibd_data[j, pc_cols])^2))
  }, pairs$i, pairs$j)

  ibd_df <- data.frame(geo_dist = geo_dist, pc_dist = pc_dist)

  jpeg(file.path(out_dir, paste0(out_name, '_IBD.jpg')), units="in", height=5, width=6, res=300)
  print(ggplot(ibd_df, aes(x = geo_dist, y = pc_dist)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(x = "Geographic distance (km)", y = "PC distance",
         title = paste(out_name, "- Isolation by Distance")) +
    theme_bw(base_size = 12))
  dev.off()
}

cat("Done. Plots written to", out_dir, "\n")
