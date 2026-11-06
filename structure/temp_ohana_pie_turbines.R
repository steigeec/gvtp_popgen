#!/usr/bin/env Rscript

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")
library(ggplot2)
library(tidyr)
library(scatterpie)
library(dplyr)
library(sf)
library(patchwork)
library(terra)
library(ggspatial)

LOCALITY<-"strauss"

#Plotting params
gvtp_theme <- theme_bw() +
  theme(
    text = element_text(size = 50), # Sets size for all text elements
    axis.text = element_text(angle = 45), # Only specify unique modifications
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    legend.box.background = element_rect(fill = "white", color = "black"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank())


proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Crop the coast off
crop_extent <- st_bbox(c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9), crs = st_crs(proj))
CA<-st_read("/private/groups/meyerlab/esteiger/env/study_area/GVTP_clip.shp")
CA<-st_make_valid(CA)
CA <- st_crop(CA, crop_extent)
CA<- st_simplify(CA)

# Read in the DEM and build contour lines for the pie-chart map
DEM <- rast("/private/groups/meyerlab/esteiger/env/tranquillon_mountain_30m.dem")
crs(DEM) <- "EPSG:26710"
DEM_transformed <- project(DEM, proj)
contour_lines <- as.contour(DEM_transformed, levels = seq(0, 800, 25))
contour_sf <- st_as_sf(contour_lines)


# Calculate Evanno's delta K
strauss_ll <- read.table(paste0(LOCALITY, "/log_likelihoods_across_runs.txt"), header=TRUE, sep="\t")

summary_ll <- strauss_ll %>%
  group_by(K) %>%
  summarise(mean_lnL = mean(loglik), sd_lnL = sd(loglik))

stats <- summary_ll %>%
  mutate(lnL_lead = lead(mean_lnL),
    lnL_lag  = lag(mean_lnL),
    delta_K  = abs(lnL_lead - 2 * mean_lnL + lnL_lag) / sd_lnL)

P_deltaK <- ggplot(stats, aes(x=K, y=delta_K)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="Delta K (Evanno method)", title="Evanno Method for Choosing K") +
  theme_minimal()

jpeg(paste0(LOCALITY, "/deltaK_", LOCALITY, ".jpg"))
print(P_deltaK)
dev.off()


# Always read in the best-likelihood run (per best_runs.txt) for each K,
# rather than hard-coding which iteration or which K values were run.
best_runs <- read.table(paste0(LOCALITY, "/best_runs.txt"), header=FALSE, sep="\t",
                         col.names=c("K", "iteration", "loglik"))

for (k in sort(unique(best_runs$K))) {
  best_iter <- best_runs$iteration[best_runs$K == k]
  q_file <- sprintf("%s/q_angsd_paralog_LD_%s_%d_%d.matrix", LOCALITY, LOCALITY, k, best_iter)
  q <- read.table(q_file, skip=1, sep="\t", row.names=NULL)
  names(q) <- paste0("C", seq_len(k))
  cat("Read K=", k, " from iteration ", best_iter, "\n")
  assign(paste0("q_", k), q)
}


my_order<-read.table(paste0("/private/groups/meyerlab/esteiger/mapping/merged/bamlist_", LOCALITY, ".txt"), header=F)
my_order$index <- seq(1, nrow(my_order))
my_order$V1 <- sub(".*/([^/]+)$", "\\1", my_order$V1)
my_order$V1 <- gsub("_rmdup.bam", "", my_order$V1)


#What's the population info for these samples?
pop<-read.table("/private/groups/meyerlab/esteiger/field_data.txt", header=T, sep="\t")
pop<-unique(pop)
pop<-merge(pop[,c("SampleID", "Year", "Lat", "Long", "Cov", "OHANA_order")], my_order, by.x="SampleID", by.y="V1", all.y=F)
my_pop<-pop[order(pop$index),]


# Deterministic ancestry-component coloring, persistent across K: the first
# plot (K=2) gets two base placeholder colors (blue, sky) by component size,
# since nothing has genuinely differentiated yet. From K=3 on, a component's
# color is tied to a *specific identity*, not to emergence order alone:
#   - "diffuse" components (no individual has membership > 0.9 in it -- i.e.
#     nobody really "belongs" to it) are always green, a single fixed
#     identity reused every time one appears at any K. This is important:
#     a component can have a locality as its plurality by argmax count and
#     still be genuinely diffuse (e.g. everyone sitting around ~30%), which
#     is not the same thing as that locality having actually broken off.
#   - Otherwise, the component's identity is whichever locality dominates
#     its high-confidence (>0.9) members; a locality's color, once assigned
#     on first appearance, is reused every time that locality resurfaces at
#     a higher K (including if it later splits into multiple components).
# This guarantees a locality's color never appears before that locality has
# genuinely broken off (real, high-confidence members), and never changes
# identity once assigned. Each locality identity maps to a fixed color;
# any real identity not in this dictionary (unanticipated) falls back to
# a small extra palette, assigned on first appearance.
identity_colors <- c(
  "L72" = "#EE8866",   # orange
  "L73" = "#99DDFF",   # sky
  "L75" = "#EEDD88",   # yellow
  "L80" = "#FFAABB",   # pink
  "L88" = "#77AADD",   # blue
  "GENERIC" = "#AAAA00" # olive/green
)
fallback_palette <- c("white", "#44BB99", "black", "#BBCC33")
STRONG_THRESHOLD <- 0.9

build_persistent_colors <- function(all_k) {
  color_state <- list()
  next_fallback_idx <- 1
  identity_color <- as.list(identity_colors)

  for (k in all_k) {
    d <- get(paste0("data_", k))
    comp_cols <- paste0("C", seq_len(k))
    mat <- as.matrix(d[, comp_cols])
    storage.mode(mat) <- "numeric"
    assignment <- apply(mat, 1, which.max)
    comp_sizes <- vapply(seq_len(k), function(c) sum(assignment == c), numeric(1))
    dominant_locality <- vapply(seq_len(k), function(c) {
      idx <- which(assignment == c)
      if (length(idx) == 0) return(NA_character_)
      tab <- table(as.character(d$OHANA_order[idx]))
      names(tab)[which.max(tab)]
    }, character(1))

    colors <- setNames(rep(NA_character_, k), comp_cols)

    if (k == all_k[1]) {
      # First K: nothing has resolved yet -- just a size-ranked base pair,
      # using the colors that L88 (by far the largest locality) and L73
      # will end up permanently owning anyway.
      ord <- order(-comp_sizes)
      colors[ord[1]] <- identity_colors[["L88"]]
      colors[ord[2]] <- identity_colors[["L73"]]
    } else {
      keys <- character(k)
      for (c in seq_len(k)) {
        strong <- which(mat[, c] > STRONG_THRESHOLD)
        if (length(strong) == 0) {
          keys[c] <- "GENERIC"
        } else {
          tab <- sort(table(as.character(d$OHANA_order[strong])), decreasing = TRUE)
          # Ties (or near-ties, common with small strong-member counts) get
          # broken toward whichever candidate is already a known identity,
          # rather than alphabetically -- an alphabetical tiebreak can hand
          # a brand-new label to what is really just noise around an
          # already-established component (e.g. the big persistent bulk
          # component occasionally has as many high-confidence members from
          # a minor locality as from its true dominant one, purely by
          # chance in a given run).
          tied <- names(tab)[tab == tab[1]]
          known <- tied[paste0("L", tied) %in% names(identity_color)]
          loc <- if (length(known) > 0) known[1] else tied[1]
          keys[c] <- paste0("L", loc)
        }
      }

      # When multiple components are simultaneously diffuse, the largest
      # gets the "minor" grey bucket and the smaller one(s) keep olive.
      generic_comps <- which(keys == "GENERIC")
      if (length(generic_comps) > 1) {
        ord <- generic_comps[order(-comp_sizes[generic_comps])]
        keys[ord[1]] <- "MINOR_GENERIC"
      }

      # When the same real locality identity claims multiple simultaneous
      # components (that locality splitting further at this K), the more
      # western one is a distinct sub-group and gets turquoise instead of
      # sharing the locality's usual color.
      real_keys <- unique(keys[!keys %in% c("GENERIC", "MINOR_GENERIC")])
      for (rk in real_keys) {
        comps <- which(keys == rk)
        if (length(comps) > 1) {
          mean_long <- vapply(comps, function(c) mean(d$Long[assignment == c], na.rm = TRUE), numeric(1))
          keys[comps[which.min(mean_long)]] <- "WEST_SPLIT"
        }
      }

      for (c in seq_len(k)) {
        key <- keys[c]
        if (key == "MINOR_GENERIC") {
          colors[c] <- "#DDDDDD"  # grey -- distinct from L80's pink
        } else if (key == "WEST_SPLIT") {
          colors[c] <- "#44BB99"
        } else {
          if (is.null(identity_color[[key]])) {
            identity_color[[key]] <- fallback_palette[next_fallback_idx]
            next_fallback_idx <- next_fallback_idx + 1
          }
          colors[c] <- identity_color[[key]]
        }
      }
    }

    color_state[[as.character(k)]] <- list(colors = colors, dominant_locality = dominant_locality, sizes = comp_sizes)
  }

  color_state
}


# Remove individuals whose coordinates I do not trust.
EXCLUDE<-c("4548", "4549", "4550", "4551", "4552", "4553", "4554", "4555", "4556", "4557", "4559", "4561", "4562")

# Merge pop data together with q matrices, sorted by locality then by
# longitude within each locality ("clumps of plants").
for (k in sort(unique(best_runs$K))) {
  qk <- get(paste0("q_", k))
  dat <- cbind(my_pop, qk)
  dat <- dat[order(dat$OHANA_order, dat$Long), ]
  dat$index <- seq(1, nrow(dat))
  dat <- dat[!(dat$SampleID %in% EXCLUDE),]
  assign(paste0("data_", k), dat)
  assign(paste0("data_long_", k), gather(dat, Cluster, Percentage, all_of(paste0("C", seq_len(k))), factor_key=T))
}


# NMI against OHANA_order sub-locality, per K -- there's no coarser
# "population" label for this within-Strauss dataset, so OHANA_order is
# the only ground truth available.
normalized_mutual_info <- function(known, assignment) {
  tab <- table(known, assignment)
  n <- sum(tab)
  p_xy <- tab / n
  p_x <- rowSums(p_xy)
  p_y <- colSums(p_xy)
  mi <- 0
  for (i in seq_len(nrow(tab))) {
    for (j in seq_len(ncol(tab))) {
      if (p_xy[i, j] > 0) {
        mi <- mi + p_xy[i, j] * log(p_xy[i, j] / (p_x[i] * p_y[j]))
      }
    }
  }
  h_x <- -sum(p_x[p_x > 0] * log(p_x[p_x > 0]))
  h_y <- -sum(p_y[p_y > 0] * log(p_y[p_y > 0]))
  if (h_x == 0 || h_y == 0) return(0)
  unname(mi / sqrt(h_x * h_y))
}

concordance <- data.frame(K = integer(), NMI = numeric())
for (k in sort(unique(best_runs$K))) {
  dat <- get(paste0("data_", k))
  comp_cols <- paste0("C", seq_len(k))
  assignment <- apply(as.matrix(dat[, comp_cols]), 1, which.max)
  known <- as.character(dat$OHANA_order)
  concordance <- rbind(concordance, data.frame(K = k, NMI = normalized_mutual_info(known, assignment)))
}

P_NMI <- ggplot(concordance, aes(x=K, y=NMI)) +
  geom_line() +
  geom_point() +
  labs(x="K", y="NMI with OHANA_order sub-locality", title="Cluster assignment vs. known sub-localities") +
  theme_minimal()

jpeg(paste0(LOCALITY, "/NMI_", LOCALITY, ".jpg"))
print(P_NMI)
dev.off()

best_k <- concordance$K[which.max(concordance$NMI)]
cat("Best K according to NMI:", best_k, "\n")


all_k <- sort(unique(best_runs$K))
color_state <- build_persistent_colors(all_k)

# Build one bar-plot panel per K (ordered by locality "clumps of plants",
# no facet strips -- this dataset is a single population, not several).
mkplot <- function(k) {
  dl <- get(paste0("data_long_", k))
  cols <- color_state[[as.character(k)]]$colors
  # At K=6 specifically, show the pink (L80) component as grey instead --
  # a one-off display choice for this row only; L80 stays pink at every
  # other K so the persistent color-per-identity scheme is otherwise
  # untouched.
  if (k == 6) {
    cols[cols == "#FFAABB"] <- "#DDDDDD"
  }
  ggplot(data=dl, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
    geom_bar(stat="identity", color="black", width=1) +
    scale_fill_manual("", values=cols) +
    scale_y_continuous(name="Fraction composition", expand = c(0,0), breaks=seq(0,1,.2)) +
    gvtp_theme +
    theme(axis.text = element_blank(),
          axis.title=element_blank(),
          legend.position = "none")
}

for (k in all_k) {
  assign(paste0("P", k), mkplot(k))
}

# Show the stack through K=7 only
stack_k <- all_k[all_k <= 7]
jpeg(paste0(LOCALITY, "/", LOCALITY, "_stack.jpg"), width=4000, height=700 * length(stack_k))
print(Reduce(`/`, mget(paste0("P", stack_k))) + plot_layout(heights = rep(1, length(stack_k))))
dev.off()

# Show K=6 on its own (fixed, independent of whichever K NMI selects as best)
display_k <- 6
jpeg(sprintf("%s/ohana_K%d.jpg", LOCALITY, display_k), width=4000, height=1000)
print(get(paste0("P", display_k)) + theme(axis.text=element_text(size=30)))
dev.off()

# Experimental variant (2026-10-22, per user request): same K=6 data, but
# split into one panel per stand (OHANA_order) instead of one continuous
# bar, with very little space between panels. Stands numbered 1-9 in the
# same order as the inbreeding/nucleotide-diversity boxplots (72,73,75,
# 78,80,85,87,89,88); 79 is folded into the same panel as 78 ("4") rather
# than getting its own, per user decision. Written to a separate file
# (not replacing ohana_K6.jpg) since this is an experiment, not a
# replacement for the original.
stand_labels_k6 <- c("72"="1", "73"="2", "75"="3", "78"="4", "79"="4",
                      "80"="5", "85"="6", "87"="7", "89"="8", "88"="9")
dl6 <- get(paste0("data_long_", display_k))
dl6$stand_label <- factor(stand_labels_k6[as.character(dl6$OHANA_order)], levels=as.character(1:9))
cols6 <- color_state[[as.character(display_k)]]$colors
cols6[cols6 == "#FFAABB"] <- "#DDDDDD"

p_ohana_k6_bystand <- ggplot(data=dl6, aes(x=as.factor(index), y=as.numeric(Percentage), fill=Cluster)) +
  geom_bar(stat="identity", color="black", width=1) +
  scale_fill_manual("", values=cols6) +
  scale_y_continuous(name="Ancestry composition", expand=c(0,0), breaks=seq(0,1,.2)) +
  # facet_grid (not facet_wrap) specifically so space="free_x" can size each
  # panel's width by its own number of individuals -- facet_wrap always
  # gives every panel equal width regardless of scales="free_x", which
  # crushed stand 88 (397 of 736 individuals) into a sliver so thin its
  # bars' black outlines blotted out the fill entirely.
  facet_grid(~stand_label, scales="free_x", space="free_x", switch="x") +
  gvtp_theme +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=55, angle=90),
        legend.position="none",
        panel.spacing=unit(4, "pt"),
        panel.border=element_rect(color="black", fill=NA, linewidth=1),
        strip.background=element_rect(fill="grey85", color="black"),
        strip.placement="outside",
        strip.text=element_text(size=55, margin=margin(t=8, b=8)))

jpeg(sprintf("%s/ohana_K%d_bystand.jpg", LOCALITY, display_k), width=2200, height=790)
print(p_ohana_k6_bystand)
dev.off()


# Pie chart on the contour map, for K=6 (fixed; see display_k above)
best_dat <- get(paste0("data_", display_k))
comp_cols <- paste0("C", seq_len(display_k))

x_range <- range(best_dat$Long, na.rm = TRUE)
y_range <- range(best_dat$Lat, na.rm = TRUE)
x_padding <- diff(x_range) * 0.05
y_padding <- diff(y_range) * 0.05
x_extent <- diff(x_range) + 2 * x_padding
y_extent <- diff(y_range) + 2 * y_padding
# coord_sf renders lon/lat with latitude foreshortening (1 degree of
# longitude covers less true distance than 1 degree of latitude away from
# the equator), so the plotting device's aspect ratio needs that same
# cos(latitude) correction -- otherwise the canvas doesn't match what
# coord_sf actually draws and the map ends up with uneven margins.
mean_lat_rad <- mean(y_range) * pi / 180
aspect_ratio <- (x_extent * cos(mean_lat_rad)) / y_extent
base_height <- 1500
plot_width <- base_height * aspect_ratio
plot_height <- base_height
base_scale <- 0.02 * (2/3) * (1.1/1.4) * 0.9 * 0.85 * 1.15 * 1.15  # was 0.20; ~10x too large, then 2/3;
                                          # then padding tightened (0.2->0.05), which
                                          # shrinks x_extent and would otherwise regrow
                                          # the pies; then a tiny bit smaller again;
                                          # then a bit bigger again (*1.15), then bigger again (*1.15)
jitter_amount <- x_extent / 80

best_dat <- best_dat %>%
  mutate(Long_jitter = jitter(Long, amount = jitter_amount),
         Lat_jitter  = jitter(Lat, amount = jitter_amount))

pie_best <- best_dat %>%
  dplyr::group_by(SampleID) %>%
  dplyr::summarise(across(all_of(comp_cols), ~ mean(as.numeric(.x))))
names(pie_best)[names(pie_best) %in% comp_cols] <- LETTERS[seq_len(display_k)]
pie_best <- merge(pie_best, best_dat[, c("SampleID", "Long_jitter", "Lat_jitter")], by = "SampleID")
pie_best <- unique(pie_best)

pie_colors <- color_state[[as.character(display_k)]]$colors
# Same K=6 display override as mkplot(): show L80's pink as grey here too.
if (display_k == 6) {
  pie_colors[pie_colors == "#FFAABB"] <- "#DDDDDD"
}
names(pie_colors) <- LETTERS[seq_len(display_k)]

# TEMP: wind turbine locations overlaid as icons (quick one-off request,
# 2026-10-23 -- this whole script is temporary, delete after use).
turbines <- read.csv("/private/groups/meyerlab/esteiger/env/Turbines.csv", stringsAsFactors = FALSE)

jpeg(sprintf("%s/%s_Map-pie_K%d_turbines_only.jpg", LOCALITY, LOCALITY, display_k), width=plot_width, height=plot_height)
print(
  ggplot2::ggplot() +
    geom_sf(data=contour_sf, color="grey40", linewidth=1.8, alpha=0.5, inherit.aes=FALSE) +
    geom_sf(data = CA, fill = NA, color = "black", linetype = "solid", linewidth = 2) +
    geom_point(data=turbines, aes(x=Lon, y=Lat), shape=8, size=6, stroke=1.8, color="black") +
    coord_sf(crs = 4326,
             xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
             ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 5.5, line_width = 2,
                      height = unit(0.8, "cm"), pad_x = unit(1, "cm"), pad_y = unit(1, "cm")) +
    theme(panel.background = element_blank()) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.ticks.length = unit(0, "pt"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
)
dev.off()
