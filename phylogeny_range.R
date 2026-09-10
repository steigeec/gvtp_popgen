#!/usr/bin/env Rscript
#
# phylogeny_range.R
#
# Two-panel figure for Deinandra increscens subsp. villosa (Gaviota tarplant):
#   left  - phylogeny of D. increscens and its close congeners, redrawn (topology
#           only, cladogram) from Simpson, M.G. (2010) Plant Systematics, 2nd ed.,
#           Fig. 19.9 (Academic Press/Elsevier, ISBN 978-0-12-374380-0), which in
#           turn draws on Baldwin's molecular (ITS) phylogenetic work on Deinandra
#           (e.g. Baldwin, B.G. 2007. Adaptive radiation of shrubby tarweeds
#           (Deinandra) in the California Islands parallels diversification of the
#           Hawaiian silversword alliance. American Journal of Botany 94:237-248).
#           Branch lengths in the source figure are parsimony steps (scale bar =
#           5 changes); this script plots the topology as a cladogram rather than
#           attempting to re-derive those lengths from a low-resolution scan, so
#           branch lengths here are NOT to scale.
#   right - range map: GBIF occurrence records for Deinandra increscens subsp.
#           villosa and any other Deinandra taxa recorded within the wider map
#           extent, zoomed out from (but still centered on) the study bounding
#           box used throughout this project (bbox_extent <- ext(-120.7, -120.1,
#           34.4, 34.9); see e.g. pheno/fire.R).
#
# Usage:
#   Rscript phylogeny_range.R <tree.nwk> <occurrences.csv> <villosa_pheno.csv> <jepson.csv> <land.shp> <output.png>
#
# occurrences.csv columns (written by phylogeny-range.sh from the GBIF API):
#   species, infraspecificEpithet, decimalLatitude, decimalLongitude,
#   coordinateUncertaintyInMeters, institutionCode
#
# D. increscens subsp. villosa is federally endangered, and essentially every
# iNaturalist record of it has been coordinate-obscured (randomized within a
# ~29 km radius) at iNaturalist's own initiative to protect the taxon -- see
# phylogeny-range.sh header. Rather than salvage what's left of GBIF's villosa
# records after dropping those, villosa's points come entirely from this
# project's own field/herbarium dataset instead (updatedvarsandpheno.csv,
# Lat._dd.dddd/Lon_dd.dddd columns -- real, unobscured coordinates for known
# populations). GBIF is still the source for every other taxon on the map;
# records with coordinateUncertaintyInMeters above UNCERTAINTY_MAX_M are
# dropped from that GBIF set so the map doesn't show precise-looking dots for
# locations that were deliberately randomized.
#
# jepson.csv is a Darwin Core Archive occurrences.csv exported by hand from
# the CCH2 portal (UC Berkeley's Jepson Herbarium, UCJEPS -- the one major
# CCH2-network collection GBIF doesn't already mirror). Its villosa records
# are dropped in favor of the project's own data for the same reason GBIF's
# are; everything else goes through the same uncertainty/water filters as the
# GBIF set.
#
# land.shp is Natural Earth's 10m land polygon (fetched/cached by
# phylogeny-range.sh) -- used both as the basemap and to drop occurrence
# records whose coordinates land in open water (a real, if modest, share of
# GBIF herbarium records have imprecise/erroneous coordinates; this is a
# simple, transparent filter for that, not a claim that every remaining point
# is error-free).

.libPaths("/private/groups/meyerlab/esteiger/software/R_libs/")

suppressMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript phylogeny_range.R <tree.nwk> <occurrences.csv> <villosa_pheno.csv> <jepson.csv> <land.shp> <output.png>")
}
tree_file   <- args[1]
occ_file    <- args[2]
villosa_csv <- args[3]
jepson_csv  <- args[4]
land_shp    <- args[5]
out_png     <- args[6]

# Above this, coordinates are treated as deliberately obscured (privacy
# protection) rather than merely imprecise, and dropped -- see header.
UNCERTAINTY_MAX_M <- 1000

# Study area (matches pheno/fire.R)
bbox_extent <- c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)

# Wider context window the range map is drawn against (south-central California
# coast: Santa Barbara Co. plus adjoining San Luis Obispo / Ventura Co. margins)
map_extent <- c(xmin = -120.8, xmax = -119.9, ymin = 34.38, ymax = 34.95)

FOCAL_TAXON <- "Deinandra increscens subsp. villosa"

gvtp_theme <- theme_bw(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# ============================================================
# Left panel: phylogeny
# ============================================================

tree <- read.tree(tree_file)
tree$tip.label <- gsub("_", " ", tree$tip.label)

is_focal <- tree$tip.label == FOCAL_TAXON
tip_df <- data.frame(
  label = tree$tip.label,
  focal = is_focal
)

p_tree <- ggtree(tree, branch.length = "none", size = 0.8) %<+% tip_df +
  geom_tiplab(aes(fontface = ifelse(focal, "bold.italic", "italic")),
              size = 4.8, offset = 0.15, color = "black") +
  geom_tippoint(size = 3, color = "black") +
  coord_cartesian(clip = "off") +
  xlim(0, 14) +
  theme(plot.margin = margin(5, 40, 5, 5))

# ============================================================
# Right panel: range map
# ============================================================

occ <- read.csv(occ_file, stringsAsFactors = FALSE)
occ <- occ %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude),
         !is.na(species), species != "") %>%
  mutate(
    taxon = ifelse(
      !is.na(infraspecificEpithet) & infraspecificEpithet != "",
      paste(species, "subsp.", infraspecificEpithet),
      species
    )
  ) %>%
  filter(
    decimalLongitude >= map_extent["xmin"], decimalLongitude <= map_extent["xmax"],
    decimalLatitude  >= map_extent["ymin"], decimalLatitude  <= map_extent["ymax"]
  )

# D. increscens records with no subspecies given are genuinely ambiguous
# between increscens and villosa in general -- but villosa is a narrow
# endemic confined to the study bbox (that's the whole reason bbox_extent is
# what it is), so a bare "Deinandra increscens" record OUTSIDE the bbox can't
# be villosa and is safely relabeled ssp. increscens. Only bare records
# INSIDE the bbox stay genuinely ambiguous and are dropped.
bare_increscens <- occ$species == "Deinandra increscens" & !grepl("subsp\\.", occ$taxon)
outside_bbox <- !(
  occ$decimalLongitude >= bbox_extent["xmin"] & occ$decimalLongitude <= bbox_extent["xmax"] &
  occ$decimalLatitude  >= bbox_extent["ymin"] & occ$decimalLatitude  <= bbox_extent["ymax"]
)
n_resolved <- sum(bare_increscens & outside_bbox)
occ$taxon[bare_increscens & outside_bbox] <- "Deinandra increscens subsp. increscens"
cat("Resolved", n_resolved, "bare Deinandra increscens records outside the study",
    "bbox to subsp. increscens (villosa can't occur there).\n")

n_dropped <- sum(occ$species == "Deinandra increscens" & !grepl("subsp\\.", occ$taxon))
occ <- occ %>%
  filter(species != "Deinandra increscens" | grepl("subsp\\.", taxon))
cat("Dropped", n_dropped, "remaining bare Deinandra increscens records inside",
    "the bbox (still genuinely ambiguous increscens vs. villosa).\n")

# --- Replace GBIF's villosa records with this project's own field data ------
# See header: villosa's GBIF records are almost all privacy-obscured. Drop
# them here and splice in the project's own precise, unobscured coordinates
# for this taxon instead (loaded below).
n_gbif_villosa <- sum(occ$taxon == FOCAL_TAXON)
occ <- occ %>% filter(taxon != FOCAL_TAXON)
cat("Dropped", n_gbif_villosa, "GBIF", FOCAL_TAXON,
    "records (replaced with project field/herbarium data below).\n")

# --- Drop coordinate-obscured records (privacy protection, not noise) --------
# See header: iNaturalist randomizes coordinates within ~29 km for threatened
# taxa. This mainly matters for taxa other than villosa now, but is kept as a
# general check in case any other Deinandra species on the map is similarly
# flagged.
n_obscured <- sum(occ$coordinateUncertaintyInMeters > UNCERTAINTY_MAX_M, na.rm = TRUE)
cat("\nDropping", n_obscured, "of", nrow(occ),
    "records with coordinate uncertainty >", UNCERTAINTY_MAX_M, "m (privacy-obscured), by taxon:\n")
print(table(occ$taxon[which(occ$coordinateUncertaintyInMeters > UNCERTAINTY_MAX_M)]))
occ <- occ %>%
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= UNCERTAINTY_MAX_M)

# --- Add the Jepson Herbarium (CCH2) export ----------------------------------
# See header: the one major CCH2-network collection not already mirrored into
# GBIF. Same taxon/uncertainty/villosa-exclusion treatment as the GBIF set.
jepson <- read.csv(jepson_csv, stringsAsFactors = FALSE) %>%
  filter(genus == "Deinandra", !is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  transmute(
    species = paste(genus, specificEpithet),
    infraspecificEpithet,
    decimalLatitude, decimalLongitude,
    coordinateUncertaintyInMeters,
    institutionCode,
    taxon = ifelse(
      !is.na(infraspecificEpithet) & infraspecificEpithet != "",
      paste(species, "subsp.", infraspecificEpithet),
      species
    )
  ) %>%
  filter(
    decimalLongitude >= map_extent["xmin"], decimalLongitude <= map_extent["xmax"],
    decimalLatitude  >= map_extent["ymin"], decimalLatitude  <= map_extent["ymax"],
    is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= UNCERTAINTY_MAX_M
  )
# Same bbox-based increscens/villosa resolution as the GBIF set above.
jepson_bare <- jepson$species == "Deinandra increscens" & !grepl("subsp\\.", jepson$taxon)
jepson_outside_bbox <- !(
  jepson$decimalLongitude >= bbox_extent["xmin"] & jepson$decimalLongitude <= bbox_extent["xmax"] &
  jepson$decimalLatitude  >= bbox_extent["ymin"] & jepson$decimalLatitude  <= bbox_extent["ymax"]
)
jepson$taxon[jepson_bare & jepson_outside_bbox] <- "Deinandra increscens subsp. increscens"
jepson <- jepson %>% filter(taxon != FOCAL_TAXON, species != "Deinandra increscens" | grepl("subsp\\.", taxon))
cat("\nAdding", nrow(jepson), "Jepson Herbarium (CCH2) records from", jepson_csv, ", by taxon:\n")
print(table(jepson$taxon))
occ <- bind_rows(occ, jepson)

# --- Splice in villosa's real coordinates from the project's own data -------
villosa_occ <- read.csv(villosa_csv, stringsAsFactors = FALSE) %>%
  filter(!is.na(Lat._dd.dddd), !is.na(Lon_dd.dddd)) %>%
  transmute(
    species = "Deinandra increscens",
    infraspecificEpithet = "villosa",
    decimalLatitude = Lat._dd.dddd,
    decimalLongitude = Lon_dd.dddd,
    coordinateUncertaintyInMeters = coordinateUncertaintyInMeters,
    institutionCode = "GVTP project field/herbarium data",
    taxon = FOCAL_TAXON
  ) %>%
  filter(
    decimalLongitude >= map_extent["xmin"], decimalLongitude <= map_extent["xmax"],
    decimalLatitude  >= map_extent["ymin"], decimalLatitude  <= map_extent["ymax"]
  )
cat("\nAdding", nrow(villosa_occ), FOCAL_TAXON, "records from", villosa_csv, "\n")
occ <- bind_rows(occ, villosa_occ)

# --- Drop occurrence records that fall in open water -------------------------
# GBIF's own hasGeospatialIssue flag doesn't catch these (they're not wrong
# enough to trip it -- just imprecise: rounded/centroid coordinates that miss
# the coast). A point-in-land-polygon test against Natural Earth's 10m
# coastline is a simple, transparent way to catch them instead.
#
# Crop the land polygon to the local region BEFORE unioning/buffering it --
# st_union()/st_buffer() on the full global polygon silently produces a
# corrupted (much too large) result here, which would let offshore points
# through as false "on land" hits.
land <- st_read(land_shp, quiet = TRUE)
land_crop <- suppressWarnings(st_crop(
  land, xmin = unname(map_extent["xmin"]) - 0.5, xmax = unname(map_extent["xmax"]) + 0.5,
  ymin = unname(map_extent["ymin"]) - 0.5, ymax = unname(map_extent["ymax"]) + 0.5
))
land_buffered <- st_buffer(st_union(land_crop), dist = 0.001)  # ~100 m, avoids
  # false negatives for real coastal-bluff records sitting right at the coast

occ_sf <- st_as_sf(occ, coords = c("decimalLongitude", "decimalLatitude"),
                    crs = 4326, remove = FALSE)
on_land <- lengths(st_intersects(occ_sf, land_buffered)) > 0
cat("Dropping", sum(!on_land), "of", nrow(occ),
    "occurrence records whose coordinates fall in open water.\n")
occ <- occ[on_land, ]

cat("\nDeinandra occurrence records plotted, by taxon:\n")
print(table(occ$taxon))

in_bbox <- occ %>%
  filter(
    decimalLongitude >= bbox_extent["xmin"], decimalLongitude <= bbox_extent["xmax"],
    decimalLatitude  >= bbox_extent["ymin"], decimalLatitude  <= bbox_extent["ymax"]
  )
cat("\nTaxa recorded within the study bbox itself:\n")
print(sort(unique(in_bbox$taxon)))

taxon_colors <- c(
  "Deinandra increscens subsp. villosa"    = "#DF536B",
  "Deinandra increscens subsp. increscens" = "black",
  "Deinandra fasciculata"                  = "gray62",
  "Deinandra paniculata"                   = "#61D04F",
  "Deinandra simplex"                      = "#28E2E5",
  "Deinandra pentactis"                    = "#F5C710"
)

# Factor with the two increscens taxa listed first so a 2-column, column-major
# legend (guide_legend's default fill order) puts them together in the left
# column, rather than alphabetical order splitting them across columns.
occ$taxon <- factor(occ$taxon, levels = names(taxon_colors))

p_map <- ggplot() +
  geom_sf(data = land_crop, fill = "grey96", color = "grey60", linewidth = 0.3) +
  geom_rect(aes(xmin = bbox_extent["xmin"], xmax = bbox_extent["xmax"],
                ymin = bbox_extent["ymin"], ymax = bbox_extent["ymax"]),
            fill = NA, color = "black", linetype = "dashed", linewidth = 0.6) +
  geom_point(data = occ, aes(decimalLongitude, decimalLatitude, color = taxon),
             shape = 16, size = 2.6, alpha = 0.8) +
  scale_color_manual(values = taxon_colors, name = NULL, drop = TRUE) +
  guides(color = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 3.5))) +
  coord_sf(xlim = c(map_extent["xmin"], map_extent["xmax"]),
           ylim = c(map_extent["ymin"], map_extent["ymax"]), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  gvtp_theme +
  theme(axis.text = element_text(size = 12, angle = 45, hjust = 1),
        legend.text = element_text(face = "italic", size = 15),
        legend.position = "bottom",
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.15, "cm"))

# ============================================================
# Combine + save
# ============================================================

combined <- p_tree + p_map + plot_layout(widths = c(1, 1.1))

ggsave(out_png, combined, width = 15, height = 7.5, dpi = 300, bg = "white")
cat("\nWrote", out_png, "\n")
