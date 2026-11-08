#!/usr/bin/env Rscript
#
# Two-panel figure for Deinandra increscens subsp. villosa:
#   left  - cladogram of D. increscens and close congeners, topology redrawn
#           from Simpson (2010) Plant Systematics, 2nd ed., Fig. 19.9 (after
#           Baldwin 2007, Am. J. Bot. 94:237-248); branch lengths not to scale
#   right - occurrence map of Deinandra taxa around the study area
#
# Usage:
#   Rscript phylogeny_range.R <tree.nwk> <occurrences.csv> <villosa_pheno.csv> <jepson.csv> <land.shp> <output.png>
#
# occurrences.csv (GBIF, from phylogeny-range.sh) columns:
#   species, infraspecificEpithet, decimalLatitude, decimalLongitude,
#   coordinateUncertaintyInMeters, institutionCode
# villosa_pheno.csv: project field/herbarium data, used for all villosa points
#   (GBIF villosa records are almost all coordinate-obscured by iNaturalist)
# jepson.csv: manual CCH2 Darwin Core export for UCJEPS, which GBIF does not mirror
# land.shp: Natural Earth 10m land polygon, used as basemap and to drop
#   records that fall in open water

BASE_DIR <- "/private/groups/meyerlab/esteiger"

.libPaths(file.path(BASE_DIR, "software/R_libs/"))

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

# records above this uncertainty are treated as privacy-obscured and dropped
UNCERTAINTY_MAX_M <- 1000

# Study area
bbox_extent <- c(xmin = -120.7, xmax = -120.1, ymin = 34.4, ymax = 34.9)

# Map extent
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

# villosa is confined to the study bbox, so D. increscens records without a
# subspecies are subsp. increscens outside the bbox and ambiguous (dropped) inside it
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

# --- Drop GBIF villosa records (replaced by project data below) ---------------
n_gbif_villosa <- sum(occ$taxon == FOCAL_TAXON)
occ <- occ %>% filter(taxon != FOCAL_TAXON)
cat("Dropped", n_gbif_villosa, "GBIF", FOCAL_TAXON,
    "records (replaced with project field/herbarium data below).\n")

# --- Drop coordinate-obscured records -----------------------------------------
n_obscured <- sum(occ$coordinateUncertaintyInMeters > UNCERTAINTY_MAX_M, na.rm = TRUE)
cat("\nDropping", n_obscured, "of", nrow(occ),
    "records with coordinate uncertainty >", UNCERTAINTY_MAX_M, "m (privacy-obscured), by taxon:\n")
print(table(occ$taxon[which(occ$coordinateUncertaintyInMeters > UNCERTAINTY_MAX_M)]))
occ <- occ %>%
  filter(is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters <= UNCERTAINTY_MAX_M)

# --- Add the Jepson Herbarium (CCH2) export ----------------------------------
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
# Crop before union/buffer: on the global polygon st_union()/st_buffer()
# returns a corrupted, oversized result
land <- st_read(land_shp, quiet = TRUE)
land_crop <- suppressWarnings(st_crop(
  land, xmin = unname(map_extent["xmin"]) - 0.5, xmax = unname(map_extent["xmax"]) + 0.5,
  ymin = unname(map_extent["ymin"]) - 0.5, ymax = unname(map_extent["ymax"]) + 0.5
))
land_buffered <- st_buffer(st_union(land_crop), dist = 0.001)  # ~100 m, keeps coastal-bluff records

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

# increscens taxa first so they share the left legend column
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
