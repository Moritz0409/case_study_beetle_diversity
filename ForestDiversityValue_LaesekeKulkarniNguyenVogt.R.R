# ============================================================
# MASTER ANALYSIS SCRIPT — BEETLE BIODIVERSITY IN THE LFF
#
# This script combines all three analysis stages into one
# reproducible pipeline:
#
#   PART 1 — Rarefaction & Bootstrapping (beetle species richness)
#   BRIDGE  — Assembly of beetle input table (beetle_data_mean.csv)
#   PART 2 — Tree Biodiversity Index (TBI) & Forest Biodiversity Index (FBI)
#   PART 3 — Sensitivity Analysis (Spearman's rho, Kendall's W)
#
# Inputs required in working directory:
#   - data_2017_species.csv
#   - data_2018_species.csv
#   - species_taxonomy_hahn.csv
#   - Results_combined_census_*.csv  (9 PPA simulation files)
#
# All outputs are written to the working directory.
#
# Fixes relative to original Rarefaction___Bootstrapping.R:
#   1. fun.aggregate = sum added to all dcast calls (duplicate rows)
#   2. Missing *_summary objects added before each bootstrapping _rl block
#   3. Rote.Liste filter corrected (was == "n", should be != "n") in _rl blocks
#   4. unique_2017/2018: Species column preserved via left_join
#   5. count_unique_rl_2trap() added for 2-trap species (Acer, Quercus rubra, Ulmus)
#   6. BRIDGE section assembles beetle_data_mean.csv automatically
# ============================================================


# ============================================================
# 0. PACKAGES & WORKING DIRECTORY
# ============================================================

setwd("/Users/Jonas/Desktop/Master/3rd Semester/Integration Module/BEETLES/test_run")

packages <- c(
  "iNEXT", "ggpubr", "ggplot2", "reshape2",
  "dplyr", "tidyr", "purrr", "tibble", "readr", "forcats"
)

install.packages(setdiff(packages, rownames(installed.packages())))
for (pkg in packages) library(pkg, character.only = TRUE)


# ============================================================
# PART 1 — RAREFACTION & BOOTSTRAPPING
# ============================================================

# ------------------------------------------------------------
# 1.1  Load & merge raw data
# ------------------------------------------------------------

raw_2017              <- read.csv("data_2017_species.csv")
data_2018_species     <- read.csv("data_2018_species.csv")
species_taxonomy_hahn <- read.csv("species_taxonomy_hahn.csv")

taxonomy_2017 <- merge(
  x = raw_2017,
  y = species_taxonomy_hahn[, c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")],
  by.x = "FHL", by.y = "Colkat"
)

taxonomy_2018 <- merge(
  x = data_2018_species,
  y = species_taxonomy_hahn[, c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")],
  by.x = "FHL", by.y = "Colkat"
)
taxonomy_2018 <- taxonomy_2018 %>%
  mutate(Species = paste(Family, Genus, Art, sep = " "))

beetles_2017_red          <- taxonomy_2017[which(taxonomy_2017$Rote.Liste != "n" & taxonomy_2017$Number != 0), ]
beetles_2017_urwaldrelikt <- taxonomy_2017[which(taxonomy_2017$Urwaldrelikt != "n"), ]


# ------------------------------------------------------------
# 1.2  Helper functions
# ------------------------------------------------------------

# Build presence/absence incidence matrix
# Fix: fun.aggregate = sum resolves duplicate Species/Trap rows
make_inc <- function(df, row_var, trap_var = "Trap", value_var = "Number") {
  mat <- reshape2::dcast(
    df,
    as.formula(paste(row_var, "~", trap_var)),
    value.var     = value_var,
    fill          = 0,
    fun.aggregate = sum
  )
  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat[mat >= 1] <- 1
  as.matrix(mat)
}

# For 2-trap species: count unique red list species directly
count_unique_rl_2trap <- function(taxonomy, sp) {
  unique_sp <- taxonomy %>%
    filter(TreeSp != "Ground") %>%
    group_by(FHL) %>%
    mutate(n_treesp = n_distinct(TreeSp)) %>%
    filter(n_treesp == 1, TreeSp == sp) %>%
    pull(FHL) %>%
    unique()
  taxonomy %>%
    filter(TreeSp == sp, FHL %in% unique_sp, Rote.Liste != "n") %>%
    summarise(n = n_distinct(FHL)) %>%
    pull(n)
}


# ------------------------------------------------------------
# 1.3  Species richness 2017 — total
# ------------------------------------------------------------

sp_per_trap_Rubra_total    <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Quercus rubra", ], "Species")
sp_per_trap_Fraxinus_total <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Fraxinus", ],      "Species")
sp_per_trap_Quercus_total  <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Quercus", ],       "Species")
sp_per_trap_Acer_total     <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Acer", ],          "Species")
sp_per_trap_Tilia_total    <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Tilia", ],         "Species")
sp_per_trap_Ulmus_total    <- make_inc(taxonomy_2017[taxonomy_2017$TreeSp == "Ulmus", ],         "Species")

Forest_treesp_total_2017 <- list(
  "Fraxinus"      = sp_per_trap_Fraxinus_total,
  "Quercus"       = sp_per_trap_Quercus_total,
  "Tilia"         = sp_per_trap_Tilia_total,
  "Quercus rubra" = sp_per_trap_Rubra_total,
  "Acer"          = sp_per_trap_Acer_total,
  "Ulmus"         = sp_per_trap_Ulmus_total
)

tree_total_2017 <- iNEXT(Forest_treesp_total_2017, q = 1,
                         datatype = "incidence_raw", endpoint = 2)

plot5  <- ggiNEXT(tree_total_2017, type = 1) +
  ggtitle("Total number species with rarefaction = 2 traps (2017)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot6  <- ggiNEXT(tree_total_2017, type = 2) +
  ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot13 <- ggiNEXT(tree_total_2017, type = 3) +
  geom_vline(xintercept = 0.65, linetype = "dashed") +
  ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
ggarrange(plot5, plot6, plot13, common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1)
ggsave("Tree_Plot_total_2017.png", width = 20, height = 10)

total2017_sc <- bind_rows(tree_total_2017$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD     = approx(SC, qD,      xout = 0.60)$y,
    qD.LCL = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL = approx(SC, qD.UCL, xout = 0.60)$y,
    SC     = 0.60
  ) %>%
  ungroup()

total2017_trap <- bind_rows(tree_total_2017$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.4  Species richness 2017 — red list
# ------------------------------------------------------------

sp_per_trap_Rubra_red    <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Quercus rubra", ], "Species")
sp_per_trap_Fraxinus_red <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Fraxinus", ],      "Species")
sp_per_trap_Quercus_red  <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Quercus", ],       "Species")
sp_per_trap_Acer_red     <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Acer", ],          "Species")
sp_per_trap_Tilia_red    <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Tilia", ],         "Species")
sp_per_trap_Ulmus_red    <- make_inc(beetles_2017_red[beetles_2017_red$TreeSp == "Ulmus", ],         "Species")

Forest_treesp_red_2017 <- list(
  "Fraxinus"      = sp_per_trap_Fraxinus_red,
  "Quercus"       = sp_per_trap_Quercus_red,
  "Tilia"         = sp_per_trap_Tilia_red,
  "Quercus rubra" = sp_per_trap_Rubra_red,
  "Acer"          = sp_per_trap_Acer_red,
  "Ulmus"         = sp_per_trap_Ulmus_red
)

tree_red_2017 <- iNEXT(Forest_treesp_red_2017, q = 1,
                       datatype = "incidence_raw", endpoint = 2)

plot3  <- ggiNEXT(tree_red_2017, type = 1) +
  ggtitle("Total number red-list species with rarefaction = 2 traps (2017)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot4  <- ggiNEXT(tree_red_2017, type = 2) +
  ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot21 <- ggiNEXT(tree_red_2017, type = 3) +
  ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
ggarrange(plot3, plot4, plot21, common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1)
ggsave("Tree_Plot_red_2017.png", width = 20, height = 10)

redlist2017_sc <- bind_rows(tree_red_2017$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD     = approx(SC, qD,      xout = 0.50)$y,
    qD.LCL = approx(SC, qD.LCL, xout = 0.50)$y,
    qD.UCL = approx(SC, qD.UCL, xout = 0.50)$y,
    SC     = 0.50
  ) %>%
  ungroup()

redlist2017_trap <- bind_rows(tree_red_2017$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.5  Unique species 2017
# ------------------------------------------------------------

beetles_2017_treesp <- taxonomy_2017[taxonomy_2017$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number), .groups = "drop")

distinct_2017 <- beetles_2017_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n())

# Fix: keep Species and Number via left_join (original script lost these in summarise)
unique_2017 <- distinct_2017 %>%
  filter(n_appearance == 1) %>%
  left_join(
    taxonomy_2017 %>% select(TreeSp, FHL, Trap, Number, Species),
    by = c("TreeSp", "FHL")
  )

unique_per_trap_Fraxinus_2017 <- make_inc(unique_2017[unique_2017$TreeSp == "Fraxinus", ], "Species")
unique_per_trap_Quercus_2017  <- make_inc(unique_2017[unique_2017$TreeSp == "Quercus", ],  "Species")
unique_per_trap_Tilia_2017    <- make_inc(unique_2017[unique_2017$TreeSp == "Tilia", ],    "Species")

unique_treesp_2017 <- list(
  "Fraxinus" = unique_per_trap_Fraxinus_2017,
  "Quercus"  = unique_per_trap_Quercus_2017,
  "Tilia"    = unique_per_trap_Tilia_2017
)

unique_inext_2017 <- iNEXT(unique_treesp_2017, q = 1,
                           datatype = "incidence_raw", endpoint = 2)

unique_inext_2017_out <- bind_rows(unique_inext_2017$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.6  Species richness 2018 — total
# ------------------------------------------------------------

sp_per_trap_Rubra_total_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Quercus rubra", ], "Species")
sp_per_trap_Fraxinus_total_2018 <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Fraxinus", ],      "Species")
sp_per_trap_Quercus_total_2018  <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Quercus", ],       "Species")
sp_per_trap_Acer_total_2018     <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Acer", ],          "Species")
sp_per_trap_Tilia_total_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Tilia", ],         "Species")
sp_per_trap_Ulmus_total_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Ulmus", ],         "Species")

Forest_treesp_total_2018 <- list(
  "Fraxinus"      = sp_per_trap_Fraxinus_total_2018,
  "Quercus"       = sp_per_trap_Quercus_total_2018,
  "Tilia"         = sp_per_trap_Tilia_total_2018,
  "Quercus rubra" = sp_per_trap_Rubra_total_2018,
  "Acer"          = sp_per_trap_Acer_total_2018,
  "Ulmus"         = sp_per_trap_Ulmus_total_2018
)

tree_total_2018 <- iNEXT(Forest_treesp_total_2018, q = 1,
                         datatype = "incidence_raw", endpoint = 2)

plot7  <- ggiNEXT(tree_total_2018, type = 1) +
  ggtitle("Total number species with rarefaction = 2 traps (2018)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot8  <- ggiNEXT(tree_total_2018, type = 2) +
  ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot14 <- ggiNEXT(tree_total_2018, type = 3) +
  geom_vline(xintercept = 0.65, linetype = "dashed") +
  ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
ggarrange(plot7, plot8, plot14, common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1)
ggsave("Tree_Plot_total_2018.png", width = 20, height = 10)

total2018_sc <- bind_rows(tree_total_2018$iNextEst) %>%
  group_by(Assemblage) %>%
  arrange(SC) %>%
  summarise(
    qD     = approx(SC, qD,      xout = 0.60)$y,
    qD.LCL = approx(SC, qD.LCL, xout = 0.60)$y,
    qD.UCL = approx(SC, qD.UCL, xout = 0.60)$y,
    SC     = 0.60
  ) %>%
  ungroup()

total2018_trap <- bind_rows(tree_total_2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.7  Species richness 2018 — red list
# ------------------------------------------------------------

beetles_2018_red <- taxonomy_2018[taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ]

sp_per_trap_Rubra_red_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Quercus rubra" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")
sp_per_trap_Fraxinus_red_2018 <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Fraxinus"      & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")
sp_per_trap_Quercus_red_2018  <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Quercus"       & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")
sp_per_trap_Acer_red_2018     <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Acer"          & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")
sp_per_trap_Tilia_red_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Tilia"         & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")
sp_per_trap_Ulmus_red_2018    <- make_inc(taxonomy_2018[taxonomy_2018$TreeSp == "Ulmus"         & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0, ], "Species")

Forest_treesp_red_2018 <- list(
  "Fraxinus"      = sp_per_trap_Fraxinus_red_2018,
  "Quercus"       = sp_per_trap_Quercus_red_2018,
  "Tilia"         = sp_per_trap_Tilia_red_2018,
  "Quercus rubra" = sp_per_trap_Rubra_red_2018,
  "Acer"          = sp_per_trap_Acer_red_2018,
  "Ulmus"         = sp_per_trap_Ulmus_red_2018
)

tree_red_2018 <- iNEXT(Forest_treesp_red_2018, q = 1,
                       datatype = "incidence_raw", endpoint = 2)

plot11 <- ggiNEXT(tree_red_2018, type = 1) +
  ggtitle("Total number red-list species with rarefaction = 2 traps (2018)") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot12 <- ggiNEXT(tree_red_2018, type = 2) +
  ggtitle("Coverage estimatation with rarefaction = 2 traps") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
plot23 <- ggiNEXT(tree_red_2018, type = 3) +
  ggtitle("Species estimation based on sample coverage") +
  scale_color_manual(values = c("purple", "blue", "green", "orange", "red", "yellow")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 10))
ggarrange(plot11, plot12, plot23, common.legend = TRUE, legend = "bottom", ncol = 3, nrow = 1)
ggsave("Tree_Plot_red_2018.png", width = 20, height = 10)

redlist2018_sc <- bind_rows(tree_red_2018$iNextEst) %>% ## 2018 not enough data for all tree species
  filter(t == 2, !is.na(SC.LCL))

redlist2018_trap <- bind_rows(tree_red_2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.8  Unique species 2018
# ------------------------------------------------------------

beetles_2018_treesp <- taxonomy_2018[taxonomy_2018$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number), .groups = "drop")

distinct_2018 <- beetles_2018_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n())

unique_2018 <- distinct_2018 %>%
  filter(n_appearance == 1) %>%
  left_join(
    taxonomy_2018 %>% select(TreeSp, FHL, Trap, Number),
    by = c("TreeSp", "FHL")
  )

unique_per_trap_Fraxinus_2018 <- make_inc(unique_2018[unique_2018$TreeSp == "Fraxinus", ], "FHL")
unique_per_trap_Quercus_2018  <- make_inc(unique_2018[unique_2018$TreeSp == "Quercus", ],  "FHL")
unique_per_trap_Tilia_2018    <- make_inc(unique_2018[unique_2018$TreeSp == "Tilia", ],    "FHL")

unique_treesp_2018 <- list(
  "Fraxinus" = unique_per_trap_Fraxinus_2018,
  "Quercus"  = unique_per_trap_Quercus_2018,
  "Tilia"    = unique_per_trap_Tilia_2018
)

unique_inext_2018 <- iNEXT(unique_treesp_2018, q = 1,
                           datatype = "incidence_raw", endpoint = 2)

unique_inext_2018_out <- bind_rows(unique_inext_2018$iNextEst) %>%
  filter(t == 2, !is.na(SC.LCL))


# ------------------------------------------------------------
# 1.9  Bootstrapping — unique red list species (2017 & 2018)
#
# Method: exact replication of original Bootstrapping.R (Ferdinand).
# 25 hardcoded trap combinations, identical for 2017 and 2018.
# Each combination contains 2 Tilia + 2 Quercus + 2 Fraxinus traps
# plus all fixed traps for Acer (27,28), Quercus rubra (21,22),
# Ulmus (29,30). Uniqueness is evaluated across ALL tree species
# simultaneously per combination. Mean across 25 combinations
# is the bootstrapping estimate.
# ------------------------------------------------------------

trap_combinations <- list(
  c( 1,  2, 11, 12, 17, 18, 21, 22, 27, 28, 29, 30),
  c( 1,  2, 13, 14, 23, 24, 21, 22, 27, 28, 29, 30),
  c( 1,  2, 15, 16, 17, 18, 21, 22, 27, 28, 29, 30),
  c( 1,  2, 37, 38, 33, 34, 21, 22, 27, 28, 29, 30),
  c( 3,  4, 11, 12, 25, 26, 21, 22, 27, 28, 29, 30),
  c( 3,  4, 13, 14, 23, 24, 21, 22, 27, 28, 29, 30),
  c( 3,  4, 15, 16, 17, 18, 21, 22, 27, 28, 29, 30),
  c( 3,  4, 37, 38, 19, 20, 21, 22, 27, 28, 29, 30),
  c( 5,  6, 11, 12, 33, 34, 21, 22, 27, 28, 29, 30),
  c( 5,  6, 13, 14, 23, 24, 21, 22, 27, 28, 29, 30),
  c( 5,  6, 15, 16, 25, 26, 21, 22, 27, 28, 29, 30),
  c( 5,  6, 37, 38, 19, 20, 21, 22, 27, 28, 29, 30),
  c( 7,  8,  9, 10, 17, 18, 21, 22, 27, 28, 29, 30),
  c( 7,  8,  9, 10, 25, 26, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 11, 12, 19, 20, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 13, 14, 25, 26, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 13, 14, 33, 34, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 15, 16, 17, 18, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 15, 16, 23, 24, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 15, 16, 25, 26, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 15, 16, 33, 34, 21, 22, 27, 28, 29, 30),
  c( 7,  8, 37, 38, 19, 20, 21, 22, 27, 28, 29, 30),
  c(35, 36, 11, 12, 19, 20, 21, 22, 27, 28, 29, 30),
  c(35, 36, 15, 16, 33, 34, 21, 22, 27, 28, 29, 30),
  c(35, 36, 37, 38, 33, 34, 21, 22, 27, 28, 29, 30)
)

run_bootstrapping_combinations <- function(taxonomy, year_label) {

  # Pre-filter: only Red-List beetles, exclude Ground
  rl_data <- taxonomy %>%
    filter(TreeSp != "Ground", Rote.Liste != "n") %>%
    select(Trap, FHL, TreeSp) %>%
    distinct()

  tree_species <- c("Fraxinus", "Quercus", "Tilia", "Acer", "Quercus rubra", "Ulmus")

  results <- matrix(0,
                    nrow = length(trap_combinations),
                    ncol = length(tree_species),
                    dimnames = list(NULL, tree_species))

  for (i in seq_along(trap_combinations)) {
    traps <- trap_combinations[[i]]

    rl_subset <- rl_data %>%
      filter(Trap %in% traps) %>%
      select(FHL, TreeSp) %>%
      distinct() %>%
      group_by(FHL) %>%
      mutate(n_treesp = n_distinct(TreeSp)) %>%
      ungroup() %>%
      filter(n_treesp == 1)

    for (sp in tree_species) {
      results[i, sp] <- sum(rl_subset$TreeSp == sp)
    }
  }

  means <- colMeans(results)
  for (sp in c("Fraxinus", "Quercus", "Tilia")) {
    cat(sprintf("%s %s — mean unique red-list: %.2f\n",
                sp, year_label, means[sp]))
  }
  as.list(means)
}

bootstrap_2017 <- run_bootstrapping_combinations(taxonomy_2017, "2017")
bootstrap_2018 <- run_bootstrapping_combinations(taxonomy_2018, "2018")


# ============================================================
# BRIDGE — Assemble beetle input table
#
# Sources:
#   total species       -> iNEXT at t=2 (trap-based rarefaction)
#   total red list      -> iNEXT at t=2 (trap-based rarefaction)
#   unique species      -> iNEXT at t=2 (Fraxinus/Quercus/Tilia only)
#   unique red list     -> Bootstrapping mean (10-trap species)
#                      -> count_unique_rl_2trap() (2-trap species)
#
# Disjoint categories:
#   unique red list        = bootstrapping / direct count
#   non-unique red list    = total red list - unique red list
#   unique non-red list    = total unique   - unique red list
#   common beetles         = total - unique non-rl - non-unique rl - unique rl
#
# Mean values rounded to integers (consistent with original methodology)
# ============================================================

tree_order <- c("Acer", "Fraxinus", "Quercus", "Quercus rubra", "Tilia", "Ulmus")

# Helper: extract qD at t=2 by assemblage name
extract_t2 <- function(inext_obj, assemblage_order) {
  bind_rows(inext_obj$iNextEst) %>%
    filter(t == 2, !is.na(SC.LCL)) %>%
    select(Assemblage, qD) %>%
    right_join(tibble(Assemblage = assemblage_order), by = "Assemblage") %>%
    arrange(match(Assemblage, assemblage_order)) %>%
    pull(qD)
}

# Helper: extract unique species qD (only Fraxinus/Quercus/Tilia)
extract_unique_sp <- function(unique_out, tree_order) {
  partial <- bind_rows(unique_out$iNextEst) %>%
    filter(t == 2, !is.na(SC.LCL)) %>%
    select(Assemblage, qD)
  sapply(tree_order, function(sp) {
    val <- partial$qD[partial$Assemblage == sp]
    if (length(val) == 0) NA_real_ else val
  })
}

# Unique red list values — alle Baumarten aus Bootstrapping
unique_rl_2017 <- c(
  Acer            = bootstrap_2017[["Acer"]],
  Fraxinus        = bootstrap_2017[["Fraxinus"]],
  Quercus         = bootstrap_2017[["Quercus"]],
  `Quercus rubra` = bootstrap_2017[["Quercus rubra"]],
  Tilia           = bootstrap_2017[["Tilia"]],
  Ulmus           = bootstrap_2017[["Ulmus"]]
)

unique_rl_2018 <- c(
  Acer            = bootstrap_2018[["Acer"]],
  Fraxinus        = bootstrap_2018[["Fraxinus"]],
  Quercus         = bootstrap_2018[["Quercus"]],
  `Quercus rubra` = bootstrap_2018[["Quercus rubra"]],
  Tilia           = bootstrap_2018[["Tilia"]],
  Ulmus           = bootstrap_2018[["Ulmus"]]
)

# Build full beetle table for one year
build_beetle_table <- function(total_inext, red_inext, unique_inext,
                               unique_rl_vec, tree_order) {
  total_sp  <- extract_t2(total_inext,  tree_order)
  total_rl  <- extract_t2(red_inext,    tree_order)
  unique_sp <- extract_unique_sp(unique_inext, tree_order)
  rl_vec    <- unique_rl_vec[tree_order]

  non_unique_rl <- total_rl  - rl_vec
  unique_non_rl <- unique_sp - rl_vec
  common        <- total_sp  - unique_non_rl - non_unique_rl - rl_vec

  df <- rbind(
    `total beetle species`               = total_sp,
    `unique non-red list beetle species` = unique_non_rl,
    `Non-unique red list beetle species` = non_unique_rl,
    `unique red list beetle species`     = rl_vec,
    `Common beetles`                     = common
  )
  colnames(df) <- tree_order
  as.data.frame(df)
}

beetle_table_2017 <- build_beetle_table(
  tree_total_2017, tree_red_2017,
  unique_inext_2017, unique_rl_2017, tree_order
)

beetle_table_2018 <- build_beetle_table(
  tree_total_2018, tree_red_2018,
  unique_inext_2018, unique_rl_2018, tree_order
)

# Mean rounded to integers — consistent with original methodology
beetle_table_mean <- as.data.frame(round((beetle_table_2017 + beetle_table_2018) / 2))

write.csv(beetle_table_2017, "beetle_data_2017.csv")
write.csv(beetle_table_2018, "beetle_data_2018.csv")
write.csv(beetle_table_mean, "beetle_data_mean.csv")

ces_input <- beetle_table_mean


# ============================================================
# PART 2 — TREE BIODIVERSITY INDEX (TBI) & FOREST INDEX (FBI)
# ============================================================

# ------------------------------------------------------------
# 2.1  Select disjoint beetle categories & rename rows
# ------------------------------------------------------------

ces_data <- ces_input[c(
  "unique red list beetle species",
  "Non-unique red list beetle species",
  "unique non-red list beetle species",
  "Common beetles"
), ]

rownames(ces_data) <- c(
  "unique_red_list",
  "non_unique_red_list",
  "unique_non_red_list",
  "common"
)


# ------------------------------------------------------------
# 2.2  Index functions
# ------------------------------------------------------------

linear_index <- function(x, a) sum(a * x)

ces_general <- function(x, a, sigma) {
  rho <- (sigma - 1) / sigma
  (sum(a * (x ^ rho))) ^ (1 / rho)
}

ces_cobb_douglas <- function(x, a) exp(sum(a * log(x)))

forest_aggregate <- function(x, sigma) {
  a <- rep(1 / length(x), length(x))
  if (abs(sigma - 1) < 1e-8) ces_cobb_douglas(x, a) else ces_general(x, a, sigma)
}


# ------------------------------------------------------------
# 2.3  Weighting schemes W1-W10
# ------------------------------------------------------------

weights_list <- list(
  W1  = c(0.25, 0.25, 0.25, 0.25),
  W2  = c(0.50, 0.30, 0.15, 0.05),
  W3  = c(0.35, 0.25, 0.25, 0.15),
  W4  = c(0.40, 0.10, 0.40, 0.10),
  W5  = c(0.45, 0.05, 0.40, 0.10),
  W6  = c(0.40, 0.30, 0.25, 0.05),
  W7  = c(0.30, 0.25, 0.25, 0.20),
  W8  = c(0.60, 0.35, 0.03, 0.02),
  W9  = c(0.50, 0.05, 0.40, 0.05),
  W10 = c(0.15, 0.15, 0.20, 0.50)
)


# ------------------------------------------------------------
# 2.4  Compute TBI (linear, all weights)
# ------------------------------------------------------------

tbi_results <- list()
for (w in names(weights_list)) {
  tbi_results[[paste0(w, "_linear")]] <- sapply(
    colnames(ces_data),
    function(tree) linear_index(ces_data[, tree], weights_list[[w]])
  )
}

ces_results_df <- as.data.frame(tbi_results)
rownames(ces_results_df) <- colnames(ces_data)

write.csv(round(ces_results_df, 3), "CES_results_only.csv")


# ------------------------------------------------------------
# 2.5  Load PPA forest simulation data
# ------------------------------------------------------------

species_lookup <- tibble(
  sp = 1:8,
  species = c(
    "Acer pseudoplatanus", "Acer campestre",
    "Fraxinus Excelsior",  "Carpinus betulus",
    "Acer platanoides",    "Quercus robur",
    "Ulmus spp.",          "Tilia spp."
  )
)

scenarios <- tibble(
  file = c(
    "Results_combined_census_dry_init_dry.csv",
    "Results_combined_census_middle_init_dry.csv",
    "Results_combined_census_moist_init_dry.csv",
    "Results_combined_census_middle_init_intermediate.csv",
    "Results_combined_census_moist_init_intermediate.csv",
    "Results_combined_census_wet_init_intermediate.csv",
    "Results_combined_census_moist_init_moist.csv",
    "Results_combined_census_wet_init_moist.csv",
    "Results_combined_census_v_wet_init_moist.csv"
  ),
  initial_state        = c("Dry","Dry","Dry",
                           "Intermediate","Intermediate","Intermediate",
                           "Moist","Moist","Moist"),
  groundwater_change_m = c(0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0)
)

compute_basal_area <- function(file, init, delta_gw) {
  read_csv(file, show_col_types = FALSE) %>%
    filter(cl == 1, time %in% c(0, 25, 50, 75, 100)) %>%
    mutate(
      dbh_m            = dbh / 100,
      basal_area_m2_ha = n * pi * (dbh_m / 2)^2
    ) %>%
    group_by(time, sp) %>%
    summarise(basal_area_m2_ha = sum(basal_area_m2_ha, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(species_lookup, by = "sp") %>%
    mutate(initial_state = init, groundwater_change_m = delta_gw)
}

all_results <- pmap_dfr(
  list(scenarios$file, scenarios$initial_state, scenarios$groundwater_change_m),
  compute_basal_area
)


# ------------------------------------------------------------
# 2.6  Compute FBI (all weights x sigma scenarios)
# ------------------------------------------------------------

species_mapping <- c(
  "Acer pseudoplatanus" = "Acer",
  "Fraxinus Excelsior"  = "Fraxinus",
  "Quercus robur"       = "Quercus",
  "Tilia spp."          = "Tilia",
  "Ulmus spp."          = "Ulmus"
)

valid_forest_species <- names(species_mapping)
valid_beetle_species <- unname(species_mapping)

sigma_scenarios <- tibble(
  aggregation = c("CES_sigma_0.4", "CES_sigma_0.7", "Cobb_Douglas_sigma_1"),
  sigma       = c(0.4, 0.7, 1.0)
)

forest_files <- tibble(
  groundwater_scenario = c("no_increase", "0.5m_increase", "1m_increase"),
  gw_value             = c(0.0, 0.5, 1.0)
)

forest_results <- list()

for (gw in forest_files$gw_value) {

  forest_mean <- all_results %>%
    filter(
      groundwater_change_m == gw,
      initial_state %in% c("Dry", "Intermediate", "Moist"),
      species %in% valid_forest_species
    ) %>%
    group_by(species, time) %>%
    summarise(basal_area_m2_ha = mean(basal_area_m2_ha), .groups = "drop") %>%
    pivot_wider(names_from = time, values_from = basal_area_m2_ha,
                names_prefix = "year_") %>%
    column_to_rownames("species")

  rownames(forest_mean) <- species_mapping[rownames(forest_mean)]
  forest_mean           <- forest_mean[valid_beetle_species, ]
  years                 <- as.numeric(gsub("year_", "", colnames(forest_mean)))

  for (w in paste0("W", 1:10, "_linear")) {

    BV_i <- ces_results_df[valid_beetle_species, w, drop = TRUE]

    for (i in seq_len(nrow(sigma_scenarios))) {

      agg       <- sigma_scenarios$aggregation[i]
      sigma_val <- sigma_scenarios$sigma[i]

      vals <- sapply(forest_mean, function(BA_vec) {
        x_i <- BV_i * as.numeric(BA_vec)
        forest_aggregate(x_i, sigma_val)
      })

      forest_results[[length(forest_results) + 1]] <- tibble(
        groundwater_scenario      = forest_files$groundwater_scenario[
          forest_files$gw_value == gw],
        year                      = years,
        weight_scenario           = w,
        aggregation               = agg,
        forest_biodiversity_index = vals
      )
    }
  }
}

forest_biodiversity_all <- bind_rows(forest_results)

write_csv(forest_biodiversity_all,
          "Forest_Biodiversity_Index_All_Scenarios_AllWeights_AllSigmas.csv")


# ------------------------------------------------------------
# 2.7  Wide-format CSVs per sigma scenario
# ------------------------------------------------------------

for (agg_name in unique(forest_biodiversity_all$aggregation)) {
  out <- forest_biodiversity_all %>%
    filter(aggregation == agg_name) %>%
    mutate(groundwater_scenario = recode(
      groundwater_scenario,
      "no_increase"   = "no_increase",
      "0.5m_increase" = "gw_0.5m",
      "1m_increase"   = "gw_1m"
    )) %>%
    pivot_wider(names_from  = groundwater_scenario,
                values_from = forest_biodiversity_index) %>%
    arrange(weight_scenario, year)
  write_csv(out, paste0("Forest_Biodiversity_", agg_name, ".csv"))
}


# ------------------------------------------------------------
# 2.8  Plots
# ------------------------------------------------------------

# Graph 2: FBI W3, all sigma
plot_forest_W3 <- forest_biodiversity_all %>%
  filter(weight_scenario == "W3_linear") %>%
  mutate(
    groundwater_scenario = recode(groundwater_scenario,
      "no_increase"   = "No increase",
      "0.5m_increase" = "0.5 m increase",
      "1m_increase"   = "1 m increase"),
    aggregation = recode(aggregation,
      "CES_sigma_0.4"        = "σ = 0.4",
      "CES_sigma_0.7"        = "σ = 0.7",
      "Cobb_Douglas_sigma_1" = "σ = 1")
  )

ggplot(plot_forest_W3,
       aes(x = year, y = forest_biodiversity_index,
           color = groundwater_scenario, group = groundwater_scenario)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~aggregation, nrow = 1) +
  labs(x = "Time (years)", y = "Forest Biodiversity Index (W3)", color = NULL) +
  theme_classic(base_size = 16) +
  theme(strip.text       = element_text(face = "bold", size = 14),
        legend.position  = "bottom",
        legend.text      = element_text(size = 20),
        legend.key.width = unit(2.2, "cm"))
ggsave("Forest_Biodiversity_W3_AllSigmas.png", width = 15, height = 5, dpi = 300)


# Graph 3: Mean basal area per species and GW scenario
plot_df <- all_results %>%
  group_by(groundwater_change_m, species, time) %>%
  summarise(basal_area_m2_ha = mean(basal_area_m2_ha), .groups = "drop") %>%
  mutate(
    scenario = case_when(
      groundwater_change_m == 0   ~ "no increase",
      groundwater_change_m == 0.5 ~ "0.5 m increase",
      groundwater_change_m == 1   ~ "1 m increase"
    ),
    scenario = factor(scenario,
                      levels = c("no increase", "0.5 m increase", "1 m increase"))
  )

species_cols <- c(
  "Acer campestre"      = "#9E9E9E",
  "Acer platanoides"    = "#9E9E9E",
  "Carpinus betulus"    = "#9E9E9E",
  "Acer pseudoplatanus" = "#D7191C",
  "Fraxinus Excelsior"  = "#D7191C",
  "Ulmus spp."          = "#D7191C",
  "Tilia spp."          = "#D7191C",
  "Quercus robur"       = "#D7191C"
)

ggplot(plot_df, aes(x = time, y = basal_area_m2_ha,
                    color = species, group = species)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~scenario, nrow = 1) +
  scale_color_manual(values = species_cols) +
  labs(x = "Time (years)",
       y = expression("Basal area ("*m^2~ha^{-1}*")"),
       color = NULL) +
  theme_classic(base_size = 20) +
  theme(strip.text       = element_text(face = "bold", size = 16),
        legend.position  = "bottom",
        legend.text      = element_text(size = 20),
        legend.key.width = unit(2.2, "cm"))
ggsave("BasalArea_Mean_Scenarios.png", width = 14, height = 5, dpi = 300)


# ============================================================
# PART 3 — SENSITIVITY ANALYSIS
# ============================================================

# ------------------------------------------------------------
# 3.1  Beetle TBI — Spearman rank correlation (W1-W10)
# ------------------------------------------------------------

beetle_df <- read_csv("CES_results_only.csv", show_col_types = FALSE) %>%
  column_to_rownames(var = colnames(.)[1])

weight_cols <- grep("^W[0-9]+_linear$", colnames(beetle_df), value = TRUE)

ranks <- apply(beetle_df[, weight_cols], 2, rank)

spearman_mat <- cor(ranks, method = "spearman")

print("Spearman rank correlation matrix:")
print(round(spearman_mat, 3))

beetle_spearman_df <- as.data.frame(round(spearman_mat, 3)) %>%
  rownames_to_column(var = "weight_scenario")
write_csv(beetle_spearman_df, "Sensitivity_Beetle_Rank_Spearman_Matrix.csv")


# ------------------------------------------------------------
# 3.2  Kendall's W — overall ranking consistency
# ------------------------------------------------------------

m   <- ncol(ranks)
n   <- nrow(ranks)
R_i <- rowSums(ranks)
W_kendall <- 12 * sum((R_i - mean(R_i))^2) / (m^2 * (n^3 - n))

cat(sprintf("\nKendall's W (ranking consistency): %.3f\n", W_kendall))
if (W_kendall < 0.3)                     cat("-> Low ranking stability\n")
if (W_kendall >= 0.3 & W_kendall < 0.7) cat("-> Moderate ranking stability\n")
if (W_kendall >= 0.7)                    cat("-> High ranking stability\n")


# ------------------------------------------------------------
# 3.3  FBI — rank stability within weighting schemes
# ------------------------------------------------------------

df04 <- read_csv("Forest_Biodiversity_CES_sigma_0.4.csv",        show_col_types = FALSE)
df07 <- read_csv("Forest_Biodiversity_CES_sigma_0.7.csv",        show_col_types = FALSE)
df1  <- read_csv("Forest_Biodiversity_Cobb_Douglas_sigma_1.csv", show_col_types = FALSE)

rank_stability_within <- function(df) {
  results <- list()
  for (w in unique(df$weight_scenario)) {
    sub   <- df %>% filter(weight_scenario == w)
    ranks <- sub %>%
      rowwise() %>%
      mutate(r_no = rank(-no_increase),
             r_05 = rank(-gw_0.5m),
             r_1  = rank(-gw_1m)) %>%
      ungroup() %>%
      select(r_no, r_05, r_1)
    results[[w]] <- cor(ranks, method = "spearman")
  }
  results
}

stab04 <- rank_stability_within(df04)
stab07 <- rank_stability_within(df07)
stab1  <- rank_stability_within(df1)

print("Rank stability within W (sigma = 0.4):")
print(stab04)

compare_sigma <- function(dfA, dfB) {
  out <- list()
  for (w in unique(dfA$weight_scenario)) {
    a <- dfA %>% filter(weight_scenario == w)
    b <- dfB %>% filter(weight_scenario == w)
    out[[w]] <- cor(
      c(a$no_increase, a$gw_0.5m, a$gw_1m),
      c(b$no_increase, b$gw_0.5m, b$gw_1m),
      method = "spearman"
    )
  }
  out
}

rho_04_07 <- compare_sigma(df04, df07)
rho_04_1  <- compare_sigma(df04, df1)
rho_07_1  <- compare_sigma(df07, df1)

print("Spearman rho (sigma 0.4 vs 0.7):"); print(rho_04_07)
print("Spearman rho (sigma 0.4 vs 1.0):"); print(rho_04_1)
print("Spearman rho (sigma 0.7 vs 1.0):"); print(rho_07_1)


# ------------------------------------------------------------
# 3.4  FBI — rank files per sigma scenario
# ------------------------------------------------------------

rank_forest_csv <- function(input_file, output_file) {
  df <- read_csv(input_file, show_col_types = FALSE)
  df_ranked <- df %>%
    rowwise() %>%
    mutate(
      ranks       = list(rank(-c(no_increase, gw_0.5m, gw_1m), ties.method = "average")),
      no_increase = ranks[[1]],
      gw_0.5m     = ranks[[2]],
      gw_1m       = ranks[[3]]
    ) %>%
    select(-ranks) %>%
    ungroup()
  write_csv(df_ranked, output_file)
}

rank_forest_csv("Forest_Biodiversity_CES_sigma_0.4.csv",
                "Forest_Biodiversity_CES_sigma_0.4_RANKS.csv")
rank_forest_csv("Forest_Biodiversity_CES_sigma_0.7.csv",
                "Forest_Biodiversity_CES_sigma_0.7_RANKS.csv")
rank_forest_csv("Forest_Biodiversity_Cobb_Douglas_sigma_1.csv",
                "Forest_Biodiversity_Cobb_Douglas_sigma_1_RANKS.csv")


# ============================================================
# END OF SCRIPT
# ============================================================
