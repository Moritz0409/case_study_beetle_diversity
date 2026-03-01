###############################################################################
# ForestDiversityValue_LaesekeKulkarniNguyenVogt
#
# Merged script for beetle biodiversity / forest biodiversity analysis
#
# PART 1: Beetle species richness analysis (iNEXT, bootstrapping, uniqueness)
# PART 2: Tree Biodiversity Index (TBI) + Forest Simulation + Forest
#          Biodiversity Index (FBI)
# PART 3: Sensitivity analysis (Spearman rank correlation)
###############################################################################


#############################################################
# GLOBAL SETUP
#############################################################

setwd("/Users/Jonas/Desktop/Master/3rd Semester/Integration Module/BEETLES/test_run")

library(dplyr)
library(reshape2)
library(iNEXT)
library(readr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(forcats)
library(ggpubr)


#############################################################
# PART 1 — BEETLE SPECIES RICHNESS ANALYSIS
#############################################################

# -----------------------------------------------------------
# 1.1  Load and merge data
# -----------------------------------------------------------

raw_2017 <- read.csv("data_2017_species.csv")
data_2018_species <- read.csv("data_2018_species.csv")
species_taxonomy_hahn <- read.csv("species_taxonomy_hahn.csv")

taxonomy_2017 <- merge(x = raw_2017, y = species_taxonomy_hahn[ , c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")], by.x ="FHL", by.y = "Colkat")
taxonomy_2018 <- merge(x = data_2018_species, y = species_taxonomy_hahn[ , c("Colkat", "xylobiont", "Rote.Liste", "Urwaldrelikt")], by.x ="FHL", by.y = "Colkat")
taxonomy_2018 <- taxonomy_2018 %>%
  mutate(Species = paste(Family, Genus, Art, sep = " "))

# Red list subsets
beetles_2017_red <- taxonomy_2017[which(taxonomy_2017$Rote.Liste != "n" & taxonomy_2017$Number != 0), ]
beetles_2018_red <- taxonomy_2018[which(taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0), ]


# -----------------------------------------------------------
# 1.2  Total beetle species 2017 (iNEXT q=1, NO xylobiont filter)
# -----------------------------------------------------------

cat("=== Computing total beetle species ===\n")

# Incidence matrices per tree species 2017
sp_per_trap_Rubra_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Rubra_total) <- sp_per_trap_Rubra_total[,1]
sp_per_trap_Rubra_total <- sp_per_trap_Rubra_total[,-1]
sp_per_trap_Rubra_total[sp_per_trap_Rubra_total>=1] <- 1

sp_per_trap_Fraxinus_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Fraxinus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Fraxinus_total) <- sp_per_trap_Fraxinus_total[,1]
sp_per_trap_Fraxinus_total <- sp_per_trap_Fraxinus_total[,-1]
sp_per_trap_Fraxinus_total[sp_per_trap_Fraxinus_total>=1] <- 1

sp_per_trap_Quercus_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Quercus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Quercus_total) <- sp_per_trap_Quercus_total[,1]
sp_per_trap_Quercus_total <- sp_per_trap_Quercus_total[,-1]
sp_per_trap_Quercus_total[sp_per_trap_Quercus_total>=1] <- 1

sp_per_trap_Acer_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Acer",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Acer_total) <- sp_per_trap_Acer_total[,1]
sp_per_trap_Acer_total <- sp_per_trap_Acer_total[,-1]
sp_per_trap_Acer_total[sp_per_trap_Acer_total>=1] <- 1

sp_per_trap_Tilia_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Tilia",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Tilia_total) <- sp_per_trap_Tilia_total[,1]
sp_per_trap_Tilia_total <- sp_per_trap_Tilia_total[,-1]
sp_per_trap_Tilia_total[sp_per_trap_Tilia_total>=1] <- 1

sp_per_trap_Ulmus_total <- dcast(taxonomy_2017[taxonomy_2017$TreeSp=="Ulmus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Ulmus_total) <- sp_per_trap_Ulmus_total[,1]
sp_per_trap_Ulmus_total <- sp_per_trap_Ulmus_total[,-1]
sp_per_trap_Ulmus_total[sp_per_trap_Ulmus_total>=1] <- 1

# iNEXT per species individually (avoids batch error for Q.rubra 2018)
total_2017_Acer     <- iNEXT(list(as.matrix(sp_per_trap_Acer_total)),     q=1, datatype="incidence_raw", endpoint=2)
total_2017_Fraxinus <- iNEXT(list(as.matrix(sp_per_trap_Fraxinus_total)), q=1, datatype="incidence_raw", endpoint=2)
total_2017_Quercus  <- iNEXT(list(as.matrix(sp_per_trap_Quercus_total)),  q=1, datatype="incidence_raw", endpoint=2)
total_2017_Rubra    <- iNEXT(list(as.matrix(sp_per_trap_Rubra_total)),    q=1, datatype="incidence_raw", endpoint=2)
total_2017_Tilia    <- iNEXT(list(as.matrix(sp_per_trap_Tilia_total)),    q=1, datatype="incidence_raw", endpoint=2)
total_2017_Ulmus    <- iNEXT(list(as.matrix(sp_per_trap_Ulmus_total)),    q=1, datatype="incidence_raw", endpoint=2)

# Extract qD at t=2
extract_qD <- function(inext_out) {
  t2 <- bind_rows(inext_out$iNextEst) %>% filter(t == 2, !is.na(SC.LCL))
  if (nrow(t2) == 0) t2 <- bind_rows(inext_out$iNextEst) %>% filter(t == 2)
  t2$qD[1]
}

# Safe iNEXT wrapper: falls back to observed Shannon diversity (exp(H'))
# when iNEXT fails (e.g. Q.rubra 2018 red list: 14 species, 2 traps)
safe_iNEXT <- function(mat, q=1, datatype="incidence_raw", endpoint=2) {
  tryCatch(
    iNEXT(list(as.matrix(mat)), q=q, datatype=datatype, endpoint=endpoint),
    error = function(e) {
      cat("  iNEXT failed, computing observed diversity directly\n")
      mat_bin <- as.matrix(mat)
      mat_bin[mat_bin >= 1] <- 1
      freq <- rowSums(mat_bin)
      freq <- freq[freq > 0]
      p <- freq / sum(freq)
      list(
        iNextEst = list(data.frame(
          t = endpoint, qD = exp(-sum(p * log(p))),
          SC.LCL = NA, Assemblage = "fallback"
        )),
        fallback = TRUE
      )
    }
  )
}

total_2017 <- c(
  "Acer"         = extract_qD(total_2017_Acer),
  "Fraxinus"     = extract_qD(total_2017_Fraxinus),
  "Quercus"      = extract_qD(total_2017_Quercus),
  "Quercus rubra"= extract_qD(total_2017_Rubra),
  "Tilia"        = extract_qD(total_2017_Tilia),
  "Ulmus"        = extract_qD(total_2017_Ulmus)
)


# -----------------------------------------------------------
# 1.3  Total beetle species 2018
# -----------------------------------------------------------

sp_per_trap_Rubra_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Rubra_total_2018) <- sp_per_trap_Rubra_total_2018[,1]
sp_per_trap_Rubra_total_2018 <- sp_per_trap_Rubra_total_2018[,-1]
sp_per_trap_Rubra_total_2018[sp_per_trap_Rubra_total_2018>=1] <- 1

sp_per_trap_Fraxinus_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Fraxinus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Fraxinus_total_2018) <- sp_per_trap_Fraxinus_total_2018[,1]
sp_per_trap_Fraxinus_total_2018 <- sp_per_trap_Fraxinus_total_2018[,-1]
sp_per_trap_Fraxinus_total_2018[sp_per_trap_Fraxinus_total_2018>=1] <- 1

sp_per_trap_Quercus_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Quercus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Quercus_total_2018) <- sp_per_trap_Quercus_total_2018[,1]
sp_per_trap_Quercus_total_2018 <- sp_per_trap_Quercus_total_2018[,-1]
sp_per_trap_Quercus_total_2018[sp_per_trap_Quercus_total_2018>=1] <- 1

sp_per_trap_Acer_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Acer",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Acer_total_2018) <- sp_per_trap_Acer_total_2018[,1]
sp_per_trap_Acer_total_2018 <- sp_per_trap_Acer_total_2018[,-1]
sp_per_trap_Acer_total_2018[sp_per_trap_Acer_total_2018>=1] <- 1

sp_per_trap_Tilia_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Tilia",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Tilia_total_2018) <- sp_per_trap_Tilia_total_2018[,1]
sp_per_trap_Tilia_total_2018 <- sp_per_trap_Tilia_total_2018[,-1]
sp_per_trap_Tilia_total_2018[sp_per_trap_Tilia_total_2018>=1] <- 1

sp_per_trap_Ulmus_total_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Ulmus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Ulmus_total_2018) <- sp_per_trap_Ulmus_total_2018[,1]
sp_per_trap_Ulmus_total_2018 <- sp_per_trap_Ulmus_total_2018[,-1]
sp_per_trap_Ulmus_total_2018[sp_per_trap_Ulmus_total_2018>=1] <- 1

total_2018_Acer     <- iNEXT(list(as.matrix(sp_per_trap_Acer_total_2018)),     q=1, datatype="incidence_raw", endpoint=2)
total_2018_Fraxinus <- iNEXT(list(as.matrix(sp_per_trap_Fraxinus_total_2018)), q=1, datatype="incidence_raw", endpoint=2)
total_2018_Quercus  <- iNEXT(list(as.matrix(sp_per_trap_Quercus_total_2018)),  q=1, datatype="incidence_raw", endpoint=2)
total_2018_Rubra    <- iNEXT(list(as.matrix(sp_per_trap_Rubra_total_2018)),    q=1, datatype="incidence_raw", endpoint=2)
total_2018_Tilia    <- iNEXT(list(as.matrix(sp_per_trap_Tilia_total_2018)),    q=1, datatype="incidence_raw", endpoint=2)
total_2018_Ulmus    <- iNEXT(list(as.matrix(sp_per_trap_Ulmus_total_2018)),    q=1, datatype="incidence_raw", endpoint=2)

total_2018 <- c(
  "Acer"         = extract_qD(total_2018_Acer),
  "Fraxinus"     = extract_qD(total_2018_Fraxinus),
  "Quercus"      = extract_qD(total_2018_Quercus),
  "Quercus rubra"= extract_qD(total_2018_Rubra),
  "Tilia"        = extract_qD(total_2018_Tilia),
  "Ulmus"        = extract_qD(total_2018_Ulmus)
)

order <- c("Acer", "Fraxinus", "Quercus", "Quercus rubra", "Tilia", "Ulmus")
total_mean <- round((total_2017[order] + total_2018[order]) / 2)

cat("Total beetle species (mean 2017/2018):\n")
cat(sprintf("  %s: %g\n", order, total_mean))


# -----------------------------------------------------------
# 1.4  Total red list species 2017 (iNEXT q=1, per species individually)
# -----------------------------------------------------------

cat("\n=== Computing total red list species ===\n")

sp_per_trap_Rubra_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Quercus rubra",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Rubra_red) <- sp_per_trap_Rubra_red[,1]
sp_per_trap_Rubra_red <- sp_per_trap_Rubra_red[,-1]
sp_per_trap_Rubra_red[sp_per_trap_Rubra_red>=1] <- 1

sp_per_trap_Fraxinus_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Fraxinus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Fraxinus_red) <- sp_per_trap_Fraxinus_red[,1]
sp_per_trap_Fraxinus_red <- sp_per_trap_Fraxinus_red[,-1]
sp_per_trap_Fraxinus_red[sp_per_trap_Fraxinus_red>=1] <- 1

sp_per_trap_Quercus_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Quercus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Quercus_red) <- sp_per_trap_Quercus_red[,1]
sp_per_trap_Quercus_red <- sp_per_trap_Quercus_red[,-1]
sp_per_trap_Quercus_red[sp_per_trap_Quercus_red>=1] <- 1

sp_per_trap_Acer_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Acer",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Acer_red) <- sp_per_trap_Acer_red[,1]
sp_per_trap_Acer_red <- sp_per_trap_Acer_red[,-1]
sp_per_trap_Acer_red[sp_per_trap_Acer_red>=1] <- 1

sp_per_trap_Tilia_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Tilia",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Tilia_red) <- sp_per_trap_Tilia_red[,1]
sp_per_trap_Tilia_red <- sp_per_trap_Tilia_red[,-1]
sp_per_trap_Tilia_red[sp_per_trap_Tilia_red>=1] <- 1

sp_per_trap_Ulmus_red <- dcast(beetles_2017_red[beetles_2017_red$TreeSp=="Ulmus",], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Ulmus_red) <- sp_per_trap_Ulmus_red[,1]
sp_per_trap_Ulmus_red <- sp_per_trap_Ulmus_red[,-1]
sp_per_trap_Ulmus_red[sp_per_trap_Ulmus_red>=1] <- 1

red_2017_Acer     <- iNEXT(list(as.matrix(sp_per_trap_Acer_red)),     q=1, datatype="incidence_raw", endpoint=2)
red_2017_Fraxinus <- iNEXT(list(as.matrix(sp_per_trap_Fraxinus_red)), q=1, datatype="incidence_raw", endpoint=2)
red_2017_Quercus  <- iNEXT(list(as.matrix(sp_per_trap_Quercus_red)),  q=1, datatype="incidence_raw", endpoint=2)
red_2017_Rubra    <- iNEXT(list(as.matrix(sp_per_trap_Rubra_red)),    q=1, datatype="incidence_raw", endpoint=2)
red_2017_Tilia    <- iNEXT(list(as.matrix(sp_per_trap_Tilia_red)),    q=1, datatype="incidence_raw", endpoint=2)
red_2017_Ulmus    <- iNEXT(list(as.matrix(sp_per_trap_Ulmus_red)),    q=1, datatype="incidence_raw", endpoint=2)

red_2017 <- c(
  "Acer"         = extract_qD(red_2017_Acer),
  "Fraxinus"     = extract_qD(red_2017_Fraxinus),
  "Quercus"      = extract_qD(red_2017_Quercus),
  "Quercus rubra"= extract_qD(red_2017_Rubra),
  "Tilia"        = extract_qD(red_2017_Tilia),
  "Ulmus"        = extract_qD(red_2017_Ulmus)
)


# -----------------------------------------------------------
# 1.5  Total red list species 2018
# -----------------------------------------------------------

sp_per_trap_Rubra_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Quercus rubra" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Rubra_red_2018) <- sp_per_trap_Rubra_red_2018[,1]
sp_per_trap_Rubra_red_2018 <- sp_per_trap_Rubra_red_2018[,-1]
sp_per_trap_Rubra_red_2018[sp_per_trap_Rubra_red_2018>=1] <- 1

sp_per_trap_Fraxinus_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Fraxinus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Fraxinus_red_2018) <- sp_per_trap_Fraxinus_red_2018[,1]
sp_per_trap_Fraxinus_red_2018 <- sp_per_trap_Fraxinus_red_2018[,-1]
sp_per_trap_Fraxinus_red_2018[sp_per_trap_Fraxinus_red_2018>=1] <- 1

sp_per_trap_Quercus_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Quercus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Quercus_red_2018) <- sp_per_trap_Quercus_red_2018[,1]
sp_per_trap_Quercus_red_2018 <- sp_per_trap_Quercus_red_2018[,-1]
sp_per_trap_Quercus_red_2018[sp_per_trap_Quercus_red_2018>=1] <- 1

sp_per_trap_Acer_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Acer" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Acer_red_2018) <- sp_per_trap_Acer_red_2018[,1]
sp_per_trap_Acer_red_2018 <- sp_per_trap_Acer_red_2018[,-1]
sp_per_trap_Acer_red_2018[sp_per_trap_Acer_red_2018>=1] <- 1

sp_per_trap_Tilia_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Tilia" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Tilia_red_2018) <- sp_per_trap_Tilia_red_2018[,1]
sp_per_trap_Tilia_red_2018 <- sp_per_trap_Tilia_red_2018[,-1]
sp_per_trap_Tilia_red_2018[sp_per_trap_Tilia_red_2018>=1] <- 1

sp_per_trap_Ulmus_red_2018 <- dcast(taxonomy_2018[taxonomy_2018$TreeSp=="Ulmus" & taxonomy_2018$Rote.Liste != "n" & taxonomy_2018$Number != 0,], Species~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(sp_per_trap_Ulmus_red_2018) <- sp_per_trap_Ulmus_red_2018[,1]
sp_per_trap_Ulmus_red_2018 <- sp_per_trap_Ulmus_red_2018[,-1]
sp_per_trap_Ulmus_red_2018[sp_per_trap_Ulmus_red_2018>=1] <- 1

red_2018_Acer     <- safe_iNEXT(sp_per_trap_Acer_red_2018)
red_2018_Fraxinus <- safe_iNEXT(sp_per_trap_Fraxinus_red_2018)
red_2018_Quercus  <- safe_iNEXT(sp_per_trap_Quercus_red_2018)
red_2018_Rubra    <- safe_iNEXT(sp_per_trap_Rubra_red_2018)
red_2018_Tilia    <- safe_iNEXT(sp_per_trap_Tilia_red_2018)
red_2018_Ulmus    <- safe_iNEXT(sp_per_trap_Ulmus_red_2018)

red_2018 <- c(
  "Acer"         = extract_qD(red_2018_Acer),
  "Fraxinus"     = extract_qD(red_2018_Fraxinus),
  "Quercus"      = extract_qD(red_2018_Quercus),
  "Quercus rubra"= extract_qD(red_2018_Rubra),
  "Tilia"        = extract_qD(red_2018_Tilia),
  "Ulmus"        = extract_qD(red_2018_Ulmus)
)

total_rl_mean <- round((red_2017[order] + red_2018[order]) / 2)

cat("Total red list (mean 2017/2018):\n")
cat(sprintf("  %s: %g\n", order, total_rl_mean))


# -----------------------------------------------------------
# 1.6  Unique total species
#       2-trap species (Acer, Q.rubra, Ulmus): plain count of n_appearance==1
#       10-trap species (Fraxinus, Quercus, Tilia): iNEXT on unique species at t=2
# -----------------------------------------------------------

cat("\n=== Computing unique total species ===\n")

### Beetles uniqueness total 2017

beetles_2017_treesp <- taxonomy_2017[taxonomy_2017$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number), .groups = "drop")

distinct_2017 <- beetles_2017_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n())

beetles_2017_b <- beetles_2017_treesp %>%
  left_join(distinct_2017, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL), .groups = "drop")

# Plain count for all species (n_appearance==1)
unique_plain_2017 <- beetles_2017_b %>%
  filter(n_appearance == 1) %>%
  select(TreeSp, n_beetlesp) %>%
  tibble::deframe()

# Unique species for iNEXT (10-trap species only)
unique_2017 <- distinct_2017 %>%
  filter(n_appearance == 1) %>%
  left_join(
    taxonomy_2017 %>% select(TreeSp, FHL, Trap, Number),
    by = c("TreeSp", "FHL")
  )

# iNEXT for uniques total 2017 (Fraxinus, Quercus, Tilia)
unique_per_trap_Fraxinus_2017 <- dcast(unique_2017[unique_2017$TreeSp=="Fraxinus",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Fraxinus_2017) <- unique_per_trap_Fraxinus_2017[,1]
unique_per_trap_Fraxinus_2017 <- unique_per_trap_Fraxinus_2017[,-1]
unique_per_trap_Fraxinus_2017[unique_per_trap_Fraxinus_2017>=1] <- 1

unique_per_trap_Quercus_2017 <- dcast(unique_2017[unique_2017$TreeSp=="Quercus",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Quercus_2017) <- unique_per_trap_Quercus_2017[,1]
unique_per_trap_Quercus_2017 <- unique_per_trap_Quercus_2017[,-1]
unique_per_trap_Quercus_2017[unique_per_trap_Quercus_2017>=1] <- 1

unique_per_trap_Tilia_2017 <- dcast(unique_2017[unique_2017$TreeSp=="Tilia",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Tilia_2017) <- unique_per_trap_Tilia_2017[,1]
unique_per_trap_Tilia_2017 <- unique_per_trap_Tilia_2017[,-1]
unique_per_trap_Tilia_2017[unique_per_trap_Tilia_2017>=1] <- 1

unique_treesp_2017 <- list("Fraxinus" = unique_per_trap_Fraxinus_2017, "Quercus" = unique_per_trap_Quercus_2017, "Tilia" = unique_per_trap_Tilia_2017)
unique_inext_2017 <- iNEXT(unique_treesp_2017, q=1, datatype="incidence_raw", endpoint=2)
unique_inext_2017_out <- bind_rows(unique_inext_2017$iNextEst) %>% filter(t == 2, !is.na(SC.LCL))

### Beetles uniqueness total 2018

beetles_2018_treesp <- taxonomy_2018[taxonomy_2018$TreeSp != "Ground", ] %>%
  select(TreeSp, FHL, Number) %>%
  group_by(TreeSp, FHL) %>%
  summarise(n_total = sum(Number), .groups = "drop")

distinct_2018 <- beetles_2018_treesp %>%
  distinct(TreeSp, FHL) %>%
  group_by(FHL) %>%
  mutate(n_appearance = n())

beetles_2018_b <- beetles_2018_treesp %>%
  left_join(distinct_2018, by = c("TreeSp", "FHL")) %>%
  group_by(TreeSp, n_appearance) %>%
  summarise(n_beetlesp = n_distinct(FHL), .groups = "drop")

unique_plain_2018 <- beetles_2018_b %>%
  filter(n_appearance == 1) %>%
  select(TreeSp, n_beetlesp) %>%
  tibble::deframe()

unique_2018 <- distinct_2018 %>%
  filter(n_appearance == 1) %>%
  left_join(
    taxonomy_2018 %>% select(TreeSp, FHL, Trap, Number),
    by = c("TreeSp", "FHL")
  )

# iNEXT for uniques total 2018 (uses FHL~Trap for 2018)
unique_per_trap_Fraxinus_2018 <- dcast(unique_2018[unique_2018$TreeSp=="Fraxinus",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Fraxinus_2018) <- unique_per_trap_Fraxinus_2018[,1]
unique_per_trap_Fraxinus_2018 <- unique_per_trap_Fraxinus_2018[,-1]
unique_per_trap_Fraxinus_2018[unique_per_trap_Fraxinus_2018>=1] <- 1

unique_per_trap_Quercus_2018 <- dcast(unique_2018[unique_2018$TreeSp=="Quercus",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Quercus_2018) <- unique_per_trap_Quercus_2018[,1]
unique_per_trap_Quercus_2018 <- unique_per_trap_Quercus_2018[,-1]
unique_per_trap_Quercus_2018[unique_per_trap_Quercus_2018>=1] <- 1

unique_per_trap_Tilia_2018 <- dcast(unique_2018[unique_2018$TreeSp=="Tilia",], FHL~Trap, value.var="Number", fill=0, fun.aggregate=sum)
rownames(unique_per_trap_Tilia_2018) <- unique_per_trap_Tilia_2018[,1]
unique_per_trap_Tilia_2018 <- unique_per_trap_Tilia_2018[,-1]
unique_per_trap_Tilia_2018[unique_per_trap_Tilia_2018>=1] <- 1

unique_treesp_2018 <- list("Fraxinus" = unique_per_trap_Fraxinus_2018, "Quercus" = unique_per_trap_Quercus_2018, "Tilia" = unique_per_trap_Tilia_2018)
unique_inext_2018 <- iNEXT(unique_treesp_2018, q=1, datatype="incidence_raw", endpoint=2)
unique_inext_2018_out <- bind_rows(unique_inext_2018$iNextEst) %>% filter(t == 2, !is.na(SC.LCL))

# Assemble unique_total: plain count for 2-trap, iNEXT for 10-trap
unique_total_2017 <- sapply(order, function(sp) {
  v <- unique_plain_2017[sp]
  if (is.na(v)) 0L else as.integer(v)
})
unique_total_2017["Fraxinus"] <- unique_inext_2017_out$qD[unique_inext_2017_out$Assemblage == "Fraxinus"]
unique_total_2017["Quercus"]  <- unique_inext_2017_out$qD[unique_inext_2017_out$Assemblage == "Quercus"]
unique_total_2017["Tilia"]    <- unique_inext_2017_out$qD[unique_inext_2017_out$Assemblage == "Tilia"]

unique_total_2018 <- sapply(order, function(sp) {
  v <- unique_plain_2018[sp]
  if (is.na(v)) 0L else as.integer(v)
})
unique_total_2018["Fraxinus"] <- unique_inext_2018_out$qD[unique_inext_2018_out$Assemblage == "Fraxinus"]
unique_total_2018["Quercus"]  <- unique_inext_2018_out$qD[unique_inext_2018_out$Assemblage == "Quercus"]
unique_total_2018["Tilia"]    <- unique_inext_2018_out$qD[unique_inext_2018_out$Assemblage == "Tilia"]

unique_total_mean <- (unique_total_2017 + unique_total_2018) / 2

cat("Unique total (mean 2017/2018):\n")
cat(sprintf("  %s: %.2f\n", order, unique_total_mean))


# -----------------------------------------------------------
# 1.7  Unique red list species (bootstrapping for 10-trap, direct for 2-trap)
# -----------------------------------------------------------

cat("\n=== Computing unique red list species ===\n")

# 25 trap combinations (original, including 11-trap combos #4 and #14)
trap_combinations <- list(
  c(1,2,11,12,17,18,21,22,27,28,29,30), c(1,2,13,14,23,24,21,22,27,28,29,30),
  c(1,2,15,16,17,18,21,22,27,28,29,30), c(1,2,37,33,34,21,22,27,28,29,30),
  c(3,4,11,12,25,26,21,22,27,28,29,30), c(3,4,13,14,23,24,21,22,27,28,29,30),
  c(3,4,15,16,17,18,21,22,27,28,29,30), c(3,4,37,38,19,20,21,22,27,28,29,30),
  c(5,6,11,12,33,34,21,22,27,28,29,30), c(5,6,13,14,23,24,21,22,27,28,29,30),
  c(5,6,15,16,25,26,21,22,27,28,29,30), c(5,6,37,38,19,20,21,22,27,28,29,30),
  c(7,8,9,10,17,18,21,22,27,28,29,30), c(7,8,9,10,26,21,22,27,28,29,30),
  c(7,8,11,12,19,20,21,22,27,28,29,30), c(7,8,13,14,25,26,21,22,27,28,29,30),
  c(7,8,13,14,33,34,21,22,27,28,29,30), c(7,8,15,16,17,18,21,22,27,28,29,30),
  c(7,8,15,16,23,24,21,22,27,28,29,30), c(7,8,15,16,25,26,21,22,27,28,29,30),
  c(7,8,15,16,33,34,21,22,27,28,29,30), c(7,8,37,38,19,20,21,22,27,28,29,30),
  c(35,36,11,12,19,20,21,22,27,28,29,30), c(35,36,15,16,33,34,21,22,27,28,29,30),
  c(35,36,37,38,33,34,21,22,27,28,29,30)
)

run_bs <- function(taxonomy) {
  rl_data <- taxonomy %>%
    filter(TreeSp != "Ground", Rote.Liste != "n", Number != 0) %>%
    select(Trap, FHL, TreeSp) %>%
    distinct()
  tree_species <- c("Acer", "Fraxinus", "Quercus", "Quercus rubra", "Tilia", "Ulmus")
  results <- matrix(0, nrow = 25, ncol = 6, dimnames = list(NULL, tree_species))
  for (i in 1:25) {
    traps <- trap_combinations[[i]]
    rl_subset <- rl_data %>%
      filter(Trap %in% traps) %>%
      select(FHL, TreeSp) %>%
      distinct() %>%
      group_by(FHL) %>%
      mutate(n_treesp = n_distinct(TreeSp)) %>%
      ungroup() %>%
      filter(n_treesp == 1)
    for (sp in tree_species) results[i, sp] <- sum(rl_subset$TreeSp == sp)
  }
  colMeans(results)
}

count_unique_rl_direct <- function(taxonomy, sp) {
  taxonomy %>%
    filter(TreeSp != "Ground", Rote.Liste != "n", Number != 0) %>%
    select(Trap, FHL, TreeSp) %>%
    distinct() %>%
    group_by(FHL) %>%
    mutate(n_treesp = n_distinct(TreeSp)) %>%
    ungroup() %>%
    filter(n_treesp == 1, TreeSp == sp) %>%
    summarise(n = n_distinct(FHL)) %>%
    pull(n) %>% { if (length(.) == 0) 0L else . }
}

bs_2017 <- run_bs(taxonomy_2017)
bs_2018 <- run_bs(taxonomy_2018)

unique_rl_2017 <- c(
  count_unique_rl_direct(taxonomy_2017, "Acer"),
  bs_2017["Fraxinus"], bs_2017["Quercus"],
  count_unique_rl_direct(taxonomy_2017, "Quercus rubra"),
  bs_2017["Tilia"],
  count_unique_rl_direct(taxonomy_2017, "Ulmus")
)
unique_rl_2018 <- c(
  count_unique_rl_direct(taxonomy_2018, "Acer"),
  bs_2018["Fraxinus"], bs_2018["Quercus"],
  count_unique_rl_direct(taxonomy_2018, "Quercus rubra"),
  bs_2018["Tilia"],
  count_unique_rl_direct(taxonomy_2018, "Ulmus")
)
names(unique_rl_2017) <- order
names(unique_rl_2018) <- order

unique_rl_mean <- (unique_rl_2017 + unique_rl_2018) / 2

cat("Unique red list (mean 2017/2018):\n")
cat(sprintf("  %s: %.2f\n", order, unique_rl_mean))


# -----------------------------------------------------------
# 1.8  Derived rows
# -----------------------------------------------------------

non_unique_rl <- total_rl_mean - unique_rl_mean
unique_non_rl <- unique_total_mean - unique_rl_mean
common <- total_mean - unique_non_rl - non_unique_rl - unique_rl_mean


# -----------------------------------------------------------
# 1.9  Assemble and write Index Total Mean CES
# -----------------------------------------------------------

index_df <- data.frame(
  row.names = c(
    "total beetle species ",
    "unique non-red list beetle species",
    "Non-unique red list beetle species",
    "unique red list beetle species",
    "Common beetles "
  ),
  Acer = c(total_mean["Acer"], unique_non_rl["Acer"], non_unique_rl["Acer"],
           unique_rl_mean["Acer"], common["Acer"]),
  Fraxinus = c(total_mean["Fraxinus"], unique_non_rl["Fraxinus"], non_unique_rl["Fraxinus"],
               unique_rl_mean["Fraxinus"], common["Fraxinus"]),
  Quercus = c(total_mean["Quercus"], unique_non_rl["Quercus"], non_unique_rl["Quercus"],
              unique_rl_mean["Quercus"], common["Quercus"]),
  Quercus.r. = c(total_mean["Quercus rubra"], unique_non_rl["Quercus rubra"], non_unique_rl["Quercus rubra"],
                 unique_rl_mean["Quercus rubra"], common["Quercus rubra"]),
  Tilia = c(total_mean["Tilia"], unique_non_rl["Tilia"], non_unique_rl["Tilia"],
            unique_rl_mean["Tilia"], common["Tilia"]),
  Ulmus = c(total_mean["Ulmus"], unique_non_rl["Ulmus"], non_unique_rl["Ulmus"],
            unique_rl_mean["Ulmus"], common["Ulmus"])
)

cat("\n=== Index Total Mean CES (computed) ===\n")
print(index_df)
cat("Part 1 complete.\n")


#############################################################
# PART 2 — TREE BIODIVERSITY INDEX + FOREST SIMULATION + FBI
#############################################################

cat("\n\n=== PART 2: Tree Biodiversity Index & Forest Biodiversity Index ===\n")

# -----------------------------------------------------------
# 2.1  Read beetle input data (uses the computed index_df)
# -----------------------------------------------------------

beetle_data <- index_df

ces_data <- beetle_data[c(
  "unique red list beetle species",
  "Non-unique red list beetle species",
  "unique non-red list beetle species",
  "Common beetles "
), ]

rownames(ces_data) <- c(
  "unique_red_list",
  "non_unique_red_list",
  "unique_non_red_list",
  "common"
)

# -----------------------------------------------------------
# 2.2  Aggregation functions
# -----------------------------------------------------------

ces_general <- function(x, a, sigma) {
  rho <- (sigma - 1) / sigma
  (sum(a * (x ^ rho)))^(1 / rho)
}

ces_cobb_douglas <- function(x, a) {
  exp(sum(a * log(x)))
}

ces_function <- function(x, a, sigma) {
  if (abs(sigma - 1) < 1e-8) {
    ces_cobb_douglas(x, a)
  } else {
    ces_general(x, a, sigma)
  }
}

linear_index <- function(x, a) {
  sum(a * x)
}

# -----------------------------------------------------------
# 2.3  Weight scenarios
# -----------------------------------------------------------

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

# -----------------------------------------------------------
# 2.4  Compute BV_i (all weights, linear)
# -----------------------------------------------------------

results <- list()
for (w in names(weights_list)) {
  results[[paste0(w, "_linear")]] <- sapply(
    colnames(ces_data),
    function(tree) linear_index(ces_data[, tree], weights_list[[w]])
  )
}

ces_results_df <- as.data.frame(results)
rownames(ces_results_df) <- colnames(ces_data)

write.csv(
  round(ces_results_df, 3),
  "CES_results_only.csv"
)

cat("CES results written to CES_results_only.csv\n")
cat("\nCES results:\n")
print(round(ces_results_df, 3))

# -----------------------------------------------------------
# 2.5  Basal area from PPA simulations
# -----------------------------------------------------------

species_lookup <- tibble(
  sp = 1:8,
  species = c(
    "Acer pseudoplatanus", "Acer campestre", "Fraxinus Excelsior",
    "Carpinus betulus", "Acer platanoides", "Quercus robur",
    "Ulmus spp.", "Tilia spp."
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
  initial_state = c(
    "Dry","Dry","Dry",
    "Intermediate","Intermediate","Intermediate",
    "Moist","Moist","Moist"
  ),
  groundwater_change_m = c(
    0.0, 0.5, 1.0,
    0.0, 0.5, 1.0,
    0.0, 0.5, 1.0
  )
)

compute_basal_area <- function(file, init, delta_gw) {
  read_csv(file, show_col_types = FALSE) %>%
    filter(cl == 1, time %in% c(0, 25, 50, 75, 100)) %>%
    mutate(
      dbh_m = dbh / 100,
      basal_area_m2_ha = n * pi * (dbh_m / 2)^2
    ) %>%
    group_by(time, sp) %>%
    summarise(
      basal_area_m2_ha = sum(basal_area_m2_ha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(species_lookup, by = "sp") %>%
    mutate(
      initial_state = init,
      groundwater_change_m = delta_gw
    )
}

all_results <- pmap_dfr(
  list(scenarios$file, scenarios$initial_state, scenarios$groundwater_change_m),
  compute_basal_area
)

cat("Basal area computed from PPA simulations.\n")

# -----------------------------------------------------------
# 2.6  Forest Biodiversity Index (all weights x sigmas)
# -----------------------------------------------------------

species_mapping <- c(
  "Acer pseudoplatanus" = "Acer",
  "Fraxinus Excelsior"  = "Fraxinus",
  "Quercus robur"       = "Quercus",
  "Tilia spp."          = "Tilia",
  "Ulmus spp."          = "Ulmus"
)

valid_forest_species <- names(species_mapping)
valid_beetle_species <- unname(species_mapping)

forest_aggregate <- function(x, sigma) {
  a <- rep(1 / length(x), length(x))
  if (abs(sigma - 1) < 1e-8) {
    exp(sum(a * log(x)))
  } else {
    rho <- (sigma - 1) / sigma
    (sum(a * (x ^ rho)))^(1 / rho)
  }
}

sigma_scenarios <- tibble(
  aggregation = c("CES_sigma_0.4", "CES_sigma_0.7", "Cobb_Douglas_sigma_1"),
  sigma = c(0.4, 0.7, 1.0)
)

forest_files <- tibble(
  groundwater_scenario = c("no_increase", "0.5m_increase", "1m_increase"),
  gw_value = c(0.0, 0.5, 1.0)
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
    pivot_wider(names_from = time, values_from = basal_area_m2_ha, names_prefix = "year_") %>%
    column_to_rownames("species")

  rownames(forest_mean) <- species_mapping[rownames(forest_mean)]
  forest_mean <- forest_mean[valid_beetle_species, ]

  years <- as.numeric(gsub("year_", "", colnames(forest_mean)))

  for (w in paste0("W", 1:10, "_linear")) {
    BV_i <- ces_results_df[valid_beetle_species, w, drop = TRUE]

    for (i in 1:nrow(sigma_scenarios)) {
      agg <- sigma_scenarios$aggregation[i]
      sigma_val <- sigma_scenarios$sigma[i]

      vals <- sapply(forest_mean, function(BA_vec) {
        x_i <- BV_i * as.numeric(BA_vec)
        forest_aggregate(x_i, sigma_val)
      })

      forest_results[[length(forest_results) + 1]] <- tibble(
        groundwater_scenario = forest_files$groundwater_scenario[forest_files$gw_value == gw],
        year = years,
        weight_scenario = w,
        aggregation = agg,
        forest_biodiversity_index = vals
      )
    }
  }
}

forest_biodiversity_all <- bind_rows(forest_results)

write_csv(
  forest_biodiversity_all,
  "Forest_Biodiversity_Index_All_Scenarios_AllWeights_AllSigmas.csv"
)

cat("Forest Biodiversity Index computed and saved.\n")

# -----------------------------------------------------------
# 2.7  Human-readable CSV per sigma
# -----------------------------------------------------------

for (agg_name in unique(forest_biodiversity_all$aggregation)) {
  out <- forest_biodiversity_all %>%
    filter(aggregation == agg_name) %>%
    mutate(
      groundwater_scenario = recode(
        groundwater_scenario,
        "no_increase"   = "no_increase",
        "0.5m_increase" = "gw_0.5m",
        "1m_increase"   = "gw_1m"
      )
    ) %>%
    pivot_wider(names_from = groundwater_scenario, values_from = forest_biodiversity_index) %>%
    arrange(weight_scenario, year)

  write_csv(out, paste0("Forest_Biodiversity_", agg_name, ".csv"))
}

cat("Per-sigma CSV files written.\n")

# -----------------------------------------------------------
# 2.8  Basal area summary tables
# -----------------------------------------------------------

for (gw_label in c("no_increase", "0.5m_increase", "1m_increase")) {
  gw_val <- forest_files$gw_value[forest_files$groundwater_scenario == gw_label]

  ba_summary <- all_results %>%
    filter(groundwater_change_m == gw_val) %>%
    group_by(species, initial_state) %>%
    summarise(
      year_0 = sum(basal_area_m2_ha[time == 0]),
      year_25 = sum(basal_area_m2_ha[time == 25]),
      year_50 = sum(basal_area_m2_ha[time == 50]),
      year_75 = sum(basal_area_m2_ha[time == 75]),
      year_100 = sum(basal_area_m2_ha[time == 100]),
      .groups = "drop"
    )

  ba_mean <- ba_summary %>%
    group_by(species) %>%
    summarise(
      initial_state = "Mean",
      year_0 = mean(year_0),
      year_25 = mean(year_25),
      year_50 = mean(year_50),
      year_75 = mean(year_75),
      year_100 = mean(year_100),
      .groups = "drop"
    )

  ba_full <- bind_rows(ba_summary, ba_mean) %>% arrange(species, initial_state)
  write_csv(ba_full, paste0("basal_area_overstory_", gw_label, ".csv"))
}

cat("Basal area summary tables written.\n")


#############################################################
# PART 3 — SENSITIVITY ANALYSIS (Spearman rank correlation)
#############################################################

cat("\n\n=== PART 3: Sensitivity Analysis ===\n")

# -----------------------------------------------------------
# 3.1  Beetle index rank correlation
# -----------------------------------------------------------

beetle_ranks <- ces_results_df %>%
  mutate(across(everything(), rank, ties.method = "average"))

beetle_spearman <- cor(beetle_ranks, method = "spearman")

beetle_spearman_df <- beetle_spearman %>%
  round(3) %>%
  as.data.frame() %>%
  rownames_to_column(var = "weight_scenario")

write_csv(beetle_spearman_df, "Sensitivity_Beetle_Rank_Spearman_Matrix.csv")

cat("Beetle Spearman matrix written.\n")

# -----------------------------------------------------------
# 3.2  Forest index ranking
# -----------------------------------------------------------

rank_forest_csv <- function(input_file, output_file) {
  df <- read_csv(input_file, show_col_types = FALSE)
  df_ranked <- df %>%
    rowwise() %>%
    mutate(
      ranks = list(rank(-c(no_increase, gw_0.5m, gw_1m), ties.method = "average")),
      no_increase = ranks[[1]],
      gw_0.5m     = ranks[[2]],
      gw_1m       = ranks[[3]]
    ) %>%
    select(-ranks) %>%
    ungroup()
  write_csv(df_ranked, output_file)
}

rank_forest_csv("Forest_Biodiversity_CES_sigma_0.4.csv", "Forest_Biodiversity_CES_sigma_0.4_RANKS.csv")
rank_forest_csv("Forest_Biodiversity_CES_sigma_0.7.csv", "Forest_Biodiversity_CES_sigma_0.7_RANKS.csv")
rank_forest_csv("Forest_Biodiversity_Cobb_Douglas_sigma_1.csv", "Forest_Biodiversity_Cobb_Douglas_sigma_1_RANKS.csv")

cat("Forest rank CSVs written.\n")

# -----------------------------------------------------------
# 3.3  Forest rank stability & sigma dependence
# -----------------------------------------------------------

df04 <- read_csv("Forest_Biodiversity_CES_sigma_0.4.csv", show_col_types = FALSE)
df07 <- read_csv("Forest_Biodiversity_CES_sigma_0.7.csv", show_col_types = FALSE)
df1  <- read_csv("Forest_Biodiversity_Cobb_Douglas_sigma_1.csv", show_col_types = FALSE)

compare_sigma <- function(dfA, dfB) {
  out <- list()
  for (w in unique(dfA$weight_scenario)) {
    a <- dfA %>% filter(weight_scenario == w)
    b <- dfB %>% filter(weight_scenario == w)
    rho <- cor(
      c(a$no_increase, a$gw_0.5m, a$gw_1m),
      c(b$no_increase, b$gw_0.5m, b$gw_1m),
      method = "spearman"
    )
    out[[w]] <- rho
  }
  return(out)
}

rho_04_07 <- compare_sigma(df04, df07)
rho_04_1  <- compare_sigma(df04, df1)
rho_07_1  <- compare_sigma(df07, df1)

cat("\nSpearman rho (sigma 0.4 vs 0.7):\n")
print(unlist(rho_04_07))
cat("\nSpearman rho (sigma 0.4 vs 1):\n")
print(unlist(rho_04_1))
cat("\nSpearman rho (sigma 0.7 vs 1):\n")
print(unlist(rho_07_1))


#############################################################
# DONE
#############################################################

cat("\n\n========================================\n")
cat("ForestDiversityValue_LaesekeKulkarniNguyenVogt complete!\n")
cat("========================================\n")

