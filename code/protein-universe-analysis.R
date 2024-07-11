source("code/protein-universe-utils.R")

###########################################################################
##### Load data (will need to change paths based on whever data end up)#####
###########################################################################
# Load all files
dat <- readRDS("data/afdb_clusters.RDS")
afdb_taxonomy <- readRDS("data/afdb_cluster_taxonomy.RDS")
cluster_stats <- readRDS("data/afdb_cluster_stats.RDS")
tree <- readRDS("data/timetree_phylogeny_cleaned.RDS")
tree_taxonomy <- readRDS("data/timetree_taxonomy.RDS")
pdb <- readRDS("data/pdb_metadata.RDS")
pdb_taxonomy <- readRDS("data/pdb_taxonomy.RDS")
genome_stats <- readRDS("data/afdb_genome_size_stats.RDS")

####################################################
##### Figure 1: PDB + AFDB species distributions#####
####################################################
# Get overall AFDB species distribution
afdb_species <- table(dat$taxonomy_ID)

# Pull out ncbi IDs
afdb_names <- names(afdb_species)

# Convert to numeric
afdb_species <- as.numeric(afdb_species)

# Get overall PDB species distribution (need to replace commas in totals)
pdb_species <- as.numeric(gsub(",", "", pdb[, 2]))

# Get PDB names
pdb_names <- name2taxid(pdb$V1, out_type = "summary")$id

# PDB: circle packing plot
png("~/Desktop/pdb_circle_packing.png", width = 2400, height = 2400)


data <- data.frame(
  group = pdb_names,
  value = pdb_species
)
data$kingdom <- pdb_taxonomy$superkingdom[match(
  data$group,
  pdb_taxonomy$ncbi_id
)]
data <- data[data$kingdom %in% c("Archaea", "Bacteria", "Eukaryota"), ]
data <- data[order(data$kingdom), ]

packing <- circleProgressiveLayout(data$value,
  sizetype = "area"
)
data <- cbind(data, packing)
dat.gg <- circleLayoutVertices(packing,
  npoints = 50
)
dat.gg$value <- rep(data$value, each = 51)
dat.gg$kingdom <- rep(data$kingdom, each = 51)

cols <- all_colors[5:7]
names(cols) <- c("Archaea", "Bacteria", "Eukaryota")
cols <- cols[match(
  data$kingdom,
  names(cols)
)]

ggplot() +

  # Make the bubbles
  geom_polygon(
    data = dat.gg, aes(x,
      y,
      group = id,
      fill = as.factor(kingdom)
    ),
    linewidth = 0.2,
    colour = "black"
  ) +
  scale_fill_manual(values = cols) +
  scale_size_continuous(range = c(1, 4)) +

  # General theme:
  theme_void() +
  theme(legend.position = "none") +
  coord_equal()

dev.off()

# AFDB: circle packing plot
png("~/Desktop/afdb_circle_packing.png", width = 2400, height = 2400)
data <- data.frame(
  group = afdb_names,
  value = afdb_species
)
data$kingdom <- afdb_taxonomy$superkingdom[match(
  data$group,
  afdb_taxonomy$ncbi_id
)]
data <- data[data$kingdom %in% c("Archaea", "Bacteria", "Eukaryota"), ]
data <- data[order(data$kingdom), ]

packing <- circleProgressiveLayout(data$value,
  sizetype = "area"
)
data <- cbind(data, packing)
dat.gg <- circleLayoutVertices(packing,
  npoints = 50
)
dat.gg$value <- rep(data$value, each = 51)
dat.gg$kingdom <- rep(data$kingdom, each = 51)

cols <- all_colors[5:7]
names(cols) <- c("Archaea", "Bacteria", "Eukaryota")
cols <- cols[match(
  data$kingdom,
  names(cols)
)]

ggplot() +

  # Make the bubbles
  geom_polygon(
    data = dat.gg, aes(x,
      y,
      group = id,
      fill = as.factor(kingdom)
    ),
    linewidth = 0.2,
    colour = "black"
  ) +
  scale_fill_manual(values = cols) +
  scale_size_continuous(range = c(1, 4)) +

  # General theme:
  theme_void() +
  theme(legend.position = "none") +
  coord_equal()

dev.off()

## Pie charts
# PDB pie chart
par(mfrow = c(1, 2))
data <- data.frame(
  group = pdb_names,
  value = pdb_species
)
data$kingdom <- pdb_taxonomy$superkingdom[match(data$group, pdb_taxonomy$ncbi_id)]
data <- data[data$kingdom %in% c("Archaea", "Bacteria", "Eukaryota"), ]
data <- lapply(split(data, data$kingdom), function(x) sum(x$value))
pie(unlist(data),
  names(data),
  col = all_colors[5:7],
  border = FALSE
)

# AFDB pie chart
data <- data.frame(
  group = afdb_names,
  value = afdb_species
)
data$kingdom <- afdb_taxonomy$superkingdom[match(data$group, afdb_taxonomy$ncbi_id)]
data <- data[data$kingdom %in% c("Archaea", "Bacteria", "Eukaryota"), ]
data <- lapply(split(data, data$kingdom), function(x) sum(x$value))
pie(unlist(data),
  names(data),
  col = all_colors[5:7],
  border = FALSE
)

## AFDB and PDB cumulative distributions
plot(ecdf(afdb_species),
  verticals = TRUE,
  main = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  lwd = 2,
  xlab = "n proteins",
  ylab = "Cumulative proportion",
  col = all_colors[1],
  yaxt = "n"
)
axis(2,
  at = seq(0.0, 1.0, 0.2),
  labels = signif(seq(0.0, 1.0, 0.2), 2),
  las = 2,
  cex.axis = 1.5
)
lines(ecdf(pdb_species),
  verticals = TRUE,
  lwd = 2,
  col = all_colors[2],
  do.points = FALSE
)
text(rep(280000, 2),
  c(0.2, 0.1),
  cex = 1.5,
  adj = 1,
  c("AFDB", "PDB"),
  col = c(
    all_colors[1],
    all_colors[2]
  )
)

###############################################
##### Figure 2: AFDB taxonomic completeness#####
###############################################
# Calculate completeness via phylogenetic distance with 'clade_PD'
pds_fam <- clade_PD(tree_taxonomy,
  afdb_taxonomy,
  tree,
  taxa1 = "phylum",
  taxa2 = "family",
  verbose = TRUE,
  show.tip.label = FALSE
)

# Extract weighted phylogenetic distances (i.e. 'taxonomic completenes')
weighted_pds <- sapply(
  pds_fam$phylogenetic_distance,
  function(x) x$weighted_PD
))

# Split by kingdom
weighted_pds <- split(
  weighted_pds,
  afdb_taxonomy$superkingdom[match(
    names(weighted_pds),
    afdb_taxonomy$phylum
  )]
)

# Test if phylogenetic distance distributions differ across kingdoms via
# Kruskal-Wallis test
kruskal.test(weighted_pds)

# ID which kingdoms differ via Dunn's test
dunn.test::dunn.test(weighted_pds)

# Plot contmap and taxonomic completeness together
par(
  mfrow = c(1, 2),
  mar = c(8.1, 4.1, 6.1, 2.1)
)

# Plot contmap
plot(pds_fam$contmap,
  ftype = "off",
  mar = c(1.1, 4.1, 4.1, 2.1)
)

# Plot
par(mar = c(8.1, 4.1, 6.1, 2.1))
vioplot::vioplot(weighted_pds,
  col = all_colors[5:7],
  side = "right",
  ylab = "",
  xlab = "",
  yaxt = "n",
  font.main = 1,
  cex.main = 1.5,
  las = 2,
  ylim = c(0, 1),
  cex.axis = 1.5,
  cex.lab = 1.5
)
axis(2,
  seq(0, 1, 0.2),
  seq(0, 1, 0.2),
  cex.axis = 1.5,
  las = 2
)
title(
  ylab = "Taxonomic completeness",
  cex.lab = 1.5
)

stripchart(weighted_pds,
  col = all_colors[5:7],
  at = seq(0.8, (length(weighted_pds) - 1) + 0.8, 1),
  jitter = 0.1,
  method = "jitter",
  vertical = TRUE,
  cex = 1,
  pch = 20,
  add = TRUE
)

## Linear model predicting phylogenetic distance
# Generate vector of weighted phylogenetic distances (outcome variable)
w_pds <- unlist(lapply(
  pds_fam$phylogenetic_distance,
  function(x) x$weighted_PD
))

# Vector of total phylogenetic distance per phylum (predictor)
total_pds <- unlist(lapply(
  pds_fam$phylogenetic_distance,
  function(x) x$total_phylogenetic_distance
))

# Match weight phylogenetic distance order to total distance vector
w_pds <- w_pds[match(
  names(total_pds),
  names(w_pds)
)]

# Vector of n families within each phylum (predictor)
n_fams <- unlist(lapply(names(total_pds), function(x) {
  length(na.omit(unique(tree_taxonomy$family[tree_taxonomy$phylum == x])))
}))

# Model
mod <- lm(w_pds ~ n_fams + total_pds)

# Model summary
summary(mod)

################################################################################
##### Figure 3: Foldseek "representative proteins" phylogenetic distribution#####
################################################################################
# Get overal taxonomy of Foldseek representative proteins
representatives <- unique(dat$taxonomy_ID[match(
  cluster_stats$cluster_ID,
  dat$member_ID
)])

# Get NCBI species IDs for each representative protein
representatives_taxonomy <- afdb_taxonomy[match(
  representatives,
  afdb_taxonomy$ncbi_id
), ]

# Calculate the distribution across phyla of representative proteins
reps_phyla <- table(representatives_taxonomy$phylum)

# Calculate phyla distribution for all proteins
all_phyla <- table(afdb_taxonomy$phylum)

# Match orders
reps_phyla <- reps_phyla[match(
  names(all_phyla),
  names(reps_phyla)
)]
names(reps_phyla) <- names(all_phyla)

# Replace NAs with 0
reps_phyla[is.na(reps_phyla)] <- 0

# Calculate completeness via phylogenetic distance with 'clade_PD'
pds_rep <- clade_PD(tree_taxonomy,
  representatives_taxonomy,
  tree,
  taxa1 = "phylum",
  taxa2 = "family",
  verbose = TRUE,
  show.tip.label = FALSE
)

# Extract weighted phylogenetic distances
weighted_pds_all <- unlist(lapply(
  pds_fam$phylogenetic_distance,
  function(x) x$weighted_PD
))

weighted_pds_rep <- unlist(lapply(
  pds_rep$phylogenetic_distance,
  function(x) x$weighted_PD
))

# Generate plot
par(mfrow = c(2, 2))
plot(pds_fam$contmap,
  ftype = "off",
  mar = c(1.1, 4.1, 4.1, 2.1)
)
plot(pds_rep$contmap,
  ftype = "off",
  direction = "leftwards",
  mar = c(1.1, 4.1, 4.1, 2.1)
)

par(mar = c(8.1, 4.1, 6.1, 2.1))

# Plot joint distribution of phyla totals
plot(log(as.numeric(all_phyla) + 1),
  log(as.numeric(reps_phyla) + 1),
  pch = 20,
  col = "gray60",
  cex = 2,
  xlab = "AFDB",
  ylab = "Foldseek",
  cex.axis = 1.5,
  cex.lab = 1.5
)

# Add regression line
abline(
  lm(log(as.numeric(reps_phyla) + 1) ~
    log(as.numeric(all_phyla) + 1)),
  lty = "dashed",
  lwd = 1.5
)

# Add title
title(
  main = "n proteins per phylum (log)",
  font.main = 1,
  cex.main = 1.5
)

# Plot joint distribution of taxonomic completeness
plot(as.numeric(weighted_pds_all),
  as.numeric(weighted_pds_rep),
  pch = 20,
  col = "gray60",
  cex = 2,
  xlab = "AFDB",
  ylab = "Foldseek",
  cex.axis = 1.5,
  cex.lab = 1.5
)

# Make trait matrix
pgls_traits <- data.frame(
  Species = names(weighted_pds_all),
  all = as.numeric(weighted_pds_all),
  reps = as.numeric(weighted_pds_rep)
)

# Make comparative data
plgs_traits <- caper::comparative.data(pds_fam$new_tree,
  pgls_traits,
  Species,
  vcv = TRUE,
  vcv.dim = 2
)

# PGLS
mod1 <- caper::pgls(reps ~ all,
  plgs_traits,
  lambda = "ML"
)

# Add regression line
abline(
  lm(as.numeric(weighted_pds_rep) ~
    as.numeric(weighted_pds_all)),
  lty = "dashed",
  lwd = 1.5
)

# Add PGLS line
abline(mod1,
  lty = "dashed",
  lwd = 1.5,
  col = "darkred"
)

# Add title
title(
  main = "Taxonomic completeness",
  font.main = 1,
  cex.main = 1.5
)

########################################################################
##### Figure 4: Foldseek representative protein pLDDT distributions######
########################################################################
# Add ncbi id to cluster_stats (to be used to split and index on species)
cluster_stats$ncbi <- dat$taxonomy_ID[match(
  cluster_stats$cluster_ID,
  dat$member_ID
)]

# Get taxonomies associated with representative proteins
representatives <- dat$taxonomy_ID[match(
  cluster_stats$cluster_ID,
  dat$member_ID
)]

# Get taxonomic IDs for each
reps_taxonomy <- afdb_taxonomy[match(
  representatives,
  afdb_taxonomy$ncbi_id
), ]

# Add phylum and species names to cluster stats
cluster_stats$phylum <- reps_taxonomy$phylum[match(
  cluster_stats$ncbi,
  reps_taxonomy$ncbi_id
)]
cluster_stats$species <- reps_taxonomy$species[match(
  cluster_stats$ncbi,
  reps_taxonomy$ncbi_id
)]

# Calculate n proteins per species
species_n <- lapply(
  split(
    dat,
    dat$taxonomy_ID
  ),
  function(x) nrow(x)
)

# Split
s <- split(
  cluster_stats,
  cluster_stats$ncbi
)

# Match 'species_n' order to cluster statistics
species_n <- species_n[match(
  names(s),
  names(species_n)
)]

## Sweep protein n cutoffs and calculate correlation
# Counter
pb <- txtProgressBar(
  min = 1,
  max = 2500,
  style = 3,
  width = 100,
  char = "."
)

# Create empty vector for results
cors <- c()
for (i in 1:2500) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Filter species based on protein number
  s_filter <- s[unlist(lapply(
    s,
    function(x) {
      nrow(x)
    }
  )) >= i]

  # Calculate mean pLDDT for each species
  p <- unlist(lapply(
    s_filter,
    function(x) mean(x$repPlddt)
  ))

  # Create vector of protein n per species
  n <- unlist(species_n[match(
    names(p),
    names(species_n)
  )])

  # Correlate protein number and mean pLDDT
  cors <- c(
    cors,
    cor(n, p, method = "spearman")
  )
}

# Plot correlations
par(mfrow = c(1, 3))
plot(cors,
  ylab = "Spearman R",
  xlab = "n representative proteins per species",
  cex.lab = 1.5,
  cex.axis = 1.5,
  type = "l",
  lwd = 1.5,
  col = all_colors[1]
)

# Plot example
# Filter
s_filter <- s[unlist(lapply(
  s,
  function(x) nrow(x)
)) >= 1000]

# Calculate mean pLDDT for each species
p <- unlist(lapply(
  s_filter,
  function(x) mean(x$repPlddt)
))

# Create vector of protein n per species
n <- unlist(species_n[match(
  names(p),
  names(species_n)
)])

# Plot
cols <- all_colors[5:7]
names(cols) <- c("Archaea", "Bacteria", "Eukaryota")
cols <- cols[match(
  afdb_taxonomy$superkingdom[match(
    names(n),
    afdb_taxonomy$ncbi_id
  )],
  names(cols)
)]
plot(n,
  p,
  ylab = "pLDDT",
  xlab = "n representative proteins per species",
  cex.axis = 1.5,
  cex.lab = 1.5,
  pch = 20,
  col = cols
)

# Correlate protein number and mean pLDDT
cor(n,
  p,
  method = "spearman"
)

# Splt pLDDT vectors by kingdom
plddts <- split(
  p,
  afdb_taxonomy$superkingdom[match(
    names(n),
    afdb_taxonomy$ncbi_id
  )]
)

# Plot
vioplot::vioplot(plddts,
  col = all_colors[5:7],
  side = "right",
  ylab = "",
  xlab = "",
  font.main = 1,
  cex.main = 1.5,
  las = 2,
  cex.axis = 1.5,
  cex.lab = 1.5
)
title(
  ylab = "pLDDT (representative protein)",
  cex.lab = 1.5
)

stripchart(plddts,
  col = all_colors[5:7],
  at = seq(0.8, (length(plddts) - 1) + 0.8, 1),
  jitter = 0.1,
  method = "jitter",
  vertical = TRUE,
  cex = 1,
  pch = 20,
  add = TRUE
)

# Compare kingdom pLDDTs with Dunn's test
dunn.test::dunn.test(plddts)

##########################################################################
##### Figure 5: Effects of data balancing on phylogenetic distribution#####
##########################################################################
# Split afdb on species
dat_species <- split(
  dat,
  dat$taxonomy_ID
)

# Calculate species n
species_n <- unlist(lapply(
  dat_species,
  function(x) nrow(x)
))

## Calculate phylogenetic distance statistics over a range of protein n cutoffs
# Create empty list to save results
pds <- list()

# Set range of cutoffs to test
toTest <- c(1, seq(100, 20000, 100))

# Counter
pb <- txtProgressBar(
  min = 1,
  max = length(toTest),
  style = 3,
  width = 100,
  char = "."
)

# Loop over cutoffs and calculate stats
for (i in 1:length(toTest)) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Run function
  pds[[as.character(toTest[i])]] <- pd_per_partition(
    species_n,
    tree,
    afdb_taxonomy,
    toTest[i]
  )
}

# Extract partition size
partitions <- unlist(lapply(
  pds,
  function(x) x$PD[1, 1]
))

# Plot percentage of phylogenetic diversity captured
par(mfrow = c(1, 3))

# Set up color gradient to match partition size
cols <- colorRampPalette(unlist(
  arcadia_magma$color_dict
))(length(seq(0, 100, 0.1)))
names(cols) <- seq(0, 100, 0.1)
cols <- cols[match(
  round(partitions / max(partitions) * 100, 1),
  names(cols)
)]

# Plot
plot(names(pds),
  partitions / max(partitions) * 100,
  xlab = "Partition size (n proteins)",
  ylab = "Phylogenetic diversity (%)",
  cex.axis = 1.5,
  cex.lab = 1.5,
  col = cols,
  cex = 1.5,
  pch = 20,
  yaxt = "n",
  ylim = c(0, 100)
)
axis(2,
  seq(0, 100, 20),
  seq(0, 100, 20),
  cex.axis = 1.5,
  las = 2
)

## Plot species diversity curves for each taxonomic level
# Extract n unique species per taxonomic level
taxa_dist <- unlist(lapply(
  pds,
  function(x) length(unique(x$taxonomy$species))
))

# Get a unique color for each taxa
cols <- all_colors[(length(all_colors) - 7):length(all_colors)]

# Plot
plot(names(pds),
  (taxa_dist / max(taxa_dist)) * 100,
  xlab = "Partition size (n proteins)",
  ylab = "Taxa diveristy (%)",
  cex.axis = 1.5,
  cex.lab = 1.5,
  lwd = 2,
  type = "l",
  col = cols[1],
  yaxt = "n",
  ylim = c(0, 100)
)
axis(2, seq(0, 100, 20), seq(0, 100, 20), cex.axis = 1.5, las = 2)

# Add taxa
for (i in 2:6) {
  taxa_dist <- unlist(lapply(
    pds,
    function(x) length(unique(x$taxonomy[, i]))
  ))
  lines(names(pds),
    (taxa_dist / max(taxa_dist)) * 100,
    lwd = 2,
    col = cols[i]
  )
}

# Add taxa names
text(20000, 95, colnames(tree_taxonomy)[1], col = cols[1], adj = 1)
text(20000, 90, colnames(tree_taxonomy)[2], col = cols[2], adj = 1)
text(20000, 85, colnames(tree_taxonomy)[3], col = cols[3], adj = 1)
text(20000, 80, colnames(tree_taxonomy)[4], col = cols[4], adj = 1)
text(20000, 75, colnames(tree_taxonomy)[5], col = cols[5], adj = 1)
text(20000, 70, colnames(tree_taxonomy)[6], col = cols[6], adj = 1)

## Assess the effect of data balancing on cluster n
# Calculate total n of clusters within data set
clusters_total <- length(unique(dat$` cluster_ID`[dat$taxonomy_ID %in%
  pds[[1]]$taxonomy$ncbi_id]))

# Calculate  cluster n per partition
# Create empty vector to save results
cluster_space_partitions <- c()

# Counter
pb <- txtProgressBar(
  min = 1,
  max = length(pds),
  style = 3,
  width = 100,
  char = "."
)

# Loop over and calculate
for (i in 1:length(pds)) {
  # Update counter
  setTxtProgressBar(pb, i)

  # Calculate n clusters within partition
  n_clusters <- length(unique(dat$` cluster_ID`[dat$taxonomy_ID %in%
    pds[[i]]$taxonomy$ncbi_id]))

  # Normalize by total number of clusters ('clusters_total')
  n_clusters_norm <- n_clusters / clusters_total

  # Add to vector
  cluster_space_partitions <- c(
    cluster_space_partitions,
    n_clusters_norm
  )
}

# Convert to percentage
cluster_space_partitions <- cluster_space_partitions * 100

# Get colors
cols <- colorRampPalette(unlist(
  arcadia_magma$color_dict
))(length(seq(0, 100, 0.1)))
names(cols) <- seq(0, 100, 0.1)
cols <- cols[match(
  round(cluster_space_partitions, 1),
  names(cols)
)]

# Plot
plot(names(pds),
  cluster_space_partitions,
  xlab = "Partition size (n proteins)",
  ylab = "% of cluster space",
  cex.axis = 1.5,
  cex.lab = 1.5,
  col = cols,
  cex = 1.5,
  pch = 20,
  yaxt = "n",
  ylim = c(0, 100)
)
axis(2,
  seq(0, 100, 20),
  seq(0, 100, 20),
  cex.axis = 1.5,
  las = 2
)

###################################################################
##### Figure 6: Effects of data balancing on pLDDT distribution#####
###################################################################
# Add ncbi id to cluster_stats
cluster_stats$ncbi <- dat$taxonomy_ID[match(
  cluster_stats$cluster_ID,
  dat$member_ID
)]

# Get n proteins per species
species_n <- lapply(
  split(
    dat,
    dat$taxonomy_ID
  ),
  function(x) nrow(x)
)

# Match order to genome stats matrix
species_n <- unlist(species_n[match(
  genome_stats$ncbi_id,
  names(species_n)
)])

# Calculate ratio of n proteins in AFDB to protein coding genes for each species
protein_ratio <- species_n /
  genome_stats$proteins

# Sort
protein_ratio <- sort(protein_ratio)

# Calculate mean pLDDT per cluster
mean_pLDDT <- unlist(lapply(
  split(
    cluster_stats,
    cluster_stats$ncbi
  ),
  function(x) mean(x$repPlddt)
))

# Calculate n proteins per cluster
s <- split(
  cluster_stats,
  cluster_stats$ncbi
)

# Plot probability density functions at different protein n cutoffs
par(mfrow = c(1, 6))
for (i in c(1, 10, 100, 500, 1000, 2000)) {
  # Filter
  s_filter <- s[unlist(lapply(s, function(x) nrow(x))) >= i]

  # Compare pLDDT and protein ratio
  g <- Reduce(intersect, list(
    names(protein_ratio),
    names(mean_pLDDT),
    names(s_filter)
  ))

  # Get kingdoms
  kingdoms <- afdb_taxonomy$superkingdom[match(g, afdb_taxonomy$ncbi_id)]
  kingdoms[kingdoms == "Bacteria"] <- "Prokaryotes"
  kingdoms[kingdoms == "Archaea"] <- "Prokaryotes"

  # Make data matrix
  plddt_ratio <- data.frame(
    ratio = log(protein_ratio[match(
      g,
      names(protein_ratio)
    )]),
    plddt = mean_pLDDT[match(g, names(mean_pLDDT))],
    kingdom = kingdoms
  )

  # Filter out NAs
  plddt_ratio <- plddt_ratio[complete.cases(plddt_ratio), ]

  # Calculate pdfs
  euk_pdf <- kde2d(plddt_ratio$ratio[plddt_ratio$kingdom == "Eukaryota"],
    plddt_ratio$plddt[plddt_ratio$kingdom == "Eukaryota"],
    n = 200,
    lims = c(
      c(-10, 10),
      c(30, 100)
    )
  )

  pro_pdf <- kde2d(plddt_ratio$ratio[plddt_ratio$kingdom == "Prokaryotes"],
    plddt_ratio$plddt[plddt_ratio$kingdom == "Prokaryotes"],
    n = 200,
    lims = c(
      c(-10, 10),
      c(30, 100)
    )
  )

  # Plot as contours
  euk_pdf$z <- euk_pdf$z / max(euk_pdf$z)
  pro_pdf$z <- pro_pdf$z / max(pro_pdf$z)

  contour(euk_pdf$x,
    euk_pdf$y,
    euk_pdf$z,
    levels = seq(0, 1, 0.025),
    # levels = pretty(euk_pdf$z, 20),
    col = "#3B9886",
    drawlabels = FALSE,
    cex.lab = 1.5,
    cex.axis = 1.5,
    ylab = "pLDDT",
    xlab = "Protein ratio (log)",
    yaxt = "n"
  )
  contour(pro_pdf$x,
    pro_pdf$y,
    pro_pdf$z,
    levels = seq(0, 1, 0.025),
    add = TRUE,
    # levels = pretty(pro_pdf$z, 20),
    col = "#f898AE",
    drawlabels = FALSE
  )
  axis(2,
    seq(30, 90, 20),
    seq(30, 90, 20),
    cex.axis = 1.5,
    las = 2
  )
  abline(
    v = 1,
    lty = "dashed",
    lwd = 1.5
  )
  title(
    main = paste("min n proteins =", i, sep = " "),
    font.main = 1,
    cex.main = 1.5
  )
}