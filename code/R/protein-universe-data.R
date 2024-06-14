source("code/R/protein-universe-utils.R")

#######################
##### Download data#####
#######################
# Download data from Foldseek
system("wget https://afdb-cluster.steineggerlab.workers.dev/1-AFDBClusters-entryId_repId_taxId.tsv.gz")
system("wget https://afdb-cluster.steineggerlab.workers.dev/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz")

# Prepare NCBI taxonomy accession database
prepareDatabase(getAccessions = FALSE)

###################################################
##### Organize data + collect relevant metadata#####
###################################################
##### AFDB IDs + cluster information#####
# Load AFDB IDs + foldseek cluster information
dat <- read.delim("1-AFDBClusters-entryId_repId_taxId.tsv",
  header = FALSE
)

# Add column names
colnames(dat) <- c(
  "member_ID",
  " cluster_ID",
  "taxonomy_ID"
)

## Taxonomy
# ID all unique taxa in AFDB
taxa <- unique(dat$taxonomy_ID)

# Collect taxonomy
afdb_taxonomy <- getRawTaxonomy(taxa)

# Change names
names(afdb_taxonomy) <- as.character(taxa)

# Simplify taxonomy
afdb_taxonomy <- simplify_ncbi(afdb_taxonomy)

##### Foldseek cluster statistics#####
# Load cluster statistics
cluster_stats <- read.delim("2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv",
  header = FALSE
)

# Add column names
colnames(cluster_stats) <- c(
  "cluster_ID",
  "isDark", "
                            n_members",
  "repLength",
  "avgLength",
  "repPlddt",
  "avgPlddt",
  "LCAtaxID"
)

##### Timetree phylogeny#####
# Load timetree
tree <- read.tree("00_data/TimeTree v5 Final.nwk")

# Clean tip labes
tree$tip.label <- gsub("\\'", "", tree$tip.label)

# Root tree on longest branch
tree <- root(tree,
  tree$edge[which.max(tree$edge.length), 2],
  resolve.root = TRUE
)

## Taxonomy
# Get taxonomic ids
ids <- name2taxid(tree$tip.label,
  out_type = "summary"
)

# Collect taxonomy
tree_taxonomy <- getRawTaxonomy(ids$id)

# Change names
names(tree_taxonomy) <- as.character(ids$id)

# Simplify taxonomy
tree_taxonomy <- simplify_ncbi(tree_taxonomy)

##### PDB species information#####
# Load PDB data (www.rcsb.org/stats/explore/scientific_name_of_source_organism)
pdb <- read.csv("00_data/pdb_taxa_distribution.csv",
  header = FALSE
)

# Get taxonomic ids
ids <- name2taxid(pdb$V1, out_type = "summary")

# Collect taxonomy
pdb_taxonomy <- getRawTaxonomy(ids$id)

# Change names
names(pdb_taxonomy) <- as.character(ids$id)

# Simplify taxonomy
pdb_taxonomy <- simplify_ncbi(pdb_taxonomy)

##### NCBI genome statistics#####
# Download genome statistics from NCBI for Eukaryotes
euk_stats <- genomes::reports("eukaryotes.txt")

# Download genome statistics from NCBI for Prokaryotes
prok_stats <- genomes::reports("prokaryotes.txt")

# Simplify
euk_stats <- data.frame(
  organism = euk_stats$Organism,
  ncbi_id = euk_stats$TaxID,
  size = euk_stats$`Size (Mb)`,
  gc = euk_stats$`GC%`,
  genes = euk_stats$Genes,
  proteins = euk_stats$Proteins
)

prok_stats <- data.frame(
  organism = prok_stats$Organism,
  ncbi_id = prok_stats$TaxID,
  size = prok_stats$`Size (Mb)`,
  gc = prok_stats$`GC%`,
  genes = prok_stats$Genes,
  proteins = prok_stats$Proteins
)

# Combine
genome_stats <- rbind(euk_stats, prok_stats)

# Split on species
g <- split(
  genome_stats,
  genome_stats$TaxID
)

# Calculate mean statistics
for (i in 1:length(g)) {
  g[[i]] <- data.frame(
    organism = g[[i]][1, 1],
    ncbi_id = g[[i]][1, 2],
    size = mean(g[[i]]$size, na.rm = TRUE),
    gc = mean(g[[i]]$gc, na.rm = TRUE),
    genes = round(mean(g[[i]]$genes, na.rm = TRUE)),
    proteins = round(mean(g[[i]]$proteins, na.rm = TRUE))
  )
}

# Filter to afdb
genome_stats <- g[names(g) %in% as.numeric(afdb_taxonomy$ncbi_id)]

# Recombine
genome_stats <- do.call(rbind, genome_stats)

###########################
##### Save cleaned data#####
###########################
saveRDS(dat, "afdb_clusters.RDS")
saveRDS(afdb_taxonomy, "afdb_cluster_taxonomy.RDS")
saveRDS(cluster_stats, "afdb_cluster_stats.RDS")
saveRDS(tree, "timetree_phylogeny_cleaned.RDS")
saveRDS(tree_taxonomy, "timetree_taxonomy.RDS")
saveRDS(pdb, "pdb_metadata.RDS")
saveRDS(pdb_taxonomy, "pdb_taxonomy.RDS")
saveRDS(genome_stats, "afdb_genome_size_stats.RDS")
