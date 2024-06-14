library(ape)
library(phytools)
library(picante)
library(RColorBrewer)

# Function to simplify an NCBI taxonomy record per species
# Inputs:
#' taxa': NCBI IDs for a set of species
simplify_ncbi <- function(taxa) {
  # Create empty list to save results in
  z <- list()

  # Loop over and clean taxonomic data for each species
  for (i in 1:length(taxa)) {
    # Collect relevant taxa
    species <- taxa[[i]][grep(
      paste("^", "species", "$", sep = ""),
      names(taxa[[i]])
    )]
    genus <- taxa[[i]][grep(
      paste("^", "genus", "$", sep = ""),
      names(taxa[[i]])
    )]
    family <- taxa[[i]][grep(
      paste("^", "family", "$", sep = ""),
      names(taxa[[i]])
    )]
    order <- taxa[[i]][grep(
      paste("^", "order", "$", sep = ""),
      names(taxa[[i]])
    )]
    class <- taxa[[i]][grep(
      paste("^", "class", "$", sep = ""),
      names(taxa[[i]])
    )]
    phylum <- taxa[[i]][grep(
      paste("^", "phylum", "$", sep = ""),
      names(taxa[[i]])
    )]
    kingdom <- taxa[[i]][grep(
      paste("^", "kingdom", "$", sep = ""),
      names(taxa[[i]])
    )]
    superkingdom <- taxa[[i]][grep(
      paste("^", "superkingdom", "$", sep = ""),
      names(taxa[[i]])
    )]
    ncbi_id <- names(taxa)[i]

    # Combine into list
    l <- list(
      species = species,
      genus = genus,
      family = family,
      order = order,
      class = class,
      phylum = phylum,
      kingdom = kingdom,
      superkingdom = superkingdom,
      ncbi_id = ncbi_id
    )

    # Replace empty elements with NA
    l[unlist(lapply(l, function(x) length(x))) == 0] <- NA

    # Replace names
    names(l) <- NULL
    l <- unlist(l)
    names(l) <- c(
      "species",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "kingdom",
      "superkingdom",
      "ncbi_id"
    )

    # Add to results lists
    z[[i]] <- l
  }

  # Combine list into data frame
  z <- as.data.frame(do.call(rbind, z))

  # Replace space with underscore in species names
  # (to keep consistent with other data structure)
  z$species <- gsub(" ", "_", z$species)

  # Add column names to data frame
  colnames(z) <- c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "superkingdom",
    "ncbi_id"
  )

  # Return data frame
  z
}

# Function to calculate per-clade phylogenetic distance statistics
# Inputs:
#' tree_taxonomy': taxonomic IDs of species within provided phylogeny
#' focal_taxonomy': taxonomic IDs of species within database (e.g. AFDB, etc.)
#' tree': phylogeny of species of interest
#' taxa1' = lower-resolution taxonomic grouping
# to compare within (e.g. phylum or class)
#' taxa2' = higher-resolution taxonomic grouping
# to measure phylogenetic distance of (e.g. genus or family)
clade_PD <- function(tree_taxonomy,
                     focal_taxonomy,
                     tree,
                     taxa1 = "class",
                     taxa2 = "family",
                     verbose = FALSE,
                     plot = FALSE,
                     ...) {
  # Get unique entries for taxonomic unit
  col1 <- grep(
    taxa1,
    colnames(tree_taxonomy)
  )
  col2 <- grep(
    taxa2,
    colnames(tree_taxonomy)
  )
  taxa <- na.omit(unique(tree_taxonomy[, col1]))

  # Loop through and calculate PD
  pds <- list()

  # Counter if desired
  if (verbose == TRUE) {
    pb <- txtProgressBar(
      min = 1,
      max = length(taxa),
      style = 3,
      width = 50,
      char = "."
    )
  }

  # Loop over and calculate taxonomic completeness for each taxa
  for (i in 1:length(taxa)) {
    # Update counter if desired
    if (verbose == TRUE) {
      setTxtProgressBar(pb, i)
    }

    # ID tree tips associated with taxa
    tips <- tree_taxonomy[grep(
      taxa[i],
      tree_taxonomy[, col1]
    ), ]

    # Remove NAs
    tips <- tips[!is.na(tips$species), ]
    tips <- tips[!is.na(tips[, col2]), ]

    # Match to tree tip order
    tips <- tips[tips$species %in% tree$tip.label, ]

    # Remove duplicated entries (i.e. keep only 1 tip per sub-taxa)
    tips <- tips[!duplicated(tips[, col2]), ]

    # Calculate taxonomic completeness if n>2 taxa
    if (nrow(tips) > 2) {
      # ID most recent common ancestor (MRCA) of taxa in tree
      mrca <- getMRCA(tree, tips$species)

      # Use the MRCA to ID species in the relevant clade
      clade <- extract.clade(tree, mrca)

      # Filter tree to be just clade
      clade <- keep.tip(clade, tips$species)

      # Replace species names with taxa
      clade$tip.label <- tips[, col2][match(
        clade$tip.label,
        tips$species
      )]

      # Reclass
      tips <- tips[match(clade$tip.label, tips[, col2]), ]

      # Get afdb in taxa
      a_tips <- focal_taxonomy[match(
        tips[, col2],
        focal_taxonomy[, col2]
      ), ]

      # Create community vector(s)
      comm <- a_tips[, col2]
      comm[!is.na(comm)] <- 1
      comm[is.na(comm)] <- 0
      comm <- as.data.frame(t(comm))
      comm <- rbind(comm, rep(1, ncol(comm)))
      colnames(comm) <- clade$tip.label
      rownames(comm) <- 1:nrow(comm)

      # Phylogenetic alpha diversity
      PD <- pd(
        samp = comm,
        tree = clade,
        include.root = TRUE
      )

      # Create results list
      l <- list(
        phylogenetic_distance = PD$PD[1],
        total_phylogenetic_distance = PD$PD[2],
        proportion_taxa = sum(comm[1, ] == 1) / ncol(comm),
        weighted_PD = PD$PD[1] / PD$PD[2]
      )
      pds[[taxa[[i]]]] <- l
    } else {
      # Create results list
      l <- list(
        phylogenetic_distance = 0,
        proportion_taxa = 0,
        weighted_PD = 0
      )
      pds[[taxa[[i]]]] <- l
    }
  }

  # Generate new tree that matches pd data
  tips <- tree_taxonomy
  tips <- tips[tips$species %in% tree$tip.label, ]
  tips <- tips[tips[, col1] %in% names(pds), ]
  tips <- tips[!duplicated(tips[, col1]), ]

  phyla_tree <- keep.tip(
    tree,
    tips$species
  )
  phyla_tree$tip.label <- tips[, col1][match(
    phyla_tree$tip.label,
    tips$species
  )]

  pds <- pds[match(phyla_tree$tip.label, names(pds))]

  # Calculate contmap
  contmap <- phytools::contMap(phyla_tree,
    unlist(lapply(pds, function(x) x$weighted_PD)),
    plot = FALSE
  )
  contmap$cols[1:length(contmap$cols)] <- colorRampPalette(brewer.pal(11, "RdYlBu"))(length(contmap$cols))

  # Plot if desired
  if (plot == TRUE) {
    plot(
      contmap,
      ...
    )

    title(main = paste(taxa1, "->", taxa2), font.main = 1, cex.main = 1)
  }

  # Return results
  l <- list(
    phylogenetic_distance = pds,
    new_tree = phyla_tree,
    contmap = contmap
  )
  l
}

# Set up magma gradient for plots
arcadia_magma <- list(
  color_dict = list(
    "arcadia:black" = "#09090A",
    "arcadia:grape" = "#5A4596",
    "arcadia:taffy" = "#E87485",
    "arcadia:orange" = "#FFB984",
    "arcadia:oat" = "#F5E4BE"
  ),
  values = c(0, 0.38, 0.72, 0.9, 1)
)

# Function to calculate PD for species with at least n proteins in AFDB
# Inputs:
# 'species_numbers': vector of protein n per species
# 'tree': phylogeny corresponding to species of interest
# 'taxonomy': taxonomic IDs associated with species of interest
# 'parition_size': minimum n per species to filter on
pd_per_partition <- function(species_numbers,
                             tree,
                             taxonomy,
                             partition_size = 1) {
  # Filter species
  species <- species_numbers[species_numbers >= partition_size]

  # Filter to just species in tree
  species <- taxonomy$species[match(
    names(species),
    taxonomy$ncbi_id
  )]
  species <- species[species %in% tree$tip.label]

  # ID MRCA
  mrca <- getMRCA(
    tree,
    species
  )

  # Extract clade
  clade <- extract.clade(
    tree,
    mrca
  )

  # Filter
  clade <- keep.tip(
    clade,
    species
  )

  # Get afdb in taxa
  a_tips <- clade$tip.label %in% species

  # Create community vector(s)
  comm <- a_tips
  comm[comm == TRUE] <- 1
  comm[comm == FALSE] <- 0
  comm <- as.data.frame(t(comm))
  colnames(comm) <- clade$tip.label
  rownames(comm) <- 1:nrow(comm)

  # Calculate phylogenetic diversity
  PD <- pd(
    samp = comm,
    tree = clade,
    include.root = TRUE
  )

  # Return
  return(list(
    PD = PD,
    taxonomy = taxonomy[match(
      species,
      taxonomy$species
    ), ]
  ))
}

# Set up color dicitionary for plots
all_colors <- c(
  "#5088C5", "#F28360", "#F7B846", "#97CD78",
  "#7A77AB", "#f898AE", "#3B9886", "#c85152",
  "#73B5E3", "#BAB0A8", "#8A99AD", "#FFB984",
  "#C6E7F4", "#F8C5C1", "#F5E4BE", "#B5BEA4",
  "#DCBFFC", "#B6C8D4", "#DAD3C7", "#DA9085"
)
