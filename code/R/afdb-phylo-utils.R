library(taxonomizr)
library(taxizedb)
library(ape)
library(picante)
library(vegan)
library(ArcadiaColorBrewer)
library(scales)
library(alluvial)
library(phytools)
library(ggiraph)
library(packcircles)
library(ggplot2)
library(pbapply)

#Function to simplify an NCBI taxonomy record per species
simplify_ncbi = function(taxa){
  z = list()
  for(i in 1:length(taxa)){
    
    species = taxa[[i]][grep(paste("^", 'species', "$", sep = ""), names(taxa[[i]]))]
    genus = taxa[[i]][grep(paste("^", 'genus', "$", sep = ""), names(taxa[[i]]))]
    family = taxa[[i]][grep(paste("^", 'family', "$", sep = ""), names(taxa[[i]]))]
    order = taxa[[i]][grep(paste("^", 'order', "$", sep = ""), names(taxa[[i]]))]
    class = taxa[[i]][grep(paste("^", 'class', "$", sep = ""), names(taxa[[i]]))]
    phylum = taxa[[i]][grep(paste("^", 'phylum', "$", sep = ""), names(taxa[[i]]))]
    kingdom = taxa[[i]][grep(paste("^", 'kingdom', "$", sep = ""), names(taxa[[i]]))]
    superkingdom = taxa[[i]][grep(paste("^", 'superkingdom', "$", sep = ""), names(taxa[[i]]))]
    ncbi_id = names(taxa)[i]
    
    l = list(species = species,
             genus = genus,
             family = family,
             order = order,
             class = class,
             phylum = phylum,
             kingdom = kingdom,
             superkingdom = superkingdom,
             ncbi_id = ncbi_id)
    l[unlist(lapply(l, function(x) length(x))) == 0] = NA
    names(l) = NULL
    l = unlist(l)
    names(l) = c('species',
                 'genus',
                 'family',
                 'order',
                 'class',
                 'phylum',
                 'kingdom',
                 'superkingdom',
                 'ncbi_id')

    z[[i]] = l
    
  }
  z = as.data.frame(do.call(rbind, z))
  z$species = gsub(' ', '_', z$species)

  colnames(z) = c('species',
                  'genus',
                  'family',
                  'order',
                  'class',
                  'phylum',
                  'kingdom',
                  'superkingdom',
                  'ncbi_id')
  
  z
}

#Function to calculate per-clade phylogenetic distance stats
#'taxa1' = lower-resolution taxonomic grouping to compare within (e.g. phylum or class)
#'taxa2' = higher-resolution taxonomic grouping to measure phylogenetic distance of (e.g. genus or family)
clade_PD = function(tree_taxonomy,
                    focal_taxonomy,
                    tree,
                    taxa1 = 'class',
                    taxa2 = 'family',
                    verbose = FALSE,
                    plot = FALSE,
                    ...){
  
  #Get unique entries for taxonomic unit
  col1 = grep(taxa1, colnames(tree_taxonomy))
  col2 = grep(taxa2, colnames(tree_taxonomy))
  taxa = na.omit(unique(tree_taxonomy[,col1]))
  
  #Loop through and calculate PD
  pds = list()
  
  #Counter if desired
  if(verbose == TRUE){
    pb <- txtProgressBar(
      min = 1,
      max = length(taxa),
      style = 3,
      width = 50,
      char = "."
    )
  }
  for(i in 1:length(taxa)){
    
    if(verbose == TRUE){
      setTxtProgressBar(pb, i)
    }

    #MRCA
    tips = tree_taxonomy[grep(taxa[i], tree_taxonomy[,col1]),]
    tips = tips[!is.na(tips$species),]
    tips = tips[!is.na(tips[,col2]),]
    tips = tips[tips$species%in%tree$tip.label,]
    tips = tips[!duplicated(tips[,col2]),]
    
    if(nrow(tips)>2){
      mrca = getMRCA(tree, tips$species)
      
      #Extract clade
      clade = extract.clade(tree, mrca)
      
      #Simplify
      clade = keep.tip(clade, tips$species)
      clade$tip.label = tips$family[match(clade$tip.label, tips$species)]
      
      #Reclass
      tips = tips[match(clade$tip.label, tips[,col2]),]
      
      #Get afdb in taxa
      a_tips = focal_taxonomy[match(tips[,col2], focal_taxonomy[,col2]),]
      
      #Create community vector(s)
      comm = a_tips[,col2]
      comm[!is.na(comm)] = 1
      comm[is.na(comm)] = 0
      comm = as.data.frame(t(comm))
      comm = rbind(comm, rep(1, ncol(comm)))
      colnames(comm) = clade$tip.label
      rownames(comm) = 1:nrow(comm)
      
      #Phylogenetic alpha diversity
      PD = pd(samp = comm, 
              tree = clade,
              include.root = TRUE)
      
      #Create results list
      l = list(phylogenetic_distance = PD$PD[1],
               total_phylogenetic_distance = PD$PD[2],
               proportion_taxa = sum(comm[1,] == 1)/ncol(comm),
               weighted_PD = PD$PD[1]/PD$PD[2])
      pds[[taxa[[i]]]] = l
    }else{
      #Create results list
      l = list(phylogenetic_distance = 0,
               proportion_taxa = 0,
               weighted_PD = 0)
      pds[[taxa[[i]]]] = l
    }
  }
  
  #Generate new tree that matches pd data
  tips = tree_taxonomy
  tips = tips[tips$species%in%tree$tip.label,]
  tips = tips[tips[,col1]%in%names(pds),]
  tips = tips[!duplicated(tips[,col1]),]
  
  phyla_tree = keep.tip(tree, tips$species)
  phyla_tree$tip.label = tips[,col1][match(phyla_tree$tip.label, tips$species)]
  
  pds = pds[match(phyla_tree$tip.label, names(pds))]
  
  #Calculate contmap
  contmap = phytools::contMap(phyla_tree, 
                              unlist(lapply(pds, function(x) x$weighted_PD)),
                              plot = FALSE)
  contmap$cols[1:length(contmap$cols)] = colorRampPalette(RColorBrewer::brewer.pal(11, 'RdYlBu'))(length(contmap$cols))
  
  #Plot if desired
  if(plot == TRUE){
    plot(contmap,
         ...)
    
    title(main = paste(taxa1, '->', taxa2), font.main = 1, cex.main = 1)
  }

  #Return results
  l = list(phylogenetic_distance = pds,
           new_tree = phyla_tree,
           contmap = contmap)
  l
}
  
# test = clade_PD(tree_taxonomy,
#                 afdb_taxonomy,
#                 tree,
#                 taxa1 = 'phylum',
#                 verbose = TRUE)

# Function to filter a tree to a given taxonomic resolution
# Function to filter a tree to a given taxonomic resolution
filter_tree <- function(input_tree,
                        resolution = c("species", "genus", "family", "order", "class", "phylum", "kingdom"),
                        verbose = FALSE,
                        plot = FALSE,
                        change_tree_labels = FALSE,
                        return_taxonomy = FALSE,
                        root_tree = FALSE,
                        ...) {
  # Get vector of species names
  species <- input_tree$tip.label
  
  # Get taxonomic ids
  if (verbose == TRUE) {
    print("grabbing taxonomic ids")
  }
  ids <- name2taxid(species, out_type = "summary")
  
  # Prepare accession database
  prepareDatabase(getAccessions = FALSE)
  
  # Get taxocomic groups for each species
  if (verbose == TRUE) {
    print("grabbing taxonomic categories")
  }
  taxonomy <- getRawTaxonomy(ids$id)
  names(taxonomy) <- as.character(ids$name)
  
  # Get unique entries for taxa of interest
  u <- unlist(lapply(taxonomy, function(x) x[grep(paste("^", resolution, "$", sep = ""), names(x))]))
  names(u) <- gsub(paste(".", resolution, sep = ""), "", names(u))
  
  # Get list of species within in taxa
  u <- split(names(u), u)
  
  # Randomly choose species from each taxa
  if (verbose == TRUE) {
    print("filtering tree")
  }
  u <- lapply(u, function(x) x[sample(1:length(x), 1)])
  
  # Subset tree
  tree_new <- drop.tip(input_tree, 
                       input_tree$tip.label[!input_tree$tip.label %in% unlist(u)],
                       rooted = root_tree)
  
  # Change tip labels
  if (change_tree_labels == TRUE) {
    tree_new$tip.label <- names(u)[match(unlist(u), tree_new$tip.label)]
  }
  
  # Plot if desired
  if (plot == TRUE) {
    plot(
      tree_new,
      ...
    )
  }
  
  # Return result
  if(return_taxonomy == TRUE){
    l = list(tree = tree_new, 
             taxonomy = taxonomy)
    return(l)
  }else{
    return(tree_new)
  }
}


# Function to calculate dispersion/distance (spread) of species as a function of phylogeny
# Inputs:
# A phylogeny, preferably rooted ('input_tree')
# A binary vector of species to test the dispersion of ('included_species')
# Function to calculate dispersion/distance (spread) of species as a function of phylogeny
# Inputs:
# A phylogeny, preferably rooted ('input_tree')
# A binary vector of species to test the dispersion of ('included_species')
# Function to calculate dispersion/distance (spread) of species as a function of phylogeny
# Inputs:
# A phylogeny, preferably rooted ('input_tree')
# A binary vector of species to test the dispersion of ('included_species')
tree_spread <- function(input_tree,
                        inlcuded_species) {
  # Root tree if necessary
  if (is.rooted(input_tree) == FALSE) {
    print("Tree not rooted; rooting on first species")
    input_tree <- root(input_tree, 1, resolve.root = TRUE)
  }
  
  # Calculate cophenic distance
  cophen <- cophenetic.phylo(input_tree)
  
  # Calculate phylogenetic alpha diversity
  PD <- pd(
    samp = comm,
    tree = input_tree,
    include.root = TRUE
  )
  
  # Calculate standard effect sizes
  ses <- ses.mpd(comm,
                 cophen,
                 null.model = "taxa.labels",
                 abundance.weighted = FALSE,
                 runs = 100
  )
  
  # Return
  l <- list(cophen, PD, ses)
  names(l) <- c("cophenic_distance", "phylogenetic_diversity", "ses.mpd")
  return(l)
}

#Function to load time tree data structures
load_tt_data = function(file){
  
  #Load
  tmp = read.delim(file, 
                   header = FALSE)
  
  #Split on '='
  tmp = data.frame(time = as.numeric(unlist(lapply(strsplit(tmp[,1], '='), function(x) x[1]))),
                   value = as.numeric(unlist(lapply(strsplit(tmp[,1], '='), function(x) x[2]))))
  
  #Return
  return(tmp)
  
}

