library(reticulate)
library(geiger)
library(caper)
pd <- import("pandas")

#Load timetree full tree
tree = read.tree('~/Desktop/protein_space_size_estimation/00_data/TimeTree v5 Final.nwk')
tree$tip.label = gsub("\\'", '', tree$tip.label)

#Root tree
tree = root(tree, tree$edge[which.max(tree$edge.length),2], resolve.root = TRUE)

# Get vector of species names
species <- tree$tip.label

# Get taxonomic ids
ids <- name2taxid(species, out_type = "summary")

#Generate tree taxonomy
tree_taxonomy <- getRawTaxonomy(ids$id)
names(tree_taxonomy) <- as.character(ids$name)

#Simplify tree taxonomy
v = list()
for(i in 1:length(tree_taxonomy)){
  
  species = tree_taxonomy[[i]][grep(paste("^", 'species', "$", sep = ""), names(tree_taxonomy[[i]]))]
  genus = tree_taxonomy[[i]][grep(paste("^", 'genus', "$", sep = ""), names(tree_taxonomy[[i]]))]
  family = tree_taxonomy[[i]][grep(paste("^", 'family', "$", sep = ""), names(tree_taxonomy[[i]]))]
  order = tree_taxonomy[[i]][grep(paste("^", 'order', "$", sep = ""), names(tree_taxonomy[[i]]))]
  class = tree_taxonomy[[i]][grep(paste("^", 'class', "$", sep = ""), names(tree_taxonomy[[i]]))]
  phylum = tree_taxonomy[[i]][grep(paste("^", 'phylum', "$", sep = ""), names(tree_taxonomy[[i]]))]

  l = list(species = species,
           genus = genus,
           family = family,
           order = order,
           class = class,
           phylum = phylum)
  l[unlist(lapply(l, function(x) length(x))) == 0] = NA
  names(l) = NULL
  
  v[[i]] = unlist(l)
  
}
v = as.data.frame(do.call(rbind, v))
v$species = gsub(' ', '_', v$species)

#Simplify afdb taxonomy
u = list()
for(i in 1:length(taxonomy)){
  
  species = taxonomy[[i]][grep(paste("^", 'species', "$", sep = ""), names(taxonomy[[i]]))]
  genus = taxonomy[[i]][grep(paste("^", 'genus', "$", sep = ""), names(taxonomy[[i]]))]
  family = taxonomy[[i]][grep(paste("^", 'family', "$", sep = ""), names(taxonomy[[i]]))]
  order = taxonomy[[i]][grep(paste("^", 'order', "$", sep = ""), names(taxonomy[[i]]))]
  class = taxonomy[[i]][grep(paste("^", 'class', "$", sep = ""), names(taxonomy[[i]]))]
  phylum = taxonomy[[i]][grep(paste("^", 'phylum', "$", sep = ""), names(taxonomy[[i]]))]
  ncbi_id = names(taxonomy)[i]
  
  l = list(species = species,
           genus = genus,
           family = family,
           order = order,
           class = class,
           phylum = phylum,
           ncbi_id = ncbi_id)
  l[unlist(lapply(l, function(x) length(x))) == 0] = NA
  names(l) = NULL
  
  u[[i]] = unlist(l)
  
}
colnames(u)[7] = 'ncbi_id'
u = as.data.frame(do.call(rbind, u))
u$species = gsub(' ', '_', u$species)

#Get unique phyla
phyla = na.omit(unique(v$class))

#Loop through and calculate PD
pds = list()
for(i in 1:length(phyla)){
  
  print(paste(i, 'out of', length(phyla)))
  #MRCA
  tips = v[grep(phyla[i], v$class),]
  tips = tips[!is.na(tips$species),]
  tips = tips[!is.na(tips$family),]
  tips = tips[tips$species%in%tree$tip.label,]
  tips = tips[!duplicated(tips$family),]
  
  if(nrow(tips)>2){
    mrca = getMRCA(tree, tips$species)

    #Extract clade
    clade = extract.clade(tree, mrca)
    
    #Simplify
    clade = keep.tip(clade, tips$species)
    clade$tip.label = tips$family[match(clade$tip.label, tips$species)]
    
    #Reclass
    tips = tips[match(clade$tip.label, tips$family),]
    
    #Get afdb in phyla
    a_tips = u[match(tips$family, u$family),]
    
    #Create community vector(s)
    comm = a_tips$family
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
             proportion_taxa = sum(comm[1,] == 1)/ncol(comm),
             weighted_PD = PD$PD[1]/PD$PD[2])
    pds[[phyla[[i]]]] = l
  }else{
    #Create results list
    l = list(phylogenetic_distance = 0,
             proportion_taxa = 0,
             weighted_PD = 0)
    pdb_pds[[phyla[[i]]]] = l
  }
}

#Plot
plot(unlist(lapply(pds, function(x) x$phylogenetic_distance$PD)))
plot(unlist(lapply(pds, function(x) x$proportion_taxa)))

plot(unlist(lapply(pds, function(x) x$phylogenetic_distance$PD)),
     unlist(lapply(pds, function(x) x$proportion_taxa)))

#Get phylum tree
tips = v
tips = tips[tips$species%in%tree$tip.label,]
tips = tips[!duplicated(tips$class),]
tips = tips[!is.na(tips$class),]
tips = tips[tips$class%in%names(pds),]
phyla_tree = keep.tip(tree, tips$species)
phyla_tree$tip.label = tips$class[match(phyla_tree$tip.label, tips$species)]

#Contmaps
par(mfrow = c(1,2))
phytools::contMap(phyla_tree, unlist(lapply(pds, function(x) x$proportion_taxa)))
phytools::contMap(phyla_tree, unlist(lapply(pds, function(x) x$weighted_PD)))

#Load cluster metadata (to get plddt)
cluster_stats = read.delim('foldseek_clusters/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv', header = FALSE)
colnames(cluster_stats) = c('cluster_ID', 'isDark', 'n_members', 'repLength', 'avgLength',
                            'repPlddt', 'avgPlddt', 'LCAtaxID')

#Get unique taxa
phyla = na.omit(unique(u$class))

#Loop through and calculate plddt for each
plddts = list()
for(i in 1:length(phyla)){
  
  print(paste(i, 'out of', length(phyla)))
  
  #Get species in taxa
  x = u[grep(phyla[i], u$class),]
  
  #Get all clusters associated with taxa
  clusters = dat$` cluster_ID`[dat$taxonomy_ID%in%x$ncbi_id]
  
  #Get plddts
  avg_plddt = cluster_stats$avgPlddt[cluster_stats$cluster_ID%in%clusters]
  rep_plddt = cluster_stats$repPlddt[cluster_stats$cluster_ID%in%clusters]
  avg_length = cluster_stats$avgLength[cluster_stats$cluster_ID%in%clusters]
  rep_length = cluster_stats$repPlddt[cluster_stats$cluster_ID%in%clusters]
  
  #Add to list
  plddts[[phyla[i]]] = list(avg_plddt = avg_plddt,
                            rep_plddt = rep_plddt,
                            avg_length = avg_length,
                            rep_length = rep_length)
  
}

#Compare to tree
trait = unlist(lapply(plddts, function(x) mean(x$rep_plddt)))
trait = trait[match(phyla_tree$tip.label, names(trait))]
trait = na.omit(trait)

#Contmap
plddt_tree = keep.tip(phyla_tree, names(trait))
phytools::contMap(plddt_tree, trait)

#Get taxa for LCA of each cluster
lca_taxonomy = getRawTaxonomy(cluster_stats$LCAtaxID)
lca_id = lapply(lca_taxonomy, function(x) x[1])
names(lca_id) = gsub(' ', '', names(lca_id))
lca_id = lca_id[!duplicated(names(lca_id))]
lca_id = unlist(lca_id)
names(lca_id) = unlist(lapply(strsplit(names(lca_id), '\\.'), function(x) x[1]))
lca_id = lca_id[match(cluster_stats$LCAtaxID, names(lca_id))]

cluster_stats$lca_name = lca_id

#Calculate mean plddt
cluster_stats_taxa = split(cluster_stats, cluster_stats$lca_name)
taxa_avg_plddt = unlist(lapply(cluster_stats_taxa, function(x) mean(x$avgPlddt)))
taxa_rep_plddt = unlist(lapply(cluster_stats_taxa, function(x) mean(x$repPlddt)))
n_lcas = unlist(lapply(cluster_stats_taxa, function(x) nrow(x)))

plot(log(n_lcas), taxa_avg_plddt)
plot(log(n_lcas), taxa_rep_plddt)

#Load PDB taxa (from https://www.rcsb.org/stats/explore/scientific_name_of_source_organism)
pdb = read.csv('~/Desktop/protein_space_size_estimation/00_data/pdb_taxa_distribution.csv', header = FALSE)

#Get NCBI ids
ids <- name2taxid(pdb$V1, out_type = "summary")

#Get taxonomy
pdb_taxonomy = getRawTaxonomy(ids$id)
names(pdb_taxonomy) = as.character(ids$id)

#Simplify pdb taxonomy
z = list()
for(i in 1:length(pdb_taxonomy)){
  
  species = pdb_taxonomy[[i]][grep(paste("^", 'species', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  genus = pdb_taxonomy[[i]][grep(paste("^", 'genus', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  family = pdb_taxonomy[[i]][grep(paste("^", 'family', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  order = pdb_taxonomy[[i]][grep(paste("^", 'order', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  class = pdb_taxonomy[[i]][grep(paste("^", 'class', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  phylum = pdb_taxonomy[[i]][grep(paste("^", 'phylum', "$", sep = ""), names(pdb_taxonomy[[i]]))]
  
  l = list(species = species,
           genus = genus,
           family = family,
           order = order,
           class = class,
           phylum = phylum)
  l[unlist(lapply(l, function(x) length(x))) == 0] = NA
  names(l) = NULL
  
  z[[i]] = unlist(l)
  
}
z = as.data.frame(do.call(rbind, z))
z$species = gsub(' ', '_', z$species)
z = na.omit(z)

#Calculate representation on tree
pdb_pds = list()
for(i in 1:length(phyla)){
  
  print(paste(i, 'out of', length(phyla)))
  #MRCA
  tips = u[grep(phyla[i], u$class),]
  tips = tips[!is.na(tips$species),]
  tips = tips[!is.na(tips$family),]
  tips = tips[tips$species%in%tree$tip.label,]
  tips = tips[!duplicated(tips$family),]
  
  if(nrow(tips)>2){
    mrca = getMRCA(tree, tips$species)
    
    #Extract clade
    clade = extract.clade(tree, mrca)
    
    #Simplify
    clade = keep.tip(clade, tips$species)
    clade$tip.label = tips$family[match(clade$tip.label, tips$species)]
    
    #Reclass
    tips = tips[match(clade$tip.label, tips$family),]
    
    #Get afdb in phyla
    a_tips = z[match(tips$family, z$family),]
    
    #Create community vector(s)
    comm = a_tips$family
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
             proportion_taxa = sum(comm[1,] == 1)/ncol(comm),
             weighted_PD = PD$PD[1]/PD$PD[2])
    pdb_pds[[phyla[[i]]]] = l
  }else{
    #Create results list
    l = list(phylogenetic_distance = 0,
             proportion_taxa = 0,
             weighted_PD = 0)
    pdb_pds[[phyla[[i]]]] = l
  }
}

#Get phylum tree
tips = v
tips = tips[tips$species%in%tree$tip.label,]
tips = tips[!duplicated(tips$class),]
tips = tips[!is.na(tips$class),]
tips = tips[tips$class%in%names(pdb_pds),]
phyla_tree = keep.tip(tree, tips$species)
phyla_tree$tip.label = tips$class[match(phyla_tree$tip.label, tips$species)]

w = pdb_pds[match(phyla_tree$tip.label, names(pdb_pds))]

#Contmaps
par(mfrow = c(1,2))
phytools::contMap(phyla_tree, unlist(lapply(w, function(x) x$proportion_taxa)))
phytools::contMap(phyla_tree, unlist(lapply(w, function(x) x$weighted_PD)))

##############
#####PGLS#####
##############
#Get common taxa
x = unlist(lapply(plddts, function(x) mean(x$avg_plddt)))
len = unlist(lapply(plddts, function(x) mean(x$avg_length, na.rm = TRUE)))
len = len[!len == 'NaN']
y = unlist(lapply(w, function(x) x$weighted_PD))
y[y>0] = 1
int = Reduce(intersect, list(names(x), names(y), names(len)))

#Make tree match data
test_tree = keep.tip(phyla_tree, int)
int = int[match(test_tree$tip.label, int)]
x = x[match(int, names(x))]
y = y[match(int, names(y))]
len = len[match(int, names(len))]

#Create data frame
datos = data.frame(plddt = x,
                   PD = y,
                   length = len,
                   Species = names(x))

#Create comparative data object
comp.data<-comparative.data(test_tree, 
                            datos, 
                            names.col="Species", 
                            vcv.dim=2, 
                            warn.dropped=TRUE)

#PGLS
mod <- pgls(plddt~PD+length, data=comp.data)

#Compare distributions
vioplot::vioplot(x[y == 0],
                 x[y == 1])




