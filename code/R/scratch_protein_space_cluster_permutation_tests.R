#Simplify tree to just phyla
tips = tips[tips$species%in%tree$tip.label,]
tips = tips[tips$phylum%in%names(effects),]
tips = tips[!duplicated(tips$phylum),]
phyla_tree = keep.tip(tree, tips$species)
phyla_tree$tip.label = tips$phylum[match(phyla_tree$tip.label, tips$species)]

#Calculate cophenetic distance
cophen = cophenetic(phyla_tree)

#Convert to three column matrix
cophen = as.data.frame(as.table(cophen))
cophen$n_clusters = rep(NA, nrow(cophen))

#Split on phyla
p = afdb_taxonomy$phylum[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id)]
s = split(dat, p)
s = s[names(s)%in%phyla_tree$tip.label]

#Calculate average cluster n per phylum via permutation
n_perms = 100
n_proteins = 100
pb <- txtProgressBar(min = 1,
                     max = nrow(cophen),
                     style = 3,
                     width = 100,
                     char = ".")
for(i in 1:nrow(cophen)){
  setTxtProgressBar(pb, i)
  
  m = c(s[grep(cophen[i,1], names(s))][[1]][,3],
        s[grep(cophen[i,2], names(s))][[1]][,3])
  if(length(m)>n_proteins){
    perms = c()
    for(k in 1:n_perms){
      perms = c(perms, length(unique(m[sample(1:length(m), n_proteins)])))
    }
    cophen$n_clusters[i] = mean(perms)
  }else{
    cophen$n_clusters[i] = NA
  }
}

#Calculate mean
cluster_means = split(cophen, cophen$Freq)
cluster_means = lapply(cluster_means, function(x) mean(x$n_clusters))







##Calculate per-species cluster diversity (n clusters/n proteins)
#Bin by taxa
x = na.omit(unique(afdb_taxonomy$phylum))
divs = list()
for(i in 1:length(x)){
  
  #Get species in taxon
  species = unique(afdb_taxonomy$ncbi_id[afdb_taxonomy$phylum == x[i]])
  
  #Get associated cluster diversities
  divs[[x[i]]] = cluster_diversity[names(cluster_diversity)%in%species]
  
}

#Remove empty elements
divs = divs[unlist(lapply(divs, function(x) length(x)))>0]

#Match species in tree and dat 
tips = tree_taxonomy
tips = tips[tips$species%in%tree$tip.label,]
int = intersect(tips$ncbi_id, dat$taxonomy_ID)
f = dat[dat$taxonomy_ID%in%int,]

#Split cluster data by species
s = split(f, afdb_taxonomy$phylum[match(f$taxonomy_ID, afdb_taxonomy$ncbi_id)])

#Filter?
#s = s[unlist(lapply(s, function(x) nrow(x)))>=100]

#Calculate cluster diversity for each
cluster_diversity = unlist(lapply(s, function(x) length(unique(x[,2]))/nrow(x)))

#Permutation test
n = 100
permutation_results = list()

pb <- txtProgressBar(min = 1,
                     max = length(s),
                     style = 3,
                     width = 100,
                     char = ".")

for(i in 1:length(s)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  #Calculate observed clusters per taxonomic distance
  n_proteins = nrow(f[f$taxonomy_ID%in%s[[i]]$taxonomy_ID,])
  n_clusters = length(unique(f$` cluster_ID`[f$taxonomy_ID%in%s[[i]]$taxonomy_ID]))
  
  #ID species in random selection
  ncbi = unique(s[[i]]$taxonomy_ID)
  
  #Calculate MRCA
  species = tree_taxonomy$species[tree_taxonomy$ncbi_id%in%ncbi]
  mrca = getMRCA(tree, species)
  
  #Extract clade
  clade = extract.clade(tree, mrca)
  clade = keep.tip(clade, species)
  
  #Create community vector
  comm = clade$tip.label
  comm = comm%in%species
  comm[comm == TRUE] = 1
  comm[comm == FALSE] = 0
  comm = as.data.frame(t(comm))
  comm = rbind(comm, rep(1, ncol(comm)))
  colnames(comm) = clade$tip.label
  rownames(comm) = 1:nrow(comm)
  
  #Phylogenetic alpha diversity
  PD = pd(samp = comm, 
          tree = clade,
          include.root = TRUE)
  
  #Calculate normalized cluster n
  obs = n_clusters/PD$PD[1]
  
  perms = c()
  for(j in 1:n){

    #Get random proteins
    random = f[sample(1:nrow(f), n_proteins),]
    
    #Calculate n clusters
    n_clusters = length(unique(random$` cluster_ID`))
    
    #ID species in random selection
    ncbi = unique(random$taxonomy_ID)
    
    #Calculate MRCA
    species = tree_taxonomy$species[tree_taxonomy$ncbi_id%in%ncbi]
    mrca = getMRCA(tree, species)
    
    #Extract clade
    clade = extract.clade(tree, mrca)
    clade = keep.tip(clade, species)
    
    #Create community vector
    comm = clade$tip.label
    comm = comm%in%species
    comm[comm == TRUE] = 1
    comm[comm == FALSE] = 0
    comm = as.data.frame(t(comm))
    comm = rbind(comm, rep(1, ncol(comm)))
    colnames(comm) = clade$tip.label
    rownames(comm) = 1:nrow(comm)
    
    #Phylogenetic alpha diversity
    PD = pd(samp = comm, 
            tree = clade,
            include.root = TRUE)
    
    #Add to results
    perms = c(perms, n_clusters/PD$PD[1])
  }
  
  p = sum(obs<perms)/n
  effect_size = obs-mean(perms)
  effect_size_norm = effect_size/sd(perms)
  permutation_results[[names(s)[i]]] = list(p_value = p,
                                            effect_size = effect_size,
                                            effect_size_norm = effect_size_norm)
  
}


