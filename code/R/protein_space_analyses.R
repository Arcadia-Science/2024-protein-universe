source('~/Documents/Research/github/afdb-phylo/01_code/R/afdb-phylo-utils.R')

#Model cluster forecasting as a function of taxa, cluster diversity, mean protein length, plddt, etc.
#cluster_forecasting ~ taxonomy + cluster_diversity + mean_protein_length + plddt + phylogenetic_uniqueness
#GLMM?

#To do: Calculate cluster number by number of SPECIES included in space
#To do: calculate taxa (phyla?) specific rates of novel protein acquisition, use this as basis of estimating how big protein structural space might be...
#Add in estimated sizes of taxa vs. amount measured in afdb? 
#Add in phylogenetic distances of phyla? 

#Get slope function

##########################################################
#####Clean data and get taxonomy for each protein set#####
##########################################################
#Set working directory
setwd('~/Desktop/protein_space_size_estimation/')
    
#####AFDB data#####
#Load data
dat = read.delim('00_data/foldseek_clusters/1-AFDBClusters-entryId_repId_taxId.tsv', header = FALSE)
colnames(dat) = c('member_ID', ' cluster_ID', 'taxonomy_ID')

#Get all unique taxa
taxa = unique(dat$taxonomy_ID)

#Prepare accession database
prepareDatabase(getAccessions=FALSE)

#Collecting taxonomy for each taxa
afdb_taxonomy = getRawTaxonomy(taxa)
names(afdb_taxonomy) = as.character(taxa)

#Simplify taxonomy
afdb_taxonomy = simplify_ncbi(afdb_taxonomy)

#####AFDB foldseek cluster metadata#####
#Load cluster metadata (to get plddt)
cluster_stats = read.delim('00_data/foldseek_clusters/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv', header = FALSE)
colnames(cluster_stats) = c('cluster_ID', 
                            'isDark', '
                            n_members', 
                            'repLength', 
                            'avgLength',
                            'repPlddt', 
                            'avgPlddt', 
                            'LCAtaxID')

#####Timetree phylogeny#####
#Load timetree full tree
tree = read.tree('00_data/TimeTree v5 Final.nwk')
tree$tip.label = gsub("\\'", '', tree$tip.label)

#Root tree
tree = root(tree, tree$edge[which.max(tree$edge.length),2], resolve.root = TRUE)

# Get taxonomic ids
ids <- name2taxid(tree$tip.label, out_type = "summary")

#Generate tree taxonomy
tree_taxonomy <- getRawTaxonomy(ids$id)
names(tree_taxonomy) <- as.character(ids$id)

#Simplify taxonomy
tree_taxonomy = simplify_ncbi(tree_taxonomy)

######PDB database metadata#####
#Load PDB taxa (from https://www.rcsb.org/stats/explore/scientific_name_of_source_organism)
pdb = read.csv('00_data/pdb_taxa_distribution.csv', header = FALSE)

#Get NCBI ids
ids <- name2taxid(pdb$V1, out_type = "summary")

#Get taxonomy
pdb_taxonomy = getRawTaxonomy(ids$id)
names(pdb_taxonomy) = as.character(ids$id)

#Simplify taxonomy
pdb_taxonomy = simplify_ncbi(pdb_taxonomy)

#####Genome statistics from NCBI#####
#Download genome statistics from NCBI
euk_stats = genomes::reports("eukaryotes.txt")
prok_stats = genomes::reports("prokaryotes.txt")

#Simplify
euk_stats = data.frame(organism = euk_stats$Organism,
                       ncbi_id = euk_stats$TaxID,
                       size = euk_stats$`Size (Mb)`,
                       gc = euk_stats$`GC%`,
                       genes = euk_stats$Genes,
                       proteins = euk_stats$Proteins)

prok_stats = data.frame(organism = prok_stats$Organism,
                        ncbi_id = prok_stats$TaxID,
                        size = prok_stats$`Size (Mb)`,
                        gc = prok_stats$`GC%`,
                        genes = prok_stats$Genes,
                        proteins = prok_stats$Proteins)

#Combine
genome_stats = rbind(euk_stats, prok_stats)

#Split on species
g = split(genome_stats, genome_stats$TaxID)

#Calculate mean statistics
for(i in 1:length(g)){
  g[[i]] = data.frame(organism = g[[i]][1,1],
                      ncbi_id = g[[i]][1,2],
                      size = mean(g[[i]]$size, na.rm = TRUE),
                      gc = mean(g[[i]]$gc, na.rm = TRUE),
                      genes = round(mean(g[[i]]$genes, na.rm = TRUE)),
                      proteins = round(mean(g[[i]]$proteins, na.rm = TRUE)))}

#Filter to afdb
genome_stats = g[names(g)%in%as.numeric(afdb_taxonomy$ncbi_id)]
genome_stats = do.call(rbind, genome_stats)

#####Save#####
# saveRDS(list(afdb_data = dat,
#              afdb_taxonomy = afdb_taxonomy,
#              cluster_stats = cluster_stats,
#              tree = tree,
#              tree_taxonomy = tree_taxonomy,
#              pdb = pdb,
#              pdb_taxonomy = pdb_taxonomy),
#         '02_output/protein_space_workspace_data.RDS')
# saveRDS(dat, '02_output/afdb_clusters.RDS')
# saveRDS(afdb_taxonomy, '02_output/afdb_cluster_taxonomy.RDS')
# saveRDS(cluster_stats, '02_output/afdb_cluster_stats.RDS')
# saveRDS(tree, '02_output/timetree_phylogeny_cleaned.RDS')
# saveRDS(tree_taxonomy, '02_output/timetree_taxonomy.RDS')
# saveRDS(pdb, '02_output/pdb_metadata.RDS')
# saveRDS(pdb_taxonomy, '02_output/pdb_taxonomy.RDS')
# saveRDS(genome_stats, '02_output/afdb_genome_size_stats.RDS')

###################
#####Load data#####
###################
#Set working directory
setwd('~/Desktop/protein_space_size_estimation/')

#Load all .RDS files
dat = readRDS('02_output/afdb_clusters.RDS')
afdb_taxonomy = readRDS('02_output/afdb_cluster_taxonomy.RDS')
cluster_stats = readRDS('02_output/afdb_cluster_stats.RDS')
tree = readRDS('02_output/timetree_phylogeny_cleaned.RDS')
tree_taxonomy = readRDS('02_output/timetree_taxonomy.RDS')
pdb = readRDS('02_output/pdb_metadata.RDS')
pdb_taxonomy = readRDS('02_output/pdb_taxonomy.RDS')
genome_stats = readRDS('02_output/afdb_genome_size_stats.RDS')

##################################################
#####Sankey diagrams of cluster stats by taxa#####
##################################################
#Get combinations of taxonomic categories
tit = data.frame(superkingdom = as.factor(afdb_taxonomy$superkingdom), 
                 phylum = as.factor(afdb_taxonomy$phylum))
tit = tit[!duplicated(tit$phylum),]
tit = tit[-is.na(tit$phylum),]

p = split(dat, afdb_taxonomy$phylum[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id)])
p = p[match(tit$phylum, names(p))]
tit$Freq = unlist(lapply(p, function(x) length(unique(x[,2]))))

tit = tit[tit$Freq>0,]

#Plot
alluvial(tit[,1:2], freq=tit$Freq, blocks = TRUE)

alluvial( tit[,1:2], freq=tit$Freq, xw=0.0, alpha=0.8,
          gap.width=0.1, col= "steelblue", border="white")

########################################################################################
#####Proportion of phylogenetic distance per taxonomic unit covered by superkingdom#####
########################################################################################
##Key result: Bacteria are significantly more completely sampled than Eukaryotes and Archaea in afdb
#Calculate phylogenetic distance as a function taxonomic resolution
toTest = c('superkingdom', 'phylum', 'class')
par(mfrow = c(1, length(toTest)))

pds_fam = list()
for(i in 1:length(toTest)){
  print(toTest[i])
  pds_fam[[toTest[i]]] = clade_PD(tree_taxonomy,
                                  afdb_taxonomy,
                                  tree,
                                  taxa1 = toTest[i],
                                  taxa2 = 'family',
                                  verbose = TRUE,
                                  show.tip.label = FALSE)
  plot(pds_fam[[i]]$contmap, ftype = 'off')
  title(main = toTest[i], font.main = 1, cex.main = 1.5)
}

#Save
saveRDS(pds_fam, '02_output/afdb_phylogenetic_distance_by_taxon.RDS')

#Compare proportion of PD covered by phylum and class
toTest = c('phylum', 'class')
par(mfrow = c(1,2))
for(i in 1:length(toTest)){
  
  #Extract weighted phylogenetic distances
  weighted_pds = unlist(lapply(pds_fam[[grep(toTest[i], names(pds_fam))]]$phylogenetic_distance, function(x) x$weighted_PD))
  
  #Split on kingdom
  weighted_pds = split(weighted_pds, afdb_taxonomy$superkingdom[match(names(weighted_pds), afdb_taxonomy[,grep(toTest[i], colnames(afdb_taxonomy))])])
  
  #Test significance
  kruskal.test(weighted_pds)
  dunn.test::dunn.test(weighted_pds)
  
  #Plot
  vioplot::vioplot(weighted_pds,
                   col = arcadia.pal(n = 3, name = 'Accent'),
                   side = "right",
                   ylab = '% PD represented per taxa',
                   xlab = '', 
                   font.main = 1,
                   cex.main = 1.5,
                   las = 2, 
                   ylim = c(0, 1),
                   cex.axis = 1.5,
                   cex.lab = 1.5)
  
  stripchart(weighted_pds,
             col = arcadia.pal(n = 3, name = 'Accent'),
             at = seq(0.8, (length(weighted_pds)-1)+0.8, 1), 
             jitter = 0.1,
             method = "jitter", 
             vertical = TRUE, 
             cex = 1,
             pch = 20, 
             add = TRUE)
  
  title(main = toTest[i], font.main = 1, cex.main = 1.5)
  
}

###############################################################
#####Amount of protein space covered by taxonomic distance#####
###############################################################
#Set up taxa to measure space coverage over
toTest = c("superkingdom", "phylum", 'class', 'order', 'family', 'genus', 'species')

#Calculate total number of clusters
total_clusters = length(unique(dat$` cluster_ID`))

#Match afdb_taxonomy to dat
afdb_long = afdb_taxonomy[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id),]

#Remove viruses from data?
dat2 = dat[-grep("Viruses", afdb_long$superkingdom),]
afdb_long = afdb_long[-grep("Viruses", afdb_long$superkingdom),]

#Calculate coverage based on taxonomic group
taxon_coverage = list()
for(i in 1:length(toTest)){
  
  print(paste(i, 'out of', length(toTest)))
  #Split data on taxonomic group
  d = split(dat2, afdb_taxonomy[,grep(toTest[i], colnames(afdb_taxonomy))][match(dat2$taxonomy_ID, afdb_taxonomy$ncbi_id)])
  
  #Count n clusters per taxa
  cl = unlist(lapply(d, function(x) length(unique(x[,2]))/total_clusters))
  
  #Add to list
  taxon_coverage[[toTest[i]]] = cl
}

#Extract superkingdoms
king = split(afdb_taxonomy, afdb_taxonomy$superkingdom)
king = king[1:3]
taxon_coverage_kingdom = list()
for(i in 1:length(king)){
  
  res = list()
  for(j in 1:length(toTest)){
    #ID taxa 
    taxa = unique(king[[i]][,grep(toTest[j], colnames(king[[i]]))])
    
    #Get coverage
    cov = taxon_coverage[[j]][names(taxon_coverage[[j]])%in%taxa]
    
    #Add to list
    res[[toTest[j]]] = cov
  }
  #Add to list
  taxon_coverage_kingdom[[names(king)[i]]] = res
}

#Compile coverage into a list
all = list(all = taxon_coverage,
           archaea = taxon_coverage_kingdom$Archaea,
           bacteria = taxon_coverage_kingdom$Bacteria,
           eukaryota = taxon_coverage_kingdom$Eukaryota)

#Normalize by total percentage
for(i in 2:4){
  all[[i]] = lapply(all[[i]], function(x) x/all[[i]][[1]])
}

#Plot
par(mfrow = c(1,4))
for(i in 1:length(all)){
  
  ses = unlist(lapply(all[[i]], function(x) sd(x)/sqrt(length(x))))
  means = unlist(lapply(all[[i]], function(x) mean(x)))
  
  plot(1:length(means),
       means,
       type = 'l',
       ylim = c(0, 1),
       ylab = '% of protein space',
       xlab = 'Taxonomic resolution',
       cex.lab = 1.5,
       cex.axis = 1.5)
  polygon(c(1:length(means), rev(1:length(means))),
          c(means+ses, rev(means-ses)),
          col = 'grey90',
          border = NA)
  lines(1:length(means),
        means,
        lwd = 1.5,
        col = 'grey50')
  title(main = names(all)[i],
        font.main = 1,
        cex.main = 1.5)
}

#Half violin plots
par(mfrow = c(1,4))
for(i in 1:length(all)){
  vioplot::vioplot(all[[i]][2:7],
                   col = 'gray80',
                   side = "right",
                   ylab = '% of protein space',
                   xlab = '', 
                   font.main = 1,
                   cex.main = 1.5,
                   las = 2, 
                   ylim = c(0, 0.3),
                   cex.axis = 1.5,
                   cex.lab = 1.5)
  
  stripchart(all[[i]][2:7],
             col = alpha('gray80', 0.5),
             at = seq(0.8, (length(all[[i]][2:7])-1)+0.8, 1), 
             jitter = 0.1,
             method = "jitter", 
             vertical = TRUE, 
             cex = 1,
             pch = 20, 
             add = TRUE)
}

#Set up taxa to measure space coverage over
toTest = c('genus', 'family', 'order', 'class', 'phylum', 'superkingdom')

#Cumulative coverage
clusters_kingdom = split(dat2, afdb_taxonomy[,grep('superkingdom', colnames(afdb_taxonomy))][match(dat2$taxonomy_ID, afdb_taxonomy$ncbi_id)])
afdb_long_kingdom = split(afdb_long, afdb_taxonomy[,grep('superkingdom', colnames(afdb_taxonomy))][match(dat2$taxonomy_ID, afdb_taxonomy$ncbi_id)])
cumulative_coverage = list()
for(i in 1:length(clusters_kingdom)){
  
  print(paste(i, 'out of', length(clusters_kingdom)))
  
  phyla_perms = list()
  for(j in 1:(length(toTest)-1)){
    
    print(toTest[j])
    
    tmp1 = split(clusters_kingdom[[i]], afdb_long_kingdom[[i]][,grep(toTest[j], colnames(afdb_long_kingdom[[i]]))])
    tmp1 = lapply(tmp1, function(x) unique(x[,2]))
    
    #Filter
    #tmp1 = tmp1[unlist(lapply(tmp1, function(x) length(x)))>1000]
    #if(length(tmp1)>1000){
    #  tmp1 = tmp1[sample(1:length(tmp1), 1000)]
    #}
    
    pb <- txtProgressBar(min = 1,
                         max = length(tmp1),
                         style = 3,
                         width = 100,
                         char = ".")
    
    perms = c()
    for(k in 1:length(tmp1)){
      setTxtProgressBar(pb, k)
      perms = c(perms, sum(tmp1[[k]]%in%unique(unlist(tmp1[-k]))==FALSE))
    }
    names(perms) = names(tmp1)
      
  phyla_perms[[toTest[j]]] = perms
    
  }
  cumulative_coverage[[names(clusters_kingdom)[i]]] = phyla_perms
}

#Plot
par(mfrow = c(1,3))
cols = arcadia.pal(n = 3, name = 'Accent')
for(i in 1:length(cumulative_coverage)){
  
  ses = unlist(lapply(cumulative_coverage[[i]], function(x) sd(x)/sqrt(length(x))))
  means = unlist(lapply(cumulative_coverage[[i]], function(x) mean(x)))
  
  plot(1:length(means),
       means,
       type = 'l',
       ylim = c(0, 12000),
       ylab = 'n unique clusters per taxon',
       xlab = '',
       col = cols[i],
       cex.lab = 1.5,
       xaxt = 'n',
       cex.axis = 1.5,
       las = 2)
  axis(1, 
       1:length(means), 
       names(means), 
       cex.axis = 1.5,
       las = 2) 
  polygon(c(1:length(means), rev(1:length(means))),
          c(means+ses, rev(means-ses)),
          col = alpha(cols[i], 0.5),
          border = NA)
  lines(1:length(means),
        means,
        lwd = 1.5,
        col = cols[i])
  title(main = names(cumulative_coverage)[i],
        font.main = 1,
        cex.main = 1.5)
}

#Compare to n training structures in each
toTest = c('genus', 'family', 'order', 'class', 'phylum')
res = list()
for(i in 1:length(toTest)){
  print(i)
  taxa = toTest[i]
  
  p = table(pdb_taxonomy[,grep(taxa, colnames(pdb_taxonomy))])
  z = c(cumulative_coverage$Archaea[[grep(taxa, names(cumulative_coverage$Archaea))]],
        cumulative_coverage$Bacteria[[grep(taxa, names(cumulative_coverage$Bacteria))]],
        cumulative_coverage$Eukaryota[[grep(taxa, names(cumulative_coverage$Eukaryota))]])
  p = p[match(names(z), names(p))]
  p[is.na(p)] = 0
  
  res[[toTest[i]]] = list(pdb_n = p,
                          cluster_n = z,
                          cor = cor(as.numeric(p), as.numeric(z)))
}

#Plot
par(mfrow = c(1,2))
cols = arcadia.pal(n = length(res), 'Accent')
plot(unlist(lapply(res, function(x) x$cor)),
     ylim = c(0,1),
     xlab = '',
     xaxt = 'n',
     col = cols,
     ylab = 'Correlation',
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 2,
     pch = 20,
     cex = 2)
axis(1, 1:length(res), names(res), cex.axis = 1.5, las = 2)

plot(as.numeric(res$phylum$pdb_n),
     as.numeric(res$phylum$cluster_n),
     pch = 20,
     col = alpha(cols[length(res)], 0.5),
     ylab = 'n private clusters',
     xlab = 'n members in pdb',
     cex.lab = 1.5,
     las = 2,
     cex = 2,
     cex.axis = 1.5)
abline(lm(as.numeric(res$phylum$cluster_n)~as.numeric(res$phylum$pdb_n)),
       lty = 'dashed',
       col = 'gray80')
title(main = 'Phylum', font.main = 1, cex.main = 1.5)

#Compare to distance from a training structure
toTest = c('genus', 'family', 'order', 'class', 'phylum')
for(i in 1:length(toTest)){
  print(i)
  taxa = toTest[i]
  
  p = table(pdb_taxonomy[,grep(taxa, colnames(pdb_taxonomy))])
  z = c(cumulative_coverage$Archaea[[grep(taxa, names(cumulative_coverage$Archaea))]],
        cumulative_coverage$Bacteria[[grep(taxa, names(cumulative_coverage$Bacteria))]],
        cumulative_coverage$Eukaryota[[grep(taxa, names(cumulative_coverage$Eukaryota))]])
  g = tree_taxonomy[tree_taxonomy[,grep(taxa, colnames(tree_taxonomy))]%in%names(z),]
  g = split(g, g[,grep(taxa, colnames(g))])
  g = lapply(g, function(x) x[!is.na(x$species),])
  g = lapply(g, function(x) x[x$species%in%tree$tip.label,])
  g = lapply(g, function(x) x$species[1])
  s = unlist(g)
  s = na.omit(s)
  phyla_tree = keep.tip(tree, s)
  phyla_tree$tip.label = names(s)[match(phyla_tree$tip.label, s)]
  cophen = cophenetic(phyla_tree)
  cophen = cophen[,colnames(cophen)%in%names(p)]
  cophen = apply(cophen, 1, function(x) min(x))
  
  plot(cophen,
       z[match(names(cophen), names(z))])
  print(cor(cophen,
      z[match(names(cophen), names(z))]))
  
}

#Contmap
cum = c(cumulative_coverage$Archaea$family,
        cumulative_coverage$Bacteria$family,
        cumulative_coverage$Eukaryota$family)
cum = log(cum)
cum[cum == '-Inf'] = 0

g = tree_taxonomy[tree_taxonomy$family%in%names(cum),]
g = split(g, g$family)
g = lapply(g, function(x) x[!is.na(x$species),])
g = lapply(g, function(x) x[x$species%in%tree$tip.label,])
g = lapply(g, function(x) x$species[1])
s = unlist(g)
s = na.omit(s)
phyla_tree = keep.tip(tree, s)
phyla_tree$tip.label = names(s)[match(phyla_tree$tip.label, s)]

cum = cum[match(names(s), names(cum))]
cum = cum[match(names(cum), phyla_tree$tip.label)]

contmap = phytools::contMap(phyla_tree, 
                            cum, 
                            plot = FALSE)

contmap$cols[1:length(contmap$cols)] = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'RdYlBu'))(length(contmap$cols)))
plot(contmap)

##############################################################################
#####Unique proteins/phyla protein diversity as a function of genome size#####
##############################################################################
#Split afdb data on species
dat_species = split(dat, dat$taxonomy_ID)

#Filter
dat_species_f = dat_species[unlist(lapply(dat_species, function(x) nrow(x)))>=5]

#Get intersection of genome size data species and afdb species
int = intersect(genome_stats$ncbi_id[!is.na(genome_stats$proteins)], names(cluster_protein_ratio))

#Match afdb to common proteins
dat_species_f = dat_species_f[match(int, names(dat_species_f))]

#Get n proteins per species genome
g_size = genome_stats$proteins[match(int, genome_stats$ncbi_id)]

#Calculate difference between protein n in afdb and n in genome
x = unlist(lapply(dat_species_f, function(x) nrow(x)))
y = abs(g_size-x)

#Filter afdb to just species with 'close' to complete proteomes
n = 5000
species = names(y[y<n])
g_size_f = g_size[y<n]
n_clusters = unlist(lapply(dat_species_f[match(names(y[y<n]), names(dat_species_f))], 
                           function(x) length(unique(x$` cluster_ID`))))
cluster_protein_ratio = unlist(lapply(dat_species_f[match(names(y[y<n]), names(dat_species_f))], 
                                      function(x) length(unique(x$` cluster_ID`))/nrow(x)))

max = 30000
n_clusters = n_clusters[g_size_f<max]
cluster_protein_ratio = cluster_protein_ratio[g_size_f<max]
g_size_f = g_size_f[g_size_f<max]

#Plot relationship
par(mfrow = c(1,2))
plot(g_size_f, 
     n_clusters)

plot(g_size_f, 
     cluster_protein_ratio)

#Compare by taxa
s_taxonomy = afdb_taxonomy[match(species, afdb_taxonomy$ncbi_id),]

#Plot
u = unique(s_taxonomy$superkingdom)
cols = arcadia.pal(n = 3, 'Accent')

par(mfrow = c(1,2))
plot(g_size_f, 
     n_clusters,
     type = 'n',
     ylab = 'n clusters',
     xlab = 'n proteins in genome',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 1:length(u)){
  points(g_size_f[s_taxonomy$superkingdom == u[i]],
         n_clusters[s_taxonomy$superkingdom == u[i]],
         pch = 20,
         col = cols[i])
  abline(lm(n_clusters[s_taxonomy$superkingdom == u[i]]~
            g_size_f[s_taxonomy$superkingdom == u[i]]),
         col = cols[i],
         lwd = 1.5)
}

plot(g_size_f, 
     cluster_protein_ratio,
     type = 'n',
     ylab = 'n clusters/n proteins',
     xlab = 'n proteins in genome',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 1:length(u)){
  points(g_size_f[s_taxonomy$superkingdom == u[i]],
         cluster_protein_ratio[s_taxonomy$superkingdom == u[i]],
         pch = 20,
         col = cols[i])
  abline(lm(cluster_protein_ratio[s_taxonomy$superkingdom == u[i]]~
              g_size_f[s_taxonomy$superkingdom == u[i]]),
         col = cols[i])
}

##################################################################################################
#####Prediction quality as a function of distance from solved structure, protein length, etc.#####
##################################################################################################
#Calculate distance of all afdb species from pdb species
tree_afdb_species = intersect(tree_taxonomy$ncbi_id, afdb_taxonomy$ncbi_id)

#Filter tree
x = tree_taxonomy$species[tree_taxonomy$ncbi_id%in%tree_afdb_species]
x = x[x%in%tree$tip.label]
tree_afdb = keep.tip(tree, x)

#Calculate cophenetic distances
cophen = cophenetic(tree_afdb)

#Filter columns to just pdb species
cophen = cophen[,colnames(cophen)%in%pdb_taxonomy$species]

#Split on species
cophen = split(cophen, rownames(cophen))

#Get minimum value
cophen = lapply(cophen, function(x) min(x))

#Change names
names(cophen) = afdb_taxonomy$ncbi_id[match(names(cophen), afdb_taxonomy$species)]

#Split afdb_data into clusters
clusters = split(dat, dat$` cluster_ID`)

#Get phylo distances from pdb organisms associated with each cluster
cluster_pdb_distances = list()
random = clusters[sample(1:length(clusters), 100000)]
pb <- txtProgressBar(min = 1,
                     max = length(random),
                     style = 3,
                     width = 100,
                     char = ".")

for(i in 1:length(random)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Get distances
  cluster_pdb_distances[[names(random)[i]]] = unlist(cophen[names(cophen)%in%random[[i]]$taxonomy_ID])
  
}

#Compute mean distances
mean_dists = unlist(lapply(cluster_pdb_distances, function(x) mean(x)))
min_dists = unlist(lapply(cluster_pdb_distances, function(x) min(x)))
max_dists = unlist(lapply(cluster_pdb_distances, function(x) max(x)))

#Compare plddts
cor(min_dists, cluster_stats$avgPlddt[match(names(min_dists), cluster_stats$cluster_ID)])
cor(mean_dists, cluster_stats$avgPlddt[match(names(mean_dists), cluster_stats$cluster_ID)])
cor(max_dists, cluster_stats$avgPlddt[match(names(max_dists), cluster_stats$cluster_ID)])

#Get superkingdoms for LCAs
lcas = cluster_stats$LCAtaxID[match(names(mean_dists), cluster_stats$cluster_ID)]
kings = tree_taxonomy$superkingdom[match(lcas, tree_taxonomy$ncbi_id)]
kings[is.na(kings)] = "cellular_organisms"

#Compare by lca
par(mfrow = c(1, 4))
for(i in 1:length(unique(kings))){
  m = unlist(lapply(cluster_pdb_distances[kings == unique(kings)[i]], function(x) max(x)))
  print(cor(m, cluster_stats$avgPlddt[match(names(m), cluster_stats$cluster_ID)]))
  plot(m, cluster_stats$avgPlddt[match(names(m), cluster_stats$cluster_ID)])
}

#Compare
plddt = cluster_stats$avgPlddt[match(names(mean_dists), cluster_stats$cluster_ID)]
length = cluster_stats$avgLength[match(names(mean_dists), cluster_stats$cluster_ID)]
lca = cluster_stats$LCAtaxID[match(names(mean_dists), cluster_stats$cluster_ID)]
n_members = cluster_stats[,8][match(names(mean_dists), cluster_stats$cluster_ID)]

mod = lm(plddt~mean_dists+length+lca+n_members)
mod = lm(plddt~mean_dists+min_dists+max_dists+length+lca+n_members)





#Get unique phyla
phyla = na.omit(unique(tree_taxonomy$phylum))

#Loop through and calculate PD
pds = list()
for(i in 1:length(phyla)){
  
  print(paste(i, 'out of', length(phyla)))
  #MRCA
  tips = tree_taxonomy[grep(phyla[i], tree_taxonomy$phylum),]
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
    a_tips = afdb_taxonomy[match(tips$family, afdb_taxonomy$family),]
    
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
plot(unlist(lapply(pds, function(x) x$phylogenetic_distance)))
plot(unlist(lapply(pds, function(x) x$proportion_taxa)))

plot(unlist(lapply(pds, function(x) x$phylogenetic_distance)),
     unlist(lapply(pds, function(x) x$proportion_taxa)))

#Get unique taxa
#phyla = na.omit(unique(afdb_taxonomy$phylum))
phyla_dat = split(dat, afdb_taxonomy$phylum[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id)])
phyla = na.omit(names(phyla_dat))

#Get phylum tree
tips = tree_taxonomy
tips = tips[tips$species%in%tree$tip.label,]
tips = tips[!duplicated(tips$phylum),]
tips = tips[!is.na(tips$phylum),]
tips = tips[tips$phylum%in%names(pds),]
phyla_tree = keep.tip(tree, tips$species)
phyla_tree$tip.label = tips$phylum[match(phyla_tree$tip.label, tips$species)]

#Contmaps
par(mfrow = c(1,2))
phytools::contMap(phyla_tree, unlist(lapply(pds, function(x) x$proportion_taxa)))
phytools::contMap(phyla_tree, unlist(lapply(pds, function(x) x$weighted_PD)))

#Loop through and calculate plddt for each
plddts = list()
for(i in 1:length(phyla)){
  
  print(paste(i, 'out of', length(phyla)))
  
  #Get species in taxa
  #x = afdb_taxonomy[grep(phyla[i], afdb_taxonomy$class),]
  
  #Get all clusters associated with taxa
  clusters = phyla_dat[grep(phyla[i], names(phyla_dat))]
  clusters = clusters[[1]]$` cluster_ID`

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

###################################
#####Cluster diversity by taxa#####
###################################
##Calculate per-species cluster diversity (n clusters/n proteins)
#Split cluster data by species
s = split(dat, dat$taxonomy_ID)

#Filter?
s = s[unlist(lapply(s, function(x) nrow(x)))>=100]

#Calculate cluster diversity for each
cluster_diversity = unlist(lapply(s, function(x) length(unique(x[,2]))/nrow(x)))

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

#Permutation test
n = 1000
permutation_results = list()
f = dat[dat$taxonomy_ID%in%tree_taxonomy$ncbi_id,]

pb <- txtProgressBar(min = 1,
                     max = length(s),
                     style = 3,
                     width = 100,
                     char = ".")

for(i in 1:length(s)){
  
  #Update progress bar
  setTxtProgressBar(pb, i)
  
  obs = length(unique(s[[i]]$` cluster_ID`))
  
  perms = c()
  for(j in 1:n){
    
    perms = c(perms, length(unique(dat$` cluster_ID`[sample(1:nrow(dat), nrow(s[[i]]))])))
  }
  
  p = sum(obs<perms)/n
  effect_size = obs-mean(perms)
  effect_size_norm = effect_size/sd(perms)
  permutation_results[[names(s)[i]]] = list(p_value = p,
                                            effect_size = effect_size,
                                            effect_size_norm = effect_size_norm)
  
}

#Save
#saveRDS(permutation_results, '02_output/taxa_cluster_diversity_permutation.RDS')

#Compare results
hist(unlist(lapply(permutation_results, function(x) x$effect_size_norm)))
plot(sort(unlist(lapply(permutation_results, function(x) x$p_value))))

#Compare by taxon
x = na.omit(unique(afdb_taxonomy$phylum))
effect_sizes = unlist(lapply(permutation_results, function(x) x$effect_size_norm))
effects = list()
for(i in 1:length(x)){
  
  #Get species in taxon
  species = unique(afdb_taxonomy$ncbi_id[afdb_taxonomy$phylum == x[i]])
  
  #Get associated cluster diversities
  effects[[x[i]]] = effect_sizes[names(effect_sizes)%in%species]
  
}

#Remove empty elements
effects = effects[unlist(lapply(effects, function(x) length(x)))>5]

#Compare to taxon size
e = unlist(lapply(effects, function(x) median(x)))
v = unlist(lapply(effects, function(x) sd(x)/mean(x)))

plot(e, 
     log(unlist(lapply(effects, function(x) length(x)))))

#Contmap
tips = tree_taxonomy
tips = tips[tips$species%in%tree$tip.label,]
tips = tips[tips$phylum%in%names(effects),]
tips = tips[!duplicated(tips$phylum),]

phyla_tree = keep.tip(tree, tips$species)
phyla_tree$tip.label = tips$phylum[match(phyla_tree$tip.label, tips$species)]

tmp = effects[match(phyla_tree$tip.label, names(effects))]

#Calculate contmap
e = unlist(lapply(tmp, function(x) median(x)))
e[e<0] = e[e<0]/abs(min(e))
e[e>0] = e[e>0]/max(e)

contmap = phytools::contMap(phyla_tree, 
                            e,
                            plot = FALSE)
contmap$cols[1:length(contmap$cols)] = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'RdYlBu'))(length(contmap$cols)))
plot(contmap)

z = names(tmp$Haptophyta)
lapply(z, function(x) sum(dat$taxonomy_ID%in%x))

#####################################
#####Protein accumulation curves#####
#####################################
#Get members of taxonomic category
p = split(afdb_taxonomy$ncbi_id, afdb_taxonomy$phylum)

#Filter p
p = p[unlist(lapply(p, function(x) length(x)))>=5]
p = p[order(unlist(lapply(p, function(x) length(x))))]

#Calculate for all
n_perm = 100
props = seq(0.1, 1, 0.05)
all_res = list()

#Progress bar
pb <- txtProgressBar(min = 1,
                     max = length(p),
                     style = 3,
                     width = 100,
                     char = ".")

for(h in 1:length(p)){
  
  #Update progress bar
  setTxtProgressBar(pb, h)
  
  x = p[[h]]
  toTest = unique(c(round(length(x)*props)))
  toTest = toTest[!toTest == 0]
  
  d = dat[dat$taxonomy_ID%in%x,]
  
  #Empty list to save results
  res = list()
  
  #Loop through and calculate
  for(i in 1:length(toTest)){
    
    #Create empty vector to save results
    tmp = c()
    
    #Loop through and calculate
    for(j in 1:n_perm){
      y = sample(x, toTest[i], replace = FALSE)
      tmp = c(tmp, length(unique(d$` cluster_ID`[d$taxonomy_ID%in%y])))
    }
    
    res[[as.character(i)]] = tmp
  }
  
  all_res[[names(p)[h]]] = res
}

#Save
saveRDS(all_res, '02_output/per_phylum_permutation_diversity_for_forecasting.RDS')

#Interpolate means for plotting
m = lapply(all_res, function(y) unlist(lapply(y, function(x) mean(x))))
m = lapply(m, function(x) approx(x, n = 19)$y)

#Plot distributions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(m))))

plot(m[[1]]/max(m[[1]]),
     type = 'l',
     ylim = c(0, 1),
     xlab = '% of data',
     col = cols[1],
     lwd = 1.5,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(m)){
  lines(m[[i]]/max(m[[i]]),
        col = cols[i],
        lwd = 2)
}

#Calculate polynomial fits
mod_fits = list()
for(i in 1:length(m)){
  x = m[[i]]
  y = 1:length(m[[i]])
  mod = lm(scale(x) ~ poly(y, 2))
  mod_fits[[names(m)[i]]] = predict(mod)
}

#Plot distributions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(m))))

plot(seq(0.1, 1, 0.05),
     mod_fits[[1]]/max(mod_fits[[1]]),
     type = 'l',
     ylim = c(-1.5, 1.5),
     xlab = '% of data',
     col = cols[1],
     lwd = 1,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(mod_fits)){
  lines(seq(0.1, 1, 0.05),
        mod_fits[[i]]/max(mod_fits[[i]]),
        col = cols[i],
        lwd = 1)
}

#Calculate minimum slope
min_slope = lapply(mod_fits, function(x){
  tmp = x/max(x)
  y = 1:length(x)
  
  slope = c()
  for(i in 1:(length(tmp)-2)){
    
    x = lm(tmp[i:(i+2)]~y[i:(i+2)])
    slope = c(slope, coef(x)[2])
    
  }
  return(slope)
})

#Compare min slope to taxa size
plot(log(unlist(lapply(p, function(x) length(x)))),
     unlist(lapply(min_slope, function(x) x[length(x)])),
     xlab = 'log (Species n)',
     ylab = 'Minimum slope',
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 20,
     cex = 1.5,
     col = 'gray60')

#Forecast with polynomial regression
mod_fits = list()
for(i in 1:length(m)){
  x = m[[i]]
  y = 1:length(m[[i]])
  mod = lm(scale(x) ~ poly(y, 2))
  
  newData = data.frame(y = 1:(length(tmp)+1000))
  pred = predict(mod, newData)
  
  if(!which.max(pred) == length(pred)){
    max_val = max(pred)
    pred[which.max(pred):length(pred)][pred[which.max(pred):length(pred)]<max_val] = max_val
  }
  
  mod_fits[[names(m)[i]]] = pred
}

#Plot predictions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(m))))
plot(mod_fits[[1]]/max(mod_fits[[1]]),
     type = 'l',
     xlim = c(0, 100),
     ylim = c(-1.5, 1),
     xlab = '% of data',
     col = cols[1],
     lwd = 1,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n')
axis(1, seq(0, 100, 20), seq(0, 100, 20)*5, cex.axis = 1.5)
for(i in 2:length(mod_fits)){
  lines(mod_fits[[i]]/max(mod_fits[[i]]),
        col = cols[i],
        lwd = 1)
}
abline(v = 20, lwd = 1.5, lty = 'dashed')

#Compare saturated values
max_val = unlist(lapply(mod_fits, function(x) which.max(x)))
names(max_val) = unlist(lapply(strsplit(names(max_val), '\\.'), function(x) x[1]))
plot(log(unlist(lapply(p, function(x) length(x)))),
     log(max_val))

#Predict species n per taxa that would be needed to saturate space
species_n_saturate = c()
for(i in 1:length(p)){
  x = ceiling(length(p[[i]])*0.05)
  species_n_saturate = c(species_n_saturate, max_val[i]*x)
}

#####################################################
#####Phylogenetic distance as a function of taxa#####
#####################################################
#Downsample tree to families
new_tree = filter_tree(tree,
                       resolution = 'family',
                       verbose = TRUE,
                       change_tree_labels = TRUE, 
                       return_taxonomy = TRUE)

#Get phyla associated with each family for afdb
u <- lapply(taxonomy, function(x) cbind(x[grep(paste("^", 'family', "$", sep = ""), names(x))],
                                        x[grep(paste("^", 'phylum', "$", sep = ""), names(x))]))
u = u[unlist(lapply(u, function(x) ncol(x)))==2]
u = do.call(rbind, u)
u = u[!duplicated(u[,1]),]

#Get phyla associated with each family for timetree
v <- lapply(new_tree$taxonomy, function(x) cbind(x[grep(paste("^", 'family', "$", sep = ""), names(x))],
                                                 x[grep(paste("^", 'phylum', "$", sep = ""), names(x))]))
v = v[unlist(lapply(v, function(x) ncol(x)))==2]
v = do.call(rbind, v)
#v = v[v[,1]%in%new_tree$tree$tip.label,]
v = v[!duplicated(v[,1]),]

#Calculate phylogenetic distance for each phylum
p = intersect(unique(u[,2]), unique(v[,2]))
dists = list()

for(i in 1:length(p)){
  
  #afdb families
  x = u[u[,2]%in%p[i],]
  
  #timetree families
  y = v[v[,2]%in%p[i],]
  
  #Filter tree to just phylum
  tip = new_tree$tree$tip.label[!new_tree$tree$tip.label %in% unique(y[,1])]
                                
  tree_f = drop.tip(new_tree$tree, 
                    tip)
  tree_f = root(tree_f, 1, resolve.root = TRUE)

  #Get binary community vector
  y = y[match(tree_f$tip.label, y[,1]),]
  z = intersect(x[,1], y[,1])
  comm = y[,1]%in%z
  comm[comm == TRUE] = 1
  comm[comm == FALSE] = 0
  comm = as.matrix(t(comm))
  colnames(comm) = tree_f$tip.label
  rownames(comm) = 1:nrow(comm)
  
  #Get phylogenetic distance
  tree_spread(tree_f, 
              comm)
  
  # Calculate cophenic distance
  cophen <- cophenetic.phylo(tree_f)
  
  # Calculate phylogenetic alpha diversity
  PD <- pd(
    samp = comm,
    tree = tree_f,
    include.root = TRUE
  )
  
  # Calculate standard effect sizes
  ses <- ses.mpd(comm,
                 cophen,
                 null.model = "taxa.labels",
                 abundance.weighted = FALSE,
                 runs = 100
  )
  
  
}


########################################################
#####Loading and cleaning catalogue of life species#####
########################################################
#Load catalogue of life species
col_species = read.delim('~/Desktop/protein_space_size_estimation/00_data/col-texttree-dataset-9923.txt', quote = "")

#Convert to dataframe
col_species = data.frame(name = unlist(lapply(strsplit(col_species[,1], '\\['), function(x) x[1])),
                         taxonomic_id = gsub('\\]', '', unlist(lapply(strsplit(col_species[,1], '\\['), function(x) x[2]))))

#ID which entries are phyla
x = which(col_species[,2] == 'phylum')
x = c(x, nrow(col_species))

#Loop through and extract entries following each phylum (will be constituent members)
phyla_totals = list()
for(i in 1:(length(x)-1)){
  phyla_totals[[col_species[x[i],1]]] = col_species[x[i]:(x[i+1]-1),]
}

#Select just species
phyla_totals = lapply(phyla_totals, function(x) x[grep('species', x[,2]),])

#Trim whitespace from names
names(phyla_totals) = lapply(names(phyla_totals), function(x) trimws(x))

#Select just first element of names
for(i in 1:length(phyla_totals)){
  names(phyla_totals)[i]= strsplit(names(phyla_totals)[i], ' ')[[1]][1]
}

#Collapse redundant entries (denoted by *)
x = grep('\\*', names(phyla_totals))
for(i in 1:length(x)){
  if(nrow(phyla_totals[[(x[i]-1)]])>0){
    phyla_totals[[(x[i]-1)]] = rbind(phyla_totals[[(x[i]-1)]], phyla_totals[[x[i]]])
  }else{
    phyla_totals[[(x[i]-1)]] = phyla_totals[[x[i]]]
  }
}

#Remove redundant entries
phyla_totals = phyla_totals[-grep('\\*', names(phyla_totals))]

##################################
#####Distribution across taxa#####
##################################
#Get member of taxonomic category
taxa = unlist(lapply(taxonomy, function(x) x[grep('^phylum', names(x))]))
n = unlist(lapply(strsplit(names(taxa), '\\.'), function(x) x[1]))
p = split(names(taxonomy)[match(n, names(taxonomy))], taxa)

#Compare species number
species_n = lapply(p, function(x) length(x))
barplot(unlist(species_n), 
        las = 2)

#Calculate summary stats (protein #, cluster #, etc.)
summary_stats = list()

pb <- txtProgressBar(min = 1,
                     max = length(p),
                     style = 3,
                     width = 100,
                     char = ".")

for(i in 1:length(p)){
  
  setTxtProgressBar(pb, i)
  
  d = dat[dat$taxonomy_ID%in%p[[i]],]
  n_proteins = nrow(d)
  n_clusters = length(unique(d$` cluster_ID`))
  n_species = length(unique(d$taxonomy_ID))
  l = list(n_proteins = n_proteins,
           n_clusters = n_clusters,
           n_species = n_species)
  summary_stats[[names(p)[i]]] = l
}

#Compare
np = unlist(lapply(summary_stats, function(x) x$n_proteins))
nc = unlist(lapply(summary_stats, function(x) x$n_clusters))
ns = unlist(lapply(summary_stats, function(x) x$n_species))
plot(np, nc)
plot(np/nc, ns)

###############################################
#####Test species effect on space coverage#####
###############################################
#Function to calculate cluster saturation
cluster_saturation = function(input_list,
                              sampling_interval = 500){
  
  #Distribution of cluster n per species
  clusters = unlist(lapply(input_list, function(x) length(unique(x$` cluster_ID`))))
  
  #Order data on n clusters per taxa
  input_list = input_list[order(clusters)]
  
  #Create vector to save results
  clust_cdf = c()
  
  #Set up progress bar
  pb <- txtProgressBar(min = 1,
                       max = length(input_list),
                       style = 3,
                       width = 50,
                       char = ".")
  
  #Calculate for each
  for(i in c(1, seq(sampling_interval, length(input_list), sampling_interval))){
    setTxtProgressBar(pb, i)
    clust_cdf = c(clust_cdf, length(unique(unlist(lapply(input_list[1:i], function(x) x$` cluster_ID`)))))
  }
  
  #Add names
  names(clust_df) = c(1, seq(sampling_interval, length(input_list), sampling_interval))
  
  #Return
  return(clust_cdf)
}


dat$genus = unlist(n)

n = lapply(taxonomy, 
           function(x) x[grep('genus', names(x))])[match(dat$taxonomy_ID, names(taxonomy))]

p = split(as.data.frame(dat), 
          n)

tmp = cluster_saturation(p)

#Distribution of cluster n per species
p = split(as.data.frame(dat), dat$taxonomy_ID)
species_clusters = unlist(lapply(p, function(x) length(unique(x$` cluster_ID`))))

#Cumulative new clusters
p = p[order(species_clusters)]

clust_cdf = c()
pb <- txtProgressBar(min = 1,
                     max = length(p),
                     style = 3,
                     width = 50,
                     char = ".")

for(i in c(1, seq(500, length(p), 500))){
  setTxtProgressBar(pb, i)
  clust_cdf = c(clust_cdf, length(unique(unlist(lapply(p[1:i], function(x) x$` cluster_ID`)))))
}

#Plot
plot(c(1, seq(500, length(p), 500)),
     clust_cdf,
     xlab = 'n species',
     ylab = 'n clusters',
     cex.axis = 1.5,
     cex.lab = 1.5,
     type = 'l',
     lwd = 1.5,
     col = arcadia.pal(3, name = 'Accent')[1])

#Predict n clusters from more species
y = c(1, seq(500, length(p), 500))
x = clust_cdf
mod = lm(x ~ poly(y, 4))

newdata <- data.frame('y'= c(1, seq(500, length(p), 500), seq(100000, 6000000, 100000)))

plot(c(1, seq(500, length(p), 500), seq(100000, 6000000, 100000)),
     predict(mod, newdata),
     xlab = 'n species',
     ylab = 'n clusters',
     cex.axis = 1.5,
     cex.lab = 1.5,
     type = 'l',
     lwd = 1.5,
     col = arcadia.pal(3, name = 'Accent')[1])
lines(c(1, seq(500, length(p), 500)),
      clust_cdf)

#Cumulative distribution function
cdf = ecdf(dat$` cluster_ID`)

#By species
#toTest = c(round(length(taxa)*seq(0.1, 0.9, 0.1)))
toTest = c(1, seq(5, 1000, 5))

#Loop through and calculate permuted completeness of space per proportion
all_res = list()
n_perms = 100

for(i in 1:length(toTest)){
  print(toTest[i])
  tmp = c()
  
  #Progress bar
  pb <- txtProgressBar(min = 1,
                       max = n_perms,
                       style = 3,
                       width = 50,
                       char = ".")
  
  for(j in 1:n_perms){
    setTxtProgressBar(pb, j)
    
    y = taxa[sample(1:length(taxa), toTest[i], replace = FALSE)]
    
    x = dat$` cluster_ID`[dat$taxonomy_ID%in%y]
    x = length(unique(x))
    tmp = c(tmp, x)
  }
  
  all_res[[as.character(toTest[i])]] = tmp
}

#########################################################################
#####Test coverage of space as a function of data % for all proteins#####
#########################################################################
#Number of proteins to test
toTest = c(round(nrow(dat)*seq(0.1, 0.9, 0.1)))

#Loop through and calculate permuted completeness of space per proportion
all_res = list()
d = length(unique(dat$` cluster_ID`))

for(i in 1:length(toTest)){
  print(seq(0.1, 0.9, 0.1)[i])
  tmp = c()
  
  #Progress bar
  pb <- txtProgressBar(min = 1,
                       max = 10,
                       style = 3,
                       width = 50,
                       char = ".")
  
  for(j in 1:10){
    setTxtProgressBar(pb, j)
    x = dat$` cluster_ID`[sample(1:nrow(dat), toTest[i], replace = FALSE)]
    x = length(unique(x))
    tmp = c(tmp, x)
  }
  
  all_res[[as.character(seq(0.1, 0.9, 0.1)[i])]] = tmp
}

#Plot coverage function
plot(seq(0.1, 0.9, 0.1),
     unlist(lapply(all_res, function(x) mean(x))),
     #ylim = c(0, 1),
     pch = 20,
     col = 'black',
     type = 'b',
     cex = 2,
     xlab = '% data',
     ylab = '#of clusters',
     cex.axis = 1.5,
     cex.lab = 1.5)

#Slope
s = coef(lm(unlist(lapply(res, function(x) mean(x)))~seq(0.1, 0.9, 0.1)))[2]

######################################################
#####Test coverage of space as a function of taxa#####
######################################################
#Get member of taxonomic category
#taxa = unlist(lapply(taxonomy, function(x) x[grep('^phylum', names(x))]))
#n = unlist(lapply(strsplit(names(taxa), '\\.'), function(x) x[1]))
#p = split(names(taxonomy)[match(n, names(taxonomy))], taxa)
p = split(afdb_taxonomy$ncbi_id, afdb_taxonomy$phylum)
  
#Filter p
p = p[unlist(lapply(p, function(x) length(x)))>=5]
p = p[order(unlist(lapply(p, function(x) length(x))))]

#Split dat
dat_phyla = split(dat, afdb_taxonomy$phylum[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id)])

#Calculate for all
n_perm = 100
props = seq(0.1, 1, 0.05)
all_res = list()

#Progress bar
pb <- txtProgressBar(min = 1,
                     max = length(p),
                     style = 3,
                     width = 100,
                     char = ".")

for(h in 1:length(p)){
  
  #Update progress bar
  setTxtProgressBar(pb, h)
  
  x = p[[h]]
  toTest = unique(c(round(length(x)*props)))
  toTest = toTest[!toTest == 0]

  d = dat_phyla[[grep(names(p)[h], names(dat_phyla))]]
  d = split(d, d$taxonomy_ID)

  #Empty list to save results
  res = list()
  
  #Loop through and calculate
  for(i in 1:length(toTest)){
    
    #print(i)
    #Create empty vector to save results
    tmp = c()
    
    #Loop through and calculate
    for(j in 1:n_perm){
      y = sample(names(d), toTest[i], replace = FALSE)
      tmp = c(tmp, length(unique(unlist(lapply(d[y], function(x) x$` cluster_ID`)))))
      #tmp = c(tmp, length(unique(d$` cluster_ID`[d$taxonomy_ID%in%y])))
    }
    
    res[[as.character(i)]] = tmp
  }
  
  all_res[[names(p)[h]]] = res
}

#SAVE
#saveRDS(all_res, '02_output/per_phyla_cluster_n_permutations_for_forecasting.RDS')

#Interpolate means for plotting
m = lapply(all_res, function(y) unlist(lapply(y, function(x) mean(x))))
m = lapply(m, function(x) approx(x, n = 19)$y)

#Plot distributions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(m))))

plot(m[[1]]/max(m[[1]]),
     type = 'l',
     ylim = c(0, 1),
     xlab = '% of data',
     col = cols[1],
     lwd = 1.5,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(m)){
  lines(m[[i]]/max(m[[i]]),
        col = cols[i],
        lwd = 2)
}

#Calculate polynomial fits
mod_fits = list()
for(i in 1:length(m)){
  x = m[[i]]
  y = 1:length(m[[i]])
  mod = lm(scale(x) ~ poly(y, 2))
  mod_fits[[names(m)[i]]] = predict(mod)
}

#Plot distributions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))(length(m))))

plot(seq(0.1, 1, 0.05),
     mod_fits[[1]]/max(mod_fits[[1]]),
     type = 'l',
     ylim = c(-1.5, 1.5),
     xlab = '% of data',
     col = cols[1],
     lwd = 1,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5)
for(i in 2:length(mod_fits)){
  lines(seq(0.1, 1, 0.05),
        mod_fits[[i]]/max(mod_fits[[i]]),
        col = cols[i],
        lwd = 1)
}

#Calculate minimum slope
min_slope = lapply(mod_fits, function(x){
  tmp = x/max(x)
  y = 1:length(x)
  
  slope = c()
  for(i in 1:(length(tmp)-2)){
    
    x = lm(tmp[i:(i+2)]~y[i:(i+2)])
    slope = c(slope, coef(x)[2])
    
  }
  return(slope)
})

#Compare min slope to taxa size
par(mfrow = c(1,2))
plot(log(unlist(lapply(p, function(x) length(x)))),
     unlist(lapply(min_slope, function(x) x[length(x)])),
     xlab = 'log (Species n)',
     ylab = 'Minimum slope',
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 20,
     cex = 1.5,
     col = 'gray60')

plot(log(unlist(lapply(p, function(x) length(x)))),
     unlist(lapply(min_slope, function(x) x[1])),
     xlab = 'log (Species n)',
     ylab = 'Maximum slope',
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 20,
     cex = 1.5,
     col = 'gray60')

#Forecast with polynomial regression
mod_fits = list()
slopes = c()
for(i in 1:length(m)){
  x = m[[i]]
  y = 1:length(m[[i]])
  mod = lm(scale(x) ~ poly(y, 2))
  slopes = c(slopes, coef(mod)[1])
  
  newData = data.frame(y = 1:(length(m[[i]])+1000))
  pred = predict(mod, newData)
  
  if(!which.max(pred) == length(pred)){
    max_val = max(pred)
    pred[which.max(pred):length(pred)][pred[which.max(pred):length(pred)]<max_val] = max_val
  }

  mod_fits[[names(m)[i]]] = pred
}
names(slopes) = names(mod_fits)

#Plot predictions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(m))))
plot(mod_fits[[1]]/max(mod_fits[[1]]),
     type = 'l',
     xlim = c(0, 100),
     ylim = c(-1.5, 1),
     xlab = '% of data',
     col = cols[1],
     lwd = 1,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n')
axis(1, seq(0, 100, 20), seq(0, 100, 20)*5, cex.axis = 1.5)
for(i in 2:length(mod_fits)){
  lines(mod_fits[[i]]/max(mod_fits[[i]]),
        col = cols[i],
        lwd = 1)
}
abline(v = 20, lwd = 1.5, lty = 'dashed')

#Compare saturated values
max_val = unlist(lapply(mod_fits, function(x) which.max(x)))
names(max_val) = unlist(lapply(strsplit(names(max_val), '\\.'), function(x) x[1]))
plot(log(unlist(lapply(p, function(x) length(x)))),
     log(max_val))

#Predict species n per taxa that would be needed to saturate space
species_n_saturate = c()
for(i in 1:length(p)){
  x = ceiling(length(p[[i]])*0.05)
  species_n_saturate = c(species_n_saturate, max_val[i]*x)
}

#Calculate total n of species in each phyla
phyla_totals = lapply(names(mod_fits), function(x) length(getDescendants(getId(x))))
names(phyla_totals) = names(mod_fits)
phyla_totals = phyla_totals[match(names(m), names(phyla_totals))]

#Forecast using species n
mod_fits = list()
slopes = c()
for(i in 1:length(m)){
  x = m[[i]]
  y = 1:length(m[[i]])
  obs_n = length(unique(dat_phyla[[grep(names(m)[i], names(dat_phyla))]][,3]))
  #mod = lm(scale(x) ~ poly(y, 2))
  mod = lm(x ~ poly(y, 2))
  slopes = c(slopes, coef(mod)[1])
  
  newData = data.frame(y = 1:(length(m[[i]])+((phyla_totals[[i]]-obs_n)/(obs_n/20))))
  pred = predict(mod, newData)
  
  if(!which.max(pred) == length(pred)){
    max_val = max(pred)
    pred[which.max(pred):length(pred)][pred[which.max(pred):length(pred)]<max_val] = max_val
  }
  
  mod_fits[[names(m)[i]]] = pred
}
names(slopes) = names(mod_fits)

#Plot predictions
cols = darken_color(rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(length(m))))
mod_fits = mod_fits[order(slopes)]
plot(log(1:length(mod_fits[[1]])),
     mod_fits[[1]]/max(mod_fits[[1]]),
     #exp(mod_fits[[1]])/max(exp(mod_fits[[1]])),
     type = 'l',
     xlim = c(0, log(max(unlist(lapply(mod_fits, function(x) length(x)))))),
     #ylim = c(-1.5, 1),
     ylim = c(0, 1),
     xlab = 'log(n Species)',
     col = cols[1],
     lwd = 1,
     ylab = 'n clusters (normalized)',
     cex.axis = 1.5,
     cex.lab = 1.5)
#axis(1, seq(0, 100, 20), seq(0, 100, 20)*5, cex.axis = 1.5)
for(i in 2:length(mod_fits)){
  lines(log(1:length(mod_fits[[i]])),
        mod_fits[[i]]/max(mod_fits[[i]]),
        #exp(mod_fits[[i]])/max(exp(mod_fits[[i]])),
        col = cols[i],
        lwd = 1)
}
#abline(v = 20, lwd = 1.5, lty = 'dashed')

#Compare saturated values
max_val = unlist(lapply(mod_fits, function(x) which.max(x)))
names(max_val) = unlist(lapply(strsplit(names(max_val), '\\.'), function(x) x[1]))
plot(log(unlist(lapply(p, function(x) length(x)))),
     log(max_val))

#Predict species n per taxa that would be needed to saturate space
species_n_saturate = c()
for(i in 1:length(p)){
  x = ceiling(length(p[[i]])*0.05)
  species_n_saturate = c(species_n_saturate, max_val[i]*x)
}

#Predict size of protein space by multiplying species by average n proteins per taxa
taxa_n = unlist(lapply(dat_phyla, function(x) mean(unlist(lapply(split(x, x$taxonomy_ID), function(y) length(unique(y[,2])))))))
taxa_n = taxa_n[match(names(species_n_saturate), names(taxa_n))]
total_clusters = species_n_saturate*taxa_n
sum(total_clusters)
#[1] 65804984

#Compare total clusters by kingdom
n = afdb_taxonomy$superkingdom[match(names(total_clusters), afdb_taxonomy$phylum)]
n = split(total_clusters, n)
total_clusters_king = lapply(n, function(x) sum(x))
obs_clusters_king = lapply(split(dat, afdb_taxonomy$superkingdom[match(dat$taxonomy_ID, afdb_taxonomy$ncbi_id)]), function(x) length(unique(x$` cluster_ID`)))

#Visualize


###################################################################################################################################
#####Test coverage of space as a function of taxa (adjusting incomplete to complete proteomes for each species via estimation)#####
###################################################################################################################################
##Assess proteome completeness
#Split afdb stats on species
dat_species = split(dat, dat$taxonomy_ID)

#Calculate n proteins for each
n_proteins = lapply(dat_species, function(x) c(length(x$member_ID), length(unique(x$` cluster_ID`))))

#Intersect with genome stats
int = intersect(names(n_proteins), genome_stats$ncbi_id)
n_proteins = n_proteins[match(int, names(n_proteins))]
g = genome_stats[match(int, genome_stats$ncbi_id),]
g$n_afdb_proteins = unlist(lapply(n_proteins, function(x) x[1]))
g$n_afdb_clusters = unlist(lapply(n_proteins, function(x) x[2]))

#Filter
g = g[g$n_afdb_proteins>5,]
g = g[!is.na(g$proteins),]

#Calculate estimated n of clusters for each
g$expected_clusters = rep(NA, nrow(g))
for(i in 1:nrow(g)){
  if(g$n_afdb_proteins[i]<g$proteins[i]){
    r = g$n_afdb_clusters[i]/g$n_afdb_proteins[i]
    g$expected_clusters[i] = g$proteins[i]*r
  }else{
    g$expected_clusters[i] = g$n_afdb_clusters[i]
  }
}

#Compare genome size and clusters
plot(g$n_afdb_clusters/g$n_afdb_proteins, g$size)
plot(g$n_afdb_clusters, g$proteins)



################################################################
#####Compare phyla saturation to protein n, species n, etc.#####
################################################################
#Create metadata matrix
n = names(max_val)
meta = data.frame(max_val = max_val,
                  n_clusters = nc[match(names(max_val), names(nc))],
                  n_proteins = np[match(names(max_val), names(np))],
                  n_species = ns[match(names(max_val), names(ns))],
                  king = rep(NA, length(max_val)),
                  clade = rep(NA, length(max_val)),
                  superking = rep(NA, length(max_val)))
for(i in 1:length(max_val)){
  print(i)
  tmp = taxonomy[grep(names(max_val)[i], taxonomy)][1]
  
  if(length(tmp[[1]][grep("^kingdom$", names(tmp[[1]]))])>0){
    meta$king[i] = tmp[[1]][grep("^kingdom$", names(tmp[[1]]))]
  }else{
    meta$king[i] = tmp[[1]][grep("^superkingdom$", names(tmp[[1]]))]
  }
  meta$superking[i] = tmp[[1]][grep("^superkingdom$", names(tmp[[1]]))]
  
  if(length(tmp[[1]][grep("^clade$", names(tmp[[1]]))])>0){
    meta$clade[i] = tmp[[1]][grep("^clade$", names(tmp[[1]]))]
  }else{
    meta$clade[i] = NA
  }
}

#Compare predictors of saturation
mod = lm(max_val ~ n_clusters + n_proteins + n_species + superking, data = meta)
mod = lm(n_clusters ~ n_proteins + n_species + king, data = meta)

#Get phyla associated with super kingdoms
king = unlist(lapply(taxonomy, function(x) x[grep('^superkingdom', names(x))]))
n = unlist(lapply(strsplit(names(king), '\\.'), function(x) x[1]))
king = split(taxonomy[match(n, names(taxonomy))], king)
king = lapply(king, function(x) unique(unlist(lapply(x, function(y) y[grep('^phylum', names(y))]))))

#Compare saturation values as a function of superkingdoms
king_max = list()
for(i in 1:length(king)){
  king_max[[names(king)[i]]] = max_val[unlist(lapply(strsplit(names(max_val), '\\.'), function(x) x[1]))%in%king[[i]]]
}

vioplot::vioplot(king_max[1:3],
                 ylim = c(0, 500))

##############################
#####Compare to phylogeny#####
##############################
#Get species and unique phyla represented in data set
s = unlist(lapply(taxonomy, function(x) x[grep('species', names(x))]))
s = lapply(s, function(x){
  x = paste(strsplit(x, ' ')[[1]][1],
            strsplit(x, ' ')[[1]][2],
            sep = '_')})
s = unlist(s)

#Intersect with tree
inTree = which(s%in%tree$tip.label)

#Get phyla associated with each
phyla = taxonomy[inTree]
p = unique(unlist(lapply(phyla, function(x) x[grep('^phylum', names(x))])))
phyla = phyla[match(p, unlist(lapply(phyla, function(x) x[grep('^phylum', names(x))])))]

#Subset tree
s = unlist(lapply(phyla, function(x) x[grep('species', names(x))]))
s = lapply(s, function(x){
  x = paste(strsplit(x, ' ')[[1]][1],
            strsplit(x, ' ')[[1]][2],
            sep = '_')})
s = unlist(s)
tree_phyla = drop.tip(tree, tree$tip.label[!tree$tip.label%in%s])



