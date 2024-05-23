source('~/Documents/Research/github/2024-protein-universe/code/R/afdb-phylo-utils.R')
source('~/Documents/Research/Rfiles/arcadia_gradients.R')

##TO DO: Split representative protein metadata by phylum, calculate mean plddt, 
##test is there's a correlation between representation in training data
##and mean plddt/distance from humans or solved structures
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

#Load human clusters
human_clusters = read.delim('00_data/foldseek_clusters/3-sapId_sapGO_repId_cluFlag_LCAtaxId.tsv',
                            header = FALSE)

#Load timetree abiotic data (donwloaded from https://github.com/timetree-of-life/MEGA-TT/blob/master)
co2 = load_tt_data('00_data/timetree_data/CO2.txt')
o2 = load_tt_data('00_data/timetree_data/O2.txt')
luminosity = load_tt_data('00_data/timetree_data/luminosity.txt')
impacts = read.csv('00_data/timetree_data/earth_impacts.csv',
                   header = FALSE)

#Uniprot reference proteome list
uniprot_reference_proteomes = read.delim('00_data/proteomes_AND_proteome_type_1_2024_02_28.tsv') 

#Major evolutionary events list (from https://www.liebertpub.com/doi/10.1089/ast.2021.0119)
evo_events = c(3.5, 3.4, 2.7, 2.4, 1.86, 1, 0.8, 0.635, 0.55)
names(evo_events) = c('LUCA', 
                      'photosynthesis',
                      'cyanobacteria',
                      'oxidation_event',
                      'LECA',
                      'multicellular_red_algae',
                      'heterotrophic_euks',
                      'metazoa',
                      'noe')

#Set up colors
all_colors <- c(
  "#5088C5", "#F28360", "#F7B846", "#97CD78",
           "#7A77AB", "#f898AE", "#3B9886", "#c85152",
           "#73B5E3", "#BAB0A8", "#8A99AD", "#FFB984",
           "#C6E7F4", "#F8C5C1", "#F5E4BE", "#B5BEA4",
           "#DCBFFC", "#B6C8D4", "#DAD3C7", "#DA9085"
)

###########################################################
#####Figure 2: Distribution of species in the PDB/AFDB#####
###########################################################
#Get overall AFDB species distribution
afdb_species = table(dat$taxonomy_ID)

#Pull out ncbi IDs
afdb_names = names(afdb_species)

#Convert to numeric
afdb_species = as.numeric(afdb_species)

#Get overall PDB species distribution (need to replace commas in totals)
pdb_species = as.numeric(gsub(',', '', pdb[,2]))

#Get PDB names
pdb_names = name2taxid(pdb$V1, out_type = "summary")$id

#Generate simulated normal distributions for statistical comparison
afdb_sim = rnorm(n = length(afdb_species),
                 mean = mean(afdb_species))
pdb_sim = rnorm(n = length(pdb_species),
                mean = mean(pdb_species))

#Compare simulated and observed with ks test
ks.test(afdb_sim, 
        afdb_species)$p.value

ks.test(pdb_sim, 
        pdb_species)$p.value

#Plot together
plot(ecdf(afdb_species),
     verticals = TRUE,
     main = '',
     cex.lab = 1.5,
     cex.axis = 1.5,
     lwd = 2,
     xlab = 'n proteins',
     ylab = 'Cumulative proportion',
     col = all_colors[1],
     yaxt = 'n')
axis(2,
     at = seq(0.0, 1.0, 0.2),
     labels = signif(seq(0.0, 1.0, 0.2), 2),
     las = 2,
     cex.axis = 1.5)
lines(ecdf(pdb_species),
      verticals = TRUE,
      lwd = 2,
      col = all_colors[2],
      do.points = FALSE)
text(rep(280000, 2),
     c(0.2, 0.1),
     cex = 1.5,
     adj = 1,
     c('AFDB', 'PDB'),
     col = c(all_colors[1],
             all_colors[2]))

#AFDB: circle packing plot
png('~/Desktop/protein_universe_pub/figures/afdb_circle_packing.png', width = 2400, height = 2400)
data = data.frame(group = afdb_names,
                  value = afdb_species)
data$kingdom = afdb_taxonomy$superkingdom[match(data$group, afdb_taxonomy$ncbi_id)]
data = data[data$kingdom%in%c('Archaea', 'Bacteria', 'Eukaryota'),]
data = data[order(data$kingdom),]

packing <- circleProgressiveLayout(data$value, 
                                   sizetype='area')
data <- cbind(data, packing)
dat.gg <- circleLayoutVertices(packing, 
                               npoints=50)
dat.gg$value <- rep(data$value, each=51)
dat.gg$kingdom <- rep(data$kingdom, each=51)

cols = all_colors[5:7]
names(cols) = c('Archaea', 'Bacteria', 'Eukaryota')
cols = cols[match(data$kingdom,
                  names(cols))]

ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill = as.factor(kingdom)), linewidth = 0.2, colour = "black") +
  
  scale_fill_manual(values = cols) +

  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void()  + 
  theme(legend.position="none") + 
  coord_equal()

dev.off()

#PDB: circle packing plot
png('~/Desktop/protein_universe_pub/figures/pdb_circle_packing.png', width = 2400, height = 2400)

data = data.frame(group = pdb_names,
                  value = pdb_species)
data$kingdom = pdb_taxonomy$superkingdom[match(data$group, pdb_taxonomy$ncbi_id)]
data = data[data$kingdom%in%c('Archaea', 'Bacteria', 'Eukaryota'),]
data = data[order(data$kingdom),]

packing <- circleProgressiveLayout(data$value, 
                                   sizetype='area')
data <- cbind(data, packing)
dat.gg <- circleLayoutVertices(packing, 
                               npoints=50)
dat.gg$value <- rep(data$value, each=51)
dat.gg$kingdom <- rep(data$kingdom, each=51)

cols = all_colors[5:7]
names(cols) = c('Archaea', 'Bacteria', 'Eukaryota')
cols = cols[match(data$kingdom,
                  names(cols))]

ggplot() + 
  
  # Make the bubbles
  geom_polygon(data = dat.gg, aes(x, y, group = id, fill = as.factor(kingdom)), linewidth = 0.2, colour = "black") +
  
  scale_fill_manual(values = cols) +
  
  scale_size_continuous(range = c(1,4)) +
  
  # General theme:
  theme_void()  + 
  theme(legend.position="none") + 
  coord_equal()

dev.off()


#PDB pie chart
par(mfrow = c(1,2))
data = data.frame(group = pdb_names,
                  value = pdb_species)
data$kingdom = pdb_taxonomy$superkingdom[match(data$group, pdb_taxonomy$ncbi_id)]
data = data[data$kingdom%in%c('Archaea', 'Bacteria', 'Eukaryota'),]
data = lapply(split(data, data$kingdom), function(x) sum(x$value))
pie(unlist(data), 
    names(data),
    col = all_colors[5:7], 
    border = FALSE)

#AFDB pie chart
data = data.frame(group = afdb_names,
                  value = afdb_species)
data$kingdom = afdb_taxonomy$superkingdom[match(data$group, afdb_taxonomy$ncbi_id)]
data = data[data$kingdom%in%c('Archaea', 'Bacteria', 'Eukaryota'),]
data = lapply(split(data, data$kingdom), function(x) sum(x$value))
pie(unlist(data), 
    names(data),
    col = all_colors[5:7], 
    border = FALSE)

##########################################################
#####Figure 3: Compare taxonomic completeness in AFDB#####
##########################################################
#Calculate completeness via phylogenetic distance with 'clade_PD'
pds_fam = clade_PD(tree_taxonomy,
                   afdb_taxonomy,
                   tree,
                   taxa1 = 'phylum',
                   taxa2 = 'family',
                   verbose = TRUE,
                   show.tip.label = FALSE)

#Extract weighted phylogenetic distances
weighted_pds = unlist(lapply(pds_fam$phylogenetic_distance, 
                             function(x) x$weighted_PD))

#Split on kingdom
weighted_pds = split(weighted_pds, 
                     afdb_taxonomy$superkingdom[match(names(weighted_pds), 
                                                      afdb_taxonomy$phylum)])

#Test significance
kruskal.test(weighted_pds)
dunn.test::dunn.test(weighted_pds)

#Plot contmap and taxonomic completeness together
par(mfrow = c(1,3),
    mar = c(8.1, 4.1, 6.1, 2.1))

#Plot distribution of taxonomic completeness
plot(sort(unlist(weighted_pds)),
     cex.main = 1.5, 
     cex.axis = 1.5,
     cex.lab = 1.5,
     cex = 1.5,
     xaxt = 'n',
     ylab = 'Taxonomic completeness',
     pch = 20,
     col = 'gray50',
     xlab = 'Phylum')

#Plot contmap
plot(pds_fam$contmap, 
     ftype = 'off',
     mar=c(1.1,4.1,4.1,2.1))

#Plot
par(mar = c(8.1, 4.1, 6.1, 2.1))
vioplot::vioplot(weighted_pds,
                 col = all_colors[5:7],
                 side = "right",
                 ylab = '',
                 xlab = '', 
                 yaxt = 'n',
                 font.main = 1,
                 cex.main = 1.5,
                 las = 2, 
                 ylim = c(0, 1),
                 cex.axis = 1.5,
                 cex.lab = 1.5)
axis(2,
     seq(0, 1, 0.2),
     seq(0, 1, 0.2),
     cex.axis = 1.5,
     las = 2)
title(ylab = 'Taxonomic completeness',
      cex.lab = 1.5)

stripchart(weighted_pds,
           col = all_colors[5:7],
           at = seq(0.8, (length(weighted_pds)-1)+0.8, 1), 
           jitter = 0.1,
           method = "jitter", 
           vertical = TRUE, 
           cex = 1,
           pch = 20, 
           add = TRUE)

#Predictors of PD
w_pds = unlist(lapply(pds_fam$phylogenetic_distance, 
                      function(x) x$weighted_PD))

total_pds = unlist(lapply(pds_fam$phylogenetic_distance, 
                          function(x) x$total_phylogenetic_distance))

w_pds = w_pds[match(names(total_pds),
                    names(w_pds))]

n_fams = unlist(lapply(names(total_pds), function(x) 
  length(na.omit(unique(tree_taxonomy$family[tree_taxonomy$phylum == x])))))

king = tree_taxonomy$superkingdom[match(names(w_pds),
                                        tree_taxonomy$phylum)]

mod = lm(w_pds ~ n_fams + total_pds)
summary(mod)

##########################################################################
#####Figure 4: Phylogenetic distribution of "representative proteins"#####
##########################################################################
#Get taxa for representative proteins
representatives = unique(dat$taxonomy_ID[match(cluster_stats$cluster_ID, dat$member_ID)])

#Get taxonomic IDs for each
representatives_taxonomy = afdb_taxonomy[match(representatives,
                                               afdb_taxonomy$ncbi_id),]

#Calculate phyla distribution for representative proteins
representatives_phyla = table(representatives_taxonomy$phylum)

#Calculate phyla distribution for all proteins
all_phyla = table(afdb_taxonomy$phylum)

#Match orders
representatives_phyla = representatives_phyla[match(names(all_phyla),
                                                    names(representatives_phyla))]
names(representatives_phyla) = names(all_phyla)

#Replace NAs with 0
representatives_phyla[is.na(representatives_phyla)] = 0

#Calculate completeness via phylogenetic distance with 'clade_PD'
pds_rep = clade_PD(tree_taxonomy,
                   representatives_taxonomy,
                   tree,
                   taxa1 = 'phylum',
                   taxa2 = 'family',
                   verbose = TRUE,
                   show.tip.label = FALSE)

#Calculate completeness via phylogenetic distance for PDB with 'clade_PD'
pds_pdb = clade_PD(tree_taxonomy,
                   pdb_taxonomy,
                   tree,
                   taxa1 = 'phylum',
                   taxa2 = 'family',
                   verbose = TRUE,
                   show.tip.label = FALSE)

plot(pds_rep$contmap, 
     ftype = 'off',
     mar=c(1.1,4.1,4.1,2.1))

#Extract weighted phylogenetic distances
weighted_pds_all = unlist(lapply(pds_fam$phylogenetic_distance, 
                                 function(x) x$weighted_PD))

weighted_pds_rep = unlist(lapply(pds_rep$phylogenetic_distance, 
                                 function(x) x$weighted_PD))

weighted_pds_pdb = unlist(lapply(pds_pdb$phylogenetic_distance, 
                                 function(x) x$weighted_PD))

#Generate plot
par(mfrow = c(2,2))
plot(pds_fam$contmap, 
     ftype = 'off',
     mar=c(1.1,4.1,4.1,2.1))
plot(pds_rep$contmap, 
     ftype = 'off',
     direction = 'leftwards',
     mar=c(1.1,4.1,4.1,2.1))

par(mar = c(8.1, 4.1, 6.1, 2.1))
#Plot joint distribution of phyla totals
plot(log(as.numeric(all_phyla)+1),
     log(as.numeric(representatives_phyla)+1),
     pch = 20,
     col = 'gray60',
     cex = 2,
     xlab = 'AFDB',
     ylab = 'Foldseek',
     cex.axis = 1.5,
     cex.lab = 1.5)
abline(lm(log(as.numeric(representatives_phyla)+1)~
            log(as.numeric(all_phyla)+1)),
       lty = 'dashed',
       lwd = 1.5)
title(main = 'n proteins per phylum (log)',
      font.main = 1,
      cex.main = 1.5)

#Correlation
cor(as.numeric(all_phyla),
    as.numeric(representatives_phyla))

#Plot joint distribution of taxonomic completeness
plot(as.numeric(weighted_pds_all),
     as.numeric(weighted_pds_rep),
     pch = 20,
     col = 'gray60',
     cex = 2,
     xlab = 'AFDB',
     ylab = 'Foldseek',
     cex.axis = 1.5,
     cex.lab = 1.5)
abline(lm(as.numeric(weighted_pds_rep)~
            as.numeric(weighted_pds_all)),
       lty = 'dashed',
       lwd = 1.5)
title(main = 'Taxonomic completeness',
      font.main = 1,
      cex.main = 1.5)

#Correlation
cor(as.numeric(weighted_pds_all),
    as.numeric(weighted_pds_rep))

###################################################################
#####Figure 5: pLDDT as a function of Foldseek representatives#####
###################################################################
#Add ncbi id to cluster_stats
cluster_stats$ncbi = dat$taxonomy_ID[match(cluster_stats$cluster_ID,
                                           dat$member_ID)]

#Get taxa for representative proteins
representatives = dat$taxonomy_ID[match(cluster_stats$cluster_ID, dat$member_ID)]

#Get taxonomic IDs for each
representatives_taxonomy = afdb_taxonomy[match(representatives,
                                               afdb_taxonomy$ncbi_id),]

#Add phylum to cluster stats
cluster_stats$phylum = representatives_taxonomy$phylum[match(cluster_stats$ncbi,
                                                             representatives_taxonomy$ncbi_id)]
cluster_stats$species = representatives_taxonomy$species[match(cluster_stats$ncbi,
                                                               representatives_taxonomy$ncbi_id)]

#Get n proteins per species
species_n = lapply(split(dat, 
                         dat$taxonomy_ID),
                   function(x) nrow(x))

#Split
s = split(cluster_stats,
          cluster_stats$ncbi)

#Match
species_n = species_n[match(names(s),
                            names(species_n))]

#Sweep protein n cutoffs and calculate correlation
cors = c()
for(i in 1:2500){
  print(i)
  #Filter
  s_filter = s[unlist(lapply(s, function(x) nrow(x)))>=i]
  
  #Recombine and simplify
  p = unlist(lapply(s_filter, function(x) mean(x$repPlddt)))
  n = unlist(species_n[match(names(p), 
                             names(species_n))])
  
  #Correlate
  cors = c(cors, cor(n, p, method = 'spearman'))
}

#Plot correlations
par(mfrow = c(1,3))
plot(cors,
     ylab = 'Spearman R',
     xlab = 'n representative proteins per species',
     cex.lab = 1.5,
     cex.axis = 1.5,
     type = 'l',
     lwd = 1.5,
     col = all_colors[1])

#Plot example
#Filter
s_filter = s[unlist(lapply(s, function(x) nrow(x)))>=1000]

#Recombine and simplify
p = unlist(lapply(s_filter, function(x) mean(x$repPlddt)))
n = unlist(species_n[match(names(p), 
                           names(species_n))])

#Plot
cols = all_colors[5:7]
names(cols) = c('Archaea', 'Bacteria', 'Eukaryota')
cols = cols[match(afdb_taxonomy$superkingdom[match(names(n), afdb_taxonomy$ncbi_id)],
                  names(cols))]
plot(n, 
     p,
     ylab = 'pLDDT',
     xlab = 'n representative proteins per species',
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 20,
     col = cols)

#Correlate
cor(n, p, method = 'spearman')

#Compare plddts by kingdom
plddts = split(p, 
               afdb_taxonomy$superkingdom[match(names(n), afdb_taxonomy$ncbi_id)])

#Plot
vioplot::vioplot(plddts,
                 col = all_colors[5:7],
                 side = "right",
                 ylab = '',
                 xlab = '', 
                 font.main = 1,
                 cex.main = 1.5,
                 las = 2, 
                 cex.axis = 1.5,
                 cex.lab = 1.5)
title(ylab = 'pLDDT (representative protein)',
      cex.lab = 1.5)

stripchart(plddts,
           col = all_colors[5:7],
           at = seq(0.8, (length(plddts)-1)+0.8, 1), 
           jitter = 0.1,
           method = "jitter", 
           vertical = TRUE, 
           cex = 1,
           pch = 20, 
           add = TRUE)

#####################################################
#####Figure 6: Protein cluster diversity by taxa#####
#####################################################
#Set up taxa to measure space coverage over
toTest = c("superkingdom", "phylum", 'class', 'order', 'family', 'genus', 'species')

#Calculate total number of clusters
total_clusters = length(unique(dat$` cluster_ID`))

#Split afdb into species
dat_species = split(dat, 
                    dat$taxonomy_ID)

#Calculate coverage based on taxonomic group
taxon_coverage = list()

pb <- txtProgressBar(min = 1,
                     max = length(toTest),
                     style = 3,
                     width = 100,
                     char = ".")
for(i in 1:length(toTest)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Get split taxonomy on taxonomic group
  taxon = split(afdb_taxonomy, afdb_taxonomy[,grep(toTest[i],
                                                   colnames(afdb_taxonomy))])
  
  #Get clusters in taxa
  taxa_clusters = lapply(taxon, function(x) unique(do.call(rbind, dat_species[names(dat_species)%in%
                                                                                x$ncbi_id])$` cluster_ID`))
  
  #Filter
  taxa_clusters = taxa_clusters[unlist(lapply(taxa_clusters, function(x) length(x)))>=5000]
  
  #Loop over and calculate n taxa that clusters occur in
  cl = list()
  for(j in 1:length(taxa_clusters)){
    
    #Extract clusters in taxon
    cl[[names(taxa_clusters)[j]]] = rowSums(do.call(cbind, lapply(taxa_clusters[-j], 
                                                                  function(x) taxa_clusters[[j]]%in%x)))/length(taxa_clusters) 
  
  }
  
  #Add to list
  taxon_coverage[[toTest[i]]] = cl
}

#Compare by kingdom
kingdoms = c('Archaea', 'Bacteria', 'Eukaryota')
kingdom_coverage = list()
for(i in 1:length(kingdoms)){
  
  k = afdb_taxonomy[afdb_taxonomy$superkingdom == kingdoms[i],]
  
  means = list()
  for(j in 1:length(taxon_coverage)){
    x = k[,toTest[j]]
    means[[toTest[j]]] = taxon_coverage[[j]][names(taxon_coverage[[j]])%in%x]
  }
  
  kingdom_coverage[[kingdoms[i]]] = means
  
}

#Plot
par(mfrow = c(1,3))
cols = all_colors[5:7]
for(i in 1:length(kingdom_coverage)){
  
  ses = c()
  means = c()
  for(j in 1:length(kingdom_coverage[[i]])){
    x = unlist(kingdom_coverage[[i]][[j]])
    ses = c(ses, sd(x)/sqrt(length(x)))
    means = c(means, mean(x))
  }
  
  plot(1:length(means),
       means,
       type = 'l',
       ylim = c(0, 0.2),
       ylab = 'Cluster uniqueness',
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
  title(main = names(kingdom_coverage)[i],
        font.main = 1,
        cex.main = 1.5)
}

#Half violin plots
par(mfrow = c(1,4))
for(i in 1:length(all)){
  vioplot::vioplot(lapply(kingdom_coverage[[i]], function(x) unlist(x)[sample(1:length(unlist(x)), 1000)]),
                   col = 'gray80',
                   side = "right",
                   ylab = '% of protein space',
                   xlab = '', 
                   font.main = 1,
                   cex.main = 1.5,
                   las = 2, 
                   ylim = c(0, 1),
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

#############################################################
#####Phylogenetic distance covered by dataset partitions#####
#############################################################
#Split afdb on species
dat_species = split(dat, dat$taxonomy_ID)

#Calculate species n
species_n = unlist(lapply(dat_species, function(x) nrow(x)))

#Function to calculate PD for species with at least n proteins in dataset
pd_per_partition = function(species_numbers,
                            tree,
                            taxonomy,
                            partition_size = 1){
  
  #Filter species
  species = species_numbers[species_numbers>=partition_size]
  
  #Filter to just species in tree
  species = taxonomy$species[match(names(species), 
                                   taxonomy$ncbi_id)]
  species = species[species%in%tree$tip.label]
  
  #ID MRCA
  mrca = getMRCA(tree, 
                 species)
  
  #Extract clade
  clade = extract.clade(tree, 
                        mrca)
  
  #Filter
  clade = keep.tip(clade, 
                   species)
  
  #Get afdb in taxa
  a_tips = clade$tip.label%in%species
  
  #Create community vector(s)
  comm = a_tips
  comm[comm == TRUE] = 1
  comm[comm == FALSE] = 0
  comm = as.data.frame(t(comm))
  #comm = rbind(comm, rep(1, ncol(comm)))
  colnames(comm) = clade$tip.label
  rownames(comm) = 1:nrow(comm)
  
  #Phylogenetic alpha diversity
  PD = pd(samp = comm, 
          tree = clade,
          include.root = TRUE)
  
  #Return
  return(list(PD = PD,
              taxonomy = taxonomy[match(species, 
                                        taxonomy$species),]))
}

#Calculate over range of protein sizes
pds = list()
test = c(1, seq(100, 20000, 100))
for(i in 1:length(test)){
  print(test[i])
  pds[[as.character(test[i])]] = pd_per_partition(species_n,
                                                  tree,
                                                  afdb_taxonomy,
                                                  test[i])
}

#Extract partition size
partitions = unlist(lapply(pds, function(x) x$PD[1,1]))

#Plot percentage of phylogenetic diversity captured
par(mfrow = c(1,3))
cols = colorRampPalette(unlist(arcadia_magma$color_dict))(length(seq(0, 100, 0.1)))
names(cols) = seq(0, 100, 0.1)
cols = cols[match(round(partitions/max(partitions)*100, 1),
                  names(cols))]
plot(names(pds),
     partitions/max(partitions)*100,
     xlab = 'Partition size (n proteins)',
     ylab = 'Phylogenetic diversity (%)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     col = cols,
     cex = 1.5,
     pch = 20,
     yaxt = 'n',
     ylim = c(0, 100))
axis(2, seq(0, 100, 20), seq(0, 100, 20), cex.axis = 1.5, las = 2)

#Plot percentage of distinct phyla
taxa_dist = unlist(lapply(pds, 
                         function(x) length(unique(x$taxonomy$species))))

cols = all_colors[(length(all_colors)-7):length(all_colors)]
plot(names(pds),
     (taxa_dist/max(taxa_dist))*100,
     xlab = 'Partition size (n proteins)',
     ylab = 'Taxa diveristy (%)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 2,
     type = 'l',
     col = cols[1],
     yaxt = 'n',
     ylim = c(0, 100))
axis(2, seq(0, 100, 20), seq(0, 100, 20), cex.axis = 1.5, las = 2)

for(i in 2:6){
  taxa_dist = unlist(lapply(pds, 
                            function(x) length(unique(x$taxonomy[,i]))))
  lines(names(pds),
        (taxa_dist/max(taxa_dist))*100,
        lwd = 2,
        col = cols[i])
}

text(20000, 95, colnames(tree_taxonomy)[1], col = cols[1], adj = 1)
text(20000, 90, colnames(tree_taxonomy)[2], col = cols[2], adj = 1)
text(20000, 85, colnames(tree_taxonomy)[3], col = cols[3], adj = 1)
text(20000, 80, colnames(tree_taxonomy)[4], col = cols[4], adj = 1)
text(20000, 75, colnames(tree_taxonomy)[5], col = cols[5], adj = 1)
text(20000, 70, colnames(tree_taxonomy)[6], col = cols[6], adj = 1)

#Amount of protein cluster space covered
all = length(unique(dat$` cluster_ID`[dat$taxonomy_ID%in%pds[[1]]$taxonomy$ncbi_id]))

cluster_space_partitions = c()
for(i in 1:length(pds)){
  print(i)
  d = length(unique(dat$` cluster_ID`[dat$taxonomy_ID%in%pds[[i]]$taxonomy$ncbi_id]))
  cluster_space_partitions = c(cluster_space_partitions, d/all)
  print(cluster_space_partitions)
}

cluster_space_partitions = cluster_space_partitions*100

cols = colorRampPalette(unlist(arcadia_magma$color_dict))(length(seq(0, 100, 0.1)))
names(cols) = seq(0, 100, 0.1)
cols = cols[match(round(cluster_space_partitions, 1),
                  names(cols))]
plot(names(pds),
     cluster_space_partitions,
     xlab = 'Partition size (n proteins)',
     ylab = '% of cluster space',
     cex.axis = 1.5,
     cex.lab = 1.5,
     col = cols,
     cex = 1.5,
     pch = 20,
     yaxt = 'n',
     ylim = c(0, 100))
axis(2, seq(0, 100, 20), seq(0, 100, 20), cex.axis = 1.5, las = 2)

#Difference distribution
plot(names(pds)[2:length(pds)],
     abs(diff(partitions))/max(partitions),
     xlab = 'n proteins in AFDB',
     ylab = 'Proportion of phylogenetic distance (derivative)',
     cex.axis = 1.5,
     cex.lab = 1.5,
     type = 'l',
     lwd = 1.5)


#Simulate phylogenetic distance given 80/20 training/test partitioning given different min protein n
#Simulate of different min protein n
train_test_simulation = list()
for(h in 1:10){
  
  print(h)
  #Extract taxonomy
  taxa = pds[[h]]$taxonomy
  
  #Rep species names n times (to match simulated partition size)
  data_all = rep(taxa$species, each = names(pds)[h])
  
  train_pds = c()
  train_phyla = c()
  for(i in 1:10){
    
    print(i)
    
    #Randomly split into 80/20
    train = sample(1:length(data_all), length(data_all)*0.8)
    train = unique(data_all[train])
    
    #ID MRCA
    mrca = getMRCA(tree, 
                   train)
    
    #Extract clade
    clade = extract.clade(tree, 
                          mrca)
    
    #Filter
    clade = keep.tip(clade, 
                     train)
    
    #Get afdb in taxa
    a_tips = clade$tip.label%in%train
    
    #Create community vector(s)
    comm = a_tips
    comm[comm == TRUE] = 1
    comm[comm == FALSE] = 0
    comm = as.data.frame(t(comm))
    #comm = rbind(comm, rep(1, ncol(comm)))
    colnames(comm) = clade$tip.label
    rownames(comm) = 1:nrow(comm)
    
    #Phylogenetic alpha diversity
    PD = pd(samp = comm, 
            tree = clade,
            include.root = TRUE)
    
    #Add to results
    train_pds = c(train_pds, PD[1,1])
    train_phyla = length(unique(tree_taxonomy$phylum[match(train, 
                                                           tree_taxonomy$species)]))
    print(train_pds)
    print(train_phyla)
    
  }
  
  train_test_simulation[[h]] = list(pds = train_pds,
                                    phyla = train_phyla)
}

#PLDDT of representative proteins by taxa
#Add ncbi id to cluster_stats
cluster_stats$ncbi = dat$taxonomy_ID[match(cluster_stats$cluster_ID,
                                           dat$member_ID)]

#Match phyla to ncbi ids
taxa = afdb_taxonomy$order[match(cluster_stats$ncbi,
                                  afdb_taxonomy$ncbi_id)]

#Split
cluster_stats_taxa = split(cluster_stats, taxa)

#Calculate mean plddt
mean_plddts = unlist(lapply(cluster_stats_taxa, function(x) mean(x$avgPlddt)))

#Compare to n proteins
plot(unlist(lapply(cluster_stats_taxa, function(x) nrow(x))),
     mean_plddts)








