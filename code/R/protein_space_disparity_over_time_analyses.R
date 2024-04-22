library(geiger)
library(phytools)
library(btw)
library(dispRity)
library(RPANDA)
library(ape)

###################
#####Functions#####
###################
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

###################
#####Load data#####
###################
#Source supporting functions
source('~/Documents/Research/github/afdb-phylo/01_code/R/afdb-phylo-utils.R')
source('~/Documents/Research/github/comparative-amoeboid-crawling/01_utils/plotting_functions.R')
source('~/Documents/Research/github/comparative-amoeboid-crawling/01_utils/treble_functions.R')
source('~/Documents/Research/Rfiles/arcadia_gradients.R')
source('~/Documents/Research/github/afdb-phylo/01_code/R/afdb-phylo-utils.R')

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

##################################################
#####ID protein clusters related to each taxa#####
##################################################
#Create filtered tree corresponding to phyla
phyla_tree = filter_tree(tree,
                         resolution = 'class',
                         verbose = TRUE,
                         change_tree_labels = TRUE, 
                         return_taxonomy = TRUE)

#Filter proteins to just species in tree
dat_proteome = dat[dat$taxonomy_ID%in%afdb_taxonomy$ncbi_id[afdb_taxonomy$class%in%phyla_tree$tree$tip.label],]

#Select just proteins from species with reference proteomes
dat_proteome = dat_proteome[dat_proteome$taxonomy_ID%in%uniprot_reference_proteomes$Organism.Id,]

#Filter to just species with "close to complete" proteomes
protein_n = table(dat_proteome$taxonomy_ID)
protein_n = protein_n/uniprot_reference_proteomes$Protein.count[
  match(names(protein_n), uniprot_reference_proteomes$Organism.Id)]
protein_n = protein_n[!is.na(protein_n)]
proteomes = protein_n[protein_n>0.25]

sort(table(afdb_taxonomy$class[match(names(proteomes), 
                                      afdb_taxonomy$ncbi_id)]))

#Refilter protein list
dat_proteome = dat_proteome[dat_proteome$taxonomy_ID%in%names(proteomes),]

#Filter clusters
clusters_filter = cluster_stats[cluster_stats$avgPlddt>80,]
clusters_filter = clusters_filter[clusters_filter$cluster_ID%in%unique(dat_proteome$` cluster_ID`),]
clusters_filter = clusters_filter[clusters_filter[,3]>100,]

#Split clusters by order
clusters_by_taxa = split(dat_proteome, 
                         afdb_taxonomy[,grep('class', colnames(afdb_taxonomy))][
                           match(dat_proteome$taxonomy_ID, afdb_taxonomy$ncbi_id)])

#Generate binary matrix of cluster presence across taxa
clusters_presence = list()
for(i in 1:length(clusters_by_taxa)){
  
  #Match taxa cluster names to overall list
  clusters_binary = clusters_filter$cluster_ID%in%
    clusters_by_taxa[[i]]$` cluster_ID`
  clusters_binary[clusters_binary == TRUE] = 1
  clusters_binary[clusters_binary == FALSE] = 0
  
  #Add to matrix
  clusters_presence[[names(clusters_by_taxa)[i]]] = clusters_binary
  
}

clusters_presence = do.call(rbind, clusters_presence)
colnames(clusters_presence) = clusters_filter$cluster_ID

#PCA
pca = prcomp(clusters_presence)

#UMAP
u = umap::umap(pca$x[,1:65], 
               verbose = TRUE)

#Plot coloring by superkingdom
s = afdb_taxonomy$superkingdom[match(rownames(u$layout),
                                     afdb_taxonomy$class)]
colors = ArcadiaColorBrewer::arcadia.pal(name = 'Accent',
                                         n = length(unique(s)))
names(colors) = unique(s)
colors = colors[match(s, names(colors))]
plot(u$layout,
     pch = 20,
     col = colors)
text(u$layout,
    rownames(u$layout),
    cex = 0.5)

#Filter tree to match taxa in final data
tree_filter = keep.tip(phyla_tree$tree, 
                       phyla_tree$tree$tip.label[phyla_tree$tree$tip.label%in%rownames(u$layout)])

##########################
#####Phylomorphospace#####
##########################
#Calculate distance
d = dist(clusters_presence)

#PCA
pca = prcomp(d)

#Plot
s = afdb_taxonomy$superkingdom[match(rownames(u$layout),
                                     afdb_taxonomy$class)]
colors = ArcadiaColorBrewer::arcadia.pal(name = 'Accent',
                                         n = length(unique(s)))
names(colors) = unique(s)
colors = colors[match(s, names(colors))]
plot(pca$x[,1:2],
     pch = 20,
     col = colors)
text(u$layout,
     rownames(u$layout),
     cex = 0.5)

#Phenogram
phenogram(tree_filter,
          pca$x[,1],
          type = "off",
          spread.labels = FALSE)

##########################################################
######Disparity through time for all protein clusters#####
##########################################################
#PCA
pca = prcomp(clusters_presence)

#LTT
all_ltt = ltt(tree_filter)

#DTT
all_dtt = dtt(tree_filter, 
              pca$x[,1])

#Regress?
mod = lm(all_dtt$dtt~log(all_ltt$ltt)[1:length(all_dtt$dtt)])
plot(all_ltt$times[1:length(all_dtt$dtt)],
     mod$residuals, 
     type = 'l')

for(i in 1:length(evo_events)){
  abline(v = max(all_ltt$times[1:length(all_dtt$dtt)])-(evo_events[i]*1000))
}

#Spiral plot
d = mod$residuals
times = all_ltt$times[1:length(all_dtt$dtt)]
times_toplot = c(0, 500, 1000, 1500, 2000, 2500, 3000, max(times))
cols = c(arcadia.pal(n = 6, "Accent"),
         arcadia.pal(n = 2, "Lighter_accents"))
ymax <- max(d)

spiral_initialize()
spiral_track(ylim=c(0, ymax*.7),
             background=FALSE, 
             background_gp = gpar(col = NA, fill = NA))

for(i in 1:length(times_toplot)){
  
  x = times>=times_toplot[i]&times<=times_toplot[i+1]
  xt = times[x]
  xd = d[x]
  
  spiral_polygon(x=c(xt, rev(xt)),
                 y=c(xd/2, -rev(xd/2)),
                 gp = gpar(col=cols[i], fill=cols[i]))
}

spiral_polygon(x=c(all_dtt$times[1:20], rev(all_dtt$times[1:20])),
               y=c(d[1:20]/2, -rev(d[1:20]/2)),
               gp = gpar(col="#d32e2b", fill="#d32e2b50"))
spiral_polygon(x=c(all_dtt$times[21:40], rev(all_dtt$times[21:40])),
               y=c(d[21:40]/2, -rev(d[21:40]/2)),
               gp = gpar(col='darkcyan', fill='darkcyan'))
#spiral_lines(x=all_dtt$times, 
#             y=0)

#############################################
#####Analyze cluster disparity over time#####
#############################################
#Calculate dtt for each cluster
pb <- txtProgressBar(min = 1,
                     max = ncol(clusters_presence),
                     style = 3,
                     width = 100,
                     char = ".")

cluster_dtts = list()
for(i in 1:ncol(clusters_presence)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Dtt
  cluster_dtts[[colnames(clusters_presence)[i]]] = dtt(tree_filter, 
                                                       clusters_presence[,i], 
                                                       plot = FALSE)
}

#Combine into matrix
dtt_matrix = do.call(rbind, lapply(cluster_dtts, function(x) x$dtt))

#Normalize
dtt_matrix = apply(dtt_matrix, 1, function(x) x/max(x))

#PCA
pca = prcomp(t(dtt_matrix))

#Leiden
d = dist(t(dtt_matrix))

#UMAP
u = umap::umap(pca$x[,1:20], verbose = TRUE)

#Plot
plot(u$layout,
     pch = 20,
     #xlim = c(-15, 15),
     #ylim = c(-15, 15),
     cex = 0.5)

#Plot PCs
layout = data.frame(x = u$layout[,1],
                    y = u$layout[,2])
layout = layout[layout$x>(-15),]
layout = layout[layout$x<(15),]
layout = layout[layout$y>(-15),]
layout = layout[layout$y<(15),]

par(mfrow = c(2,3))
for(i in 1:6){
  plot_parameter(pca$x[match(rownames(layout), rownames(pca$x)),i],
                 layout,
                 col = colorRampPalette(unlist(arcadia_poppies$color_dict))(100),
                 n_bins = 200)
}

################################################
#####Compare dtt to environmental variables#####
################################################
#Correlate
cor_dtt_metadata = function(dtt,
                            dtt_tree,
                            metadata,
                            interpolate_n = 2000){
  
  #Interplote dtt
  x = approx(dtt, n = interpolate_n)$y
  
  #Filter metadata on time
  y = metadata[metadata$time<max(dtt_tree$edge.length),]
  
  #Interpolate
  y = approx(y$value, n = interpolate_n)$y
  
  #Get time
  time = seq.int(0, max(dtt_tree$edge.length), length.out = interpolate_n)
  
  #Correlate
  corr = cor(x, y)
  corr_derivative = cor(x[1:(length(x)-1)],
                        diff(y))
  
  #Return
  return(list(correlation = corr,
              correlation_derivative = corr_derivative,
              dtt = x,
              meta = y,
              time = time))
  
}

#CO2
co2_correlations = list()
for(i in 1:ncol(dtt_matrix)){
  co2_correlations[[colnames(dtt_matrix)[i]]] = cor_dtt_metadata(dtt_matrix[,i],
                                                                 tree_filter,
                                                                 co2)
  
}

tail(sort(unlist(lapply(co2_correlations, function(x) x$correlation))))

m = which.min(unlist(lapply(co2_correlations, function(x) x$correlation)))

plot(co2_correlations[[m]]$time,
     co2_correlations[[m]]$dtt/max(co2_correlations[[m]]$dtt), 
     type = 'l')
lines(co2_correlations[[m]]$time,
      co2_correlations[[m]]$meta/max(co2_correlations[[m]]$meta), col = 'red')

#O2
o2_correlations = list()
for(i in 1:ncol(dtt_matrix)){
  o2_correlations[[colnames(dtt_matrix)[i]]] = cor_dtt_metadata(dtt_matrix[,i],
                                                                tree_filter,
                                                                o2)
  
}

tail(sort(unlist(lapply(o2_correlations, function(x) x$correlation))))

m = which.min(unlist(lapply(o2_correlations, function(x) x$correlation)))

plot(o2_correlations[[m]]$time,
     o2_correlations[[m]]$dtt/max(o2_correlations[[m]]$dtt), 
     type = 'l',
     ylab = 'Disparity',
     xlab = "Time (mya)",
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(o2_correlations[[m]]$time,
      o2_correlations[[m]]$meta/max(o2_correlations[[m]]$meta), col = 'red')

human_o2_correlations = unlist(lapply(o2_correlations, function(x) x$correlation))
human_o2_correlations = human_o2_correlations[names(human_o2_correlations)%in%human_clusters[,3]]
h = human_clusters[match(names(human_o2_correlations), human_clusters[,3]),]
h$corr = human_o2_correlations

go = list()
for(i in 1:nrow(h)){
  x = strsplit(h[i,2], ';')[[1]]
  y = rep(h$corr[i], length(x))
  names(y) = x
  go[[i]] = y
}
go = data.frame(go = names(unlist(go)),
                 corr = unlist(go))
write.csv(go, '~/Desktop/test.csv', row.names = FALSE)

library(rrvgo)
simMatrix <- calculateSimMatrix(go$go,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(go$corr, go$go)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores, 
                                threshold=0.5,
                                orgdb="org.Hs.eg.db")

scatterPlot(simMatrix, reducedTerms)

######################################
#####Compare to species diversity#####
######################################
#Load data (from https://github.com/UW-Macrostrat/PNAS_201702297)
dat = read.csv('~/Desktop/RawFossilOccurrences.csv')

#Order on max_ma
dat = dat[order(dat$max_ma, decreasing = TRUE),]

#Set bin range to compare over
bins = (round(max(dat$max_ma))*-1):0

#Calculate total n genera over bin range
bin_genera = list()
range_size = 30
for(i in 1:(length(bins)-range_size)){
  
  #Extract genera for corresponding bin
  genera = dat$genus[(dat$max_ma*-1)>=bins[i]&(dat$max_ma*-1)<=bins[(i+range_size)]]
  
  #Add to list
  bin_genera[[as.character(bins[i])]] = unique(genera)
  
}

#Count
bin_genera_n = unlist(lapply(bin_genera, function(x) length(x)))
plot(bin_genera_n, 
     type = 'l')

#Smooth
bin_genera_smooth = ksmooth(as.numeric(names(bin_genera_n)), 
                            bin_genera_n, 
                            bandwidth = 10)
plot(bin_genera_smooth, 
     type = 'l')

#Rescale dtts
fraction = round(540/3772.466, 2)
cors = c()
ccfs = list()
for(i in 1:length(cluster_dtts)){
  d = cluster_dtts[[i]]$dtt[cluster_dtts[[i]]$times>=(1-fraction)]
  d = approx(d, n = length(bin_genera_smooth$y))
  cors = c(cors, cor(d$y, 
                     bin_genera_smooth$y))
  
  ccfs[[names(cluster_dtts)[i]]] = ccf(d$y, 
                                      bin_genera_smooth$y, 
                                      lag.max = 100, 
                                      plot = FALSE)
}

unlist(lapply(ccfs, function(x) x$acf[,,1][which.max(abs(x$acf[,,1]))]))

######################################################
#####Calculate evolutionary rate for each cluster#####
######################################################
#Run
library(evorates)
fit.evorates(tree_filter, 
             pca$x[,1],
             chains = 1)

out = multirateBM(tree_filter, 
                  pca$x[,1],
                  method = 'REML',
                  optim = 'L-BFGS-B',
                  lambda = 1,
                  parallel = TRUE)

#Ancestral state reconstructions for each cluster (PC)
cluster_anc_state = list()
for(i in 1:50){
  cluster_anc_state[[colnames(clusters_presence)[i]]] = fastAnc(tree_filter, pca$x[,i])
}

###################################################
#####Accumulation of unique clusters over time#####
###################################################
#Loop over clusters and, for each, ID deepest phylogenetic split
pb <- txtProgressBar(min = 1,
                     max = ncol(clusters_presence),
                     style = 3,
                     width = 100,
                     char = ".")

cluster_dates = list()
for(i in 1:ncol(clusters_presence)){
  
  #Update counter
  setTxtProgressBar(pb, i)
  
  #Extract taxa with cluster
  s = clusters_presence[,i][clusters_presence[,i]==1]
  
  #Filter tree to just these taxa
  tmp_tree = keep.tip(tree_filter, names(s))
  
  #Get max time
  max_time = max(dispRity::tree.age(tmp_tree)[,1])
  
  #Add
  cluster_dates[[colnames(clusters_presence)[i]]] = max_time
}
cluster_dates = unlist(cluster_dates)

#Set bin range to compare over
bins = (ceiling(max(cluster_dates))*-1):0

#Calculate total n genera over bin range
bin_clusters = list()
range_size = 100
for(i in 1:(length(bins)-range_size)){
  
  #Extract genera for corresponding bin
  cl = names(cluster_dates)[(cluster_dates*-1)>=bins[i]&(cluster_dates*-1)<=bins[(i+range_size)]]
  
  #Add to list
  bin_clusters[[as.character(bins[i])]] = unique(cl)
  
}

#Add clusters in adjacent bins to each other
for(i in 2:length(bin_clusters)){
  bin_clusters[[i]] = c(bin_clusters[[i-1]], bin_clusters[[i]])
}

#Make unique
bin_clusters_n = lapply(bin_clusters, function(x) unique(x))

#Count
bin_clusters_n = unlist(lapply(bin_clusters_n, function(x) length(x)))

#Plot
plot(bin_clusters_n[20:length(bin_clusters_n)],
     lwd = 1.5,
     type = 'l',
     ylab = 'Unique clusters',
     xlab = 'Billion years ago',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     bty = 'n')




