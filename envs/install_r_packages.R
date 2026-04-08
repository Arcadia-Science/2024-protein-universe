library(devtools)

devtools::install_github("Arcadia-Science/arcadia-color-brewer")
devtools::install_version("taxonomizr", version = "0.10.2")
devtools::install_version("packcircles", version = "0.3.6")
devtools::install_version("rsqlite", version = "2.3.0")
devtools::install_version("sm", version = "2.2")
devtools::install_version("taxizedb", version = "0.3.1")
devtools::install_version("vioplot", version = "0.5.0")

# Install Bioconductor package for NCBI genome statistics (used in protein-universe-data.R)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("genomes")
