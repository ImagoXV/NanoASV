#NanoASV phylosequisation
#Arthur Cousson - 2023
#Contact : arthur.cousson@ird.fr

#Workaround conflicting host vs local libraries
.libPaths(c(file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")))

# usage: Rscript script.r INPUT_DIR OUTPUT_DIR R_CLEANING TREE METADATA_DIR
# (INPUT_DIR is not used)
args <- commandArgs(trailingOnly = TRUE)
results_dir <- file.path(args[2], "Results")
r_cleaning <- as.integer(args[3]) == 1
with_tree <- as.integer(args[4]) == 1
metadata_dir <- args[5]

# Chloroplasts and mitochondria (mainly) have more ranks than the seven
# standard ones: extra fields are kept until the cleaning step
TAXONOMY_COLUMNS <- c("Kingdom", "Phylum", "Class", "Order", "Family",
                      "Genus", "Species", paste0("other", 1:6),
                      paste0("others", c(7:17, 19)))
STANDARD_RANKS <- TAXONOMY_COLUMNS[1:7]


read_metadata <- function(directory) {
  read.csv(file.path(directory, "metadata.csv"), row.names = 1,
           header = TRUE, check.names = FALSE)
}

# vsearch cluster table, and a taxonomy of "Unknown" for each cluster;
# NULL when no cluster was kept
read_unknown_clusters <- function(results_dir) {
  file <- file.path(results_dir, "Unknown_clusters", "unknown_clusters.tsv")
  if (!file.exists(file)) {
    return(NULL)
  }
  abundances <- read.csv(file, sep = "\t", header = TRUE, row.names = 1,
                         check.names = FALSE)
  taxonomy <- as.data.frame(matrix("Unknown",
                                   nrow = nrow(abundances),
                                   ncol = length(STANDARD_RANKS),
                                   dimnames = list(rownames(abundances),
                                                   STANDARD_RANKS)))
  taxonomy[setdiff(TAXONOMY_COLUMNS, STANDARD_RANKS)] <- NA
  list(abundances = abundances, taxonomy = taxonomy)
}

# one column (the sample), one row per reference sequence, then one row
# per unknown cluster when the sample has reads in the cluster table
read_abundances <- function(results_dir, file, sample, unknown) {
  abundances <- read.csv(file.path(results_dir, "ASV", file), sep = " ",
                         row.names = 2, header = FALSE)
  colnames(abundances) <- sample
  if (!is.null(unknown) && sample %in% colnames(unknown$abundances)) {
    abundances <- rbind(abundances, unknown$abundances[sample])
  }
  abundances
}

# lines look like "ID Kingdom;Phylum;...": the first field holds both the
# reference ID (row name) and the kingdom
read_taxonomy <- function(results_dir, sample, unknown) {
  file <- file.path(results_dir, "Tax", paste0("Taxonomy_", sample, ".csv"))
  taxonomy <- read.csv2(file, sep = ";", header = FALSE, check.names = FALSE,
                        fill = TRUE, col.names = TAXONOMY_COLUMNS)
  id_and_kingdom <- strsplit(taxonomy[, 1], split = " ")
  rownames(taxonomy) <- sapply(id_and_kingdom, "[[", 1)
  taxonomy[, 1] <- sapply(id_and_kingdom, "[", 2)
  if (!is.null(unknown)) {
    taxonomy <- rbind(taxonomy, unknown$taxonomy)
  }
  taxonomy
}

sample_phyloseq <- function(results_dir, file, metadata, unknown) {
  sample <- sub("_abundance\\.tsv$", "", file)
  abundances <- read_abundances(results_dir, file, sample, unknown)
  taxonomy <- read_taxonomy(results_dir, sample, unknown)
  phyloseq::phyloseq(phyloseq::otu_table(abundances, taxa_are_rows = TRUE),
                     phyloseq::tax_table(as.matrix(taxonomy)),
                     phyloseq::sample_data(metadata))
}

remove_organelles_and_eukaryotes <- function(physeq) {
  physeq <- phyloseq::subset_taxa(physeq, Kingdom != "Eukaryota")
  physeq <- phyloseq::subset_taxa(physeq, Family != "Mitochondria")
  physeq <- phyloseq::subset_taxa(physeq, Order != "Chloroplast")
  phyloseq::tax_table(physeq) <- phyloseq::tax_table(physeq)[, STANDARD_RANKS]
  physeq
}


metadata <- read_metadata(metadata_dir)
unknown <- read_unknown_clusters(results_dir)
abundance_files <- list.files(file.path(results_dir, "ASV"), pattern = "*.tsv")

NanoASV <- Reduce(phyloseq::merge_phyloseq,
                  lapply(abundance_files, sample_phyloseq,
                         results_dir = results_dir, metadata = metadata,
                         unknown = unknown))

if (with_tree) {
  tree <- ape::read.tree(file.path(results_dir, "Phylogeny", "ASV.tree"))
  phyloseq::phy_tree(NanoASV) <- phyloseq::phy_tree(tree)
}

if (r_cleaning) {
  NanoASV <- remove_organelles_and_eukaryotes(NanoASV)
}

print("Exporting results")
save(NanoASV, file = file.path(results_dir, "Rdata", "NanoASV.rdata"))

write.csv(data.frame(phyloseq::tax_table(NanoASV), phyloseq::otu_table(NanoASV),
                     check.names = FALSE),
          file = file.path(results_dir, "CSV", "Taxonomy-Abundance_table.csv"))
