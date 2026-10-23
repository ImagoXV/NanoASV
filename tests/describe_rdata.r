# Usage: Rscript describe_rdata.r NanoASV.rdata
#
# Print a text description of the phyloseq object saved by NanoASV.
# Abundances and taxonomy are already pinned through the CSV output.

rdata <- commandArgs(trailingOnly = TRUE)[1]
env <- new.env()
load(rdata, envir = env)
NanoASV <- env$NanoASV

cat("taxa:", phyloseq::ntaxa(NanoASV), "\n")
cat("samples:", phyloseq::nsamples(NanoASV), "\n")
cat("ranks:", phyloseq::rank_names(NanoASV), "\n")
cat("sample variables:", phyloseq::sample_variables(NanoASV), "\n")
print(methods::as(phyloseq::sample_data(NanoASV), "data.frame"))
cat("tree tips match taxa:",
    setequal(phyloseq::phy_tree(NanoASV)$tip.label,
             phyloseq::taxa_names(NanoASV)), "\n")
