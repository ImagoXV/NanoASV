#NanoASV phylosequisation
#Arthur Cousson - 2023
#Contact : arthur.cousson@ird.fr
#Workaround conflicting host vs local libraries

.libPaths(c("/usr/local/lib/R/site-library","/usr/lib/R/site-library","/usr/lib/R/library"))

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]
OUTPWD <- args[2]
R_CLEANING <- args[3]
TREE <- args[4]
METADATA <- args[5]

#Phylosequization -----
library(phyloseq)
library(ape) #To handle trees

metadata <- read.csv(paste0(METADATA,"/metadata.csv"), row.names = 1, header = TRUE, check.names = FALSE)
barcodes <- rownames(metadata)

if (TREE ==1){
  ASV.tree <- ape::read.tree(paste0(OUTPWD, "/Results/Phylogeny/ASV.tree"))
}

##Unknown OTUs ----

if (file.exists(paste0(OUTPWD,"/Results/Unknown_clusters/unknown_clusters.tsv"))) {
U_OTU <- read.csv(paste0(OUTPWD,"/Results/Unknown_clusters/unknown_clusters.tsv"), sep = "\t", header = T, row.names = 1)

#U_OTU <- U_OTU[rowSums(U_OTU) > 5,] #Remove clusters with total abundance inferior to five

unknown_taxonomy <- rep("Unknown", times = 7) #Adapt taxonomy dataframe so it fits later on with bad entries for cleaning

###inpu %>% read_csv() %>% filter(rowSums()>5) -> U_OTU

###Unknown taxonomy ----
#Taking care of unknown taxonomy table
U_TAX <- data.frame(matrix(nrow = nrow(U_OTU), ncol = 7))
rownames(U_TAX) <- rownames(U_OTU)
colnames(U_TAX) <- c("Kingdom", "Phylum", "Class",
                     "Order", "Family", "Genus", "Species")
U_TAX[,] <- unknown_taxonomy

U_TAX <- data.frame(U_TAX,
                    other1 = NA, 
                    other2 = NA, 
                    other3 = NA, 
                    other4 = NA, 
                    other5 = NA, 
                    other6 = NA,
                    others7 = NA, 
                    others8 = NA, 
                    others9 = NA, 
                    others10 = NA, 
                    others11 = NA, 
                    others12 = NA, 
                    others13 = NA, 
                    others14 = NA, 
                    others15 = NA, 
                    others16 = NA, 
                    others17 = NA, 
                    others19 = NA)
}

##ASV tables ----
#Individuals ASV tables loading
temp_ASV = list.files(path = paste0(OUTPWD,"/Results/ASV/"), pattern = "*.tsv")
for (i in 1:length(temp_ASV)) {
  if (file.size(paste(OUTPWD,"/Results/ASV/", temp_ASV[i], sep = "")) == 0) {
    assign(temp_ASV[i], data.frame())
  } else {
    assign(temp_ASV[i], read.csv(file = paste(OUTPWD,"/Results/ASV/", temp_ASV[i], sep = ""), 
                                 sep = " ", 
                                 row.names = 2, header = FALSE))
  }
}


temp_ASV <- mget(temp_ASV) #To get files as objects from names

names(temp_ASV) <- lapply(strsplit(names(temp_ASV), split = "_"), function(x) x[[1]])

barcodes <- names(temp_ASV)

if(file.exists(paste0(OUTPWD,"/Results/Unknown_clusters/unknown_clusters.tsv"))){
  for (i in 1:ncol(U_OTU)){
    for (j in 1:length(temp_ASV)){
      if (colnames(U_OTU)[i] == names(temp_ASV[j])) {
        colnames(temp_ASV[[j]]) <- barcodes[j]
        temp_ASV[[j]] <- rbind(data.frame(temp_ASV[j]), U_OTU[i])
      }
    }
  }
}

#The following function will deal with empty barcodes (like blanks)
for(i in 1:length(temp_ASV)) {
  if (lapply(temp_ASV[i], function(df) nrow(df)) == 0) {
    temp_ASV[i] <- list(data.frame(V1 = rep(0, times = 1)))
    rownames(temp_ASV[[i]]) <- rownames(temp_ASV[[i-1]])[nrow(temp_ASV[[i-1]])]
  }
}

for (i in 1:length(temp_ASV)){
  if (rownames(temp_ASV[[i]])[1] == "*"){
    rn <- rownames(temp_ASV[[i]])[-1]
    temp_ASV[[i]] <- data.frame(temp_ASV[[i]][-1,], check.rows = F, check.names = F)
    rownames(temp_ASV[[i]]) <- rn
  }
}

for (i in 1:length(temp_ASV)) colnames(temp_ASV[[i]]) <- barcodes[i]

names(temp_ASV) <- barcodes[1:length(temp_ASV)]

##ASV Taxonomy ----
#Individual taxonomy tables loading
temp_TAX = list.files(path = paste0(OUTPWD,"/Results/Tax/"), pattern="*.csv")
for (i in 1:length(temp_TAX)) assign(temp_TAX[i], data.frame(read.csv2(file = paste(OUTPWD,"/Results/Tax/",temp_TAX[i], sep = ""), sep = ";", header = F, check.names = F, fill = TRUE, 
                                                                       col.names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
                                                                                     "Species", "other1", "other2", "other3", "other4", "other5", "other6",
                                                                                     "others7", "others8", "others9", "others10", "others11", "others12", "others13", "others14", "others15", "others16", "others17", "others19"))))
#The reason we allow so many fields is that Chloroplasts and Mitochondria (mainly), have non-coherent taxonomy sizes compared to Bacterial full taxonomy.
#To avoid any "Out of field" errors or alike, we allow them to exist in a first time

temp_TAX <- mget(temp_TAX)

for(i in 1: length(temp_TAX)){
  if (lapply(temp_TAX[i], function(df) nrow(df)) != 0) {#This step allows to correct for Kingdom, which is currently merged with SILVA ID
    rownames(temp_TAX[[i]]) <- sapply(strsplit(temp_TAX[[i]][,1], split = " "), "[[", 1)
    temp_TAX[[i]][,1] <- sapply(strsplit(temp_TAX[[i]][,1], split = " "), "[", 2)
    if(file.exists(paste0(OUTPWD,"/Results/Unknown_clusters/unknown_clusters.tsv"))){
      temp_TAX[[i]] <- rbind(temp_TAX[[i]], U_TAX)
    }
  }
}

for(i in 1: length(temp_TAX)) {
  if (nrow(temp_TAX[[i]]) == 0) {
    temp_TAX[[i]][1,] <- temp_TAX[[i-1]][1,]
    rownames(temp_TAX[[i]]) <- rownames(temp_ASV[[i-1]])[length(rownames(temp_ASV[[i-1]]))]
  }
}

temp_phyloseq <- barcodes

##Phyloseq objects ----

physeq_list <- list()

for (i in 1:length(temp_ASV)) {
  # Create the phyloseq object
  physeq_object <- phyloseq::phyloseq(phyloseq::otu_table(temp_ASV[[i]], taxa_are_rows = TRUE), 
                            phyloseq::tax_table(as.matrix(temp_TAX[[i]])), 
                            phyloseq::sample_data(metadata))
  
  # Get the name from the barcodes vector
  barcode_name <- barcodes[i]
  # Assign the phyloseq object to the list with the dynamic name
  physeq_list[[barcode_name]] <- physeq_object
}

i<-1 #Reset the incrementation
#Initialize the phyloseq object
NanoASV <- phyloseq::merge_phyloseq(physeq_list[[i]], physeq_list[[i + 1]])

# If more than 2 samples, then, adding them all together
if (length(physeq_list) > 2) {
  for (i in 3:length(physeq_list)) {
    NanoASV <- phyloseq::merge_phyloseq(NanoASV, physeq_list[[i]])
  }
}

##Phylogeny ----
if(TREE == 1){
phyloseq::phy_tree(NanoASV) <- phyloseq::phy_tree(ASV.tree)
}
##Dataset cleaning ----
#Delete bad entries such as Eukaryota, Cyanobacteria and Archea if any
if(R_CLEANING == 1){
  NanoASV <- phyloseq::subset_taxa(NanoASV, Kingdom != "Eukaryota")
  NanoASV <- phyloseq::subset_taxa(NanoASV, Family != "Mitochondria")
  NanoASV <- phyloseq::subset_taxa(NanoASV, Order != "Chloroplast")
  ##Taxonomy cleaning ----
  #After those functions, there is no more taxa with mixed up names so we can remove supp fields of taxa table
  tax_table(NanoASV) <- phyloseq::tax_table(NanoASV)[,1:7]
}

#Phyloseq export ----
print("Exporting results")
save(NanoASV, file = paste0(OUTPWD,"/Results/Rdata/NanoASV.rdata"))

#CSV export -----
#Create the TAXOTU file, encompassing taxonomy and abundance
TAXOTU <- data.frame(NanoASV@tax_table, NanoASV@otu_table)
#Exporting the file for people needing this format as well
write.csv(TAXOTU, file = paste0(OUTPWD,"/Results/CSV/Taxonomy-Abundance_table.csv"))
