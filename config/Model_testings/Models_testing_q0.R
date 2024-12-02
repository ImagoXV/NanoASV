library(ggplot2)
library(phyloseq)

setwd("~/Documents/Thesis/NanoASV/Software_dev/NanoASV/config/Tests_with_real_dataset/q0/")

load("map-ont-tree-fixed-q0.out/Results/Rdata/NanoASV.rdata")
mapont <- NanoASV
mapont@sam_data <- sample_data(data.frame(mapont@sam_data,
                                          model = "map-ont"))
load("asm10-tree-fixed-q0.out//Results/Rdata/NanoASV.rdata")
asm10 <- NanoASV
asm10@sam_data <- sample_data(data.frame(asm10@sam_data,
                                          model = "asm10"))
load("asm5-tree-fixed-q0.out//Results/Rdata/NanoASV.rdata")
asm5 <- NanoASV
asm5@sam_data <- sample_data(data.frame(asm5@sam_data,
                                         model = "asm5"))

sample_names(mapont) <- paste0(sample_names(mapont), "_map-ont")
sample_names(asm10) <- paste0(sample_names(asm10), "_asm10")
sample_names(asm5) <- paste0(sample_names(asm5), "_asm5")

mapont@phy_tree <- NULL
asm10@phy_tree <- NULL
asm5@phy_tree <- NULL

Finest <- merge_phyloseq(mapont, asm10, asm5)
metadata <- data.frame(Finest@sam_data)

metadata[rownames(metadata)%in%c("barcode33_map-ont", "barcode33_asm10", "barcode33_asm5"), "Refs"] <- "SFR24-B3"
metadata[rownames(metadata)%in%c("barcode33_map-ont", "barcode33_asm10", "barcode33_asm5"), "B"] <- "B3"

Finest@sam_data <- sample_data(metadata)

# Genus  <- tax_glom(Finest, taxrank = "Genus")
# Family  <- tax_glom(Finest, taxrank = "Family")
# Order  <- tax_glom(Finest, taxrank = "Order")
# Class  <- tax_glom(Finest, taxrank = "Class")
# Phylum  <- tax_glom(Finest, taxrank = "Phylum")
#save(Genus, Family, Order, Class, Phylum, file = "q0_Tax_glom_taxo_levels.rdata")

load("q0_Tax_glom_taxo_levels.rdata")

phy.list <- list(Finest = Finest,
                 Genus = Genus,
                 Family = Family,
                 Order = Order,
                 Class = Class,
                 Phylum = Phylum)


alpha_tab <- data.frame()
j <- 0

for(i in phy.list){
  j <- j+1
  OTU <- data.frame(i@otu_table)
  rich <- data.frame(apply(OTU>0,2,sum))
  Shan <- data.frame(vegan::diversity(t(OTU), index="shannon"))
  df <- data.frame( richness = rich,
                           shannon = Shan)
  colnames(df) <- c("richness", "shannon")
  metadata <- i@sam_data
  metadata$level <- paste0(names(phy.list)[j])
  #metadata$Barcode <- c("Barcode_01", "Barcode_02")
  df <- cbind(metadata, df)
  df$model <- paste0(metadata$model)
  df$Refs <- paste0(df$model, "_", df$level)
  rownames(df) <- paste0(names(phy.list)[j], "_",rownames(df))
  alpha_tab <- rbind(alpha_tab, df)
}

library(reshape2)

alpha.melted <- melt(alpha_tab)

alpha.melted$model <- factor(alpha.melted$model, levels = c("map-ont", "asm10", "asm5"))
alpha.melted$level <- factor(alpha.melted$level, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Finest"))

pdf("q0_Numerical_richness.pdf", he =18, wi = 18)
ggplot(data = alpha.melted[alpha.melted$variable == "richness",], aes(x = reorder(SFR, value) , y = value, color = alpha_tab$SFR)) +
  facet_grid(rows = vars(level), cols = vars(model), scales = "free_y", switch = "y") +  # Free y-axis for each row
  geom_boxplot() +
  geom_point() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        legend.position = "none") + 
  ylab("Numerical richness") + xlab("Model") + 
  labs(title = "Numerical richness\n q0 filter",
       caption = date())
dev.off()

pdf("q0_shannon_index", he =18, wi = 18)
ggplot(data = alpha.melted[alpha.melted$variable == "shannon",], aes(x =  reorder(SFR, value) , y = value, color = alpha_tab$SFR)) +
  facet_grid(rows = vars(level), cols = vars(model), scales = "free_y", switch = "y") +  # Free y-axis for each row
  geom_boxplot() +
  geom_point() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        legend.position = "none") + 
  ylab("Shannon index") + xlab("Model") + 
  labs(title = "Shannon index\n q0 filter",
       caption = date())
dev.off()


#To count raw hits and affiliation with models

smp_sums <- as.data.frame(sample_sums(Finest))
#smp_sums$Barcode <- factor(rownames(smp_sums), level = rownames(smp_sums))
colnames(smp_sums)[1] <- "Reads"

metadata <- data.frame(Finest@sam_data)

smp_sums <- merge(smp_sums, metadata, by = 0)

pdf("q0_Sample_sums_barplot.pdf", he = 6, wi = 12)
ggplot(data = smp_sums, aes(x = Refs, y = Reads)) +
  geom_bar(stat="identity") + 
  facet_wrap(~model) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=8))
dev.off()

# Taxonomical composition ----

couleurs_genus <- c("#FF0000FF","#FF9900FF","#FFCC00FF","#00FF00FF","#6699FFFF",
                    "#CC33FFFF","#99991EFF","#FF00CCFF", "#999999FF","#CC0000FF",
                    "#FFCCCCFF","#FFFF00FF","#CCFF00FF","#358000FF","#0000CCFF",
                    "#99CCFFFF","#00FFFFFF","#CCFFFFFF","#9900CCFF","#CC99FFFF",
                    "#996600FF","#79CC3DFF","#CCCCCCFF","#79CC3DFF","#CCCC99FF","black")

#Composition ----

pdf("q0_Composition_with_different_models.pdf", he = 8, wi = 16)

Genus.norm <- transform_sample_counts(Genus, function(x) x/sum(x))
Genus.melted <- psmelt(Genus.norm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR != "Bl",] #Remove blanks
sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR != "Mock",] #Remove positive control mock
sub_Genus.melted[sub_Genus.melted$Abundance<0.014,32:38] <- "Z_Others"

ggplot(sub_Genus.melted, aes(x = Refs, y = Abundance, fill = Genus)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Genus-level composition",
       #subtitle = "",
       caption = date()) 


Family.norm <- transform_sample_counts(Family, function(x) x/sum(x))
Family.melted <- psmelt(Family.norm)
sub_Family.melted <- Family.melted
sub_Family.melted <- sub_Family.melted[sub_Family.melted$SFR != "Bl",] #Remove blanks
sub_Family.melted <- sub_Family.melted[sub_Family.melted$SFR != "Mock",] #Remove positive control mock
sub_Family.melted[sub_Family.melted$Abundance<0.02,32:38] <- "Z_Others"

ggplot(sub_Family.melted, aes(x = Refs, y = Abundance, fill = Family)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Family-level composition",
       #subtitle = "",
       caption = date()) 



Order.norm <- transform_sample_counts(Order, function(x) x/sum(x))
Order.melted <- psmelt(Order.norm)
sub_Order.melted <- Order.melted
sub_Order.melted <- sub_Order.melted[sub_Order.melted$SFR != "Bl",] #Remove blanks
sub_Order.melted <- sub_Order.melted[sub_Order.melted$SFR != "Mock",] #Remove positive control mock
sub_Order.melted[sub_Order.melted$Abundance<0.022,32:38] <- "Z_Others"

ggplot(sub_Order.melted, aes(x = Refs, y = Abundance, fill = Order)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Order-level composition",
       #subtitle = "",
       caption = date()) 



Class.norm <- transform_sample_counts(Class, function(x) x/sum(x))
Class.melted <- psmelt(Class.norm)
sub_Class.melted <- Class.melted
sub_Class.melted <- sub_Class.melted[sub_Class.melted$SFR != "Bl",] #Remove blanks
sub_Class.melted <- sub_Class.melted[sub_Class.melted$SFR != "Mock",] #Remove positive control mock
sub_Class.melted[sub_Class.melted$Abundance<0.01,32:38] <- "Z_Others"

ggplot(sub_Class.melted, aes(x = Refs, y = Abundance, fill = Class)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Class-level composition",
       #subtitle = "",
       caption = date()) 



Phylum.norm <- transform_sample_counts(Phylum, function(x) x/sum(x))
Phylum.melted <- psmelt(Phylum.norm)
sub_Phylum.melted <- Phylum.melted
sub_Phylum.melted <- sub_Phylum.melted[sub_Phylum.melted$SFR != "Bl",] #Remove blanks
sub_Phylum.melted <- sub_Phylum.melted[sub_Phylum.melted$SFR != "Mock",] #Remove positive control mock
sub_Phylum.melted[sub_Phylum.melted$Abundance<0.0022,32:38] <- "Z_Others"

ggplot(sub_Phylum.melted, aes(x = Refs, y = Abundance, fill = Phylum)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Phylum-level composition",
       #subtitle = "",
       caption = date()) 

dev.off()















Genus.norm <- transform_sample_counts(Genus, function(x) x/sum(x))
Genus.melted <- psmelt(Genus.norm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR %in% c("Bl", "Mock"),] #Remove blanks
#sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR != "Mock",] #Remove positive control mock
sub_Genus.melted[sub_Genus.melted$Abundance<0.05,32:38] <- "Z_Others"

pdf("q30_Controls_relative_composition.pdf", he = 6, wi = 18)
ggplot(sub_Genus.melted, aes(x = Refs, y = Abundance, fill = Genus)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Genus-level composition",
       #subtitle = "",
       caption = date()) 
dev.off()


Genus.unnorm <-Genus
Genus.melted <- psmelt(Genus.unnorm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR %in% c("Bl", "Mock"),] #Remove blanks
#sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR != "Mock",] #Remove positive control mock
sub_Genus.melted[sub_Genus.melted$Abundance<11,32:38] <- "Z_Others"

pdf("q30_Controls_Absolute_composition.pdf", he = 6, wi = 18)
ggplot(sub_Genus.melted, aes(x = Refs, y = Abundance, fill = Genus)) +
  theme_bw() +
  facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Genus-level composition",
       #subtitle = "",
       caption = date()) 
dev.off()



Genus.norm <- transform_sample_counts(Genus, function(x) x/sum(x))
Genus.melted <- psmelt(Genus.norm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR %in% c("Mock"),] #Remove blanks
#sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR != "Mock",] #Remove positive control mock
sub_Genus.melted[sub_Genus.melted$Abundance<0.01,32:38] <- "Z_Others"
# 
# Genus_mock_merged <- merge_samples(Genus, Genus@sam_data$model)
# #Genus_mock_merged@sam_data$model <- sample_names(Genus_mock_merged)
# Genus_mock_merged <- transform_sample_counts(Genus_mock_merged, function(x) x/sum(x))
# Genus_mock_merged.melted <- psmelt(Genus_mock_merged)
# sub_Genus_mock_merged.melted <- Genus_mock_merged.melted
# sub_Genus.melted <- sub_Genus.melted[sub_Genus.melted$SFR == "Mock",] #Remove positive control mock
# sub_Genus_mock_merged.melted[sub_Genus_mock_merged.melted$Abundance<0.01,32:38] <- "Z_Others"

#pdf("Positive_Controls_relative_composition.pdf", he = 6, wi = 18)
ggplot(sub_Genus_mock_merged.melted, aes(x = model, y = Abundance, fill = Genus)) +
  theme_bw() +
  #facet_wrap(~ model) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = couleurs_genus) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Genus-level composition",
       #subtitle = "",
       caption = date()) 
#dev.off()

Mock_phy <- subset_samples(Genus, Genus@sam_data$SFR == "Mock")
View(Mock_phy)
TAXOTU <- data.frame(tax_table(Mock_phy), otu_table(Mock_phy))

Mock_phy.norm <- transform_sample_counts(Mock_phy, function(x) (x/sum(x))*100)
TAXOTU.norm <- data.frame(tax_table(Mock_phy.norm), otu_table(Mock_phy.norm))

Mock_phy <- subset_samples(Finest, Finest@sam_data$SFR == "Mock")
View(Mock_phy)
TAXOTU <- data.frame(tax_table(Mock_phy), otu_table(Mock_phy))
