library(ggplot2)
library(phyloseq)
load("../../map-ont-mock.out/Results/Rdata/NanoASV.rdata")
mapont <- NanoASV
mapont@sam_data <- sample_data(data.frame(mapont@sam_data,
                                          model = "map-ont"))
load("../../asm10-mock.out/Results/Rdata/NanoASV.rdata")
asm10 <- NanoASV
asm10@sam_data <- sample_data(data.frame(asm10@sam_data,
                                          model = "asm10"))
load("../../asm5-mock.out/Results/Rdata/NanoASV.rdata")
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

Genus  <- tax_glom(Finest, taxrank = "Genus")
Family  <- tax_glom(Finest, taxrank = "Family")
Order  <- tax_glom(Finest, taxrank = "Order")
Class  <- tax_glom(Finest, taxrank = "Class")
Phylum  <- tax_glom(Finest, taxrank = "Phylum")

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
  metadata$Barcode <- c("Barcode_01", "Barcode_02")
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


ggplot(data = alpha.melted[alpha.melted$variable == "richness",], aes(x =  Barcode , y = value), color = alpha_tab$model) +
  facet_grid(rows = vars(level), cols = vars(model), scales = "free_y", switch = "y") +  # Free y-axis for each row
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        legend.position = "none") + 
  ylab("Numerical richness") + xlab("Model") + 
  labs(title = "Numerical richness",
       caption = date())

ggplot(data = alpha.melted[alpha.melted$variable == "shannon",], aes(x =  Barcode , y = value), color = alpha_tab$model) +
  facet_grid(rows = vars(level), cols = vars(model), scales = "free_y", switch = "y") +  # Free y-axis for each row
  geom_bar(stat = "identity") + theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=14),
        legend.position = "none") + 
  ylab("Shannon index") + xlab("Model") + 
  labs(title = "SHannon index",
       caption = date())

TAX <- data.frame(asm5@tax_table)
TAX <- data.frame(asm10@tax_table)
TAX <- data.frame(mapont@tax_table)
