library(ggplot2)
library(phyloseq)

setwd("~/Documents/Thesis/NanoASV/Software_dev/NanoASV/Mock_run_OUTPUT/")

load("Results/Rdata/NanoASV.rdata")

Finest@sam_data <- sample_data(metadata)

Genus  <- tax_glom(NanoASV, taxrank = "Genus")
Family  <- tax_glom(NanoASV, taxrank = "Family")
Order  <- tax_glom(NanoASV, taxrank = "Order")
Class  <- tax_glom(NanoASV, taxrank = "Class")
Phylum  <- tax_glom(NanoASV, taxrank = "Phylum")






# Taxonomical composition ----

couleurs_genus <- c("#FF0000FF","#FF9900FF","#FFCC00FF","#00FF00FF","#6699FFFF",
                    "#CC33FFFF","#99991EFF","#FF00CCFF", "#999999FF","#CC0000FF",
                    "#FFCCCCFF","#FFFF00FF","#CCFF00FF","#358000FF","#0000CCFF",
                    "#99CCFFFF","#00FFFFFF","#CCFFFFFF","#9900CCFF","#CC99FFFF",
                    "#996600FF","#79CC3DFF","#CCCCCCFF","#79CC3DFF","#CCCC99FF","black")

#Composition ----

#pdf("Composition_with_different_models.pdf", he = 8, wi = 16)

Genus.norm <- transform_sample_counts(Genus, function(x) x/sum(x))
Genus.melted <- psmelt(Genus.norm)
sub_Genus.melted <- Genus.melted
sub_Genus.melted[sub_Genus.melted$Abundance<0.008,15:10] <- "Z_Others"

ggplot(sub_Genus.melted, aes(x = Sample, y = Abundance, fill = Genus)) +
  theme_bw() +
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
sub_Family.melted[sub_Family.melted$Abundance<0.00875,15:10] <- "Z_Others"

ggplot(sub_Family.melted, aes(x = Sample, y = Abundance, fill = Family)) +
  theme_bw() +
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
sub_Order.melted[sub_Order.melted$Abundance<0.0091,15:10] <- "Z_Others"

ggplot(sub_Order.melted, aes(x = Sample, y = Abundance, fill = Order)) +
  theme_bw() +
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
sub_Class.melted[sub_Class.melted$Abundance<0.00555,15:10] <- "Z_Others"

ggplot(sub_Class.melted, aes(x = Sample, y = Abundance, fill = Class)) +
  theme_bw() +
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
sub_Phylum.melted[sub_Phylum.melted$Abundance<0.00555,15:10] <- "Z_Others"


situseq_colors <- c("#2f4858ff", "#33658aff", "#86bbd8ff", "#830689ff", "#f5a614ff", "#f26419ff", "#bb3551ff", "#c1d7aeff" ,"#68ac5dff", "#7f7f7fff")
pdf("../../SituSeq/MOCK/NanoASV_phylum_mock.pdf", he = 10, wi = 9)
ggplot(sub_Phylum.melted, aes(x = Sample, y = Abundance, fill = Phylum)) +
  theme_bw() +
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = situseq_colors) + 
  ylab("Relative Abundance") +
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "right") +  #remove minor-grid labels
  labs(title = "Soil samples Phylum-level composition",
       #subtitle = "",
       caption = date()) 
dev.off()
