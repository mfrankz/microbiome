###downstream 16s analyses using phyloseq and ggplot2: ASV FREQUENCY
##create publication-level plots 
##convert data from class phyloseq to class dataframe for manipulation/analysis 

#load libraries
library(phyloseq)
library(ggplot2)

#if you have not yet installed BioConductor or phyloseq: 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13") 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#import ps.rds and rename variables
ps <- readRDS("path/ps.rds") #change "path" to the name of your folder containing ps.rds
sample_data(ps)$Timepoint<-sample_data(ps)$When

#transform abundance to relative counts
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))

#isolate and filter phylum data
phy <- tax_glom(ps_rel, taxrank="Phylum") 
phy <- prune_taxa(taxa_sums(phy)>=0.01, phy)

#convert phylum data to dataframe format
phy_df <- psmelt(phy)
View(phy_df)



#advanced plotting 
#set theme
my_theme<-theme(
  plot.title = element_text(size=22, face="bold"),
  axis.title.x = element_text(size=20, face="bold"),
  axis.title.y = element_text(size=20, face="bold"),
  axis.text.y = element_text(size=16, face="bold", color="black"),
  axis.text.x = element_text(size=16, angle=0, hjust = 0.4, face="bold", color="black"),
  legend.title = element_text(size = 20, face="bold"),
  legend.text = element_text(size = 16, face="bold"),
  strip.text.x = element_text(size = 16, face="bold"), 
  strip.background = element_rect(color="white", fill="white"),
  panel.background = element_rect(fill="white", colour="white"),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.key=element_blank(),
  legend.key.height = unit(1.5, "cm"),
  legend.key.width = unit(1, "cm")
) 

#plot individual samples
ggplot(phy_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+  	
  ggtitle("Phyla Abundance Across Individual Samples") +
  ylab("Relative Abundance")+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x = element_text(size=14,angle=90))

#save plot
ggsave("individual_phyla_plot.png", width = 25, height = 18, units = "cm")

#plot aggregate data

#calculate average phyla abundance across timepoints
library(tidyr)
phy_avg <- aggregate(phy_df[,"Abundance"],
                     by = list(phy_df[,"Timepoint"], 
                               phy_df[,"Phylum"]),
                     FUN = function(x) mean = mean(x))
colnames(phy_avg) <- c("Timepoint", "Phylum", "Abundance")


#plot aggregate abundance across timepoint
ggplot(phy_avg, aes(x = Phylum, y = Abundance, fill=Phylum))+
  facet_wrap(~Timepoint)+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+  	
  ggtitle("Phyla Abundance Across Time")+
  ylab("Relative Abundance")+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x = element_text(size=14,angle=-45, vjust=-0.5))

#save plot
ggsave("avg_phyla_plot.png", width = 25, height = 18, units = "cm")


