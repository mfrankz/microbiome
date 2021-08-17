###downstream 16s analyses using phyloseq and ggplot2: BETA DIVERSITY
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

#import ps.rds
ps <- readRDS("path/ps.rds")
  #change "path" to the name of your folder containing ps.rds

#create basic plot for Bray-Curtis dissimilarity index 
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
p<-plot_ordination(ps.prop, ord.bray, color="When", title="Beta Diversity (Bray-Curtis)")

#create publication-quality plot for Bray-Curtis dissimilarity index

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

#create beta diversity plot
plot_ordination(ps.prop, ord.bray, color="When", shape="When")+
  geom_point(aes(fill=When),color="black",size=9, alpha=0.9)+
  scale_fill_manual(values = c("#2DA05A", "#234664"))+
  scale_color_manual(values = c("#2DA05A", "#234664"))+
  scale_shape_manual(values=c(21, 24))+
  ggtitle("Beta Diversity (Bray-Curtis Index)")+
  stat_ellipse(type = "norm", linetype = 2, size=1.5)+ #ellipses represent a 95% confidence interval 
  my_theme

#save plot
ggsave("beta_plot.png", width = 16, height = 10, units = "cm")


#convert beta diversity scores to data frame format for easier manipulation
library(reshape2)
beta_df = phyloseq::distance(ps.prop, "bray")
beta_df = melt(as.matrix(beta_df))
beta_df <- reshape(beta_df, 
                   timevar = "Var2",
                   idvar = c("Var1"),
                   direction = "wide")

#create dataframe containing metadata for easier access and merge with beta diversity values
samdf<-data.frame(sample_data(ps))
samdf$SampleID <- rownames(samdf)
names(beta_df)[names(beta_df) == "Var1"] <- "SampleID"
beta_df<-merge(beta_df, samdf, by=c("SampleID"),all=T)

#divide into 2 dataframes (independent vs. dependent variables) as needed for PERMANOVA
IVs <- subset(beta_df, select = c("SampleID", "Subject", "Gender", "Day", "When")) 
DVs <- subset(beta_df, select = -c(SampleID, Subject, Gender, Day, When)) 

#conduct PERMANOVA
library(vegan)
permanova <- adonis(DVs ~ When, data=IVs, permutations=999)
permanova
