# Downstream plotting and analysis of 16s microbiome data in R using phyloseq and ggplot
### This tutorial will allow you to create publication-level graphs and convert phyloseq objects to dataframes for easier manipulation and analysis. Below you will find code for 3 variables: alpha diversity, beta diversity, and taxa frequency

To begin, we will load packages and import a phyloseq object. ***These steps are necessary for all variables in this tutorial.***

1. Load packages. Note: you will also need an up-to-date version of R studio and the phyloseq package (installed through BioConductor): 
```
#BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13") 

#phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#load additional packages
library(phyloseq)
library(ggplot2)
library(dplyr)
```
2. Import a phyloseq object. We will use the phyloseq object from the [dada2 tutorial for R users](https://benjjneb.github.io/dada2/tutorial.html). You can create the object yourself using the dada2 tutorial or import mine [here](https://github.com/mfrankz/microbiome/blob/main/ps.rds). 
This object contains 16s data collected at both early and late timepoints (represented by the variable "When"). 

```
#read in phyloseq object
ps <- readRDS("path/ps.rds") #change path to your directory containing ps.rds
```


# Alpha Diversity 
### Alpha diversity is a gross measurement of species abundance/richness within a sample. There are several different indices you can use to quantify alpha diversity. Here, we will use the Shannon and Simpson metrics.
An R syntax file containing the alpha diversity code can be found [here](https://github.com/mfrankz/microbiome/blob/main/phyloseq_alpha.R).


1. Create a basic alpha diversity plot. This is the type of plot you will find in basic phyloseq tutorials.
```
plot_richness(ps, x="When", measures=c("Shannon", "Simpson"), color="When")
```

2. Create a publication-quality alpha diversity plot. If you would like to change any features (e.g., colors, axes), use ggplot2 syntax to edit. If you have additional grouping variables, consider using +facet_wrap(~VAR_NAME).
```
#set theme
my_theme<-theme(
  plot.title = element_text(size=40, face="bold"),
  axis.title.x = element_text(size=30, face="bold"),
  axis.title.y = element_text(size=30, face="bold"),
  axis.text.y = element_text(size=20, face="bold", color="black"),
  axis.text.x = element_text(size=20, angle=0, hjust = 0.4, face="bold", color="black"),
  legend.title = element_text(size = 30, face="bold"),
  legend.text = element_text(size = 25, face="bold"),
  strip.text.x = element_text(size = 25, face="bold"), 
  strip.background = element_rect(color="white", fill="white"),
  panel.background = element_rect(fill="white", colour="white"),
  panel.border = element_rect(colour = "black", fill=NA, size=2),
  legend.key=element_blank(),
  legend.key.height = unit(1.5, "cm"),
  legend.key.width = unit(1, "cm")
) 

#create plot
plot_richness(ps, x="When", measures=c("Shannon", "Simpson"), 
                 color="When", shape="When")+   
  geom_point(aes(fill=When),size=9, alpha=0.9, color="black")+
  scale_fill_manual(values = c("#2DA05A", "#234664"))+
  scale_color_manual(values = c("#2DA05A", "#234664"))+
  scale_shape_manual(values=c(21, 24))+
  ylab("Diversity Score")+
  ggtitle("Alpha Diversity Across Time")+ 
  xlab("Collection Timepoint")+
  my_theme+
  theme(legend.position = "none")
```

![alpha_plot](https://user-images.githubusercontent.com/88938223/129733837-651a2a56-aa66-4244-8358-a061cb059e6e.png)




3. Convert alpha diversity data into a dataframe for easier manipulation and analyses. The dataframe alpha_df contains a row for each sample with the metadata (SampleID, Subject, Gender, Day, When) and Shannon alpha diversity score.
```
alpha_df <- estimate_richness(ps, split = TRUE, measure = "Shannon")
alpha_df$SampleID <- rownames(alpha_df) %>%
  as.factor()
alpha_df <- merge(alpha_df, sample_data(ps), by=0)
View(alpha_df)

#the dataframe format can now be used for easier manipulation/analyses. See examples below: 
#example 1: check distribution of alpha scores
hist(alpha_df$Shannon)

#example 2: linear regression
alpha_df$When<-as.factor(alpha_df$When)
model<- lm(Shannon~When, data=alpha_df) #determine effects of collection timepoint
summary(model)
```

# Beta Diversity 
### Beta diversity is a gross measurement of dissimilarity between samples. There are several methods that can be used to conduct dimensionality reduction to calculate beta diversity. Here, we will use the NMDS method to calculate Bray-Curtis dissimilarity. Other methods (e.g., PCoA), can be found [here](https://joey711.github.io/phyloseq/plot_ordination-examples.html).
An R syntax file containing the beta diversity code can be found here.

Note: you must import the phyloseq/ggplot packages as well as a phyloseq object. See instructions above if you have not already completed this step. 

1. Create basic Bray-Curtis plot
```
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
p<-plot_ordination(ps.prop, ord.bray, color="When", title="Beta Diversity (Bray-Curtis)")
```
2. Create publication-quality Bray-Curtis plot
```
#set theme
my_theme<-theme(
  plot.title = element_text(size=32, face="bold"),
  axis.title.x = element_text(size=28, face="bold"),
  axis.title.y = element_text(size=28, face="bold"),
  axis.text.y = element_text(size=20, face="bold", color="black"),
  axis.text.x = element_text(size=20, angle=0, hjust = 0.4, face="bold", color="black"),
  legend.title = element_text(size = 28, face="bold"),
  legend.text = element_text(size = 25, face="bold"),
  strip.text.x = element_text(size = 25, face="bold"), 
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
  stat_ellipse(type = "norm", linetype = 2, size=1.5)+
  my_theme
```

