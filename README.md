# Downstream plotting and analysis of 16s microbiome data in R using phyloseq and ggplot
### This beginner-friendly tutorial will allow you to create publication-level graphs and convert phyloseq objects into dataframes for easier manipulation and analysis. Below you will find R code for extracting alpha diversity, beta diversity, and taxa abundance. 

To begin, we will load packages and import a phyloseq object. ***This is necessary for all the following steps in this tutorial.***

1. Load packages. Note: you will also need an up-to-date version of R studio and the phyloseq package (installed through BioConductor): 
```
#install BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13") 

#install phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#load additional packages
library(phyloseq)
library(ggplot2)
library(dplyr)
```
2. Import a phyloseq object. We will use the phyloseq object from the [dada2 tutorial for R users](https://benjjneb.github.io/dada2/tutorial.html). You can create the object yourself using the dada2 tutorial or import mine [here](https://github.com/mfrankz/microbiome/blob/main/ps.rds). 
This object contains 16s data collected at both early and late timepoints (represented by the variable "When", which we will rename as "Timepoint"). 

```
#read in phyloseq object and rename "When" variable
ps <- readRDS("path/ps.rds") #change "path" to the name of your folder containing ps.rds
sample_data(ps)$Timepoint<-sample_data(ps)$When
```


# 1. Alpha Diversity 
### Alpha diversity is a gross measurement of species abundance/richness within a sample. There are several indices you can use to quantify alpha diversity. Here, we will use the Shannon and Simpson metrics.
Note: click [here](https://github.com/mfrankz/microbiome/blob/main/phyloseq_alpha.R) if you would prefer to view the following alpha diversity syntax in an R script file.


1. Create a basic alpha diversity plot. This is the type of plot you will find in phyloseq tutorials.
```
plot_richness(ps, x="Timepoint", measures=c("Shannon", "Simpson"), color="Timepoint")
```

2. Create a publication-quality alpha diversity plot. If you would like to change any features (e.g., colors, axes), use ggplot2 syntax to edit. If you have additional grouping variables, consider using +facet_wrap(~VAR_NAME).
```
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

#create plot
plot_richness(ps, x="Timepoint", measures=c("Shannon", "Simpson"), 
                 color="Timepoint", shape="Timepoint")+   
  geom_point(aes(fill=Timepoint),size=6, alpha=0.9, color="black")+
  scale_fill_manual(values = c("#2DA05A", "#234664"))+
  scale_color_manual(values = c("#2DA05A", "#234664"))+
  scale_shape_manual(values=c(21, 24))+
  ylab("Diversity Score")+
  ggtitle("Alpha Diversity Across Time")+ 
  xlab("Collection Timepoint")+
  my_theme
```


<img src="https://user-images.githubusercontent.com/88938223/129901354-1b5820e7-1d83-4896-a310-b962b9abdd8f.png" width="500">



3. Convert alpha diversity data into a dataframe for easier manipulation and analyses. The dataframe alpha_df contains a row for each sample with the metadata (SampleID, Subject, Gender, Day, When) and Shannon alpha diversity score.
```
alpha_df <- estimate_richness(ps, split = TRUE, measure = "Shannon")
alpha_df$SampleID <- rownames(alpha_df) %>%
  as.factor()
alpha_df <- merge(alpha_df, sample_data(ps), by=0)
names(alpha_df)[names(alpha_df) == "Row.names"] <- "SampleID"
View(alpha_df)
```
<img src="https://user-images.githubusercontent.com/88938223/130094308-848888d1-ed08-4b1b-91fa-7ad9f253b1a8.png" width="500">



4. The dataframe format can now be used for easier manipulation/analyses. See examples below: 
```
#example 1: check distribution of alpha scores
hist(alpha_df$Shannon)

#example 2: linear regression
alpha_df$Timepoint<-as.factor(alpha_df$Timepoint)
model<- lm(Shannon~Timepoint, data=alpha_df) #determine effects of collection timepoint
summary(model)
```

# 2. Beta Diversity 
### Beta diversity is a gross measurement of dissimilarity between samples. There are several methods that can be used to conduct dimensionality reduction to calculate beta diversity. Here, we will use the NMDS method to calculate Bray-Curtis dissimilarity. Other methods (e.g., PCoA), can be found [here](https://joey711.github.io/phyloseq/plot_ordination-examples.html).
Note: click [here](https://github.com/mfrankz/microbiome/blob/main/phyloseq_beta.R) if you would prefer to view the following beta diversity syntax in an R script file.

Note: you must import the phyloseq/ggplot packages and the phyloseq object ps. See instructions above if you have not already completed this step. 

1. Create basic Bray-Curtis plot
```
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.bray, color="Timepoint", title="Beta Diversity (Bray-Curtis)")
```
2. Create publication-quality Bray-Curtis plot
```
plot_ordination(ps.prop, ord.bray, color="Timepoint", shape="Timepoint")+
  geom_point(aes(fill=Timepoint),color="black",size=9, alpha=0.9)+
  scale_fill_manual(values = c("#2DA05A", "#234664"))+
  scale_color_manual(values = c("#2DA05A", "#234664"))+
  scale_shape_manual(values=c(21, 24))+
  ggtitle("Beta Diversity (Bray-Curtis Index)")+
  stat_ellipse(type = "norm", linetype = 2, size=1.5)+ #ellipses represent a 95% confidence interval 
  my_theme
```

<img src="https://user-images.githubusercontent.com/88938223/129900230-ac865a0f-fad8-45e9-b40d-fb70cc8c806a.png" width="600">

3. Convert beta diversity scores to data frame format for easier manipulation and analysis. 
```
#create data frame with beta diversity scores
library(reshape2)
beta_df = phyloseq::distance(ps.prop, "bray")
beta_df = melt(as.matrix(beta_df))
beta_df <- reshape(beta_df, 
                   timevar = "Var2",
                   idvar = c("Var1"),
                   direction = "wide")
View(beta_df)
                   
#create dataframe containing metadata for easier access and merge with beta diversity values
samdf<-data.frame(sample_data(ps))
samdf <- subset(samdf, select = -c(When)) 
samdf$SampleID <- rownames(samdf)
names(beta_df)[names(beta_df) == "Var1"] <- "SampleID"
beta_df<-merge(beta_df, samdf, by=c("SampleID"),all=T)
```
4. Analyze beta diversity scores using PERMANOVA. This type of analysis is used to determine the effects of your independent variables on beta diversity. To conduct PERMANOVA, we will need to create 2 separate dataframes: one that contains all possible independent variables (metadata) and one that contains the dependent variables (beta diversity scores).
```
#divide into 2 dataframes (independent vs. dependent variables) as needed for PERMANOVA
IVs <- subset(beta_df, select = c("SampleID", "Subject", "Gender", "Day", "Timepoint")) 
DVs <- subset(beta_df, select = -c(SampleID, Subject, Gender, Day, Timepoint)) 

#conduct PERMANOVA
library(vegan)
permanova <- adonis(DVs ~ Timepoint, data=IVs, permutations=999)
permanova
```
In the permanova output, we can see that collection timepoint has a significant effect on beta diversity.


<img src="https://user-images.githubusercontent.com/88938223/130094075-803254db-948d-4577-9a11-e182fcc43cce.png" width="600">
          
# 3. Taxa Abundance
### The next step is to determine how changes in specific levels of taxonomy may be driving these broader changes in alpha and beta diversity. Below, we will quantify taxa abundance at the level of the phylum, but please note that you can easily change this to other levels of taxonomy (e.g., Species, Family).
Note: click [here](https://github.com/mfrankz/microbiome/blob/main/phyloseq_abundance.R) if you would prefer to view the following taxa abundance syntax in an R script file.

Note: you must import the phyloseq/ggplot packages and the phyloseq object ps. See instructions above if you have not already completed this step. 

1. Prepare phylum data by (1) converting counts to relative abundance, (2) isolating phylum data, and (3) converting the phylum data to a dataframe. 
```
#transform abundance to relative counts
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))

#isolate and filter phylum data
phy <- tax_glom(ps_rel, taxrank="Phylum") 
phy <- prune_taxa(taxa_sums(phy)>=0.01, phy)

#convert phylum data to dataframe format
phy_df <- psmelt(phy)
View(phy_df)
```

<img src="https://user-images.githubusercontent.com/88938223/129935058-0beb9c71-ffb8-49d4-a572-29f9431df772.png" width="600">

2. Create publication-quality plot of individual samples. Note: if you have individual subjects data, you can add +facet_wrap(~Subject). 
```
ggplot(phy_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))+  	
  ggtitle("Phyla Abundance Across Individual Samples") +
  ylab("Relative Abundance")+
  scale_fill_brewer(palette="Paired")+
  my_theme+
  theme(axis.text.x = element_text(size=14,angle=90))
```

<img src="https://user-images.githubusercontent.com/88938223/129935966-e6a24f25-738c-4938-80b7-95d5c17fcd02.png" width="600">

3. Create publication-quality plot of aggregate data. We will manually calculate the mean of each phylum per timepoint before creating the plot. 
```
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
  theme(axis.text.x = element_text(size=14,angle=-45, hjust=0))
```

<img src="https://user-images.githubusercontent.com/88938223/129966269-eb3f4ae0-6630-4317-985f-73a337efefc6.png" width="600">
