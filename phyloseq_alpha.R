###downstream 16s analyses using phyloseq and ggplot2: ALPHA DIVERSITY
##create publication-level plots 
##convert data from class phyloseq to class dataframe for manipulation/analysis 

#load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

#import ps.rds
ps <- readRDS("path/ps.rds")
  #change "path" to the name of your folder containing ps.rds

#basic alpha diversity plot 
plot_richness(ps, x="When", measures=c("Shannon", "Simpson"), color="When")


#advanced alpha diversity plot
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

#save plot
ggsave("alpha.png", width = 30, height = 20, units = "cm")



#convert alpha diversity scores to data frame format for easier manipulation
alpha_df <- estimate_richness(ps, split = TRUE, measure = "Shannon")
alpha_df$SampleID <- rownames(alpha_df) %>%
  as.factor()
alpha_df <- merge(alpha_df, sample_data(ps), by=0)

#data frame format can be used to more easily run analyses and check model assumptions
#see examples below:

#example 1: check distribution of alpha scores
hist(alpha_df$Shannon)

#example 2: linear regression
alpha_df$When<-as.factor(alpha_df$When)
model<- lm(Shannon~When, data=alpha_df) #determine effects of collection timepoint
summary(model)

