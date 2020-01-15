#may 21, 2019

#make NMDS with different shapes for Time points and depths and pyrocosms

#Reset R's Brain
rm(list=ls())

#Packages to load 
library(tidyverse)
library(SPECIES)
library(vegan)
library(BiodiversityR)
library(scales)
library(ggpubr)
#load library for nonlinear mixed effect models
library(nlme)
library(multcomp)
library("stringr")
library("plyr")

#calculate beta diversity and plot ordinations for all fungi
setwd("~/Dropbox/StatsandProgramming/Pyrocosms")
#upload total fungi table

#source in functions
source('~/Dropbox/StatsandProgramming/source/pvaluelegend.R', chdir = TRUE)

#read in dataframes
otu_e3950 <- read.csv("otu_table_pyro_e3950.csv", row.names=1)
ncol(otu_e3950)

##upload map
map_e3950 <-read.csv("mappingfile_e3950.csv", row.names=1)
names(map_e3950)

## We will create a new data frame for the last two columns as it provides different set of information. We'll call it "taxon" data.frame
taxon_pyro <- (otu_e3950[ ,25:26])
## make community matrix; make it a data frame, transpose it, remove taxonomy
otu_trans <- data.frame(t(otu_e3950[ ,1:24]))


#now calculate richness with BioDiversityR package
#run function once to see what it does 

#run function on transformed table to switch rows and columns
OTU_richness<- t(estimateR(otu_trans))

#make a function to get estimates and make plots
#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("Figures/richnesscores_fungalrichness_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
  #perform correlation test
  cortest1 <- cor.test(estimates2[,2],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest1$estimate, cortest1$p.value)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  #perform correlation test
  cortest2 <- cor.test(estimates2[,4],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest2$estimate, cortest2$p.value)
  mtext("B",side=3,adj=0)
  dev.off()
  
}


#run function
estimates_plot_function(otu_trans,"otu_e3950")




##############################
#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
otu.H <- diversity(otu_trans) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(otu_trans, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
OTU_ITS_e3950_richness <- cbind(otu.richness,OTU_richness)

write.csv(OTU_ITS_e3950_richness , "OTU_ITS_e3950_richness.csv")

#test if data are normal distirbuted, if significant not normal
shapiro.test(OTU_ITS_e3950_richness$S.obs)
histogram(OTU_ITS_e3950_richness$S.obs)

########################################################################################
########Make some figures
##################################################################

#check that names match, if they dont match, use a matching function
row.names(OTU_ITS_e3950_richness) == map_e3950$Samplename2


fungaldata <- cbind(OTU_ITS_e3950_richness,map_e3950)

#rename and reorder time points, can rename and reorder replicates if I know what they are
fungaldata$Time <- factor(fungaldata$Time, levels=c("T1","T2","T3"), labels=c("1 week", "2 weeks", "4 weeks"))
fungaldata$Pyrocosm <- factor(fungaldata$Pyrocosm, levels=c("Pyrocosm1Neg","Pyrocosm1","Pyrocosm2"), labels=c("Control", "Pyrocosm 1", "Pyrocosm 2"))


fungibyplot <- ggplot(fungaldata, aes(x=Pyrocosm, y=S.obs)) +   
  geom_boxplot(aes(fill=Pyrocosm, alpha=Time))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal species (OTUs)") +  #change yaxis label
  theme(#legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=1),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  # stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
  #             label = c("a","a","a","b","a","b"), vjust = -0.5, col="black") +
  scale_fill_manual(labels=c("Control","Pyrocosm 1","Pyrocosm 2"), #manual labels for legend
                    values=c("burlywood4", "darkorange1", "firebrick1"))  #add in manual colors for points/lines
 #geom_jitter(shape=16, position=position_jitter(0.1))
#figure out how to get tukey letters above
#stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
#           label = c("a","a","b"), vjust = -0.5, col="black")

#now you have to figure out how to visualize the fact that Plot58 is not burned etc
fungibyplot


fungibyplot%>%
  ggexport(filename = "Figures/fungalrichness_pyrocosms_boxplot.pdf")

#################ggplot plot with mean and SE instead of boxplot ##########


A <-ggplot(fungaldata, aes(x=Pyrocosm, y=S.obs, group=Time, col=Pyrocosm,shape=Time))+
  stat_summary(fun.y=mean,geom="point", size=4)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=0.5,
               alpha=0.7,position = position_dodge(0.01))+
  scale_color_manual(labels=c("Control","Pyrocosm1","Pyrocosm2"), #manual labels for legend
                    values=c("burlywood4", "darkorange1", "firebrick1"))  +#add in manual colors for points/lines
  theme_bw() + #make black and white
  ylab("Observed number of fungal species (OTUs)") +  #change yaxis label
  #geom_jitter(shape=16, position=position_jitter(0.1))+
  theme(#legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=1),
        axis.text.x = element_text(size=18, angle=45, hjust=1))   #make x axis text larger and angled
A

A %>%
  ggexport(filename = "Figures/fungalrichness_pyrocosms_meanandSE.pdf")

#what's the % reduction from pre to post fire in the 2 plots?
#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

fungal_plots1 <- ddply(fungaldata, c("Pyrocosm","Time","Replicates"), summarise,
                   T1_mean = mean(S.obs, na.rm=TRUE),
                   T1_sd = sd(S.obs, na.rm=TRUE),
                   T1_n = sum(!is.na( S.obs)),
                   T1_se = T1_sd/sqrt(T1_n)
)

head(fungal_plots1)

write.csv(fungal_plots1,"fungalrichnessbypyrocosmandtimeandreplicates_meanandSE.csv")

fungal_plots2 <- ddply(fungaldata, c("Pyrocosm","Time"), summarise,
                      T1_mean = mean(S.obs, na.rm=TRUE),
                      T1_sd = sd(S.obs, na.rm=TRUE),
                      T1_n = sum(!is.na( S.obs)),
                      T1_se = T1_sd/sqrt(T1_n)
)

head(fungal_plots2)

write.csv(fungal_plots2,"fungalrichnessbypyrocosmandtime_meanandSE.csv")


########################################################################################################
#---------------------------------Do stats on richness: lme ------------------------------------------
########################################################################################################

#used linear mixed effect model from packaage nlme to run nested anova
#load library for nonlinear mixed effect models
library(nlme)
names(fungaldata)

#model
Fit1 <- (aov(S.obs ~  Pyrocosm*Time, data = fungaldata))
summary(Fit1 )

#if data are not normally distributed, use glmer
#generate generalized linear mixed effects model to deal with poisson distrubition, covariate, and nested subsampling with plot random effects
#install.packages("lme4")
#library(lme4)

#fit1<-glmer(S.obs~Fire*Burn + (1|Plot),family=poisson, data=bacteriadata)
#summary(fit1)

########################################################################################################
#---------------------------------calculate beta diversity value ------------------------------------------
########################################################################################################

#make dissimilarity matrix - bray-curtis, abundance based
otu_e3950.bray <- vegdist(otu_trans,method="bray",binary=FALSE,upper=TRUE)

#make dissimilarity matrix - jaccard, presence-absence
otu_e3950.jaccard <- vegdist(otu_trans,method="jaccard",binary=TRUE,upper=TRUE)


#https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation
plotnmds1 <- metaMDS(otu_e3950.bray, k=2, trymax=100)
plotnmds1 #stress=0.12
stressplot(plotnmds1)

#make dataframe of NMDS scores
scores <- as.data.frame(scores(plotnmds1))

#add scores to metadata
betadiversitydata <- cbind(fungaldata, scores)

#Do adonis to test if there are significant differences of treatments on bray curtis dissimilarity
adonis(otu_e3950.bray ~ Pyrocosm*Time, data = betadiversitydata, permutations=999)
adonis(otu_e3950.jaccard ~ Pyrocosm*Time, data = betadiversitydata, permutations=999)
#sig effects of pyrocosm, time, pyrocosm by time


#######

#separate unburned v burned
betadiversitydata

#betadiversitydata_Burnplots <- betadiversitydata[which(betadiversitydata$Burn=="Burned"), ]
#betadiversitydata_unburnplot <- betadiversitydata[which(betadiversitydata$Burn=="Unburned"), ]



#NMDS for pyrocosms
NMDS1_plot <- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2, col=Pyrocosm, shape=Time, group=Pyrocosm, alpha=Replicates1)) +
  geom_point(size=3) +
  theme_bw() +
  labs(title = "Fungal Bray-Curtis Dissimilarity")+
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(#legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=20), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Control","Pyrocosm1","Pyrocosm2"), #manual labels for legend
                     values=c("burlywood4", "darkorange1", "firebrick1"))  #add in manual colors for points/lines
  


  NMDS1_plot
  
  NMDS1_plot %>%
    ggexport(filename = "Figures/Pyrocosms_NMDS_pyrotimereplicates.pdf")
  
  
  #NMDS for burn plots
  NMDS2_plot <- ggplot(betadiversitydata, aes(x=NMDS1, y=NMDS2, col=Pyrocosm, shape=Time, group=Pyrocosm)) +
    geom_point(size=3) +
    theme_bw() +
    labs(title = "Fungal Bray-Curtis Dissimilarity")+
    #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
    theme(#legend.position = "bottom", #put legend under graph,
          legend.title = element_blank(),
          strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
          axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
          axis.title=element_text(size=18),#change size of x and y axis labels
          legend.spacing = unit(0,"cm"), #change spacing between legends
          legend.text=element_text(size=20), #change size of legend text
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
    stat_ellipse() + #include ellipse around the group
    scale_color_manual(labels=c("Control","Pyrocosm1","Pyrocosm2"), #manual labels for legend
                       values=c("burlywood4", "darkorange1", "firebrick1"))  #add in manual colors for points/lines
  
  
  
  NMDS2_plot 
  
  NMDS2_plot %>%
    ggexport(filename = "Figures/Pyrocosms_NMDS_pyrotime.pdf")
  


#install.packages("ggordiplots")
#library(ggordiplots)
  #https://stackoverflow.com/questions/47516448/how-to-get-ordispider-like-clusters-in-ggplot-with-nmds


pdf("pyrocosms_NMDS_jaccard_e3950.pdf",bg="white")

map_e3950$colors <- as.character(map_e3950$colors)
#how do I account for nestedness of samples within plots? Not sure how to do that here
plotnmds.jac <- metaMDS(as.dist(otu_e3950.jaccard))
#plot the NMDS
plot(scores(plotnmds.jac),col=map_e3950$colors,pch = 19, cex = 2)
#add in legend to color trees
legend("topright", pch = 19, col =c("burlywood4","firebrick1", "darkorange1"), legend = c("control","pyrocosm 1", "pyrocosm 2"), cex=1)

#add in ordiellipse
#ordiellipse(plotnmds, groups=map$Tree, label = T, kind = "sd", draw = "polygon", col ="lightblue")
#add in ordispider
ordispider(scores(plotnmds.jac), groups = map_e3950$Pyrocosm, label = T, col = "darkgrey")
#conduct adonis on distance matrix by Tree
dev.off()
adonis_pyro1 <- adonis(otu_e3950.jaccard~ Pyrocosm*Time, data=map_e3950, permutation=9999 )
adonis_pyro1





pdf("pyrocosms_NMDS_bray_e3950.pdf",bg="white")

#how do I account for nestedness of samples within plots? Not sure how to do that here
plotnmds.bray <- metaMDS(as.dist(otu_e3950.bray))
#plot the NMDS
plot(scores(plotnmds.bray ),col=map_e3950$colors,pch = 19, cex = 2)
#add in legend to color trees
legend("topright", pch = 19, col =c("burlywood4","firebrick1", "darkorange1"), legend = c("control","pyrocosm 1", "pyrocosm 2"), cex=1)
#add in ordispider
ordispider(scores(plotnmds.bray), groups = map_e3950$Pyrocosm, label = T, col = "darkgrey")
#conduct adonis on distance matrix by Tree
dev.off()
adonis_pyro2 <- adonis(otu_e3950.bray~ Pyrocosm*Time, data=map_e3950, permutation=9999 )
adonis_pyro2


