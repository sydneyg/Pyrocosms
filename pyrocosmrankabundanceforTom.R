
#July 24, 2019
#make rank abundance figure for tom

#Reset R's Brain
rm(list=ls())

##set working directory
setwd("~/Dropbox/StatsandProgramming/Pyrocosms/Tomstuff/")

#load libraries
library(tidyverse)

##call file to variable "sanddata"
rankdata<-read.csv("rankabundance.csv",header=TRUE)

#it looks like Tom only plotted top 50
rankdatatop50 <- rankdata[1:50, ]

names(rankdata)

#how to add points from different columns to same plot in ggplot2: https://www.sixhat.net/how-to-plot-multpile-data-series-with-ggplot.html
rankabundancefigure <- ggplot(rankdatatop50, aes(x=rank, y= value, color=variable)) +
            geom_point(aes(y=Control, col="Control")) + 
            geom_point(aes(y=Pyrocosm1, col="Pyrocosm1"))+
            geom_point(aes(y=Pyrocosm2, col="Pyrocosm2")) +
            theme_bw() +
             scale_y_continuous(trans='log10') + #make scale log
            ylab("Percent Sequence Abundance (log scale)") + xlab ("Rank of OTUs") +
            scale_color_manual(labels=c("Control","Pyrocosm 1","Pyrocosm 2"), #manual labels for legend
                             values=c("burlywood4", "darkorange1", "firebrick1")) +  #add in manual colors for points/lines
    theme(legend.position =  c(0.7, 0.9), #put legend inside graph upper right corner
                  legend.title=element_blank(), #remove legend title
                  text = element_text(size=18),
                  axis.text.y = element_text(size=18, angle=90,hjust=1),
                  axis.text.x = element_text(size=18, angle=45, hjust=1))   #make x axis text larger and angled

rankabundancefigure


pdf("pyrocosm_rankabundance_fortom_logscale.pdf", height=6, width=6 )
rankabundancefigure 
dev.off()

#try it without log scale but add jitter to see points easier
rankabundancefigurewithjitter <- ggplot(rankdatatop50, aes(x=rank, y= value, color=variable)) +
  geom_point(aes(y=Control, col="Control")) + 
  geom_point(aes(y=Pyrocosm1, col="Pyrocosm1"))+
  geom_point(aes(y=Pyrocosm2, col="Pyrocosm2")) +
  theme_bw() +
  geom_jitter(aes(y=Pyrocosm2, col="Pyrocosm2"))+ #add jitter to points so you can see them better
  geom_jitter(aes(y=Pyrocosm1, col="Pyrocosm1"))+ #add jitter to points so you can see them better
  ylim(0.1,100)+
  ylab("Percent Sequence Abundance") + xlab ("Rank of OTUs") +
  scale_color_manual(labels=c("Control","Pyrocosm 1","Pyrocosm 2"), #manual labels for legend
                     values=c("burlywood4", "darkorange1", "firebrick1")) +  #add in manual colors for points/lines
  theme(legend.position = "bottom", #put legend under graph
        legend.title=element_blank(), #remove legend title
        text = element_text(size=18),
        axis.text.y = element_text(size=18, angle=90,hjust=1),
        axis.text.x = element_text(size=18, angle=45, hjust=1))   #make x axis text larger and angled

rankabundancefigurewithjitter

  
  pdf("pyrocosm_rankabundance_fortom_jitter.pdf", height=6, width=6 )
  rankabundancefigurewithjitter
  dev.off()
  
  
  #plot without pyronema
  
  names(rankdatatop50)
  
  rankabundancefigure_nopyronema <- ggplot(rankdatatop50, aes(x=rank, y= value, color=variable)) +
    geom_point(aes(y=Control, col="Control")) + 
    geom_point(aes(y=Pyrocosm.1.no.Pyronema, col="Pyrocosm.1.no.Pyronema"))+
    geom_point(aes(y=Pyrocosm.2.no.Pyronema, col="Pyrocosm.2.no.Pyronema")) +
    theme_bw() +
    scale_y_continuous(trans='log10') + #make scale log
    ylab("Percent Sequence Abundance (log scale)") + xlab ("Rank of OTUs") +
    scale_color_manual(labels=c("Control","Pyrocosm 1 no Pyronema","Pyrocosm 2 no Pyronema"), #manual labels for legend
                       values=c("burlywood4", "darkorange1", "firebrick1")) +  #add in manual colors for points/lines
    theme(legend.position =  c(0.7, 0.9), #put legend inside graph upper right corner
          legend.title=element_blank(), #remove legend title
          text = element_text(size=18),
          axis.text.y = element_text(size=18, angle=90,hjust=1),
          axis.text.x = element_text(size=18, angle=45, hjust=1))   #make x axis text larger and angled
  
  rankabundancefigure_nopyronema
  
  pdf("pyrocosm_rankabundance_fortom_nopyronema.pdf", height=6, width=6 )
  rankabundancefigure_nopyronema
  dev.off()
