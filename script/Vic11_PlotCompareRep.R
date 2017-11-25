#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)

ReadFitnessTable <- function(filename,inputcutoff){
  Table <- read_tsv(filename) %>%
             filter(InputCount>=inputcutoff)
  mean  <- rowMeans(as.matrix(select(Table, R1Rep1Fitness, R1Rep2Fitness, R1Rep3Fitness)))
  Table <- Table %>%
             mutate(fit=mean)
  return (Table)
  }

coloring <- function(mut){
  if (grepl('_',mut)){return ('Nonsense')}
  else{return ('Missense')}
  }

PlotCompareFit_Rep <- function(Tlib,ID){
  Plot_RepvsRep <- function(p,xlab,ylab){
    p <- p +
           geom_point(size=0.1) +
           xlab(xlab) +
           ylab(ylab) +
           scale_color_manual(values=c('#A0A000AF','#7001BAAF')) +
           #guides(colour = guide_legend(override.aes = list(size=1))) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 legend.title=element_blank(),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='none')
    return (p)
    } 
  print (ID)
  print (paste('Total Variants',length(Tlib$mut),sep=": "))
  print (paste('  Nonsense Variants',length(grep('_',Tlib$mut)),sep=": "))
  print (paste('  Missense Variants',length(grep('_',Tlib$mut,invert=TRUE)),sep=": "))
  print (paste("Rep1 vs Rep2",cor(log10(Tlib$R1Rep1Fitness),y=log10(Tlib$R1Rep2Fitness)),sep=": "))
  print (paste("Rep1 vs Rep3",cor(log10(Tlib$R1Rep1Fitness),y=log10(Tlib$R1Rep3Fitness)),sep=": "))
  print (paste("Rep2 vs Rep3",cor(log10(Tlib$R1Rep2Fitness),y=log10(Tlib$R1Rep3Fitness)),sep=": "))
  Tlib <- Tlib %>% 
             mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
             arrange(desc(color))
  textsize <- 14
  p1vs2 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep1Fitness),y=log10(R1Rep2Fitness),color=color)),
                         'Replicate 1','Replicate 2')
         
  p1vs3 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep1Fitness),y=log10(R1Rep3Fitness),color=color)),
                         'Replicate 1','Replicate 3')
         
  p2vs3 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep2Fitness),y=log10(R1Rep3Fitness),color=color)),
                         'Replicate 2','Replicate 3')
  p <- grid.arrange(p1vs2,p1vs3,p2vs3,ncol=3)
  return (p)
  }

inputcutoff  <- 20
Vic11WT_Tlib    <- ReadFitnessTable('result/Vic11WT_Tlib.count',inputcutoff) %>%
                     select(mut,R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness)
Vic11D190E_Tlib <- ReadFitnessTable('result/Vic11D190E_Tlib.count',inputcutoff) %>%
                     select(mut,R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness)
HK68WT_Tlib     <- ReadFitnessTable('result/HK68WT_Tlib.count',inputcutoff) %>%
                     select(mut,R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness)
HK68E190D_Tlib  <- ReadFitnessTable('result/HK68E190D_Tlib.count',inputcutoff) %>%
                     select(mut,R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness)
p_Vic11WT <- PlotCompareFit_Rep(Vic11WT_Tlib,'Vic11WT')
p_Vic11D190E <- PlotCompareFit_Rep(Vic11D190E_Tlib,'Vic11D190E')
p_HK68WT  <- PlotCompareFit_Rep(HK68WT_Tlib,'HK68WT')
p_HK68E190D  <- PlotCompareFit_Rep(HK68E190D_Tlib,'HK68E190D')
p <- grid.arrange(p_Vic11WT,p_Vic11D190E,p_HK68WT,p_HK68E190D,nrow=4)
ggsave('graph/Tlib_RepCompare.png',p,height=8,width=6.5)
