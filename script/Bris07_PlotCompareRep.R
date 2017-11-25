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

PlotCompareFit_Rep <- function(Tlib,ID){
  Plot_RepvsRep <- function(p,xlab,ylab){
    p <- p +
           geom_point(size=0.1,color='#7001BAAF') +
           xlab(xlab) +
           ylab(ylab) +
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
  textsize <- 14
  p1vs2 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep1Fitness),y=log10(R1Rep2Fitness))),
                         'Replicate 1','Replicate 2')
         
  p1vs3 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep1Fitness),y=log10(R1Rep3Fitness))),
                         'Replicate 1','Replicate 3')
         
  p2vs3 <- Plot_RepvsRep(ggplot(Tlib,aes(x=log10(R1Rep2Fitness),y=log10(R1Rep3Fitness))),
                         'Replicate 2','Replicate 3')
  p <- grid.arrange(p1vs2,p1vs3,p2vs3,ncol=3)
  return (p)
  }

inputcutoff  <- 7
Bris07_lib   <- ReadFitnessTable('result/Bris07_MultiMutLib.tsv',inputcutoff) %>%
                     select(mut,R1Rep1Fitness,R1Rep2Fitness,R1Rep3Fitness)
p_Bris07 <- PlotCompareFit_Rep(Bris07_lib,'Bris07')
ggsave('graph/Bris07_RepCompare.png',p_Bris07,height=3,width=6.5)
