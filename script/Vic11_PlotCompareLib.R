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

PlotCompareFit_Lib <- function(Tlib1, Tlib2, xlab, ylab, graphname){
  Tlib1 <- Tlib1 %>% rename(fit1=fit)
  Tlib2 <- Tlib2 %>% rename(fit2=fit)
  t     <- inner_join(Tlib1,Tlib2) %>%
             mutate(color=factor(mapply(coloring,mut),levels=c('Nonsense','Missense'))) %>%
             arrange(desc(color))
  textsize <- 11
  p <- ggplot(t,aes(x=log10(fit1),y=log10(fit2),color=color)) +
         #geom_rect(data=NULL,aes(xmin=-3.5,xmax=-0.5,ymin=4,ymax=5), color=NA, fill='grey85', alpha=1) +
         #geom_rect(data=NULL,aes(xmin=log10(0.5),xmax=0.8,ymin=-0.2,ymax=2.3), color=NA, fill='grey85', alpha=1) +
         geom_point(size=0.1) +
         xlab(bquote(bold(Log['10']~.(xlab)))) +
         ylab(bquote(bold(Log['10']~.(ylab)))) +
         scale_color_manual(values=c('#A0A000AF','#7001BAAF')) +
         theme(axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         guides(colour = guide_legend(override.aes = list(size=1)))
  ggsave(graphname,p,height=3,width=4)
  print (graphname)
  print (cor(log10(t$fit1),log10(t$fit2)))
  print (filter(t,fit2>10000))
  }

inputcutoff  <- 20
Vic11WT_Tlib    <- ReadFitnessTable('result/Vic11WT_Tlib.count',inputcutoff) %>%
                     select(mut,fit)
Vic11D190E_Tlib <- ReadFitnessTable('result/Vic11D190E_Tlib.count',inputcutoff) %>%
                     select(mut,fit)
HK68WT_Tlib     <- ReadFitnessTable('result/HK68WT_Tlib.count',inputcutoff) %>%
                     select(mut,fit)
HK68E190D_Tlib  <- ReadFitnessTable('result/HK68E190D_Tlib.count',inputcutoff) %>%
                     select(mut,fit)
PlotCompareFit_Lib(Vic11WT_Tlib,Vic11D190E_Tlib,'RF index (Vic11 WT)','RF index (Vic11 D190E)',
                   'graph/Tlib_RFcompare_Vic11WTvsVic11D190E.png')
PlotCompareFit_Lib(HK68WT_Tlib,Vic11WT_Tlib,'RF index (HK68 WT)','RF index (Vic11 WT)',
                   'graph/Tlib_RFcompare_Vic11WTvsHK68WT.png')
PlotCompareFit_Lib(HK68WT_Tlib,HK68E190D_Tlib,'RF index (HK68 WT)','RF index (HK68 E190D)',
                   'graph/Tlib_RFcompare_HK68E190DvsHK68WT.png')

HK68WT_Tlib <- rename(HK68WT_Tlib,HK68WT_fit=fit)
Vic11WT_Tlib <- rename(Vic11WT_Tlib,Vic11WT_fit=fit)
Vic11D190E_Tlib <- rename(Vic11D190E_Tlib,Vic11D190E_fit=fit)
t <- inner_join(HK68WT_Tlib, Vic11WT_Tlib)
filter(t, HK68WT_fit>1 & Vic11WT_fit>1)
