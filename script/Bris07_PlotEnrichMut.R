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

trimD190E <- function(m,MutofInterest){
  if (m==MutofInterest){return (m)}
  else {return (str_replace(str_replace(m,paste('/',MutofInterest,sep=''),''),paste(MutofInterest,'/',sep=''),''))}
  }

labeling <- function(mut,D190Efit,WTfit){
  if (D190Efit>10){return (mut)}
  else {return ('')}
  }

fittable <- read_tsv('result/Bris07_MutfitTable.tsv') %>%
              mutate(mut=mapply(function(m){return (str_replace_all(m,'-','/'))},mut))
MutofInterest <- 'D190E'
WT_table    <- fittable %>% 
                 filter(!grepl(MutofInterest,mut)) %>%
                 select(mut,mutclass,fitness) %>%
                 rename(WTfit=fitness)
D190E_table <- fittable %>% 
                 filter(grepl(MutofInterest,mut)) %>%
                 mutate(mut=mapply(trimD190E,mut,rep(MutofInterest,length(mut)))) %>%
                 select(mut,fitness) %>%
                 rename(D190Efit=fitness)

combinetable <- inner_join(WT_table,D190E_table) %>%
                  mutate(label=mapply(labeling,mut,D190Efit,WTfit))
combinetable_label <- filter(combinetable,label!='')

textsize <- 9
colorscale  <- c(brewer.pal(9,"GnBu"))
p <- ggplot(combinetable,aes(x=log10(WTfit),y=log10(D190Efit),fill=mutclass)) +
       geom_point(pch=21,color='black') + 
       scale_fill_gradient2(limits=c(0,9),midpoint=4.5,low=colorscale[1],mid=colorscale[5],high=colorscale[9],
                            breaks=c(0,9),labels=c(0,9),
                            name='# of substitutions',
                            guide=guide_colorbar(title.position="top",
                                                 title.vjust=0,
                                                 label.position="bottom")) +
       #xlab(bquote(bold(Log['10']~RF~index~'('*without~.(MutofInterest)*')'))) +
       #ylab(bquote(bold(Log['10']~RF~index~'('*with~.(MutofInterest)*')'))) +
       xlab(bquote(bold(Log['10']~RF~index~'(with D190)'))) +
       ylab(bquote(bold(Log['10']~RF~index~'(with E190)'))) +
       ylim(-1.5,1.5) +
       theme(axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(2,'mm'),
               legend.position='top',
               legend.justification='center',
               legend.direction = "horizontal") +
       geom_text_repel(data=combinetable_label,aes(x=log10(WTfit),y=log10(D190Efit),label=label),
                         size = 2.4,
                         fontface = 'bold',
                         segment.size = 0.2,
                         point.padding = unit(0.1, "lines"),
                         force = 0.2,
                         nudge_x=-2,
                         show.legend = FALSE)
ggsave(paste('graph/Bris07_MutEnich_',MutofInterest,'.png',sep=''),p,height=3,width=3.2)
