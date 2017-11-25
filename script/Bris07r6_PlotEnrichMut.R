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

fittable2genomap <- function(fittable,mutsofinterest){
  rank_list <- c()
  mut_list  <- c()
  values    <- c()
  for (n in 1:length(fittable$mut)){
    mut  <- fittable$mut[n]
    rank <- fittable$rank[n]
    muts <- strsplit(mut,'-')[[1]]
    restofmuts <- mutsofinterest 
    if (mut=='WT'){next}
    for (m in 1:length(muts)){
      rank_list  <- c(rank_list, rank)
      mut_list   <- c(mut_list, muts[m])
      values     <- c(values, 1)
      restofmuts <- restofmuts[restofmuts!=muts[m]]
      }
    if (length(restofmuts)==0){next}
    for (m in 1:length(restofmuts)){
      rank_list  <- c(rank_list, rank)
      mut_list   <- c(mut_list, restofmuts[m])
      values     <- c(values, 0)
      }
    }
  genomap <- data.frame(cbind(rank_list,mut_list,values)) %>% 
	       rename(rank=rank_list) %>% 
	       rename(mut=mut_list) %>%
               rename(value=values) %>%
               mutate(value=as.numeric(value)) %>%
               mutate(mut=factor(mut,levels=rev(mutsofinterest)))
  return (genomap)
  }

plotrankmap <- function(fittable, graphname){
  textsize <- 9
  colorscale  <- c(brewer.pal(9,"GnBu"))
  p <- ggplot(fittable,aes(x=rank, y=log10(fitness),fill=mutclass)) +
         geom_rect(data=NULL,aes(xmin=500,xmax=520,ymin=3,ymax=4.2), color=NA, fill='grey85', alpha=1) +
	 geom_point(color='black',pch=20,size=1) +
	 scale_fill_gradient2(limits=c(0,9),midpoint=4.5,low=colorscale[1],mid=colorscale[5],high=colorscale[9],
			      breaks=c(0,9),labels=c(0,9),
			      name='# of substitutions',
			      guide=guide_colorbar(title.position="top",
						   title.vjust=0,
						   label.position="bottom")) +
	 xlab(bquote(bold(Mutant~ranked~by~RF~index))) +
	 ylab(bquote(bold(Log['10']~RF~index))) +
	 theme(axis.title=element_text(size=textsize,face="bold"),
		 axis.text.y=element_text(size=textsize,face="bold"),
		 axis.text.x=element_blank(),
		 axis.ticks.x=element_blank(),
		 legend.title=element_text(size=textsize,face="bold"),
		 legend.text=element_text(size=textsize,face="bold"),
		 legend.key.size=unit(2,'mm'),
		 legend.position='none',
		 legend.justification='center',
		 legend.direction = "horizontal")
  ggsave(paste(graphname,sep=''),p,height=2,width=3)
  }

plotheatmap <- function(genomap, graphname){
  p <- ggplot(data=genomap, aes(x=rank, y=mut, fill=value)) +
              geom_tile() +
              scale_fill_gradientn(colours=c("white", "yellow", "red"),
                         values=rescale(c(0, 0.5, 1)),
                         guide='none',
                         na.value="grey") +
              theme_classic() +
              theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
                    text = element_text(size=9,colour='black', face='bold'),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.title.x=element_blank())
  ggsave(graphname,p, height=1.2, width=3)
  }

mutsofinterest <- c('T155H','H156Q','F159Y','K160R','N189S','I192T','A196T','I202V','R222W')
fittable <- read_tsv('result/Bris07r6_MutfitTable.tsv') %>%
              arrange(fitness) %>%
              mutate(rank=1:length(.$fitness)) %>%
              filter(fitness>-100)
genomap  <- fittable2genomap(fittable,mutsofinterest)
plotrankmap(fittable, 'graph/Bris07r6_MutEnich.png')
