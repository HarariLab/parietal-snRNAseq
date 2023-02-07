library(RColorBrewer)
library(ggplot2)

setwd('~/SingleCellProjects/MyProjects/67BrainsPaper/Oligodendrocytes/unfiltered_nebula_results/')

genes = c('HNRNPA2B1','HNRNPC','HNRNPD','HNRNPH3','HNRNPM','HNRNPU','APP','CLU','MAP1B','PICALM')
tables = lapply(0:8,function(x) { return( read.table(sprintf('nebula_compareSubclusters_oligo_%s_all_03.25.21.txt',x)) ) } )
fullDF = data.frame(Cluster=c(0:8))
pvalDF = data.frame(Cluster=c(0:8))

plots  = lapply(genes[genes %in% rownames(tables[[1]])], function(gene){
	df = data.frame(Cluster = c(0:8),Estimate = unlist(lapply(tables, function(t){ return(t[gene,1]) })), pval = unlist(lapply(tables, function(t){ return(t[gene,3]) })) )
	df$log2FC = log2(exp(df$Estimate))
	fullDF[,gene] <<- df$log2FC
	pvalDF[,gene] <<- df$pval
})

library(tidyr)
df = gather(fullDF,'Gene','log2FC',2:dim(fullDF)[2])

df2 = gather(pvalDF,'Gene','pval',2:dim(pvalDF)[2])
df$pval = df2$pval

p = ggplot(df,aes(factor(as.character(Cluster)),y=factor(Gene,levels=rev(unique(Gene))),fill=log2FC)) + geom_tile() + geom_point(data = df[df$pval < 0.05,]) +
	scale_fill_distiller(palette='RdBu',limits = max(df$log2FC) * c(-1,1)) +#scale_fill_manual()+
	theme(panel.background=element_rect(fill = "white",color='white'),axis.title.x=element_blank(),axis.title.y=element_blank(),
	axis.ticks=element_blank(),axis.text = element_text(size=16,color='black'),legend.text=element_text(size=13,color='black'),legend.title=element_text(size=15,color='black')) + labs(fill='Log2(FC)')
ggsave(p,h=5.25,w=6.75,file = '~/scratch/HNRNP_oligo3_heatmap.pdf')


write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig7i.csv', sep = ', ', quote = F, row.names = F, col.names = T)
