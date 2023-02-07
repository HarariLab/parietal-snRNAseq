activated = c('ABCA1','RELB','GPNMB','CD68','C5AR1','TNFAIP3','CD83')
resting = c('TMEM119','P2RY13','CX3CR1','BIN1','MED12L','SELPLG')
genes = c(resting, activated)

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
myColors = c(brewer.pal(9, 'Set1'))
myColors[6] = '#EBEB00'

#corrected = readRDS('~/SingleCellProjects/Collabs/KipnisCollab/correctedCounts_betweenClustersDE.rds')
corrected = readRDS('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/microglia_resting_activated_genes_partialResidual.rds')
corrected = corrected[,c('clusterSet',genes)]
df = aggregate(corrected,by=list(corrected$clusterSet),FUN='mean')[-c(4:9),]
df = gather(df,'Gene','Expression',3:dim(df)[2])
df$Cluster = df$Group.1
df$Expression = df$Expression - (min(df$Expression)-0.1)

write.table(df[,-c(1,2)], '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig4b.csv', sep = ', ', quote = F, row.names = F, col.names = T) # nolint: line_length_linter.

p = ggplot(df,aes(x=factor(Gene,levels=rev(genes)), y=Expression, fill=factor(as.character(Cluster),levels=c(2,1,0)) ))+#, color=factor(as.character(Cluster)), alpha=factor(as.character(Cluster)))) +
	geom_bar(width=0.5,stat="identity", position =position_dodge())+#"identity") +
  #scale_colour_manual(values=rev(myColors[1:3])) + #values=c('green','red','blue')) + #
  scale_fill_manual(values=rev(myColors[1:3]),guide=guide_legend(reverse=T)) + #values=c('green','red','blue')) + #
  theme(axis.title.y=element_blank(), axis.text.y=element_text(size=14),axis.text.x=element_text(size=14) ) +
	labs(fill='Cluster') + coord_flip() 

ggsave(p,w=3.5,file='~/scratch/mic2_restingActivatedMarkers_boxplot.pdf')
