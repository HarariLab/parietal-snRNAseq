library(RColorBrewer)
library(ggplot2)

#o = readRDS('SingleCellProjects/dataObjects/oligo.rds')
df = data.frame(Cluster=c(0:8),Estimate=c(-0.0716154594641995,0.0559990010777891,0.0449209227658119,-0.108518980112217,0.0618656679358859,0.102760034290545,-0.103419292657633,0.0638871048734107,-0.0847450222094133))
df$log2FC = log2(exp(df$Estimate))
df$Cluster = factor(df$Cluster,levels=rev(df$Cluster))

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig4f.csv', sep = ', ', quote = F, row.names = F, col.names = T)


#myColors = c(brewer.pal(9, 'Set1'))
#myColors[6] = '#EBEB00'
#scale_fill_manual(values=rev(myColors))
#p = ggplot(df,aes(x=Cluster,y=log2FC,fill=Cluster)) + geom_bar(stat='identity') + 
#		theme(panel.background=element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray"),
#					panel.grid.major.y = element_blank(),legend.position = 'none',
#					axis.text.y= element_text(size=24),axis.text.x= element_text(size=24),
#					axis.title.x = element_text(size=28),axis.title.y = element_text(size=28)) + 
#		scale_fill_manual(values=rev(myColors)) + scale_y_continuous(breaks=c(-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15)) + 
#		coord_flip() + ylab('TFEB log2FC')
#ggsave(p,w=3.5,file='~/scratch/oligo5_TFEB_logfc.pdf')
#