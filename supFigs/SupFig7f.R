library(RColorBrewer)
library(ggplot2)

setwd('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/astro/bh_nebula_results/')

daa = c('SLC1A3','CST3','NCAM2','VIM','KCNIP4','PLCE1','CADM2','GRM5','ENOX1')
daa2 = c('APOE','CLU','MT1X','CPE','CD9','ID3','CD81','FXYD1','CKB','PRDX6','CSMD1','CTSB','DBI','HSD17B4','FTH1')
gfapHigh = c('GFAP','ID3','CKB','AQP4','SPARCL1','FBXO2','ID1','SERPINF1','FABP7')#'MYOC' not found in data
gfapLow = c('GPC5','TRPM3','LSAMP','SLC7A10')
genes = c(daa,daa2,gfapHigh,gfapLow)
genes = c('LUZP2','SLC7A10','MFGE8','GFAP','ID3','AQP4','ID1','FABP7','CTSB','VIM','OSMR','GSN')#'GGTA1','SERPINA3'
tables = lapply(0:4,function(x) { return( read.table(sprintf('nebula_compareSubclusters_astrocyte_%s_all_06.02.21.txt',x)) ) } )
fullDF = data.frame(Cluster=c(0:4))
pvalDF = data.frame(Cluster=c(0:4))

plots  = lapply(genes[genes %in% rownames(tables[[1]])], function(gene){
	df = data.frame(Cluster = c(0:4),Estimate = unlist(lapply(tables, function(t){ return(t[gene,1]) })), pval = unlist(lapply(tables, function(t){ return(t[gene,3]) })) )
	df$log2FC = log2(exp(df$Estimate))
	fullDF[,gene] <<- df$log2FC
	pvalDF[,gene] <<- df$pval
})

library(tidyr)
df = gather(fullDF,'Gene','log2FC',2:dim(fullDF)[2])
df$log2FC[df$log2FC < -2.5] = -2.5
df$log2FC[df$log2FC > 2.5] = 2.5

df2 = gather(pvalDF,'Gene','pval',2:dim(pvalDF)[2])
df$pval = df2$pval

p = ggplot(df,aes(factor(as.character(Cluster)),y=factor(Gene,levels=rev(unique(Gene))),fill=log2FC)) + geom_tile() + geom_point(data = df[df$pval < 0.05,]) +
	scale_fill_distiller(palette='RdBu',limits = max(df$log2FC) * c(-1,1)) +#scale_fill_manual()+
	theme(panel.background=element_rect(fill = "white",color='white'),axis.title.x=element_blank(),axis.title.y=element_blank(),
	axis.ticks=element_blank(),axis.text = element_text(size=16,color='black'),legend.text=element_text(size=13,color='black'),legend.title=element_text(size=15,color='black')) + labs(fill='Log2(FC)')
ggsave(p,h=5.25,w=4.25,file = '~/scratch/DAA_GFAP_heatmap.pdf')

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig7f.csv', sep = ', ', quote = F, row.names = F, col.names = T)



#myLabels = c('x<-5','-5<x<-4','-4<x<-3','-3<x<-2','-2<x<-1','-1<x<0','0<x<1','1<x<2','2<x<3','3<x<4','4<x<5','5<x')
#df$Z.Score.split <- cut(df$Z.Score,breaks = c(min(df$Z.Score),-5,-4,-3,-2,-1,0,1,2,3,4,5, max(df$Z.Score)),#c(min(df$Z.Score),-10,-5,0,5,10,15,20, max(df$Z.Score)),
#                        labels = myLabels, include.lowest = T)#labels = c("-15<x<-10", "-10<x<-5", "-5<x<0", "0<x<5", "5<x<10", "10<x<15", "15<x<20", "20<x<25"), include.lowest = T)

#p= ggplot(df,aes(x=Cluster,y=factor(Genes,levels=rev(rownames(tmp2))),fill=factor(Z.Score.split,levels=rev(myLabels)))) + geom_tile() + scale_fill_manual(values=colorRampPalette(c('blue','white','red'))(10)) +#scale_fill_manual(values=c(colorRampPalette(c('red','white'))(5),colorRampPalette(c('white','blue'))(5)))#
#	theme(panel.background=element_rect(fill = "white",color='white'),axis.title.x=element_blank(),axis.title.y=element_blank(),
#	axis.ticks=element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12)) + labs(fill='Z.Score')