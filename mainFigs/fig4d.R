library(Seurat)
library(ggplot2)
library(ggsignif)
library(lme4)
library(lmerTest)
library(tibble)
library(grDevices)
library(ggpubr)
library(RColorBrewer)

m = readRDS('~/SingleCellProjects/dataObjects/Colonna_TREM2_data/colonna_seurat_object_microglia.rds')

files = system('ls /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/micro/bh_nebula_results/bh_nebula_compareSubclusters_microglia_*_all*',intern=T)
sigGenes = lapply(files,function(f){
	tmp = read.table(f,stringsAsFactors=F,head=T,sep='\t')
	z = list(rownames(tmp[tmp$z.score > 0 & abs(tmp$logFC_cluster) > 0.25,]),rownames(tmp[tmp$z.score < 0 & abs(tmp$logFC_cluster) > 0.25,]))
	names(z) = c('upReg','downReg')
	return(z)
})

upGenes = lapply(1:length(sigGenes),function(cluster){
	up = sigGenes[[cluster]]$upReg[sigGenes[[cluster]]$upReg%in% rownames(m)]
	cat(sprintf('cluster_%s upReg genes: %s/%s\n',cluster-1,length(up),length(sigGenes[[cluster]]$upReg)))
	return(up)
})
m = AddModuleScore(m,upGenes,name='upSig')
colnames(m@meta.data)[grep('upSig',colnames(m@meta.data))] = paste0('upSig',0:(length(files)-1))

cluster = 2
df = data.frame(Cluster = as.factor(m$newClusters), upSig =  m@meta.data[,paste0('upSig',cluster)], Sample = m$orig.ident)
df2 = df
df2$Cluster = factor(df2$Cluster, levels=c(cluster,c(0:(length(files)-1))[-(cluster+1)]))
re = lmer(upSig ~ Cluster + (1|Sample), data = df2)
stats = tibble(group1=as.character(cluster),group2=levels(df2$Cluster)[-1],p.value=sprintf('β=%s, P=%s',signif(summary(re)$coefficients[-1,'Estimate'], digits=1),signif(summary(re)$coefficients[-1,'Pr(>|t|)'], digits=3)))
stats$p.value = '***'
stats$p.value[8] = '**'	
Cluster = df$Cluster

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig4d.csv', sep = ', ', quote = F, row.names = F, col.names = T)

myColors = c(brewer.pal(9, 'Set1'))
myColors[6] = '#EBEB00'
p = ggplot(df,aes(x=Cluster,y=upSig)) + geom_jitter(size=0.1) + geom_violin(draw_quantiles=c(0.25,0.5,0.75),alpha=0.8,aes(fill=Cluster)) + 
		ggtitle(paste0('ROSMAP Mic_reduced Signature')) + ylab('Up Signature') + xlab('Micro State') + scale_fill_manual(values=myColors) + theme_bw() + 
		theme(axis.title = element_text(size=16),axis.text= element_text(size=14,color='black')) +
		stat_summary(fun = "mean", geom = "crossbar", width = 0.7, colour = "red") + stat_pvalue_manual(stats,label='p.value',y.position=max(df$upSig), step.increase = 0.05, tip.length=0.01) + rremove('legend')

cairo_pdf(h=7.09,w=3.545,file='~/scratch/micro3_upSig_rs1582763_replication.pdf')
p
dev.off()