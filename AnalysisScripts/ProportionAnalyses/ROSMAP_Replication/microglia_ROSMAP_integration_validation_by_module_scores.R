library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(lme4)
library(lmerTest)
library(tibble)
library(grDevices)

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
downGenes = lapply(1:length(sigGenes),function(cluster){
	down = sigGenes[[cluster]]$downReg[sigGenes[[cluster]]$downReg%in% rownames(m)]
	cat(sprintf('cluster_%s downReg genes: %s/%s\n',cluster-1,length(down),length(sigGenes[[cluster]]$downReg)))
	return(down)
})
	
m = AddModuleScore(m,upGenes,name='upSig')
colnames(m@meta.data)[grep('upSig',colnames(m@meta.data))] = paste0('upSig',0:(length(files)-1))

m = AddModuleScore(m,downGenes,name='downSig')
colnames(m@meta.data)[grep('downSig',colnames(m@meta.data))] = paste0('downSig',0:(length(files)-1))

upSigPlots = lapply(0:(length(files)-1),function(cluster){
	df = data.frame(Cluster = as.factor(m$newClusters), upSig =  m@meta.data[,paste0('upSig',cluster)], Sample = m$orig.ident)
	df2 = df
	df2$Cluster = factor(df2$Cluster, levels=c(cluster,c(0:(length(files)-1))[-(cluster+1)]))
	re = lmer(upSig ~ Cluster + (1|Sample), data = df2)
	stats = tibble(group1=as.character(cluster),group2=levels(df2$Cluster)[-1],p.value=sprintf('β=%s, P=%s',signif(summary(re)$coefficients[-1,'Estimate'], digits=1),signif(summary(re)$coefficients[-1,'Pr(>|t|)'], digits=3)))

	p = ggviolin(df,x='Cluster',y='upSig',fill = 'Cluster',draw_quantiles=c(0.25,0.5,0.75),trim=T) + geom_jitter(size=0.1) + ggtitle(sprintf('Microglia.%s upregulated signature (ROSMAP nuclei)',cluster)) + 
		stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red") + stat_pvalue_manual(stats,label='p.value',y.position=max(df$upSig), step.increase = 0.05, tip.length=0.01) + ylab(sprintf('Microglia.%s upregulated gene signature',cluster)) + xlab('Microglia cluster') + rremove('legend')
	#ggsave(p,file = sprintf('~/scratch/micro%s_upSig_plot_logFC_0.25_genes.png',cluster),h=10)
	return(p)})
saveRDS(upSigPlots, file='~/scratch/micro_upSig_plots_logFC_0.25_genes.rds')
q = ggarrange(plotlist=upSigPlots,ncol=3,nrow=3)
ggsave(q,file='scratch/microgliaUpregulatedSignatures.pdf',h=30,w=21,device=cairo_pdf)

downSigPlots = lapply(0:(length(files)-1),function(cluster){
	df = data.frame(Cluster = as.factor(m$newClusters), downSig =  m@meta.data[,paste0('downSig',cluster)], Sample = m$orig.ident)
	df2 = df
	df2$Cluster = factor(df2$Cluster, levels=c(cluster,c(0:(length(files)-1))[-(cluster+1)]))
	re = lmer(downSig ~ Cluster + (1|Sample), data = df2)
stats = tibble(group1=as.character(cluster),group2=levels(df2$Cluster)[-1],p.value=sprintf('β=%s, P=%s',signif(summary(re)$coefficients[-1,'Estimate'], digits=1),signif(summary(re)$coefficients[-1,'Pr(>|t|)'], digits=3)))

	p = ggviolin(df,x='Cluster',y='downSig',fill = 'Cluster',draw_quantiles=c(0.25,0.5,0.75),trim=T) + geom_jitter(size=0.1) + ggtitle(sprintf('Microglia.%s downregulated signature (ROSMAP nuclei)',cluster)) + 
		stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red") + stat_pvalue_manual(stats,label='p.value',y.position=max(df$downSig), step.increase = 0.05, tip.length=0.01) + ylab(sprintf('Microglia.%s downregulated gene signature',cluster)) + xlab('Microglia cluster')+ rremove('legend')
	#ggsave(p,file = sprintf('~/scratch/micro%s_downSig_plot_logFC_0.25_genes.png',cluster),h=10)
	return(p)})
saveRDS(downSigPlots, file='~/scratch/micro_downSig_plots_logFC_0.25_genes.rds')
q = ggarrange(plotlist=downSigPlots,ncol=3,nrow=3)
ggsave(q,file='scratch/microgliaDownregulatedSignatures.pdf',h=30,w=21,device=cairo_pdf)
