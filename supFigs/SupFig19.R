library(ggplot2)
library(ggpubr)
library(grid)

s40Removed_oligo = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample40_removed_oligo_all.txt',sep='\t',head=T)
ADAD_oligo = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample40_removed_oligo_all.txt',sep='\t',head=T)
df = data.frame(s40_Removed = s40Removed_oligo[rownames(ADAD_oligo),1], s40_Retained = ADAD_oligo[,1])
#ggplot(df[df[,1] < 19,],aes(s40_Retained, s40_Removed)) + geom_point() + ggtitle('ADAD Oligo Estimates') + geom_smooth(method='lm') + geom_abline(slope=1,intercept=0,color='red')
df2 = data.frame(s40_Removed = -log10(s40Removed_oligo[rownames(ADAD_oligo),'p.value']), s40_Retained = -log10(ADAD_oligo[,'p.value']))
#ggplot(df2,aes(s40_Retained, s40_Removed)) + geom_point() + ggtitle('ADAD Oligo P-values (-log10)') + geom_smooth(method='lm') + geom_abline(slope=1,intercept=0,color='red')

comparisons = list(
	c('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample40_removed_oligo_all.txt', '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_both_ADAD_included_oligo_all.txt','ADAD_Oligo'),
	c('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample40_removed_microglia_all.txt', '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_both_ADAD_included_microglia_all.txt','ADAD_Micro'),
	c('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample59_removed_oligo_all.txt', '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_both_sAD_included_oligo_all.txt','sAD_Oligo'),
	c('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_sample59_removed_microglia_all.txt', '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_both_sAD_included_microglia_all.txt', 'sAD_Micro')
)
myplots=list()
myEst = NULL
myPs = NULL
tmp = lapply(comparisons,function(files){
	Removed_data = read.table(files[1],sep='\t',head=T)
	Retained_data = read.table(files[2],sep='\t',head=T)
	df = data.frame(Removed = Removed_data[rownames(Retained_data),1], Retained = Retained_data[,1])
	df = na.omit(df)
    tmp = df
    tmp$Plot = sprintf('%s Estimates',files[3])
    myEst <<- rbind(myEst,tmp)
	slope = sd(df$Removed)/sd(df$Retained)
	p = ggplot(df,aes(Retained, Removed)) + geom_point() + ggtitle(sprintf('%s Estimates',files[3])) + geom_smooth(method='lm') + geom_abline(slope=1,intercept=0,color='red') + annotation_custom(grobTree(textGrob(sprintf('Slope = %s',round(slope,3)),x=0.2,y=0.9, gp=gpar(fontsize=20,fontface='italic')))) + theme_minimal() + theme(title = element_text(size=24),axis.title=element_text(size=20),axis.text=element_text(size=16,color='black'))
	#ggsave(p,file=sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_%s_Estimates.pdf',files[3]))
	myplots[[length(myplots)+1]] <<- p
	df2 = data.frame(Removed = -log10(Removed_data[rownames(Retained_data),'p.value']), Retained = -log10(Retained_data[,'p.value']))
	df2 = na.omit(df2)
    tmp = df2
    tmp$Plot = sprintf('%s P-values (-log10)',files[3])
    myPs <<- rbind(myPs,tmp)
	slope = sd(df2$Removed)/sd(df2$Retained)
	p = ggplot(df2,aes(Retained, Removed)) + geom_point() + ggtitle(sprintf('%s P-values (-log10)',files[3])) + geom_smooth(method='lm') + geom_abline(slope=1,intercept=0,color='red') + annotation_custom(grobTree(textGrob(sprintf('Slope = %s',round(slope,3)),x=0.2,y=0.9, gp=gpar(fontsize=20,fontface='italic')))) + theme_minimal() + theme(title = element_text(size=24),axis.title=element_text(size=20),axis.text=element_text(size=16,color='black'))
	#ggsave(p,file=sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_%s_Pvalues.pdf',files[3]))
	myplots[[length(myplots)+1]] <<- p
})
q = ggarrange(plotlist=myplots,ncol=2,nrow=4,labels=letters[1:length(myplots)],font.label=list(size=24))
ggsave(q,file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/GeneticallyRelatedAnalyses/Genetically_Related_corPlots.pdf',h=28,w=14)


write.table(myEst, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig19est.csv', sep = ', ', quote = F, row.names = F, col.names = T)
write.table(myPs, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig19p.csv', sep = ', ', quote = F, row.names = F, col.names = T)
