library(ggplot2)
library(stringr)
library(ggpubr)

files = system('ls /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/*/*/nebula*',intern=T)
files = files[grep('all',files)]
fs = NULL

phenotypes = c('ADADvsCO','sADvsCO')
celltypes = c('astro','micro','n_ex','n_inh','oligo','opc')

allOverlaps = lapply(celltypes, function(ctype){
	fs <<- files[grep(ctype,files)]
	DEGs = lapply(phenotypes, function(pheno){
		tab = read.table(fs[grep(pheno,fs)],sep='\t',head=T,stringsAsFactor=F)
		tmp = tab[tab$p.value < 0.05,1]
		names(tmp) = rownames(tab[tab$p.value < 0.05,])
		return(tmp)
	})
	inter = intersect(names(DEGs[[1]]),names(DEGs[[2]]))
	myCor = cor(DEGs[[1]][inter],DEGs[[2]][inter])
	df = data.frame(DEGs[[1]][inter],DEGs[[2]][inter])
	colnames(df) = phenotypes
	p = ggplot(df,aes(x=sADvsCO,y=ADADvsCO)) + geom_point() + ggtitle(str_to_title(sprintf('%s (R=%s)',ctype,round(myCor,2)))) + geom_smooth(method='lm') +
		xlab('sADvsCO Effect') + ylab('ADADvsCO Effect') + geom_abline(slope=1,intercept=0,color='red') +
		theme(axis.title = element_text(size=16),title = element_text(size=20),axis.text=element_text(size=14,color='black')) +
		theme_minimal()
	#plot(DEGs[[1]][inter],DEGs[[2]][inter])
	#ggsave(p,file=sprintf('~/scratch/sAD_ADAD_effects_%s.pdf',ctype),h=6,w=6)
	return(p)
})
p = ggarrange(plotlist=allOverlaps,labels=letters[1:6],ncol=3,nrow=2)
ggsave(p,file='~/scratch/sAD_ADAD_effects_scatterplots.pdf',h=12,w=18)


df = NULL
tmp = lapply(seq_len(length(allOverlaps)), function(i){
    tmp = allOverlaps[[i]]$data
    tmp$celltype = celltypes[i]
    df <<- rbind(df,tmp)
})

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig4.csv', sep = ', ', quote = F, row.names = T, col.names = T) # nolint: line_length_linter.
