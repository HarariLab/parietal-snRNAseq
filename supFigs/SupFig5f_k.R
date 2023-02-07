library(grid)
#library(UpSetR)
library(ComplexUpset) #https://krassowski.github.io/complex-upset/articles/Examples_R.html#0-basic-usage
library(ggplot2)
library(openxlsx)
library(tidyr)

files = system('ls /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/*/*/nebula*',intern=T)
files = files[grep('all',files)]
fs = NULL

phenotypes = c('ADADvsCO','ADADvssAD','APOE','sADvsCO','TREM2vsCO')
celltypes = c('astro','micro','n_ex','n_inh','oligo','opc')
DEGlist = list()
binaryTables = list()
allOverlaps = lapply(celltypes, function(ctype){
	fs <<- files[grep(ctype,files)]
	DEGs = lapply(phenotypes, function(pheno){
		tab = read.table(fs[grep(pheno,fs)],sep='\t',head=T,stringsAsFactor=F)
		return(rownames(tab[tab$BH < 0.05,]))
	})
	names(DEGs) = phenotypes 
	DEGlist[[ctype]] <<- DEGs
	overlap = as.data.frame.array(table(stack(DEGs)))
	overlap$groups = apply(overlap,1,function(values) return(paste(colnames(overlap)[which(values == 1)],collapse='_')) )
	binaryTables[[length(binaryTables)+1]] <<- overlap	

	maxG = max(table(overlap$groups))
	groups = apply(crossing(v1 = 0:1, v2=0:1, v3=0:1,v4=0:1,v5=0:1),1,function(values) return(paste(colnames(overlap)[which(values == 1)],collapse='_')) )[-1]
	myDF = data.frame(X = rep('',maxG))
	tmp = lapply(groups[order(unlist(lapply(groups,function(g) length(unlist(strsplit(g,'_'))))),decreasing=T)], function(group) {
		genes = rownames(overlap[overlap$groups == group,])
		myDF[,group] <<- c(genes, rep('', maxG-length(genes)))
	})
	return(myDF[,-1])
})
names(allOverlaps) = paste0(celltypes,'_geneOverlaps')
names(binaryTables) = celltypes
#saveRDS(allOverlaps,file = '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/overlap_DEGs_between_phenos_within_cell_type_v2.rds')
#write.xlsx(allOverlaps,file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/overlap_DEGs_between_phenos_within_cell_type_v2.xlsx')


binaryTables = lapply(binaryTables,function(tab) {colnames(tab) = c('ADAD vs CO','ADAD vs sAD','APOEe4+ vs e4-','sAD vs CO','TREM2 vs CO','Groups'); return(tab)} )
names(binaryTables) = c('Astrocytes','Microglia','E. Neurons','I. Neurons','Oligodendrocytes','OPCs')
myPlots = lapply(names(binaryTables),function(ctype) {
	maxDeg = 1
	#if(grepl('ADAD',pheno)) maxDeg=3
	p = upset(as.data.frame(binaryTables[[ctype]]), colnames(binaryTables[[1]]),sort_intersections_by = c('degree','cardinality'),min_degree=maxDeg,width_ratio=0.1,stripes=c('white', 'grey90'),
			matrix=(intersection_matrix(outline_color=list(active='black',inactive='gray95'))+scale_color_manual(values=c('TRUE'='black','FALSE'='gray95')) ),
			set_sizes=(upset_set_size()+theme(axis.text.x=element_text(angle=90,size=14),axis.title.x=element_text(size=16,color='black'))),
			base_annotations = list('Intersection size'=(intersection_size() + ggtitle(ctype) + theme(title=element_text(size=20),axis.text.y=element_text(size=18,color='black'),axis.title.y=element_text(size=19)) )),
			themes=upset_modify_themes(list('overall_sizes'=theme(axis.text=element_text(color='black')),'intersections_matrix'=theme(legend.position='none',text=element_text(size=16),axis.text=element_text(color='black')))) 
	)
	#ggsave(p,file=sprintf('~/scratch/upsetPlot_%s.pdf',ctype),h=5,w=14)
	return(p)
})
#saveRDS(myPlots,file='~/scratch/upsetPlots_cellType.rds')

df = NULL
tmp = lapply(seq_len(length(binaryTables)),function(i) {
    tab = binaryTables[[i]]
    tab$Plot = names(binaryTables)[i]
    df <<- rbind(df,tab)
    return(NULL)
})

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig5f_k.csv', sep = ', ', quote = F, row.names = T, col.names = T)

ctype= celltypes[1]
	pdf(sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.pdf',ctype),h=5,w=14.5)#h&w are in inches
	upset(fromList(DEGlist[[ctype]]), order.by = 'degree',keep.order=T)
	grid.text(ctype,x = 0.6, y = 0.97,gp = gpar(fontsize = 20))
	dev.off()
ctype= celltypes[2]
	pdf(sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.pdf',ctype),h=5,w=14.5)#h&w are in inches
	upset(fromList(DEGlist[[ctype]]), order.by = 'degree',keep.order=T)
	grid.text(ctype,x = 0.6, y = 0.97,gp = gpar(fontsize = 20))
	dev.off()
ctype = celltypes[3]
	pdf(sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.pdf',ctype),h=5,w=14.5)#h&w are in inches
	upset(fromList(DEGlist[[ctype]]), order.by = 'degree',keep.order=T)
	grid.text(ctype,x = 0.6, y = 0.97,gp = gpar(fontsize = 20))
	dev.off()
ctype= celltypes[4]
	pdf(sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.pdf',ctype),h=5,w=14.5)#h&w are in inches
	upset(fromList(DEGlist[[ctype]]), order.by = 'degree',keep.order=T)
	grid.text(ctype,x = 0.6, y = 0.97,gp = gpar(fontsize = 20))
	dev.off()
ctype = celltypes[5]
	pdf(sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.pdf',ctype),h=5,w=14.5)#h&w are in inches
	upset(fromList(DEGlist[[ctype]]), order.by = 'degree',keep.order=T)
	grid.text(ctype,x = 0.6, y = 0.97,gp = gpar(fontsize = 20))
	dev.off()
	#write.table(myDF[,-1],file=sprintf('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/%s_geneOverlaps.csv',pheno),sep=',',quote=F,row.names=F,col.names=T)
	
	
		
	
		

myTable = NULL
tmp = lapply(files,function(file){
	tab = read.table(file,sep='\t',head=T,stringsAsFactor=F)
	z = unlist(strsplit(file,'/'))
	condition = z[9]
	z = unlist(strsplit(z[10],'_|\\.'))
	cellState= paste(z[length(z)-2],z[length(z)-1],sep='_')
	nUp = dim(tab[tab$z.score > 0 & tab$BH< 0.05,])[1]
	nDown = dim(tab[tab$z.score < 0 & tab$BH< 0.05,])[1]
	nTested = dim(tab)[1]
	myTable <<- rbind(myTable,c(condition, cellState,nUp, nDown,nTested))
})
colnames(myTable) = c('Condition', 'CellState','#UpReg','#DownReg','#Tested')

write.table(as.data.frame(myTable),file= '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/num_DEGs_by_pheno.csv',sep=',',quote=F,row.names=F,col.names=T)