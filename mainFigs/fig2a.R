library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(stringr)



files = system('ls /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/*/*/nebula*all*', intern = T)

cellTypes = c('astro', 'micro', 'n_ex', 'n_inh', 'oligo', 'opc')#c('astro', 'micro', 'n_inh', 'n_ex')#
pheno = c('sADvsCO', 'TREM2vsCO', 'ADADvsCO')
#cellType = 'astro'
chunks = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/dataObjects/Estimate_chunks.tsv', sep = '\t', stringsAsFactors = F, head = T)
myFactors = c('NS', 'Nominal', 'BH', levels(as.factor(chunks$Group)))
myFactors = factor(myFactors, levels = c("all_up", "all_down", "ADAD_up", "ADAD_down", "TREM2_up", "TREM2_down", 'TREM2_ADAD_down', "TREM2_up_ADAD_down", 'sAD_ADAD_up', "sAD_TREM2_up", "NS", "Nominal", "BH"))

sigHits_up = list()
chunkList <<-  list()
distributionPlots = list()
nGenes = 500 #500 or 25
df2 = NULL
tmp = lapply(cellTypes, function(cellType){
	f = files[grep(cellType, files)]

	myTable = list()
	tabs = lapply(pheno, function(p){
		file = f[grep(p, f)]
		tab = read.table(file, sep = '\t', head = T, stringsAsFactor = F)
		z = unlist(strsplit(file, '/'))
		condition = z[9]
		z = unlist(strsplit(z[10], '_|\\.'))
		cellState =  paste(z[length(z)-2], z[length(z)-1], sep = '_')
		myTable[[paste(condition, cellState, sep = '_')]] <<- rownames(tab[tab$BH < 0.05, ])
		return(tab)
	})

	g = union(union(myTable[[1]], myTable[[2]]), myTable[[3]])

	chart = data.frame(Genes = g, stringsAsFactors = F)
	pval = data.frame(Genes = g, stringsAsFactors = F)
	BH = data.frame(Genes = g, stringsAsFactors = F)
	tmp = lapply(tabs, function(tab){
		chart <<- cbind(chart, tab[g, 1])
		pval <<- cbind(pval, tab[g, 'p.value'])
		BH   <<- cbind(BH, tab[g, 'BH'])
	})

	colnames(chart) = c('Gene',  'sAD', 'TREM2', 'ADAD')
	colnames(pval) = c('Gene',  'sAD', 'TREM2', 'ADAD')
	rownames(pval) = pval$Gene
	pval_all = pval
	colnames(BH) = c('Gene',  'sAD', 'TREM2', 'ADAD')
	rownames(BH) = pval$Gene
	chart = na.omit(chart)
	chart_all = chart 
	rownames(chart_all) = chart_all[, 1] 

	df = tidyr::gather(chart_all[, -1], 'Pheno', 'Estimate', 1:3)
	df$Pheno = factor(df$Pheno, levels = c('sAD', 'TREM2', 'ADAD'))
	df2 = tidyr::gather(pval_all[rownames(chart_all), -1], 'Pheno', 'p_val', 1:3)
	df$p_val = df2$p_val
	myC = brewer.pal(3,  'Accent')#("white", 'purple', 'green', "red")
	p3 = ggplot(df[df$p_val < 0.05, ], aes(x = Estimate, y = Pheno, fill = Pheno)) + geom_density_ridges() + scale_fill_brewer(palette = 'Pastel1') + xlim(c(-3, 3)) + ylab(str_to_title(cellType)) + theme_bw()+ theme(axis.text = element_text(size = 16, color = 'black'), axis.title = element_text(size = 20)) + rremove('legend')
	distributionPlots[[length(distributionPlots)+1]] <<- p3
	cat(print(table(df[df$p_val < 0.05, 'Pheno'])))
	#ggsave(p3, file = sprintf('~/scratch/Estimate_all_Distributions_%s.pdf', cellType))
})

#q = ggarrange(plotlist = distributionPlots, ncol = 3, nrow = 2)
#ggsave(q, units = 'mm', h = 100, w = 400, file = '~/scratch/DEG_estimate_distributions_byCellType.pdf')


df = NULL
tmp = lapply(seq_along(distributionPlots), function(ind){
	distributionPlots[[ind]]$data$CellType = cellTypes[ind]
	df <<- rbind(df, distributionPlots[[ind]]$data)
})

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig2a.csv', sep = ', ', quote = F, row.names = F, col.names = T) # nolint: line_length_linter.
