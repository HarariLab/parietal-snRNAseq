library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggridges)
library("grid")
library("ggdendro")
library(heatmaply)
library(plotly)
library(enrichR)
library(stringr)



files = system('ls /home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/*/*/nebula*all*',intern=T)

cellTypes = c('astro','micro','n_ex','n_inh','oligo','opc')#c('astro','micro','n_inh','n_ex')#
pheno = c('sADvsCO','TREM2vsCO','ADADvsCO')
#cellType = 'astro'
chunks = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/dataObjects/Estimate_chunks.tsv',sep='\t',stringsAsFactors=F,head=T)
myFactors = c('NS','Nominal','BH',levels(as.factor(chunks$Group)))
myFactors = factor(myFactors,levels=c("all_up","all_down","ADAD_up","ADAD_down","TREM2_up","TREM2_down",'TREM2_ADAD_down',"TREM2_up_ADAD_down",'sAD_ADAD_up',"sAD_TREM2_up","NS","Nominal","BH"))

sigHits_up = list()
chunkList <<-  list()
distributionPlots = list()
nGenes = 25 #500 or 25
df2 = NULL
myFullDF = NULL
tmp = lapply(cellTypes,function(cellType){
	f = files[grep(cellType,files)]

	myTable = list()
	tabs = lapply(pheno,function(p){
		file = f[grep(p,f)]
		tab = read.table(file,sep='\t',head=T,stringsAsFactor=F)
		z = unlist(strsplit(file,'/'))
		condition = z[9]
		z = unlist(strsplit(z[10],'_|\\.'))
		cellState= paste(z[length(z)-2],z[length(z)-1],sep='_')
		myTable[[paste(condition,cellState,sep='_')]] <<- rownames(tab[tab$BH < 0.05,])
		return(tab)
	})

	g = union(union(myTable[[1]],myTable[[2]]),myTable[[3]])

	chart = data.frame(Genes = g,stringsAsFactors=F)
	pval = data.frame(Genes = g,stringsAsFactors=F)
	BH = data.frame(Genes = g,stringsAsFactors=F)
	tmp = lapply(tabs,function(tab){
		chart <<- cbind(chart,tab[g,1])
		pval <<- cbind(pval,tab[g,'p.value'])
		BH   <<- cbind(BH,tab[g,'BH'])
	})
	
	colnames(chart) = c('Gene', 'sAD','TREM2','ADAD')
	colnames(pval) = c('Gene', 'sAD','TREM2','ADAD')
	rownames(pval) = pval$Gene
	pval_all = pval
	colnames(BH) = c('Gene', 'sAD','TREM2','ADAD')
	rownames(BH) = pval$Gene
	chart = na.omit(chart)
	chart_all = chart ######################################### cut to distribution code at bottom ####################################3
	
	#print out average estimate size split by pos and neg
	cat(sprintf('%s\n',cellType))
	tmp = lapply(colnames(chart)[-1],function(x) {
		cat(sprintf('%s\t#<-1.5: %s\t#>1.5: %s\n',x,sum(chart[,x] < -1.5),sum(chart[,x] > 1.5)))
	})
	#filter to top 500 estimates
	maxEst = unlist(apply(chart[,-1],1,function(x) max(abs(x)) ))
	chart = chart[order(maxEst,decreasing=T),]
	chart = chart[1:min(c(nrow(chart),500)),] #only 50 for main figure, 500 for supplement; na.omit just incase there are less than 500 BH genes
	#chart = chart[maxEst >= 1,]
	pval = pval[chart$Gene,-1]
	BH = BH[chart$Gene,-1]

	pval = pval < 0.05
	pval = apply(pval,2,as.numeric)
	BH = BH < 0.05
	BH = apply(BH,2,function(col) as.numeric(col)*2 )
	Sig = data.frame(sAD_sig=apply(cbind(pval[,1],BH[,1]),1,max), TREM2_sig=apply(cbind(pval[,2],BH[,2]),1,max), ADAD_sig=apply(cbind(pval[,3],BH[,3]),1,max))
	Sig[Sig == 0] = 'NS'	
	Sig[Sig == 1] = 'Nominal'
	Sig[Sig ==2] = 'BH'
	Sig = apply(Sig,2,function(col) factor(col,levels=myFactors) )

	#### heatmaply method ####
	rownames(chart) = chart[,1]
	chart = chart[,-1]
	df = cbind(chart,Sig)
	p= heatmaply(df,k_row=8,Colv=NA,scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0),
		plot_method='ggplot',fontsize_row=10,subplot_widths=c(0.4,0.2,0.4),colorbar_thickness=100,side_color_colorbar_len=0.1,col_side_palette='RdYlBu')
	geneOrder = rev(unlist(lapply(p$x$data[[1]]$text,function(x) unlist(strsplit(x,' |<'))[2] )))[1:min(c(nrow(chart),500))]
	#lncRNA = sum(grepl('\\.|LINC|-',geneOrder)); 	#cat(cellType); 	#cat('\n');	#cat(sprintf('%s/%s = %s\n',lncRNA,500,lncRNA/500)); 	#lncRNA = sum(grepl('\\.|LINC|-',chart_all$Gene)); 	#cat(sprintf('%s/%s = %s\n',lncRNA,length(chart_all$Gene),lncRNA/length(chart_all$Gene)))
	chunkInd = c(1,which(geneOrder %in% chunks[chunks$CellType == cellType,'Gene']))
	chunkList <<-  list()
	delete = lapply(1:(length(chunkInd)-1), function(ind){
		condition = chunks[chunks$CellType == cellType,'Group'][ind]
		chunkList[[condition]] <<- c(chunkList[[condition]],geneOrder[chunkInd[ind]:chunkInd[ind+1]])
	})
	#modules
	Modules = rep('Error',nrow(chart))
	chunkInd[1] = 0 
	delete = lapply(1:(length(chunkInd)-1), function(ind){
		Modules[(chunkInd[ind]+1):chunkInd[ind+1]] <<- chunks[chunks$CellType == cellType,'Group'][ind]
	})
	df2 <<- cbind(df[geneOrder,],Modules)
	tmp = lapply(4:7,function(ind) df2[,ind] <<- factor(df2[,ind],levels=myFactors) )
	width = 3333/2
	height = 5000
	rowv = T
	subW = c(0.4,0.2,0.4)
	if (nGenes != 500){
		mainFigKeep = unlist(lapply(unique(df2$Modules),function(mod) {
			tmp = df2[df2$Modules == mod,]
			tmp = apply(tmp[,1:3],1,sum)
			tmp = tmp[order(tmp,decreasing=T)]
			return(names(tmp)[1:round(length(tmp)/(500/nGenes))])
		}))
		df2 <<- df2[rownames(df2) %in% mainFigKeep,]
		height = height/(500/nGenes*0.7)
		width = width*0.3
		rowv = F
		subW = c(0.67,0.33)
	}
	#p= heatmaply(df2,k_row=8,Colv=NA,Rowv=rowv, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red", midpoint = 0),
			#plot_method='ggplot',fontsize_row=10,subplot_widths=subW,colorbar_thickness=100,side_color_colorbar_len=0.1,col_side_palette='Paired')
	#plotly::orca(p,file = sprintf('Estimates_heatmap_dendrogram_sigHits_%s_%sgenes.pdf',cellType,nGenes),width=width, height=height,more_args=c('-d','/home/brasel/scratch/'))
    df2$cellType = cellType
    myFullDF <<- rbind(myFullDF, df2)
})

write.table(myFullDF, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig2b.csv', sep = ', ', quote = F, row.names = T, col.names = T) # nolint: line_length_linter.
