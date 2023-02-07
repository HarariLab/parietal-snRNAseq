#------------------------------------------------------------------------------------------
# -log10(pval) matrix
library(openxlsx)
library(enrichR)

outfile = '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly.csv'

Hits = data.frame(Term = NA)
degs = read.xlsx('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/Sup Table 4 - DEG_by_cellState_summary.xlsx')
pages = grep('all',degs[[1]][-1]) + 1
setName = degs[[1]][pages[1]]
tmp = lapply(pages,function(p){
	d = read.xlsx('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/Sup Table 4 - DEG_by_cellState_summary.xlsx',sheet=p)
	d = d[d$logFC_cluster > 0,]
	if(nrow(d) >= 2){
		#write(paste(sheets[p],paste(paste(d$X1,signif(d$logFC_nAPOE34,3),sep=','),collapse='\t'),sep='\t\t'),file=outfile,append=T)
		r = enrichr(d$X1,c('KEGG_2021_Human','GO_Biological_Process_2021'))
		r = rbind(r$KEGG_2021_Human, r$GO_Biological_Process_2021)
		r[,degs[[1]][p]] = -log10(r$P.value)
		Hits <<- merge(Hits, r[,c('Term',degs[[1]][p])], by='Term', all=T)
	}
})
Hits = Hits[!is.na(Hits$Term),]
rownames(Hits) = Hits$Term
Hits = Hits[,-1]
write.table(Hits,outfile,sep=',')

# padj matrix --------------------------------------------------------------------------------------------
library(openxlsx)
library(enrichR)

outfile = '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly_padj.csv'

Hits = data.frame(Term = NA)
degs = read.xlsx('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/Sup Table 4 - DEG_by_cellState_summary.xlsx')
pages = grep('all',degs[[1]][-1]) + 1
#setName = degs[[1]][pages[1]]
tmp = lapply(pages,function(p){
	d = read.xlsx('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/Sup Table 4 - DEG_by_cellState_summary.xlsx',sheet=p)
	d = d[d$logFC_cluster > 0,]
	if(nrow(d) >= 2){
		#write(paste(sheets[p],paste(paste(d$X1,signif(d$logFC_nAPOE34,3),sep=','),collapse='\t'),sep='\t\t'),file=outfile,append=T)
		r = enrichr(d$X1,c('KEGG_2021_Human','GO_Biological_Process_2021'))
		r = rbind(r$KEGG_2021_Human, r$GO_Biological_Process_2021)
		r[,degs[[1]][p]] = r$Adjusted.P.value
		Hits <<- merge(Hits, r[,c('Term',degs[[1]][p])], by='Term', all=T)
	}
})
Hits = Hits[!is.na(Hits$Term),]
rownames(Hits) = Hits$Term
Hits = Hits[,-1]
write.table(Hits,outfile,sep=',')

#------------------------------------------------------------------------------------------
#must open up rrvgo mamba environment
#results = Hits #read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/EnricherConsensus_APOE_Neuron_cellStates.txt')
results = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly.csv',sep=',')
padj  = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly_padj.csv',sep=',')
results = results[grep('GO:',rownames(results)),]
padj = padj[grep('GO:',rownames(padj)),]
results[is.na(results)] = 0
padj[is.na(padj)] = 1
#run using rrvgo instead of revigo.irb.hr -----------------------------------------------------------------------
library(rrvgo)#v1.6.0
SetNames = unique(unlist(lapply(1:ncol(results),function(col){
	return(rownames(results[order(results[,col],decreasing=T),])[1:10])
})))
					#SetNames = unique(c(SetNames,'cytokine-mediated signaling pathway (GO:0019221)','cytoplasmic translation (GO:0002181)','negative regulation of programmed cell death (GO:0043069)','MAPK signaling pathway','NOD-like receptor signaling pathway','Necroptosis','Estrogen signaling pathway','mRNA splicing, via spliceosome (GO:0000398)'))
simM = calculateSimMatrix(sub('\\)','',sub(".*[\\(]","",SetNames)),orgdb="org.Hs.eg.db", ont="BP", method="Rel") 
score = unlist(apply(results[SetNames,],1,mean))
names(score) = sub('\\)','',sub(".*[\\(]","",SetNames))
reducedTerms = reduceSimMatrix(simM, score, threshold=0.6,orgdb="org.Hs.eg.db") 
GOs = unique(reducedTerms[,'parent'])
terms = unique(reducedTerms[,'parentTerm'])
#re-rank and order the terms now that some have been summarized out
rankR2 = apply(results[paste(terms,'('),],2,scale)
rownames(rankR2) = rownames(results[paste(terms,'('),])
rsum2 = rowMeans(rankR2)
rsum2 = rsum2[order(rsum2,decreasing=T)]
SetNames = names(rsum2)#[1:25] #manually did extra pruning from this list to get the list below (manualFilter)

filt = results[SetNames,]
colnames(filt) = substr(colnames(filt),0,nchar(colnames(filt))-4)
filt_padj = padj[SetNames,]
colnames(filt_padj ) = substr(colnames(filt_padj ),0,nchar(colnames(filt_padj ))-4)
data_wide = data.frame(row= rownames(filt),filt)
padj_wid = data.frame(row=rownames(filt_padj),filt_padj)
dm = data.matrix(filt)



# plot the heatmap #######################################################
# Input: a melted data.frame with columns "row", "column", "val", and "padj"
library(tidyr)
library(dplyr)
library(ggplot2)

# Melt my data keeping the "padj" column unmelted
#not for full figure: manualFilter = c('cytoplasmic translation (GO:0002181)','mRNA splicing, via spliceosome (GO:0000398)','chemical synaptic transmission (GO:0007268)','regulation of neuron projection development (GO:0010975)','regulation of hydrolase activity (GO:0051336)','response to calcium ion (GO:0051592)','establishment or maintenance of apical/basal cell polarity (GO:0035088)','protein localization to membrane (GO:0072657)','mitotic cell cycle phase transition (GO:0044772)','ciliary basal body-plasma membrane docking (GO:0097711)','mitochondrial respiratory chain complex assembly (GO:0033108)','positive regulation of telomerase RNA localization to Cajal body (GO:1904874)','cholesterol biosynthetic process (GO:0006695)','response to lipopolysaccharide (GO:0032496)','antigen receptor-mediated signaling pathway (GO:0050851)','nuclear migration (GO:0007097)','receptor-mediated endocytosis (GO:0006898)','extracellular matrix organization (GO:0030198)','positive regulation of cell differentiation (GO:0045597)','chondroitin sulfate proteoglycan biosynthetic process (GO:0050650)','regulation of cellular response to stress (GO:0080135)','regulation of apoptotic process (GO:0042981)','transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)','regulation of cell migration (GO:0030334)','phosphorylation (GO:0016310)','neutrophil degranulation (GO:0043312)','cytokine-mediated signaling pathway (GO:0019221)')
mydata_long <- data_wide %>%	gather("col", "val", -row)
padj_long = padj_wid %>% gather('col','padj',-row)
mydata_long$padj = padj_long$padj
mydata_long = mydata_long[mydata_long$row %in% manualFilter,]

# Convert data in to a numeric matrix for clustering
#dm <- mydata_long %>%
#  select(row, col, val) %>%
  # optional case_when to deal with infinite values
  #mutate(val = case_when(
    #is.infinite(val) & val < 0 ~ min(val[!is.infinite(val)]),  # Make +Inf values = max observed
    #is.infinite(val) & val > 0 ~ max(val[!is.infinite(val)]),  # Make -Inf values = min observed
    #T ~ val)  # Keep the others as they are (otherwise they become NA with this function
  #) %>%
  #select(-padj) %>%  # remove the padj column to only have three columns
 # spread(col, val, fill = 0) %>%  # make data wide
  #tibble::column_to_rownames("col") %>%  # Need rownames so that there are no non-numeric columns
 # data.matrix()  # convert to numeric matrix
#tmp = lapply(c("euclidean","maximum","manhattan","canberra",'p','s'),function(colMeth) {
	#tmp = lapply(c("euclidean","maximum","manhattan","canberra",'p','s'),function(rowMeth) {
colMeth = 'p'
rowMeth = 'euclidean'
		if(colMeth %in% c("euclidean","maximum","manhattan","canberra")){
			# Cluster columns using some distance calculation
			col_ord <- dm %>% t() %>%
			  dist(method = colMeth) %>%  # look at ?dist for options ‘"euclidean"’, ‘"maximum"’, ‘"manhattan"’, ‘"canberra",‘"binary"’
			  hclust()
			p1 = plot(col_ord, hang = -1)  # generate dendrogram
			col_ord <- col_ord$labels[col_ord$order]  # extract ordered col labels
		} else{
			# Alternative: cluster columns using correlation
			col_ord <- dm %>%
			 # t() %>%  # correlation matrix transposes data, so we need to pre-transpose
			  cor(method = "p") %>%  # Spearman correlation matrix
			  dist(method = "canberra") %>%  # default euclidean distance
			  hclust()
			hclust(dist(dm))
			p1 = plot(col_ord, hang = -1)  # generate dendrogram
			col_ord <- col_ord$labels[col_ord$order]  # extract ordered col labels
		}
		if(rowMeth %in% c("euclidean","maximum","manhattan","canberra")){
			## Cluster rows
			row_ord <- hclust(dist(dm, method = "manhattan"))  # one-liner doing the same thing - note transposition to cluster by rows
			p2 = plot(row_ord, hang = -1)  # generate dendrogram
			row_ord <- row_ord$labels[row_ord$order]  # extract ordered row labels
		} else{
			# Alternative: cluster rows using correlation
			row_ord <- dm %>%  t() %>% cor(method = "s") %>%  dist(method = "canberra") %>%  hclust()
			hclust(dist(dm))
			p1 = plot(row_ord, hang = -1)  # generate dendrogram
			row_ord <- row_ord$labels[row_ord$order]  # extract ordered col labels
		}
		col_lim <- max(mydata_long$val)  # make color scale symmetric
		row_ord = row_ord[row_ord %in% manualFilter] # for man fig#############################################################################
		p <- mydata_long %>%
		  # Convert rows and columns in factor to show them ordered
		  mutate(col = factor(col, levels = col_ord)) %>%
		  mutate(row = factor(row, levels = row_ord)) %>%
		
		  ggplot(aes(x = col, y = row, fill = val)) +
		  geom_tile() +
		  # Optional: add black points on tiles with significant p-adj
		  geom_point(data = filter(mydata_long, padj < 0.05)) +
			#geom_point(data = filter(mydata_long, val > -log10(0.05))) +
		  # Color scale
		  #scale_fill_distiller(palette = "RdBu", limits = c((-log10(0.05)*2-col_lim), col_lim)) +
			scale_fill_gradientn(colors=c('#A0C5FF','white','orange','red'), values=scales::rescale(c(0,-log10(0.05),5,col_lim)),guide='colorbar',limits=c(0,col_lim)) +
		  # # Alternative color scale using log10-transformation, but showing actual values
		  #scale_fill_distiller(
		  #   palette = "RdBu", limits = c((-log10(0.05)*2-col_lim), col_lim),
		  #   trans = "log10", labels = scales::comma
		  # ) +
		  labs(x = "Cell State", y = "Biological Process Terms", fill = expression("P value (-Log"[10]*")"),
		       title = "My heatmap") +
		  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
		        panel.grid = element_blank())
		#plot(p)
		#ggsave(p,file=sprintf('~/scratch/GO_KEGG_heatmap_cellStates_col_%s_row_%s.pdf',colMeth,rowMeth),h=17,w=15) #full fig
		ggsave(p,file='~/scratch/GO_KEGG_heatmap_cellStates_v6.pdf',h=183,w=300,units='mm') #main Fig
#})})
#########################################################################
breaks = seq(0,max(mat),length.out=1000)
gradient1 = colorpanel(300,'blue','white' )[(301-sum(breaks[-1]<=log(-log10(0.05)))):300]#sum( breaks[-1]<=-log10(0.05) ), "blue", "white" )
gradient2 = colorpanel( sum( breaks[-1]>log(-log10(0.05)) ), "lightyellow", 'orange',"red" )
hm.colors = c(gradient1,gradient2)

pdf('~/scratch/GO_KEGG_heatmap_cellStates_v4.pdf',h=17,w=15)
heatmap.2(mat,margins=c(10, 27),trace='none',col=hm.colors,breaks=breaks,
	#distfun = function(x) as.dist(1-cor(t(x))),
  #hclustfun = function(x) hclust(x, method="average")
	distfun = function(x) dist(x, method="euclidean"),
  hclustfun = function(x) hclust(x, method="ward.D2")
)#rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)))  #scale='column'
dev.off()

#replace ',' with ';' in mydata_long$row
mydata_long$row = gsub(',',';',mydata_long$row)
write.table(mydata_long, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig6.csv', sep = ', ', quote = F, row.names = F, col.names = T)
