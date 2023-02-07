#------------------------------------------------------------------------------------------
#results = Hits #read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_pheno/EnricherConsensus_APOE_Neuron_cellStates.txt')
results = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly.csv',sep=',')
padj  = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/EnricherConsensus_cellStates_upOnly_padj.csv',sep=',')
results = results[grep('GO:',rownames(results)),]
rownames(results)[rownames(results) == "'de novo' posttranslational protein folding (GO:0051084)"] = "'de novo' post-translational protein folding (GO:0051084)"
padj = padj[grep('GO:',rownames(padj)),]
rownames(padj)[rownames(padj) == "'de novo' posttranslational protein folding (GO:0051084)"] = "'de novo' post-translational protein folding (GO:0051084)"
results[is.na(results)] = 0
padj[is.na(padj)] = 1
manualFilter = c('cytoplasmic translation (GO:0002181)','mRNA splicing, via spliceosome (GO:0000398)','chemical synaptic transmission (GO:0007268)','regulation of neuron projection development (GO:0010975)','regulation of hydrolase activity (GO:0051336)','response to calcium ion (GO:0051592)','establishment or maintenance of apical/basal cell polarity (GO:0035088)','protein localization to membrane (GO:0072657)','mitotic cell cycle phase transition (GO:0044772)','ciliary basal body-plasma membrane docking (GO:0097711)','mitochondrial respiratory chain complex assembly (GO:0033108)','positive regulation of telomerase RNA localization to Cajal body (GO:1904874)','cholesterol biosynthetic process (GO:0006695)','response to lipopolysaccharide (GO:0032496)','antigen receptor-mediated signaling pathway (GO:0050851)','nuclear migration (GO:0007097)','receptor-mediated endocytosis (GO:0006898)','extracellular matrix organization (GO:0030198)','positive regulation of cell differentiation (GO:0045597)','chondroitin sulfate proteoglycan biosynthetic process (GO:0050650)','regulation of cellular response to stress (GO:0080135)','regulation of apoptotic process (GO:0042981)','transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)','regulation of cell migration (GO:0030334)','phosphorylation (GO:0016310)','neutrophil degranulation (GO:0043312)','cytokine-mediated signaling pathway (GO:0019221)')
data_wide = data.frame(row=rownames(results),results,stringsAsFactors=F)
padj_wid = data.frame(row=rownames(padj),padj,stringsAsFactors=F)
mydata_long <- data_wide %>%	gather("col", "val", -row)
padj_long = padj_wid %>% gather('col','padj',-row)
mydata_long$padj = padj_long$padj
mydata_long = mydata_long[mydata_long$row %in% manualFilter,]
#sanity check
manualFilter %in% unique(mydata_long$row)
#replace ',' with ';' in mydata_long$row
mydata_long$row = gsub(',',';',mydata_long$row)

write.table(mydata_long, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig3g.csv', sep = ', ', quote = F, row.names = F, col.names = T) # nolint: line_length_linter.