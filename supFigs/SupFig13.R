myPath='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/AD_GWAS_hits_plot/'

#files to manually change
genesToPlot=            read.csv(paste0(myPath,'data/2022.08.19_GWAS_loci_genes_inData.csv') )
largestEstimateFile=    read.csv(paste0(myPath,'all_ADgenes_merged_2022.08.22.csv'))
meanExpressionFile=     read.csv(paste0(myPath,'all_avgExprByCellType_2022.08.22.csv'))

#collapse the largestEstimateFile into tidy format using gather function from tidyr
library(tidyr)
est=gather(largestEstimateFile,cellType,estimate,-Row.names)

meanExpressionFile$Row.names = rownames(meanExpressionFile)
colnames(meanExpressionFile) = colnames(largestEstimateFile)[c(2:ncol(largestEstimateFile),1)]
exp=gather(meanExpressionFile,cellType,expression,-Row.names)

#set Gene column as rownames for genesToPlot
rownames(genesToPlot) = genesToPlot$Gene
#merge est and exp by Row.names and celltype columns keeping all rows from both
estExp = merge(est,exp,by=c('Row.names','cellType'),all=TRUE)

#add Locus column to estExp by passing in the Locus column from genesToPlot
estExp$Locus = genesToPlot[estExp$Row.names,'Locus']

#change Row.names to Gene
colnames(estExp)[1] = 'Gene'
#move Locus column to first column
estExp = estExp[,c(ncol(estExp),1:(ncol(estExp)-1))]
#reorder by Locus and then Gene
estExp = estExp[order(estExp$Locus,estExp$Gene),]

write.table(estExp, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig13.csv', sep = ', ', quote = F, row.names = F, col.names = T)
