orig_plot =read.csv('~/SingleCellProjects/MyProjects/67BrainsPaper/AD_GWAS_hits_plot/mergedPriority_avgExprByCellType_9.24.21.csv')

swarup_avgExpr =read.csv(file = '/home/brasel/SingleCellProjects/MyProjects/SwarupData/GWASgenes_avgExr.csv')
swarup_maxEstTable = read.csv(file = '/home/brasel/SingleCellProjects/MyProjects/SwarupData/GWASgenes_maxEst_p_0.05.csv',row.names=1)

washu_avgExpr =read.csv('~/SingleCellProjects/MyProjects/67BrainsPaper/AD_GWAS_hits_plot/mergedPriority_avgExprByCellType_22.03.08.csv')
washu_maxEstTable = read.csv('~/SingleCellProjects/MyProjects/67BrainsPaper/AD_GWAS_hits_plot/mergedPriority_ADgenes_merged_22.03.08.csv',row.names=1)

swarup_avgExpr = swarup_avgExpr[rownames(orig_plot ),]
swarup_maxEstTable = swarup_maxEstTable[rownames(orig_plot ),]
washu_avgExpr = washu_avgExpr[rownames(orig_plot ),]
washu_maxEstTable = washu_maxEstTable[rownames(orig_plot ),]

#add Row.names column to swarup_avgExpr, swarup_maxEstTable, washu_avgExpr, washu_maxEstTable
swarup_avgExpr$Row.names = rownames(swarup_avgExpr)
swarup_maxEstTable$Row.names = rownames(swarup_maxEstTable)
washu_avgExpr$Row.names = rownames(washu_avgExpr)
washu_maxEstTable$Row.names = rownames(washu_maxEstTable)

library(tidyr)
sexp = gather(swarup_avgExpr,cellType,expression,-Row.names)
sest = gather(swarup_maxEstTable,cellType,estimate,-Row.names)
wexp = gather(washu_avgExpr,cellType,expression,-Row.names)
west = gather(washu_maxEstTable,cellType,estimate,-Row.names)

sexp$estimate = sest$estimate
sexp$Study = 'Swarup'

wexp$estimate = west$estimate
wexp$Study = 'WashU'

#change the values in sexp$cellType to match those in wexp$cellType
sexp$cellType = gsub('ASC','astro',sexp$cellType)
sexp$cellType = gsub('OPC','opc',sexp$cellType)
sexp$cellType = gsub('EX','n_ex',sexp$cellType)
sexp$cellType = gsub('INH','n_inh',sexp$cellType)
sexp$cellType = gsub('MG','microglia',sexp$cellType)
sexp$cellType = gsub('ODC','oligo',sexp$cellType)

all = rbind(sexp,wexp)

write.table(all, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig14.csv', sep = ', ', quote = F, row.names = F, col.names = T)
