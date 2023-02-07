library(Seurat)

m = readRDS('~/SingleCellProjects/dataObjects/merged_mouseHuman.rds')
m = subset(m, subset = species == 'M.musculus')

df = data.frame(Mic4Sig = m$Mic4_Sig1, Clusters = m$seurat_clusters, stringsAsFactors = F)

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig3h.csv', sep = ', ', quote = F, row.names = T, col.names = T) 
