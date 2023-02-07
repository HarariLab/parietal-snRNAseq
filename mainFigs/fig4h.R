library(tidyverse)
library(igraph)
#install.packages('ggraph')
library(ggraph)

dfc = read.csv('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/SCENIC/ColonnaNetwork_bh.csv',col.names = c('X','TF','Gene','Weight'))[,-1]
rownames(dfc) = paste(dfc$TF,dfc$Gene,sep = '_')
dfk = read.csv('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/SCENIC/KnightNetwork_bh.csv',col.names = c('X','TF','Gene','Weight'))[,-1]
rownames(dfk) = paste(dfk$TF,dfk$Gene,sep = '_')
o5DEG = read.table('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/oligo/bh_nebula_results/bh_nebula_compareSubclusters_oligo_5_all_05.30.21.txt')
ADgenes = read.table('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/AD_GWAS_hits_plot/data/2022.08.19_GWAS_loci_genes_inData.csv',sep=',',head=T)

library(RColorBrewer)
myEdgeCols = brewer.pal(n=11,name='Set1')

edges = rownames(dfk)
dfmerg = dfk
dfmerg$edgeColor = myEdgeCols[4]#'Knight'
dfmerg$edgeWidth = 2
tmp = lapply(rownames(dfc),function(r){
  if(r %in% edges){
    dfmerg[r,'edgeColor'] <<- myEdgeCols[3]#'Replicated'
    dfmerg[r,'edgeWidth'] <<- 8
  }else{
    new = c(dfc[r,],myEdgeCols[5],2)#'ROSMAP')
    names(new)[4] = 'edgeColor'
    names(new)[5] = 'edgeWidth'
    dfmerg <<- rbind(dfmerg,new)
  }
})

scale = 2
gr <- dfmerg %>% graph_from_data_frame(directed = T)
V(gr)$logFC = o5DEG[names(degree(gr)),'logFC_cluster']
#get node colors
myCols = colorRampPalette(brewer.pal(n = 11, name = 'RdBu'))
myCols = myCols(length(V(gr)$logFC)+100)
V(gr)$color = rev(myCols[rank(V(gr)$logFC)+50])

V(gr)$type = names(degree(gr))
V(gr)$type[V(gr)$type %in% c('SOX8','SREBF1','NFE2L2','NKX6-2','ZNF518A')] = 'TF'
V(gr)$type[V(gr)$type %in% df1$Gene] = 'Replicated'
V(gr)$type[!(V(gr)$type %in% c('TF','Replicated'))] = 'Target'
V(gr)$size = 2
V(gr)$size[V(gr)$type == 'TF'] = 15
V(gr)$size[V(gr)$type == 'Replicated'] = 12
V(gr)$label.cex = 0.7
V(gr)$label = names(degree(gr))
V(gr)$label[V(gr)$type == 'Target'] = ''
V(gr)$label.font = 2#bold
V(gr)$frame.color	 = 'black'
V(gr)$frame.color[names(degree(gr)) %in% ADgenes$Gene] = 'red'
#V(gr)$frame.color	[V(gr)$type == 'Replicated'] = 'gold'
V(gr)$color	 = 'white'
V(gr)$color	[V(gr)$type == 'TF'] = 'cyan'
V(gr)$color	[V(gr)$type == 'Replicated'] = 'gold'
E(gr)$color = dfmerg$edgeColor
E(gr)$width = dfmerg$edgeWidth/scale
E(gr)$arrow.size = dfmerg$edgeWidth/4/scale
plot(gr,vertex.label.family='Helvetica')

#ggraph(gr) + geom_edge_link() + geom_node_point() + geom_node_text(aes(label=name))
p = ggraph(gr) + geom_edge_link(aes(color = E(gr)$color)) + #scale_edge_color_gradient2(high='red',mid='gray25') +  #edge_width = log(Weight+1),edge_alpha = log(Weight+1), 
  geom_node_point(color = 'white', aes(size = size)) +# geom_node_point(aes(fill = type), size = 22) + 
 # geom_node_text(aes(label=name)) + 
  theme_graph()

pData = p$data

write.table(pData, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig4h.csv', sep = ', ', quote = F, row.names = F, col.names = T)

#pData$nodeSize[267:287] = c(1,3:14,16:23)
#pData$nodeSize = factor(pData$nodeSize,levels=c(1:24))
#pData$name[!(pData$type %in% c('TF','Replicated'))] = ''
p + geom_point(data=pData,aes(x,y,fill=logFC,color=type,size=size,stroke=sqrt(size)),shape=21) + scale_color_manual(values=c("TF" = 'green','Target'='black','Replicated' = 'gold')) + 
  scale_fill_gradient2(low=myCols[50],high=myCols[length(V(gr)$logFC)+50]) +
  geom_text(data=pData,aes(x,y,label=label))
