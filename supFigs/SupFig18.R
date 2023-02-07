#adapted from /home/brasel/Fenix/SingleCellProjects/Neurons/excit_inhib_gene_scores.R
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)

setwd('/home/brasel/scratch')

n = readRDS('/home/brasel/SingleCellProjects/dataObjects/neuron.rds')

layer.markers = c('SLC17A7','GAD2','LAMP5','CUX2','GLRA3','CARTPT','PRSS12','RORB','TPBG','TOX','FOXP2','ETV1','RPRM','RXFP1','TLE4','GRIK4','NTNG2','OPRK1','NR4A2','ADRA2A')
names(layer.markers) = c('n_ex','n_inh','l2_3','l2_3_4','l2','l2_3.2','l3_4','l4','l3_5_6','l5_6','l6','l5_6.2','l5_6.3','l5_6.4','l6.2','l4_6','l6.3','l6.4','l6_6b','l6b')
names(layer.markers) = sprintf('%s (%s)',names(layer.markers),layer.markers)

d = as.data.frame(t(log2(as.matrix(GetAssayData(n,slot='counts')[layer.markers,]+1))))
d$Sample = n$Sample_ID
d$Cluster = Idents(n)

z = aggregate(d,by=list(d$Sample,d$Cluster),FUN='mean')
z.1 = aggregate(z,by=list(z$Group.2),FUN='mean')
scores = z.1[,4:23]
colnames(scores) = toupper(names(layer.markers))
rownames(scores) = 0:6#c('cluster_0','cluster_1','cluster_2','cluster_3','cluster_4','cluster_5','cluster_6')
#scores = as.data.frame(scale(scores))
scores$clusters = 0:6


library(RColorBrewer)
library(tidyr)
mydata_long <- scores %>% gather("col", "val",-clusters)
mydata_long$clusters = as.factor(mydata_long$clusters)

col_lim <- max(mydata_long$val)  # make color scale symmetric

p <- mydata_long %>%
  # Convert rows and columns in factor to show them ordered
#  mutate(col = factor(col, levels = )) %>%
#  mutate(row = factor(row, levels = row_ord)) %>%

  ggplot(aes(x = clusters, y = col, fill = val)) +
  geom_tile() +
  # Optional: add black points on tiles with significant p-adj
  #geom_point(data = filter(mydata_long, padj < 0.05)) +
  # Color scale
  scale_fill_distiller(palette = "RdBu", limits = c(-col_lim, col_lim)) +
  # # Alternative color scale using log10-transformation, but showing actual values
  # scale_fill_distiller(
  #   palette = "RdBu", limits = c(-col_lim, col_lim),
  #   trans = "log10", labels = scales::comma
  # ) +
  labs(x = "Cluster", y = "Signature", fill = "Avg.\nExpression",
       title = "Cortical Layers") +
  theme(axis.text= element_text(color='black'),axis.title= element_text(color='black'),
        panel.grid = element_blank())
#plot(p)
ggsave(p,w=3,h=5,file='~/scratch/neuron_ex_inh_layers.pdf')
#ggsave(p,w=9*(length(unique(mydata_long$cluster))+3)/length(unique(mydata_long$col)),file='~/scratch/neuron_ex_inh_layers.pdf')

write.table(mydata_long, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig18.csv', sep = ', ', quote = F, row.names = F, col.names = T)
