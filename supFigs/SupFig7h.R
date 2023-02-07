library(ggplot2)
library(dplyr)
library(ggrepel)

t = read.table('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/micro/bh_nebula_results/nebula_compareSubclusters_microglia_4_all_06.01.21.txt')
celltype = 'Micro_4vAll'

log2FC = log2(exp(t$logFC_cluster))
color = rep('grey',length(log2FC))
color[(log2FC > 0) & (t$BH < 0.05)] = '#FFAAA0'
color[(log2FC < 0) & (t$BH < 0.05)] = '#A0C5FF'
df = data.frame(rownames(t),t,log2FC,color)
colnames(df)[1] = 'Genes'

df = df[abs(df$log2FC) < 15,]

label = NULL
df$label = 0
df$label[df$log2FC > 0] = -log10(df$p.value[df$log2FC > 0])/max(-log10(df$p.value[df$log2FC > 0])) + abs(df$log2FC[df$log2FC > 0])/max(abs(df$log2FC[df$log2FC > 0]))
label = c(label,rownames(df[order(df$label,decreasing=T),])[1:7])
df$label = 0
df$label[df$log2FC < 0] = -log10(df$p.value[df$log2FC < 0])/max(-log10(df$p.value[df$log2FC < 0])) + abs(df$log2FC[df$log2FC < 0])/max(abs(df$log2FC[df$log2FC < 0]))
label = c(label,rownames(df[order(df$label,decreasing=T),])[1:7])

p = ggplot(df[abs(df$log2FC) < 15,], aes(x = log2FC, y = -log10(p.value),color = color)) + #ifelse(abs(log2FC)>0.6 & BH < 0.05,"red","grey")
    geom_point() +
    xlab(expression("Fold Change, Log"[2]*"")) + xlim(c(-8,8)) + 
    ylab(expression("P value, Log"[10]*"")) +
		ggtitle('Mic-stress') +
    theme_bw() +
	  theme(legend.position = "none",axis.title=element_text(size=18),axis.text=element_text(size=16,color='black'),title=element_text(size=18)) +
    theme(legend.position = "none")+
    scale_colour_manual(values = names(table(color))) +
    geom_text_repel(data=df[c(label,'MECP2'),],
                    aes(log2FC, -log10(p.value), label = Genes),size = 5, color="black",box.padding = 0.5,min.segment.length=0.6,max.iter = 8000, max.overlaps =20)

ggsave(p,file=paste0('~/scratch/',celltype,'_volcano.pdf'),h=4,w=4)

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig7h.csv', sep = ', ', quote = F, row.names = F, col.names = T)
