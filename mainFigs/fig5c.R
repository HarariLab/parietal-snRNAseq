library(ggplot2)
library(dplyr)
library(ggrepel)

pvalue.extreme <- function(z){#,len) {
   log.pvalue <-  log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)# + log(len)
   log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
	 return(-log10.pvalue)
}


t = read.table('~/SingleCellProjects/MyProjects/67BrainsPaper/DE_by_cluster/micro/bh_nebula_results/nebula_compareSubclusters_microglia_1_3_06.01.21.txt')
celltype = 'Micro_1v3'

log2FC = log2(exp(t$logFC_cluster))
color = rep('grey',length(log2FC))
color[(log2FC > 0) & (t$BH < 0.05)] = '#A0C5FF'
color[(log2FC < 0) & (t$BH < 0.05)] = '#b572d4'
df = data.frame(rownames(t),t,log2FC,color)
colnames(df)[1] = 'Genes'

df = df[abs(df$log2FC) < 15,]
df$p.value[df$p.value == 0] = min(df$p.value[df$p.value > 0])
df$log10P = unlist(lapply(df$z.score,pvalue.extreme))
df$log10P[df$log10P > 500] = 500

label = NULL
df$label = 0
df$label[df$log2FC > 0 & df$BH < 0.05] = -log10(df$p.value[df$log2FC > 0 & df$BH < 0.05])/max(-log10(df$p.value[df$log2FC > 0 & df$BH < 0.05])) + abs(df$log2FC[df$log2FC > 0 & df$BH < 0.05])/max(abs(df$log2FC[df$log2FC > 0 & df$BH < 0.05]))
label = c(label,rownames(df[order(df$label,decreasing=T),])[1:10])
df$label = 0
df$label[df$log2FC < 0 & df$BH < 0.05] = -log10(df$p.value[df$log2FC < 0 & df$BH < 0.05])/max(-log10(df$p.value[df$log2FC < 0 & df$BH < 0.05])) + abs(df$log2FC[df$log2FC < 0 & df$BH < 0.05])/max(abs(df$log2FC[df$log2FC < 0 & df$BH < 0.05]))
label = c(label,rownames(df[order(df$label,decreasing=T),])[1:7])

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig5c.csv', sep = ', ', quote = F, row.names = F, col.names = T)

p = ggplot(df[abs(df$log2FC) < 15,], aes(x = log2FC, y = log10P,color = color)) + #ifelse(abs(log2FC)>0.6 & BH < 0.05,"red","grey")
    geom_point() +
    xlab(expression("Fold Change, Log"[2]*"")) +
    ylab(expression("P value, Log"[10]*"")) +
		ggtitle('Mic-proinflammatory\tMic-activated') +
    theme_bw() +
    theme(legend.position = "none",axis.title=element_text(size=20),axis.text=element_text(size=19,color='black'),title=element_text(size=20)) +
    scale_colour_manual(values = names(table(color))) +
    geom_text_repel(data=df[unique(c(label,'C5','C3','ACVR1','TGFBR1','TGFBR2','BMPR2')),],
                    aes(log2FC, log10P, label = Genes,segment.inflect=T),size = 7, color="black",box.padding = 0.5,min.segment.length=0.5,max.iter = 32000,
										#force = 1,#nudge_x = 0.25,
								    #direction         = "y",
								    #hjust             = 1,
								   # segment.size      = 0.2,
								   # xlim=c(1.5,NA))
										)

ggsave(p,file=paste0('~/scratch/',celltype,'_volcano.pdf'))