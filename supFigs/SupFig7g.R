library(Seurat)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(stringr)
options(stringsAsFactors=F)

m = readRDS('~/SingleCellProjects/dataObjects/microglia.rds')

tab = table(m$Sample_ID,Idents(m))
tab = tab[rowSums(tab) > 60,]
tab = data.frame(Sample = rownames(tab), Status = unlist(lapply(rownames(tab),function(sID) m$Status[m$Sample_ID == sID][1] )), as.data.frame.array(tab),check.names=F)
tab$Status[tab$Status == 'Neuro_AD'] = 'sAD'
tab$Status[tab$Status == 'Neuro_ADAD'] = 'ADAD'
tab$Status[tab$Status == 'Neuro_CO'] = 'CO'
tab$Status[tab$Status == 'Neuro_OT'] = 'OTH'
tab$Status[tab$Status == 'Neuro_Presympt'] = 'Pres'
tab = gather(tab,'CellState','NucleiCounts',-Sample, -Status)
tab$Status = factor(tab$Status,levels=c('ADAD','sAD','Pres','CO','OTH'))
tab$Status_Sample = paste(tab$Status,tab$Sample,sep='_')

p = ggplot(tab,aes(x=Status_Sample,y=NucleiCounts,fill=CellState)) + geom_bar(position="stack", stat="identity")+ # rs1582763 plots
		scale_fill_manual(values=c(brewer.pal(9, 'Set1')[1:length(unique(tab[,3]))])) +
		theme_bw()+
		theme(panel.grid.major.x = element_blank() ,panel.background=element_rect(fill='white'),#fill = "white",
		axis.text.x = element_text(color='black',size=20,angle = 90, vjust = 0.5, hjust=1), 
		axis.text.y = element_text(color='black',size=20), axis.title.y = element_text(size=20),
		axis.line = element_line(linetype='blank'),axis.ticks = element_blank(),
		legend.title = element_text(size=20),
		legend.text = element_text(size=16)) + 
		labs(x=element_blank(),fill=str_wrap('Micro Cell State',width=7)) 

ggsave(p,w=20,file='~/scratch/fig2b_barplot_absolute_counts_bySample.pdf')

write.table(tab, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig7g.csv', sep = ', ', quote = F, row.names = F, col.names = T)
