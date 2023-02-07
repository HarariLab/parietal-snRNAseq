library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

magic = read.csv('~/SingleCellProjects/dataObjects/neurons_magic.csv',row.names=1)
meta = read.csv('~/SingleCellProjects/dataObjects/neuron_metaData.csv')
fulldf = cbind(magic,meta)
fulldf$Status[fulldf$Status == 'Neuro_CO'] = 'CO'
fulldf$Status[fulldf$Status == 'Neuro_Presympt'] = 'Pres'
fulldf$Status[fulldf$Status == 'Neuro_AD'] = 'sAD'
fulldf$Status[fulldf$Status == 'Neuro_ADAD'] = 'ADAD'
fulldf$Status[fulldf$Status == 'Neuro_OT'] = 'OTH'
fulldf$Status = factor(fulldf$Status,levels=c('CO','Pres','sAD','ADAD','OTH'))
fulldf$cellState[fulldf$cellState == 0] = 'N_Ex.0'
fulldf$cellState[fulldf$cellState == 1] = 'N_Ex.1'
fulldf$cellState[fulldf$cellState == 2] = 'N_Ex.2'
fulldf$cellState[fulldf$cellState == 3] = 'N_Inh.0'
fulldf$cellState[fulldf$cellState == 4] = 'N_Inh.1'
fulldf$cellState[fulldf$cellState == 5] = 'N_Inh.2'
fulldf$cellState[fulldf$cellState == 6] = 'N_Ex.3'
fulldf$cellState = factor(fulldf$cellState)

threshold = (2 * sd(magic$APOE)) + median(magic$APOE)

#merged status
apoeHighExpression = fulldf %>% group_by(Status) %>% summarize(highAPOE = sum(APOE > threshold))
propHighMerged = apoeHighExpression$highAPOE/table(fulldf$Status)

#by subject then status
bySubj = fulldf %>% group_by(Sample_ID,Status) %>% summarize(highAPOE = sum(APOE > threshold))
bySubj$prop = bySubj$highAPOE/table(fulldf$Sample_ID)
propHighBySubj = aggregate(bySubj,by=list(bySubj$Status),FUN=mean)
propHighBySubj = propHighBySubj[,c(-2,-3)]
colnames(propHighBySubj)[1] = 'Status'
propHighBySubj$sd = aggregate(bySubj[,3:4],by=list(bySubj$Status),FUN=sd)[,'prop']
propHighBySubj$Status = factor(propHighBySubj$Status,levels=c('Neuro_CO','Neuro_Presympt','Neuro_AD','Neuro_ADAD','Neuro_OT'))
#p = ggplot(propHighBySubj,aes(x=Status,y=prop)) + geom_bar(stat='identity',fill='black') +
	#geom_errorbar(aes(ymin=prop, ymax=prop+sd), width=.2, position=position_dodge(.9)) + 
	#geom_dotplot(data=bySubj,aes(x=Status,y=prop),binaxis = "y", stackdir = "center",dotsize =0.4,fill='white') +
p = ggplot(bySubj,aes(x=Status,y=prop)) +	geom_boxplot()+ #geom_dotplot(binaxis = "y", stackdir = "center") +
	ylab('Proportion APOE-High')+ theme_minimal() + 
	theme(axis.text.y = element_text(size=20,color='black'),axis.text.x = element_text(size=20,color='black',angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24),legend.title = element_text(size=24),legend.text = element_text(size=20))
ggsave(p,units='mm',h=220,w=150,file='~/scratch/prop_APOE_high_boxplot.pdf')
write.table(bySubj, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig21a.csv', sep = ', ', quote = F, row.names = F, col.names = T)

#by subject then status but split by subcluster
myColors = c(brewer.pal(9, 'Set1'))
myColors[6] = '#EBEB00'

byCellState = fulldf %>% group_by(Sample_ID,Status,cellState) %>% summarize(highAPOE = sum(APOE > threshold)) %>% as.data.frame()
subjByCellState = as.data.frame(table(fulldf$Sample_ID,fulldf$cellState))
colnames(subjByCellState) = c('Sample_ID','cellState','Counts')
subjByCellState  = subjByCellState[order(subjByCellState$Sample_ID),]
byCellState = merge(byCellState,subjByCellState,by=c('Sample_ID','cellState'),all=T)
byCellState$Prop = byCellState$highAPOE/byCellState$Counts
byCellState$Prop[is.na(byCellState$Prop)] = 0
tmp = lapply(byCellState[is.na(byCellState$Status),'Sample_ID'],function(sid){
	byCellState[byCellState$Sample_ID == sid,'Status'] <<- byCellState[byCellState$Sample_ID == sid,'Status'][1]
})
byCellState$StatusName = byCellState$Status
byCellState$Status = as.numeric(byCellState$Status)
byCellState = aggregate(byCellState,by=list(byCellState$StatusName,byCellState$Status,byCellState$cellState),FUN=mean)
byCellState = byCellState[,c(-4,-5,-6,-10)]
colnames(byCellState) = c('StatusName','Status','CellState','highAPOE','Counts','Proportion_APOE_High')
p = ggplot(byCellState,aes(x=Status,y=Proportion_APOE_High,color=CellState,fill=CellState)) + geom_line() + geom_point()+ theme_minimal() + 
	theme(axis.text.y = element_text(size=20,color='black'),axis.text.x = element_text(size=20,color='black'),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24),legend.title = element_text(size=24),legend.text = element_text(size=20)) +
	scale_color_manual(values=myColors)
ggsave(p,units='mm',h=220,w=220,file='~/scratch/prop_APOE_high_CellState.pdf')
write.table(byCellState, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig21b.csv', sep = ', ', quote = F, row.names = F, col.names = T)

ylim = max(fulldf$MHC1)
p = ggplot(fulldf,aes(x=APOE,y=MHC1,color=cellState)) + geom_point(size=0.75)  + geom_smooth(method='lm',fill=NA,size=2)+ #+ geom_vline(xintercept=threshold, linetype="dashed", color = "black")
			theme_minimal() + theme(title = element_text(size=24),axis.text. = element_text(size=20,color='black'),axis.title = element_text(size=32),legend.title = element_text(size=24),legend.text = element_text(size=20)) +
			stat_cor(method = "pearson",label.x=0,label.y=ylim*1.1,size=10) + scale_color_manual(values=myColors) + rremove('legend')
#ggsave(p,units='mm',h=220,w=220,file='~/scratch/test.pdf')
q = facet(p, facet.by=c('cellState','Status'),panel.labs.font = list(size=28),panel.labs.background = list(color = 'lightgray', fill = 'lightgray'))
ggsave(q,units='mm',h=110*length(levels(fulldf$cellState)),w=110*length(levels(fulldf$Status)),file='~/scratch/neuron_APOE_MHC1_correlation_plots.pdf',limitsize=F)

