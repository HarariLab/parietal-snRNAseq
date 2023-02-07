library(ggplot2)
library(ggpubr)
library(RColorBrewer)

magic = read.csv('~/SingleCellProjects/dataObjects/neurons_magic_noiseAdded.csv',row.names=1)
meta = read.csv('~/SingleCellProjects/dataObjects/neuron_metaData.csv')
fulldf = cbind(magic,meta)
fulldf$Status[fulldf$Status == 'Neuro_CO'] = 'CO'
fulldf$Status[fulldf$Status == 'Neuro_Presympt'] = 'Pres'
fulldf$Status[fulldf$Status == 'Neuro_AD'] = 'sAD'
fulldf$Status[fulldf$Status == 'Neuro_ADAD'] = 'ADAD'
fulldf$Status[fulldf$Status == 'Neuro_OT'] = 'OTH'
fulldf$Status = factor(fulldf$Status,levels=c('CO','Pres','sAD','ADAD','OTH'))
fulldf$cellState[fulldf$cellState == 0] = 'EN.0'
fulldf$cellState[fulldf$cellState == 1] = 'EN.1'
fulldf$cellState[fulldf$cellState == 2] = 'EN.2'
fulldf$cellState[fulldf$cellState == 3] = 'IN.0'
fulldf$cellState[fulldf$cellState == 4] = 'IN.1'
fulldf$cellState[fulldf$cellState == 5] = 'IN.2'
fulldf$cellState[fulldf$cellState == 6] = 'EN.3'
fulldf$cellState = factor(fulldf$cellState)


tmp = lapply(unique(fulldf$Status),function(stat){
	tap = lapply(unique(fulldf$cellState),function(cState){
		fulldf[fulldf$Status == stat & fulldf$cellState == cState,'APOE'] <<- rank(fulldf[fulldf$Status == stat & fulldf$cellState == cState,'APOE'])
		fulldf[fulldf$Status == stat & fulldf$cellState == cState,'MHC1'] <<- rank(fulldf[fulldf$Status == stat & fulldf$cellState == cState,'MHC1'])
	})
})

myColors = c(brewer.pal(9, 'Set1'))
myColors[6] = '#EBEB00'

myPlots = list()						
tap = lapply(as.factor(levels(fulldf$cellState)),function(cState){
	tmp = lapply(levels(fulldf$Status),function(stat){
		#stat = fulldf$Status[1]
		ylim = max(fulldf[fulldf$Status == stat & fulldf$cellState == cState,'APOE'])
		p = ggplot(fulldf[fulldf$Status == stat & fulldf$cellState == cState,],aes(x=APOE,y=MHC1)) + geom_point(size=0.75,color=myColors[cState])  + geom_smooth(method='lm',fill=NA,size=2,color=myColors[cState])+ #+ geom_vline(xintercept=threshold, linetype="dashed", color = "black")
			theme_minimal() + theme(title = element_text(size=24),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24),legend.title = element_text(size=24),legend.text = element_text(size=20)) +
			stat_cor(method = "spearman",label.x=0,label.y=ylim*1.1,size=10) #+ scale_y_continuous(limits=c(0,ylim*1.25))#+ ggtitle(stat) 
		#ggsave(p,units='mm',h=220,w=220,file='~/scratch/test.pdf')
		myPlots[[length(myPlots)+1]] <<- p
#ggsave(p,units='mm',h=220,w=220,file=sprintf('~/scratch/%s_APOE_MHC1.pdf',stat))
	})
})
q = ggarrange(plotlist=myPlots,ncol=length(levels(fulldf$Status)),nrow=length(levels(fulldf$cellState)))
ggsave(q,units='mm',h=220*length(levels(fulldf$cellState)),w=220*length(levels(fulldf$Status)),file='~/scratch/neuron_APOE_MHC1_correlation_plots.pdf',limitsize=F)


write.table(fulldf[,-c(9,11,12,13,14,15)], '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig21c.csv', sep = ', ', quote = F, row.names = F, col.names = T)
