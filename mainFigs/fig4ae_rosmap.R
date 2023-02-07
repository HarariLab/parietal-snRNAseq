#dragon
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

myFullDF = NULL
######### mic TREM2_reduced #########################
cmic = readRDS('~/SingleCellProjects/dataObjects/Colonna_TREM2_data/colonna_seurat_object_microglia.rds')
t=unique(cmic@meta.data[,c('Sample','Status')])
status = t[order(t$Sample),'Status']
prop = as.data.frame.array(table(cmic$Sample,cmic$newClusters))
prop = prop/rowSums(prop)
#prop = cbind(prop,status)
prop = aggregate(prop,by=list(status),FUN=mean)
prop = prop[-2,]

df = prop
df$Group.1 = c('Other','TREM2 Reduced')
df = gather(df,'Cluster','Proportion',2:10)
df$CellType = 'Microglia'
myFullDF = rbind(myFullDF,df)

prop$Group.1 = c('Other',str_wrap('TREM2 Reduced',width=10))
prop = gather(prop,'Cluster','Proportion',2:10)

#colnames(prop) = c('Status','Cluster','Proportion')
#p = ggplot(prop,aes(x=factor(Group.1,levels=c('TREM2_reduced','other')),y=Proportion,fill=Cluster)) + geom_bar(position="stack", stat="identity")+ 
#	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab('TREM2') + ylab('Proportion') + labs(fill='Cell State')
p = ggplot(prop,aes(x=factor(Group.1,levels=c(str_wrap('TREM2 Reduced',width=10),'Other')),y=Proportion,fill=Cluster)) + geom_bar(position="stack", stat="identity")+ # for trem2_reduced plots
	##p = ggplot(sDF,aes(x=factor(Group.1,levels=c('TREM2','Wild Type')),y=Proportion,fill=Group.2)) + geom_bar(position="stack", stat="identity")+ # for trem2 plots
	##p = ggplot(sDF,aes(x=factor(Group.1,levels=c('GG','AG','AA')),y=Proportion,fill=Group.2)) + geom_bar(position="stack", stat="identity")+ # rs1582763 plots
		scale_fill_manual(values=c(brewer.pal(9, 'Set1')[1:length(unique(sDF[,2]))])) +
		theme(panel.background=element_rect(fill = "white",color='white'),
		axis.text.x = element_text(color='black',size=20,angle = 90, vjust = 0.5, hjust=1), 
		axis.text.y = element_text(color='black',size=20), axis.title.y = element_text(size=20),
		axis.line = element_line(linetype='blank'),axis.ticks = element_blank(),
		legend.title = element_text(size=20),
		legend.text = element_text(size=16)) + 
		labs(x=element_blank(),fill=str_wrap('Cell State',width=7))
ggsave(p,width = 3.5,file='/home/brasel/scratch/colonna_micro_trem2_reduced_merged.pdf')

######### oligo TREM2_reduced #########################
cmic = readRDS('~/SingleCellProjects/dataObjects/Colonna_TREM2_data/colonna_seurat_object_oligodendrocytes.rds')
t=unique(cmic@meta.data[,c('Sample','Status')])
status = t[order(t$Sample),'Status']
prop = as.data.frame.array(table(cmic$Sample,cmic$newClusters))
prop = prop/rowSums(prop)
#prop = cbind(prop,status)
prop = aggregate(prop,by=list(status),FUN=mean)
prop = prop[-2,]

df = prop
df$Group.1 = c('Other','TREM2 Reduced')
df = gather(df,'Cluster','Proportion',2:10)
df$CellType = 'Oligodendrocytes'
myFullDF = rbind(myFullDF,df)
write.table(myFullDF, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig4ae_rosmap.csv', sep = ', ', quote = F, row.names = T, col.names = T)

prop$Group.1 = c('Other',str_wrap('TREM2 Reduced',width=10))
prop = gather(prop,'Cluster','Proportion',2:10)

p = ggplot(prop,aes(x=factor(Group.1,levels=c(str_wrap('TREM2 Reduced',width=10),'Other')),y=Proportion,fill=Cluster)) + geom_bar(position="stack", stat="identity")+ # for trem2_reduced plots
	scale_fill_manual(values=c(brewer.pal(9, 'Set1')[1:length(unique(sDF[,2]))])) +
	theme(panel.background=element_rect(fill = "white",color='white'),
	axis.text.x = element_text(color='black',size=20,angle = 90, vjust = 0.5, hjust=1), 
	axis.text.y = element_text(color='black',size=20), axis.title.y = element_text(size=20),
	axis.line = element_line(linetype='blank'),axis.ticks = element_blank(),
	legend.title = element_text(size=20),
	legend.text = element_text(size=16)) + 
	labs(x=element_blank(),fill=str_wrap('Cell State',width=7))
ggsave(p,width=3.5,file='/home/brasel/scratch/colonna_oligo_trem2_reduced_merged.pdf')