library(ggplot2)
library(ggrepel)
library(ggpubr)
#options(stringsAsFactors=F)

sampleE = read.csv('~/SingleCellProjects/dataObjects/sampleEntropy.csv')
batchE = read.csv('~/SingleCellProjects/dataObjects/batchEntropy.csv')
#sampleE = read.csv('Dragon/SingleCellProjects/dataObjects/sampleEntropy.csv')
#batchE = read.csv('Dragon/SingleCellProjects/dataObjects/batchEntropy.csv')

sE = stack(sampleE)
sE$CellType = sampleE$Cell.Type
sE$ind = as.character(sE$ind)
colnames(sE)[1] = 'Normalized_Entropy'
sE$Overall = sE$ind == 'Overall'
sE$CellType[sE$CellType == 'all'] = 'All'
sE$Normalized_Entropy = as.numeric(sE$Normalized_Entropy)
sE$ind[grep('c',sE$ind)] = unlist(substr(sE$ind[grep('c',sE$ind)],2,4))
p= ggplot(sE[!sE$Overall,],aes(x=Normalized_Entropy,color=CellType)) + geom_density()+ ggtitle('Sample Entropy')+ theme_minimal()+theme(axis.text=element_text(color='black'))+ rremove('xlab')
p2 = ggplot(sE,aes(x=Normalized_Entropy,y=CellType,color=CellType)) + geom_point() + 
	geom_point(data=subset(sE,subset=Overall),aes(x=Normalized_Entropy,y=CellType,fill=CellType),color='black',shape=24, stroke=1.5,size=2) +
	geom_label_repel(aes(label=ifelse(Normalized_Entropy <0.5,as.character(ind),'')),box.padding = 0.15, point.padding = 0.15) +
	theme_minimal()+theme(axis.text=element_text(color='black')) + scale_y_discrete(limits=rev)

bE = stack(batchE)
bE$CellType = batchE$Cell.Type
bE$ind = as.character(bE$ind)
colnames(bE)[1] = 'Normalized_Entropy'
bE$Overall = sE$ind == 'Overall'
bE$CellType[bE$CellType == 'all'] = 'All'
bE$Normalized_Entropy = as.numeric(bE$Normalized_Entropy)
bE$ind[grep('c',bE$ind)] = unlist(substr(bE$ind[grep('c',bE$ind)],2,4))
p3= ggplot(bE[!bE$Overall,],aes(x=Normalized_Entropy,color=CellType)) + geom_density() + ggtitle('Batch Entropy')+ theme_minimal()+theme(axis.text=element_text(color='black'))+rremove('xlab')
p4 = ggplot(bE,aes(x=Normalized_Entropy,y=CellType,color=CellType)) + geom_point() + 
	geom_point(data=subset(bE,subset=Overall),aes(x=Normalized_Entropy,y=CellType,fill=CellType),color='black',shape=24, stroke=1.5,size=2) +
	geom_label_repel(aes(label=ifelse(Normalized_Entropy <0.5,as.character(ind),'')),box.padding = 0.15, point.padding = 0.15) +
	theme_minimal()+theme(axis.text=element_text(color='black')) + scale_y_discrete(limits=rev)
p5 =ggarrange(p3,p,p4,p2,labels=c('a','c','b','d'),ncol=2,nrow=2)
ggsave(p5,w=14,file='~/scratch/Batch_Sample_Entropy_Distributions.pdf')

sE$Group = 'SampleEntropy'
bE$Group = 'BatchEntropy'
df = rbind(bE,sE)

write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig1.csv', sep = ', ', quote = F, row.names = F, col.names = T)
