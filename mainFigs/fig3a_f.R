#dragon
m = read.csv('/home/brasel/SingleCellProjects/dataObjects/Proportions_merged.csv',head=T,row.names=1)

m = m[-which(rownames(m) %in% c('sample18','sample45','sample48')),]
#m = m^(1/3) ##############################################################################################################
pheno = read.csv('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/Brain_pheno_jorge.csv',head=T)

m = merge(m,pheno[,c('Sample','TREM2_type','Final_Status','SEX','APOE','rs1582763')],by.x='row.names',by.y='Sample')
rownames(m) = m$Row.names
m$TREM2_type = as.character(m$TREM2_type)
m$TREM2_type[is.na(m$TREM2_type)] = 'Wild Type'
m$TREM2_all = m$TREM2_type
m$TREM2_all[m$TREM2_all != 'Wild Type'] = 'TREM2'
m$TREM2_reduced = m$TREM2_all
m$TREM2_reduced[m$TREM2_type %in% c('R47H','R62H','H157Y')] = 'TREM2 Reduced'
m$TREM2_reduced[m$TREM2_reduced != 'TREM2 Reduced'] = 'Other'
m$SEX = as.character(m$SEX)
m$SEX[m$SEX == 1] = 'Male'
m$SEX[m$SEX == 2] = 'Female'
m$APOE = as.character(m$APOE)
m$rs1582763 = factor(m$rs1582763,levels=c('GG','AG','AA'))
m$Final_Status = as.character(m$Final_Status)
m$Final_Status[m$Final_Status== 'Neuro_ADAD_noncarrier'] = 'Neuro_CO'
m$Final_Status[m$Final_Status== 'Neuro_DLB'] = 'Neuro_OT'

m$Final_Status[m$Final_Status== 'Neuro_OT'] = 'OTH'
m$Final_Status[m$Final_Status== 'Neuro_AD'] = 'sAD'
m$Final_Status[m$Final_Status== 'Neuro_CO'] = 'CO'
m$Final_Status[m$Final_Status== 'Neuro_Presyntomatic'] = 'Pres'
m$Final_Status[m$Final_Status== 'ADAD'] = 'ADAD'

library(tidyr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(GGally)
library(RColorBrewer)
library(stringr)

myFullDF = NULL
tmp = lapply(c('astro','micro','neuron','oligo','opc','all'),function(cellType){
	myCols = grep(cellType,colnames(m))
	if (cellType == 'all'){
		myCols = grep('Astro|Micro|Neuro|Oligo|OPC|Endo',colnames(m))
	}
	assign(cellType,gather(m[,c(myCols,1,45:51)],'Cluster','Proportion',1:length(myCols)))
	tmp = lapply(c('Final_Status'),function(covar) { #colnames(m)[45:51]
		if (cellType == 'neuron' & covar == 'TREM2_all'){
			sub = m[,c(myCols[c(1:3,length(myCols))],1,45:51)]
			sub[,1:4] = sub[,1:4]/rowSums(sub[,1:4])
			assign(cellType,gather(sub,'Cluster','Proportion',1:4))
		}
		sDF = eval(as.symbol(cellType))
		sDF = sDF[!is.na(sDF$Proportion),]
		if(grepl('TREM2',covar)){
			sDF = sDF[sDF$Final_Status == 'sAD',]
		}
		if(covar == 'rs1582763'){ sDF = sDF[!is.na(sDF$rs1582763),] }
		sDF = aggregate(sDF,list(sDF[,covar],sDF$Cluster),mean)
		sDF$Group.1 = str_wrap(sDF$Group.1,width=10)
        df = data.frame(Status = sDF$Group.1, CellState = sDF$Group.2, Proportion = sDF$Proportion)
        df$CellType = cellType
        myFullDF <<- rbind(myFullDF,df)
		p = ggplot(sDF,aes(x=factor(Group.1,levels=c('ADAD','sAD','Pres','CO','OTH')),y=Proportion,fill=Group.2)) + geom_bar(position="stack", stat="identity")+ #for ADAD plots
			scale_fill_manual(values=c(brewer.pal(9, 'Set1')[1:length(unique(sDF[,2]))])) +
			theme(panel.background=element_rect(fill = "white",color='white'),
			axis.text.x = element_text(color='black',size=20,angle = 90, vjust = 0.5, hjust=1), 
			axis.text.y = element_text(color='black',size=20), axis.title.y = element_text(size=20),
			axis.line = element_line(linetype='blank'),axis.ticks = element_blank(),
			legend.title = element_text(size=20),
			legend.text = element_text(size=16)) + 
			labs(x=element_blank(),fill=str_wrap('Cell State',width=7))
			
		#ggsave(p,width=3.5,file=paste0('/home/brasel/scratch/',cellType,'_',covar,'_merged.pdf'))
	})
})

write.table(myFullDF, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig3a_f.csv', sep = ', ', quote = F, row.names = F, col.names = T) # nolint: line_length_linter.
