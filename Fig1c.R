library(qs)
library(Seurat)
library(ggplot2)

all = qread('/home/brasel/SingleCellProjects/dataObjects/allCellTypes_garnettCleaned_finalObj_294114nuclei.qs')

p = DotPlot(all,col.min=0.00,cols=c("white","red"),
features=rev(c('TNR','PDGFRA', #'CNTN1', #OPCs
'QDPR','TULP4','PIP4K2A','TMEM144','SLC44A1','CNP','SLAIN1','ANLN','MOBP','SCD','CNDP1','MBP','ERMN','CLDND1','UGT8','TTLL7','SLC24A2','ENPP2','TF','PLP1', #'OPALIN','SEPT4', #Oligos
'PCLO','GABRA1','SCG2','UCHL1','GABRB2','GAD1','VSNL1','STMN2','RTN1','SNAP25','SYT1','SCN2A','DLX6-AS1','SYNPR', #'CNR1','RELN', #Neurons
'C3','CCL4','MSR1','OLR1','GPR183','CD53','LCP1','HAVCR2','LHFPL2','PLEK','HLA-DRA','CD74','CX3CR1','C3AR1','CLEC7A','B3GNT5', #'CD83','IL1B','CH25H', #Microglia
'NET1','FN1','APOLD1','TM4SF1','ABCG2','RGS5','DCN', #Endothelial
'SLCO1C1','ATP1B2','CPE','AQP4','EDNRB','SLC4A4','GJA1','SLC1A2','CLU','ETNPPL','F3','RYR3')))  + #'DIO2','GPR37L1', #Astrocytes
theme(axis.text.x = element_text(size = 18,angle=90, hjust=1),axis.text.y = element_text(size=14),legend.position= 'none') #+ scale_y_discrete(position='right')#+ coord_flip() + scale_x_discrete(position='top')

write.table(p$data,'/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig1c.csv',sep=',',quote=F,row.names=T,col.names=T)
