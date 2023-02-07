# calculate the corrected counts for all clusters
library(Matrix)
library(Seurat)
library(glmmTMB)
library(ggplot2)

date <- format(Sys.Date(),format='%m.%d.%y')


key <- read.table('/home/brasel/SingleCellProjects/Collabs/KipnisCollab/70BrainPheno_Jan2020.csv',sep=',',head=T,stringsAsFactors=F)
colnames(key)[1] <- 'Sample'
ADstatus <- key[,c('Sample','Final_Status')]
ADADstatus <- key[,c('Sample','Final_Status')]
OTstatus <- key[,c('Sample','Final_Status')]
AODstatus <- key[,c('Sample','AOD')]
SEXstatus <- key[,c('Sample','SEX')]
MS4Astatus<- key[,c('Sample','rs1582763')]
PMIstatus <- key[,c('Sample','PMI')]
cellCounts <- key[,c('Sample','count_clusterMicro')]

#clusterSet = readRDS('/home/brasel/SingleCellProjects/KipnisCollab/new.human.clusters.rds')
#ablated = readRDS('/home/brasel/SingleCellProjects/KipnisCollab/ablation.sig.human.rds')

data = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
meta = data@meta.data
clusterSet = as.character(Idents(data))
names(clusterSet) = colnames(data)
data = GetAssayData(data,slot='counts')


#apply the subject label to each cell
subjs <- unlist(lapply(colnames(data),function(x){ return(unlist(strsplit(x,'_'))[1]) }))

#apply subject SEX status to the cells
nm = SEXstatus[,1]
SEXstatus <- SEXstatus[,-1]
names(SEXstatus) <- nm
SEXstatus[SEXstatus == 2] <- 0
SEX <- unlist(lapply(subjs,function(x){ return(SEXstatus[x]) }))

#apply subject AOD status to the cells
AODstatus <- AODstatus[,-1]
names(AODstatus) <- nm
#AODstatus <- AODstatus - (mean(AODstatus))#center around the mean
AODstatus <- scale(AODstatus)
AOD <- unlist(lapply(subjs,function(x){ return(AODstatus[x,1]) }))

subjs <- as.numeric(factor(subjs))
AOD <- as.numeric(AOD); SEX<- as.numeric(SEX); 

#pResid = readRDS('/home/brasel/SingleCellProjects/KipnisColab/ablated_partialResiduals_postMouseIntegration_02.01.21.rds')
pResid = data.frame(clusterSet)
if( dim(pResid)[2] > 1){
    tmp = data.frame(as.numeric(as.character(pResid[,'clusterSet'])),log2(pResid[,gene]+1-min(pResid[,gene])))
    colnames(tmp) = c('cluster','Partial_Residuals')
} else{
    tmp = lapply(c('ABCA1','RELB','GPNMB','CD68','C5AR1','TNFAIP3','CD83','TMEM119','P2RY13','CX3CR1','BIN1','MED12L','SELPLG'),function(gene){
        y = data[gene,]
       newExpr = data.frame(clusterSet,stringsAsFactors = F)
       for (cur in seq.int(0, 8)){
        cSet = clusterSet
        cSet[cSet != cur] = -1
        cSet[cSet == cur] = 1
        cSet[cSet == -1] = 0
        re <- glmmTMB(y~ cSet + SEX + AOD + (1|subjs),ziformula=~1,family=nbinom2)
        #re <- glmmTMB(y~ cSet + SEX + ADgroups + ADADgroups + OTgroups + (1|subjs),ziformula=~0,family=gaussian)
        #re <- glmmTMB(y~ cSet + SEX + AOD + (1|subjs),ziformula=~0,family=gaussian)

        X <- getME(re,"X")
        beta <- fixef(re)$cond
        beta_X <- sweep(X,MARGIN=2,STATS=beta,FUN="*")
        p_resid <- sweep(beta_X,MARGIN=1,STATS=residuals(re),FUN="+")
        newExpr = cbind(newExpr,p_resid[,2])
       }
        colnames(newExpr) = c('clusterSet',paste0('c',seq.int(0,8)))
        newExpr$y = y
        newExpr$y[newExpr$clusterSet == 0] = newExpr$c0[newExpr$clusterSet == 0]
        newExpr$y[newExpr$clusterSet == 1] = newExpr$c1[newExpr$clusterSet == 1]
        newExpr$y[newExpr$clusterSet == 2] = newExpr$c2[newExpr$clusterSet == 2]
        newExpr$y[newExpr$clusterSet == 3] = newExpr$c3[newExpr$clusterSet == 3]
        newExpr$y[newExpr$clusterSet == 4] = newExpr$c4[newExpr$clusterSet == 4]
        newExpr$y[newExpr$clusterSet == 5] = newExpr$c5[newExpr$clusterSet == 5]
        newExpr$y[newExpr$clusterSet == 6] = newExpr$c6[newExpr$clusterSet == 6]
        newExpr$y[newExpr$clusterSet == 7] = newExpr$c7[newExpr$clusterSet == 7]
        newExpr$y[newExpr$clusterSet == 8] = newExpr$c8[newExpr$clusterSet == 8]
        pResid[,gene] <<- newExpr$y

        #tmp = data.frame(clusterSet,log2(newExpr$y+1-min(newExpr$y)))
        #colnames(tmp) = c('cluster','Partial_Residuals')
    })
}
rownames(pResid) = colnames(data)
saveRDS(pResid,'/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/microglia_resting_activated_genes_partialResidual.rds')

#plot the partial residuals
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(RColorBrewer)
myColors = c(brewer.pal(9, 'Set1'))
myColors[6] = '#EBEB00'

activated = c('ABCA1','RELB','GPNMB','CD68','C5AR1','TNFAIP3','CD83')
resting = c('TMEM119','P2RY13','CX3CR1','BIN1','MED12L','SELPLG')
genes = c(activated,resting)

corrected = readRDS('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/microglia_resting_activated_genes_partialResidual.rds')
corrected = corrected[,c('clusterSet',genes)]
df = aggregate(corrected,by=list(corrected$clusterSet),FUN='mean')
df = gather(df,'Gene','Expression',3:dim(df)[2])
df$Cluster = df$Group.1
df$Expression = df$Expression - (min(df$Expression)-0.1)
write.table(df[,-c(1,2)], '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFig8.csv', sep = ', ', quote = F, row.names = F, col.names = T) # nolint: line_length_linter.

p = ggplot(df[df$Gene %in% activated,],aes(x=factor(Gene,levels=rev(unique(Gene))), y=Expression, fill=factor(as.character(Cluster),levels=c(0:8)) ))+#, color=factor(as.character(Cluster)), alpha=factor(as.character(Cluster)))) +
geom_bar(width=0.7,stat="identity", position =position_dodge())+#"identity") +
 #scale_colour_manual(values=rev(myColors[1:3])) + #values=c('green','red','blue')) + #
 scale_fill_manual(values=rev(myColors[9:1]),guide=guide_legend(reverse=F)) + theme_minimal() +#values=c('green','red','blue')) + #
 theme(axis.title.x=element_blank(), axis.text.y=element_text(size=14,color='black'),axis.text.x=element_text(size=14,color='black') ) +
labs(fill='Cluster') 
p2 = ggplot(df[df$Gene %in% resting,],aes(x=factor(Gene,levels=rev(unique(Gene))), y=Expression, fill=factor(as.character(Cluster),levels=c(0:8)) ))+#, color=factor(as.character(Cluster)), alpha=factor(as.character(Cluster)))) +
geom_bar(width=0.7,stat="identity", position =position_dodge())+#"identity") +
 ##scale_colour_manual(values=rev(myColors[1:3])) + #values=c('green','red','blue')) + #
 scale_fill_manual(values=rev(myColors[9:1]),guide=guide_legend(reverse=F)) + theme_minimal() +#values=c('green','red','blue')) + #
 theme(axis.title.x=element_blank(), axis.text.y=element_text(size=14,color='black'),axis.text.x=element_text(size=14,color='black') ) +
labs(fill='Cluster')
q = ggarrange(p2,p, labels = c("a", "b"), ncol = 1, nrow = 2, common.legend=T,legend="right")

#ggsave(q,units='mm',w=183,h=100,file='~/scratch/mic2_restingActivatedMarkers_boxplot_all_cellStates.pdf')
