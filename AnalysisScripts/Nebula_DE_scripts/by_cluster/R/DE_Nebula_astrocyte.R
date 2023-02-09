#!/bin/env Rscript

library(nebula)
library(Seurat)

pvalue.extreme <- function(z) {
  log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
  mantissa <- 10^(log10.pvalue %% 1)
  exponent <- log10.pvalue %/% 1
  ## or return(c(mantissa,exponent))
  #return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
  return(sprintf("%1.2fE%d",mantissa,exponent))
}

BH.extreme <- function(Zs) {
  len = length(Zs)

  log_pval <- log(2) + pnorm(abs(Zs), lower.tail = FALSE, log.p = TRUE) + log(len) - log(1:len)
  log10_pval <- log_pval/log(10)
  mantissa <- 10^(log10_pval %% 1)
  exponent <- log10_pval %/% 1
  BHs <- sprintf("%1.2fE%d",mantissa,exponent)
  
  for(i in 1:length(BHs)){ # for each BH calculated so far
    change = which(as.numeric(BHs[i]) < as.numeric(BHs[1:(i-1)])) # find if previous BH is larger 
    if(length(change) > 0){  BHs[change] = BHs[i]  } # If a previous BH is larger, replace it with the current BH
  }
  return(BHs)
}

args = commandArgs(trailingOnly=TRUE)
subclust1 = args[1]
subclust2 = args[2]

print(paste("Running", subclust1, subclust2))
t1 <- Sys.time()

date <- format(Sys.Date(),format='%m.%d.%y')

m = readRDS('/home/brasel/SingleCellProjects/dataObjects/astro.rds')

m@meta.data$Status = toupper(m@meta.data$Status)
m@meta.data$Status[m@meta.data$Status == 'NEURO_ADAD_NONCARRIER'] = 'NEURO_CO'
m = subset(m,subset=Status %in% c('NEURO_AD','NEURO_ADAD','NEURO_CO') )
# table(m@meta.data$Sample_ID)[order(table(m@meta.data$Sample_ID))]
m = subset(m,subset=Sample_ID %in% names(table(m@meta.data$Sample_ID)[table(m@meta.data$Sample_ID) > 50]))
m@meta.data$Gender = as.factor(m@meta.data$Gender)
m@meta.data$AOD = as.numeric(m@meta.data$AOD)

clusterSet = as.numeric(as.character(Idents(m)))
if (subclust2 != 'all'){
  # 1 vs 1
  clusterSet[clusterSet != subclust1 & clusterSet != subclust2] = 9
  clusterSet[clusterSet == subclust1] = -1
  clusterSet[clusterSet == subclust2] = 0
  clusterSet[clusterSet == -1] = 1
}else{
  ## 1 vs all
  clusterSet[clusterSet != subclust1] = 9
  clusterSet[clusterSet == subclust1] = 1
  clusterSet[clusterSet == 9] = 0
}
m@meta.data$cluster = clusterSet

m = subset(m,subset=cluster != 9)


pred = model.matrix(~cluster + Gender ,data=m@meta.data)


re_ln = nebula(GetAssayData(m,slot='counts'),m@meta.data$Sample_ID,pred=pred,method='LN',model='NBLMM', verbose = TRUE)


tmp = data.frame(gene=re_ln$summary$gene, logFC_cluster=re_ln$summary$logFC_cluster)
tmp$z.score = re_ln$summary$logFC_cluster/re_ln$summary$se_cluster
tmp$p.value = pvalue.extreme(tmp$z.score)
tmp <- tmp[order(-abs(tmp$z.score)),]
tmp$BH = BH.extreme(tmp$z.score)
merged = merge(tmp,re_ln$summary,by=c('gene', 'logFC_cluster'))
rownames(merged) = merged$gene
merged = merged[order(abs(merged$z.score),decreasing=T),c(2:14,1)]

dir <- '/home/bnovotny/single_nuclei/astrocyte/nebula_de/'
mainfilename <- paste('nebula_compareSubclusters_astrocyte', subclust1, subclust2,date, sep='_')
write.table(merged,file=paste0(dir, mainfilename,'.txt'),quote=F,sep='\t')
write.table(merged[as.numeric(merged$BH) < 0.05,],file=paste0(dir, 'bh_',mainfilename,'.txt'),quote=F,sep='\t')

t2 <- Sys.time()
print(t2-t1)

print(paste("Finished", subclust1, subclust2))
print("---")


