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

# Could also specify statuses to compare
args = commandArgs(trailingOnly=TRUE)
celltype = args[1]
subclust = args[2]
cutoff = args[3]

print(paste("Running", celltype, subclust))
t1 <- Sys.time()

m = readRDS(paste0('/home/brasel/SingleCellProjects/dataObjects/', celltype, '.rds'))


m@meta.data$Status = toupper(m@meta.data$Status)

# Could specify status in arguments
m = subset(m,subset=Status %in% c('NEURO_ADAD','NEURO_AD') )
m@meta.data$Status = factor(m@meta.data$Status, levels = c("NEURO_AD", "NEURO_ADAD"))

m = subset(m,subset=Sample_ID %in% names(table(m@meta.data$Sample_ID)[table(m@meta.data$Sample_ID) > as.numeric(cutoff)]))
m@meta.data$Gender = as.factor(m@meta.data$Gender)
m@meta.data$AOD = as.numeric(m@meta.data$AOD)

# Subset to cluster we want
clusterSet = as.numeric(as.character(Idents(m)))


if (subclust != 'all'){
  # 1 vs 1
  clusterSet[clusterSet != subclust] = 99
}

m@meta.data$cluster = clusterSet

m = subset(m,subset=cluster != 99)


pred = model.matrix(~Status + Gender,data=m@meta.data)

re_ln = nebula(GetAssayData(m,slot='counts'),m@meta.data$Sample_ID,pred=pred,method='LN',model='NBLMM')


tmp = data.frame(gene=re_ln$summary$gene, logFC_StatusNEURO_ADAD=re_ln$summary$logFC_StatusNEURO_ADAD)
tmp$z.score = re_ln$summary$logFC_StatusNEURO_ADAD/re_ln$summary$se_Status
tmp$p.value = pvalue.extreme(tmp$z.score)
tmp <- tmp[order(-abs(tmp$z.score)),]
tmp$BH = BH.extreme(tmp$z.score)
merged = merge(tmp,re_ln$summary,by=c('gene', "logFC_StatusNEURO_ADAD"))
rownames(merged) = merged$gene
merged = merged[order(abs(merged$z.score),decreasing=T),c(2:14,1)]


dir <- paste0('/home/bnovotny/single_nuclei/', celltype, '/nebula_de/DE_by_status/ADADvssAD/')
mainfilename <- paste('nebula_compareStatus_ADADvssAD', celltype, subclust, sep='_')


write.table(merged,file=paste0(dir, mainfilename,'.txt'),quote=F,sep='\t')
write.table(merged[as.numeric(merged$BH) < 0.05,],file=paste0(dir, 'bh_',mainfilename,'.txt'),quote=F,sep='\t')

t2 <- Sys.time()
print(t2-t1)
print(paste("Finished", celltype, subclust))
print("---")






