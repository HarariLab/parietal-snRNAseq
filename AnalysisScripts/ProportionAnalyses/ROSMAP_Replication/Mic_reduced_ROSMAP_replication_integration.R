#!/bin/env Rscript
#this script reads in the subject matrices and column/row data and combines them into a seurat object. This object is
#then SCTransformed and integrated with our data to replicated our TREM2 enriched clusters in both neurons and 
#microglia.

library(Matrix)
library(Seurat)

#---------------make a seurat object out of the colonna data ------------------------
allM = c()

for (i in 1:13){#merge all subject matrices into a single matrix
	if (file.exists(paste0('AD',i,'_matrix.mtx.gz'))){
		amat = readMM(gzfile(paste0('AD',i,'_matrix.mtx.gz')))
		abar = read.table(gzfile(paste0('AD',i,'_barcodes.tsv.gz')))$V1
		afea = read.table(gzfile(paste0('AD',i,'_features.tsv.gz')))$V2
		rownames(amat) = afea
		colnames(amat) = paste0('AD',i,'_',abar)
		allM = cbind(allM,amat)
	}
	if (file.exists(paste0('C',i,'_matrix.mtx.gz'))){
		cmat = readMM(gzfile(paste0('C',i,'_matrix.mtx.gz')))
		cbar = read.table(gzfile(paste0('C',i,'_barcodes.tsv.gz')))$V1
		cfea = read.table(gzfile(paste0('C',i,'_features.tsv.gz')))$V2
		rownames(cmat) = cfea
		colnames(cmat) = paste0('C',i,'_',cbar)
		allM = cbind(allM,cmat)
	}
	if (file.exists(paste0('P',i,'_matrix.mtx.gz'))){
		pmat = readMM(gzfile(paste0('P',i,'_matrix.mtx.gz')))
		pbar = read.table(gzfile(paste0('P',i,'_barcodes.tsv.gz')))$V1
		pfea = read.table(gzfile(paste0('P',i,'_features.tsv.gz')))$V2
		rownames(pmat) = pfea
		colnames(pmat) = paste0('P',i,'_',pbar)
		allM = cbind(allM,pmat)
	}
}
#read in the meta data file and add it the AD, trem2, control statuses
met = read.csv('clusters_cellID_all.csv',stringsAsFactors=F)
met$Status = substr(met$Sample,1,1) 
met$Status[met$Status == 'A'] = 'AD'
met$Status[met$Status == 'C'] = 'CO'
met$Status[met$Status == 'P'] = 'R62H'

#make the format for the matrix column names and the metadata rownames match
rownames(met) = met$Barcodes
tmp =data.frame(strsplit(colnames(allS),'-'))
colnames(allM) = tmp[,1]

#filter the matrix by the nuclei contained in the meta data
allM = allM[,rownames(met)]

allS = CreateSeuratObject(allM)
allS@meta.data = cbind(allS@meta.data,met)
saveRDS(allS,file='colonna_seurat_object.rds')

#isolate the microglia
mic = subset(allS,subset=Label == 'Micro')
saveRDS(mic,file='colonna_seurat_object_microglia.rds')


#-----------------Integrate the colonna data and the harari data--------------
m = readRDS('/home/brasel/SingleCellProjects/dataObjects/microglia.rds')
options(future.globals.maxSize=10*1024^3)

mic = SCTransform(mic, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = T)
so.list= list(mic,m)
names(so.list) = c('colonna','Human')

features <- SelectIntegrationFeatures(so.list, nfeatures = 3000)
so.list <- PrepSCTIntegration(so.list, anchor.features = features, verbose = T)
reference_dataset <- which(names(so.list) == "Human")
anchors <- FindIntegrationAnchors(so.list, normalization.method = "SCT", anchor.features = features, verbose = T, reference=reference_dataset)
allso <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
DefaultAssay(allso) <- "integrated"

saveRDS(allso, file='colonna_harari_merged.rds')


#----------split integrated data into super small sublusters and merge based on overlap with harari labels----------

allso = readRDS('colonna_harari_merged.rds')

allso = FindNeighbors(allso)
allso = FindClusters(allso,res=15)

tmp = table(Idents(allso),allso$SCT_snn_res.0.2)
tmp2 =  as.data.frame.array(tmp)
tmp2$id = unlist(apply(tmp2,1,function(x) which(x==max(x))))-1
allso$newClusters = unlist(lapply(Idents(allso),function(x){ return(tmp2[x,'id']) }))

saveRDS(allso, file='colonna_harari_merged.rds')

