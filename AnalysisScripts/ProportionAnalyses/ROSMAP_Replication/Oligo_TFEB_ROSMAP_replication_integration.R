#!/bin/env Rscript
#this script reads in the subject matrices and column/row data and combines them into a seurat object. This object is
#then SCTransformed and integrated with our data to replicated our TREM2 enriched clusters in both neurons and
#microglia.

cAll = readRDS('/home/brasel/SingleCellProjects/Colonna_TREM2_data/colonna_seurat_object.rds')
oligo = subset(cAll,subset=Label %in% c('Oli0','Oli1'))
saveRDS(oligo,file='/home/brasel/SingleCellProjects/Colonna_TREM2_data/colonna_seurat_object_oligodendrocytes.rds')

#---------------------------------------------------------------------------------------------------
oligo = readRDS('/home/brasel/SingleCellProjects/Colonna_TREM2_data/colonna_seurat_object_oligodendrocytes.rds')

o = readRDS('/home/brasel/SingleCellProjects/dataObjects/oligo.rds')
options(future.globals.maxSize=10*1024^3)

oligo = SCTransform(oligo, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = T)
so.list= list(oligo,o)
names(so.list) = c('colonna','Human')

features <- SelectIntegrationFeatures(so.list, nfeatures = 3000)
so.list <- PrepSCTIntegration(so.list, anchor.features = features, verbose = T)
reference_dataset <- which(names(so.list) == "Human")
anchors <- FindIntegrationAnchors(so.list, normalization.method = "SCT", anchor.features = features, verbose = T, reference=reference_dataset)
allso <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
DefaultAssay(allso) <- "integrated"

allso = RunPCA(allso)
allso <- FindNeighbors(allso,dims=1:8)
allso <- FindClusters(allso, resolution = 15)

tmp = table(Idents(allso),allso$SCT_snn_res.0.2)
tmp2 =  as.data.frame.array(tmp)
tmp2$id = unlist(apply(tmp2,1,function(x) which(x==max(x))))-1
allso$newClusters = factor(unlist(lapply(Idents(allso),function(x){ return(tmp2[x,'id']) })),levels=c(0,1,2,3,4,5,6,7,8))

saveRDS(allso, file='colonna_harari_merged_oligodendrocytes.rds')

oligo$newClusters = unlist(lapply(rownames(oligo@meta.data),function(x){ return(allso@meta.data[x,'newClusters']) }))
Idents(oligo) = oligo$newClusters
saveRDS(oligo,file='colonna_seurat_object_oligodendrocytes.rds')

#-----------------------------------------------------------------------------------------------
merged = read.csv('colonna_rs1582763_genotype.csv')

allso$rs1582763 = unlist(lapply(allso$Sample,function(x){return(merged[x,'rs1582763'])}))
allso$age_death = as.character(unlist(lapply(allso$Sample,function(x){return(merged[x,'age_death'])})))
allso$age_death[allso$age_death == '90+'] = 90
allso$msex = unlist(lapply(allso$Sample,function(x){return(merged[x,'msex'])}))

oligo$rs1582763 = unlist(lapply(oligo$Sample,function(x){return(merged[x,'rs1582763'])}))
oligo$age_death = as.character(unlist(lapply(oligo$Sample,function(x){return(merged[x,'age_death'])})))
oligo$age_death[oligo$age_death == '90+'] = 90
oligo$msex = unlist(lapply(oligo$Sample,function(x){return(merged[x,'msex'])}))

saveRDS(allso, file='colonna_harari_merged_oligodendrocytes.rds')
saveRDS(oligo,file='colonna_seurat_object_oligodendrocytes.rds')
