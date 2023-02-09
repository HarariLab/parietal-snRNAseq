library(Seurat)

cellTypes = c('astro','microglia','n_ex','n_inh','oligo','opc')
#cellTypes = c('microglia','n_ex','n_inh')

tmp_pheno = NULL
tmp = lapply(cellTypes,function(cellType){
	if (cellType %in% c('n_ex','n_inh')){
		n = readRDS('/home/brasel/SingleCellProjects/dataObjects/neuron.rds')
		n$clusters = Idents(n)
		if(cellType == 'n_ex') { n = subset(n,subset=clusters %in% c(0:2,6))
		} else{ n = subset(n,subset=clusters %in% c(3:5)) }
	} else { n = readRDS(sprintf('/home/brasel/SingleCellProjects/dataObjects/%s.rds',cellType)) }
	
	clusterCounts = table(n@meta.data$Sample_ID,Idents(n))
	totalCellType = rowSums(clusterCounts)
	clusterCounts = clusterCounts/totalCellType
	
	pheno = data.frame(rbind(clusterCounts))
	colnames(pheno) = paste0('prop_',0:(ncol(clusterCounts)-1)) #c('prop_0','prop_1','prop_2','prop_3','prop_4')
	clusters = 0:(ncol(clusterCounts)-1)
	
	phenoFile = read.csv('/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/Brain_pheno_jorge.csv',row.names=2)
	phenoFile = phenoFile[rownames(pheno),]
	
	### Add in Age of Death ###
	AOD = unique(n@meta.data[,c('Sample_ID','AOD')])
	pheno = merge(pheno,AOD,by.x='row.names',by.y='Sample_ID')
	pheno$AOD <- as.numeric(pheno$AOD)
	
	library(lme4)
	### add in one other covariate at a time to see if it removes the significance
	covarPos <- c()
	
	#Sex (additive Model)
	SEX = unique(n@meta.data[,c('Sample_ID','Gender')])
	colnames(SEX)[2] = 'SEX'
	pheno = merge(pheno,SEX,by.x='Row.names',by.y='Sample_ID')
	
	rs1582763 = unique(n@meta.data[,c('Sample_ID','MS4')])
	colnames(rs1582763)[2] = 'rs1582763'
	pheno = merge(pheno,rs1582763,by.x='Row.names',by.y='Sample_ID')
	pheno[grep('GG',pheno$rs1582763),'rs1582763'] <- 0
	pheno[grep('AG',pheno$rs1582763),'rs1582763'] <- 1
	pheno[grep('AA',pheno$rs1582763),'rs1582763'] <- 2
	pheno$rs1582763 <- as.numeric(pheno$rs1582763)
	
	#TREM2
	nTREM2 = unique(n@meta.data[,c('Sample_ID','nTREM2')])
	nTREM2$nTREM2[nTREM2$nTREM2 != 'TREM2'] = 0
	nTREM2$nTREM2[nTREM2$nTREM2 == 'TREM2'] = 1
	pheno = merge(pheno,nTREM2,by.x='Row.names',by.y='Sample_ID')
	pheno$nTREM2 = as.numeric(pheno$nTREM2)
	
	pheno$TREM2_reduced = as.numeric(as.factor((phenoFile$TREM2_type %in% c('R62H','H157Y','R47H') )))-1
	pheno$TREM2_reduced[is.na(pheno$TREM2_reduced)] = 0
	
	nAPOE = unique(n@meta.data[,c('Sample_ID','nAPOE')])
	pheno = merge(pheno,nAPOE,by.x='Row.names',by.y='Sample_ID')
	pheno[-grep('4',pheno$nAPOE),'nAPOE'] <- 0
	pheno[grep('4',pheno$nAPOE),'nAPOE'] <- 1
	pheno$nAPOE <- as.numeric(pheno$nAPOE)
	
	Final_Status = unique(n@meta.data[,c('Sample_ID','Status')])
	colnames(Final_Status)[2] = 'Final_Status'
	Final_Status$Final_Status = factor(Final_Status$Final_Status, levels=c('Neuro_CO','Neuro_Presympt','Neuro_AD','Neuro_ADAD','Neuro_OT'))
	pheno = merge(pheno,Final_Status,by.x='Row.names',by.y='Sample_ID')
	
	pheno$count_clusterCellType = totalCellType
	
	rownames(pheno) = pheno$Row.names
	pheno = pheno[,-1]
	
	minCells = 60
	filt_pheno <<- pheno[-which(pheno[,'count_clusterCellType'] < minCells),] #remove subjects with less than 60 cells in the cluster
	
	tmp = lapply(clusters,function(clust){
	#ADAD
	tmp_pheno <<- filt_pheno
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'ADAD'
	Group2 = 'nonADAD'
	Group3 = ''
	Group1Samples = sum(tmp_pheno$Final_Status == 'Neuro_ADAD')
	Group2Samples = sum(tmp_pheno$Final_Status != 'Neuro_ADAD')
	Group3Samples = ''
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$Final_Status == 'Neuro_ADAD')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$Final_Status != 'Neuro_ADAD')
	G3S_GT_0.01 = ''
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
        model <- paste0(paste0('prop_',clust), ' ~ Final_Status + SEX')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['Final_StatusNeuro_ADAD',],collapse=','),sep=','))
	cat('\n')

	#TREM2
	tmp_pheno <<- filt_pheno#[filt_pheno$Final_Status == 'Neuro_AD',]
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'TREM2'
	Group2 = 'nonTREM2'
	Group3 = ''
	Group1Samples = sum(tmp_pheno$nTREM2 == '1')
	Group2Samples = sum(tmp_pheno$nTREM2 != '1')
	Group3Samples = ''
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$nTREM2 == '1')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$nTREM2 != '1')
	G3S_GT_0.01 = ''
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
        model <- paste0(paste0('prop_',clust), ' ~ nTREM2 + SEX + AOD')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['nTREM2',],collapse=','),sep=','))
	cat('\n')

	#TREM2_reduced
	tmp_pheno <<- filt_pheno#[filt_pheno$Final_Status == 'Neuro_AD',]
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'TREM2_reduced'
	Group2 = 'other'
	Group3 = ''
	Group1Samples = sum(tmp_pheno$TREM2_reduced == '1')
	Group2Samples = sum(tmp_pheno$TREM2_reduced != '1')
	Group3Samples = ''
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$TREM2_reduced == '1')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$TREM2_reduced != '1')
	G3S_GT_0.01 = ''
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
        model <- paste0(paste0('prop_',clust), ' ~ TREM2_reduced + SEX + AOD')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['TREM2_reduced',],collapse=','),sep=','))
	cat('\n')

	#rs1582763
	tmp_pheno <<- na.omit(filt_pheno)
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'AA'
	Group2 = 'AG'
	Group3 = 'GG'
	Group1Samples = sum(tmp_pheno$rs1582763 == '2')
	Group2Samples = sum(tmp_pheno$rs1582763 == '1')
	Group3Samples = sum(tmp_pheno$rs1582763 == '0')
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$rs1582763 == '2')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$rs1582763 == '1')
	G3S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$rs1582763 == '0')
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
        model <- paste0(paste0('prop_',clust), ' ~ rs1582763 + SEX + Final_Status')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['rs1582763',],collapse=','),sep=','))
	cat('\n')

	#sAD
	tmp_pheno <<- filt_pheno
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'sAD'
	Group2 = 'non_sAD'
	Group3 = ''
	Group1Samples = sum(tmp_pheno$Final_Status == 'Neuro_AD')
	Group2Samples = sum(tmp_pheno$Final_Status != 'Neuro_AD')
	Group3Samples = ''
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$Final_Status == 'Neuro_AD')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$Final_Status != 'Neuro_AD')
	G3S_GT_0.01 = ''
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
	#tmp_pheno$ADstatus = tmp_pheno$Final_Status == 'Neuro_AD'
        model <- paste0(paste0('prop_',clust), ' ~ Final_Status + SEX')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['Final_StatusNeuro_AD',],collapse=','),sep=','))
	cat('\n')
	
	#APOE
	tmp_pheno <<- filt_pheno[filt_pheno$Final_Status == 'Neuro_AD',]
	TotalSamples = nrow(tmp_pheno)
	Group1 = 'APOEe4+'
	Group2 = 'APOEe4-'
	Group3 = ''
	Group1Samples = sum(tmp_pheno$nAPOE == '1')
	Group2Samples = sum(tmp_pheno$nAPOE != '1')
	Group3Samples = ''
	SamplesGT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01)
	G1S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$nAPOE == '1')
	G2S_GT_0.01 = sum(tmp_pheno[,paste0('prop_',clust)] > 0.01 & tmp_pheno$nAPOE != '1')
	G3S_GT_0.01 = ''
	write(paste(cellType,clust,Group1,Group2,Group3,TotalSamples,Group1Samples,Group2Samples,Group3Samples,SamplesGT_0.01,G1S_GT_0.01,G2S_GT_0.01,G3S_GT_0.01,sep=','),file='/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/PropAnalyses/SampleNumbersForEachPropAnalysis/sampleNumbersForPropAnalyses_v2.csv',append=T)
	tmp = lapply(colnames(tmp_pheno)[grep('prop',colnames(tmp_pheno))],function(col) tmp_pheno[,col] <<- tmp_pheno[,col]^(1/3) )
        model <- paste0(paste0('prop_',clust), ' ~ nAPOE + SEX + AOD')
        re <- glm(formula = model,  data=tmp_pheno)
        coef <- summary(re)$coefficient
        cat(paste(cellType,clust,Group1,paste(coef['nAPOE',],collapse=','),sep=','))
	cat('\n')

})
})
