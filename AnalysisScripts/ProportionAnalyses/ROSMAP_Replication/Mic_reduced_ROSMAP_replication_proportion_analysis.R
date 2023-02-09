#dragon
library(Seurat)
colonna = readRDS('~/SingleCellProjects/Colonna_TREM2_data/colonna_seurat_object_microglia.rds')

tmp = as.data.frame.array(table(colonna$Sample,colonna$newClusters)/rowSums(table(colonna$Sample,colonna$newClusters)))
tmp = tmp[rowSums(table(colonna$Sample,colonna$newClusters)) > 30,]
df = merge(tmp[,3,drop=F],unique(colonna@meta.data[,c('Sample','Status','msex','age_death','rs1582763')]), by.x='row.names',by.y='Sample')
colnames(df)[2] = 'c2'
df$age_death = as.numeric(as.character(df$age_death))
df$nTREM2 = as.numeric(df$Status == 'R62H')
df$Status[df$Status == 'R62H'] = 'AD'
df$cohort = 'colonna'

harari = readRDS('~/SingleCellProjects/dataObjects/microglia.rds')

tmp = as.data.frame.array(table(harari$Sample_ID,harari$SCT_snn_res.0.2)/rowSums(table(harari$Sample_ID,harari$SCT_snn_res.0.2)))
tmp = tmp[rowSums(table(harari$Sample_ID,harari$SCT_snn_res.0.2)) > 60,]
pheno = read.csv('~/SingleCellProjects/Brain_pheno_jorge.csv',row.names=2)
df2 = merge(tmp[,3,drop=F],unique(harari@meta.data[,c('Sample_ID','Status','nTREM2','Gender','AOD','MS4')]), by.x='row.names',by.y='Sample_ID')
colnames(df2)[2:7] = c('c2','Status','nTREM2','msex','age_death','rs1582763')
rownames(df2) = df2$Row.names
df2$age_death = as.numeric(as.character(df2$age_death))
df2$nTREM2[df2$Row.names %in% rownames(pheno)[which(pheno$TREM2_type %in% c('R62H','H157Y','R47H'))]] = 'TREM2_reduced'
df2$nTREM2 = as.numeric(df2$nTREM2 == 'TREM2_reduced')
df2$Status[df2$Status %in% c('Neuro_ADAD_noncarrier','neuro_CO','Neuro_CO')] = 'CO'
df2$Status[df2$Status == 'Neuro_AD'] = 'AD'
df2$Status[df2$Status %in% c('Neuro_DLB','Neuro_OT')] = 'OTH'
df2$Status[df2$Status %in% c('Neuro_presyntomatic','Neuro_Presyntomatic')] = 'presympt'
df2$msex[df2$msex == 2] = 0
df2$cohort = 'harari'

d = rbind(df,df2)

library(glmmTMB)
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$Status + d$msex + d$cohort,ziformula=~0,family='gaussian'))

mydata = d
d = mydata[mydata$cohort == 'colonna',]
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$Status + d$msex + d$age_death,ziformula=~0,family='gaussian'))
d = mydata[mydata$cohort == 'harari',]
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$Status + d$msex ,ziformula=~0,family='gaussian'))

##meta
#Estimate Std. Error z value Pr(>|z|)
#(Intercept)         0.292779   0.064585   4.533 5.81e-06 ***
#d$nTREM2            0.155900   0.068372   2.280   0.0226 *
#d$StatusCO         -0.023993   0.074136  -0.324   0.7462
#d$StatusNeuro_ADAD -0.097521   0.078803  -1.238   0.2159
#d$StatusOTH         0.009048   0.098846   0.092   0.9271
#d$Statuspresympt   -0.091308   0.143874  -0.635   0.5257
#d$msex1             0.031197   0.052002   0.600   0.5486
#d$cohortharari      0.080503   0.061948   1.300   0.1938
#
##Colona
#Estimate Std. Error z value Pr(>|z|)
#(Intercept) -1.08976    1.09270  -0.997   0.3186
#d$nTREM2     0.06804    0.09832   0.692   0.4889
#d$StatusCO  -0.19551    0.09768  -2.002   0.0453 *
#d$msex       0.08743    0.08520   1.026   0.3048
#d$age_death  0.01658    0.01236   1.341   0.1800
#
##harari
#Estimate Std. Error z value Pr(>|z|)
#(Intercept)         0.34671    0.06680   5.191  2.1e-07 ***
#d$nTREM2            0.19514    0.08751   2.230   0.0258 *
#d$StatusCO          0.13958    0.10367   1.346   0.1782
#d$StatusNeuro_ADAD -0.04836    0.08175  -0.592   0.5542
#d$StatusOTH         0.03946    0.09961   0.396   0.6920
#d$Statuspresympt   -0.04583    0.14447  -0.317   0.7511
#d$msex             -0.01276    0.06323  -0.202   0.8401