#dragon
library(Seurat)
colonna = readRDS('~/SingleCellProjects/dataObjects/Colonna_TREM2_data/colonna_seurat_object_oligodendrocytes.rds')

tmp = as.data.frame.array(table(colonna$Sample,colonna$newClusters)/rowSums(table(colonna$Sample,colonna$newClusters)))
tmp = tmp[rowSums(table(colonna$Sample,colonna$newClusters)) > 30,]
df = merge(tmp[,6,drop=F],unique(colonna@meta.data[,c('Sample','Status','msex','age_death','rs1582763')]), by.x='row.names',by.y='Sample')
colnames(df)[2] = 'c5'
df$age_death = as.numeric(as.character(df$age_death))
df$nTREM2 = as.numeric(df$Status == 'R62H')
df$Status[df$Status == 'R62H'] = 'AD'
df$cohort = 'colonna'

harari = readRDS('~/SingleCellProjects/dataObjects/oligo.rds')

tmp = as.data.frame.array(table(harari$Sample_ID,Idents(harari))/rowSums(table(harari$Sample_ID,Idents(harari))))
tmp = tmp[rowSums(table(harari$Sample_ID,Idents(harari))) > 60,]
pheno = read.csv('~/SingleCellProjects/MyProjects/67BrainsPaper/Brain_pheno_jorge.csv',row.names=2)
df2 = merge(tmp[,6,drop=F],unique(harari@meta.data[,c('Sample_ID','Status','nTREM2','Gender','AOD','MS4')]), by.x='row.names',by.y='Sample_ID')
colnames(df2)[2:7] = c('c5','Status','nTREM2','msex','age_death','rs1582763')
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
d = d[d$Status == 'AD',]

library(glmmTMB)
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$age_death + d$msex + d$cohort,ziformula=~0,family='gaussian'))

mydata = d
d = mydata[mydata$cohort == 'colonna',]
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$msex + d$age_death,ziformula=~0,family='gaussian'))
d = mydata[mydata$cohort == 'harari',]
summary(glmmTMB(d[,2]^(1/3) ~ d$nTREM2 + d$age_death + d$msex ,ziformula=~0,family='gaussian'))

##meta
#Estimate Std. Error z value Pr(>|z|)
#(Intercept)     0.166781   0.288976   0.577  0.56384
#d$nTREM2        0.096754   0.035285   2.742  0.00611 **
#d$age_death     0.001116   0.003264   0.342  0.73248
#d$msex1        -0.172643   0.034675  -4.979  6.4e-07 ***
#d$cohortharari -0.018728   0.039541  -0.474  0.63576
#
#
##colonna
#Estimate Std. Error z value Pr(>|z|)
#(Intercept)  0.564897   0.292567   1.931   0.0535 .
#d$nTREM2     0.052384   0.020585   2.545   0.0109 *
#d$msex      -0.148000   0.022988  -6.438 1.21e-10 ***
#d$age_death -0.003342   0.003297  -1.014   0.3108
#
##harari
#Estimate Std. Error z value Pr(>|z|)
#(Intercept) -0.095915   0.349436  -0.274 0.783712
#d$nTREM2     0.132260   0.058915   2.245 0.024773 *
#d$age_death  0.001688   0.004375   0.386 0.699616
#d$msex       0.197710   0.054845   3.605 0.000312 ***
#