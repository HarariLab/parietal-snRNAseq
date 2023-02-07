#oligo.5 (oligo-TFEB)
library(ggplot2)
effect = c(0.132260,0.052384,0.096754)#c(discovery,replication,meta)
sEr = c(0.058915,0.020585,0.035285)
range = sEr*3.92
ci_l = effect - (range/2)
ci_u = effect + (range/2)
df = data.frame(effect,ci_l,ci_u)
df$index = c(3,2,1)

p = ggplot(df,aes(y=index,x=effect,xmin=ci_l,xmax=ci_u)) + geom_point() + geom_errorbarh(height=.1) +
geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
scale_y_continuous(name = "", breaks=1:3, labels = c('META','ROSMAP','Discovery'))+ #, trans="reverse")
scale_x_continuous(limits=c(-0.2,0.4))+#, breaks = c(-2.5:1)
theme_minimal()

ggsave(p,h=1.5,file='~/scratch/oligo5_replication_treeplot.pdf')

#mic.2 (mic-reduced)
library(ggplot2)
effect = c(0.229137,0.07293,0.152622)#c(discovery,replication,meta)
sEr = c(0.091269,0.10168,0.071014)
range = sEr*3.92
ci_l = effect - (range/2)
ci_u = effect + (range/2)
df2 = data.frame(effect,ci_l,ci_u)
df2$index = c(3,2,1)

p = ggplot(df2,aes(y=index,x=effect,xmin=ci_l,xmax=ci_u)) + geom_point() + geom_errorbarh(height=.1) +
geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
scale_y_continuous(name = "", breaks=1:3, labels = c('META','ROSMAP','Discovery'))+ #, trans="reverse")
scale_x_continuous(limits=c(-0.2,0.6))+#, breaks = c(-2.5:1)
theme_minimal()

ggsave(p,h=1.5,file='~/scratch/mic2_replication_treeplot.pdf')

df$CellType = 'oligo-TFEB'
df2$CellType = 'mic-reduced'
df$Group = df2$Group = c('Discovery','Replication','Meta')
df3 = rbind(df2,df)
write.table(df3, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/SupFigs/SupFig9gh.csv', sep = ', ', quote = F, row.names = F, col.names = T)
