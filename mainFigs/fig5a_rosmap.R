######### mic rs1582763 #########################
cmic = readRDS('~/SingleCellProjects/dataObjects/Colonna_TREM2_data/colonna_seurat_object_microglia.rds')
t=unique(cmic@meta.data[,c('Sample','rs1582763')])
rs1582763 = t[order(t$Sample),'rs1582763']
prop = as.data.frame.array(table(cmic$Sample,cmic$newClusters))
prop = prop/rowSums(prop)
prop = aggregate(prop,by=list(rs1582763),FUN=mean)
prop = gather(prop,'Cluster','Proportion',2:10)

df = prop
write.table(df, '/home/brasel/SingleCellProjects/MyProjects/67BrainsPaper/createSourceFiles_NatComm/fig5a_rosmap.csv', sep = ', ', quote = F, row.names = T, col.names = T)

p = ggplot(prop,aes(x=factor(Group.1,levels=c('GG','AG','AA')),y=Proportion,fill=Cluster)) + geom_bar(position="stack", stat="identity")+ # rs1582763 plots
	scale_fill_manual(values=c(brewer.pal(9, 'Set1')[1:length(unique(sDF[,2]))])) +
	theme(panel.background=element_rect(fill = "white",color='white'),
	axis.text.x = element_text(color='black',size=20,angle = 90, vjust = 0.5, hjust=1), 
	axis.text.y = element_text(color='black',size=20), axis.title.y = element_text(size=20),
	axis.line = element_line(linetype='blank'),axis.ticks = element_blank(),
	legend.title = element_text(size=20),
	legend.text = element_text(size=16)) + 
	labs(x=element_blank(),fill=str_wrap('Cell State',width=7))
ggsave(p,width = 3.5,file='/home/brasel/scratch/colonna_micro_rs1582763_merged.pdf')