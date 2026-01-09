a <- AverageExpression(HCC, assays = "RNA",slot = "data",group.by =c("tissue_source","cell_type1"))
a$RNA["SAMSN1",]
# Adjacent liver_T Cells               Adjacent liver_NK          Adjacent liver_B Cells 
# 3.8078647                       2.9788686                       1.9731207 
# Adjacent liver_Mye   Adjacent liver_Parenchyma Cells            Adjacent liver_Tumor 
# 3.0190963                       0.9124489                       0.4297335 
# Tumor_T Cells   Tumor_NK     Tumor_B Cells 
#4.5858959      3.7964177      2.7326824 
# Tumor_Mye          Tumor_Parenchyma Cells                     Tumor_Tumor 
#2.8797813                       0.8912194                       0.2741002 

Stacked_VlnPlot(HCC,features = "SAMSN1",group.by = "cell_type1",split.by = "tissue_source",
                split.plot = F, log = F,colors_use =c("red","blue"))



plotdata <- read.xlsx("./2206/CK5fc.xlsx",sheet = 1)


pdf("./2206/CK5FC.pdf",height = 5,width = 5)
ggplot(plotdata,aes(x = log(FoldChange),y = Celltypes,fill = log(FoldChange))) +
  geom_bar(stat = "identity",position = "stack",color = "#e0e0e0",size = 1,width = 0.3) +
  theme_par() + #背景为空
  theme(axis.title.x = element_text(color = "Black", size = 15, face = "plain"),
        axis.title.y.left  = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(colour="Black", size=15 ,face="plain"),
        axis.text = element_text(colour="Black", size=15 ,face="plain"),
        axis.ticks.y.left = element_blank()
  )+
  labs(fill = "Cluster",title = " ")+
 guides(fill="none")+
  geom_bar(position = position_stack(), stat ="identity", width = .7)+
  geom_text(aes(label = " "), position = position_stack(vjust = 0.5), size = 5)+
  scale_y_discrete(limits=rev(c('B Cells','NK ','T Cells','Parenchyma Cells','Myeloid',"Malignant Cells")))+
  scale_fill_gradient2(low = "Blue",high = "Red",mid = "Purple",midpoint = 0)
dev.off()