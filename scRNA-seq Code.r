library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(patchwork)
#library(FlexDot)
library(scales)#调色板包
library(ggthemes)
library(viridis)
library(openxlsx)
library(scCustomize)

HCC.data <- read.table(file = "./Data/HCC_log_tpm_expression_matrix.txt.gz",header = T)

metadata <- read.xlsx("./Data/metadata.xlsx", rowNames = T)
#数据修理
a <- HCC.data$gene
rownames(HCC.data)<- a
HCC.data <- HCC.data[,-1]
rm(a)

HCC_18 <- CreateSeuratObject(counts = HCC.data, meta.data = metadata)
rm(HCC.data,metadata)

HCC_18@meta.data$cell_type <- as.factor(HCC_18@meta.data$cell_type)
HCC_18@meta.data$orig.ident <- as.factor(HCC_18@meta.data$orig.ident)

sam.name <- "220330"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

##QC

HCC_18 <- Add_Mito_Ribo_Seurat(seurat_object = HCC_18, species = "Human")

# All functions contain
pdf(paste0("./220330/QC.pdf"),width = 12,height = 6)
p1 <- QC_Plot_UMIvsGene(seurat_object = HCC_18, 
                        low_cutoff_gene = 800, high_cutoff_gene = 12000, 
                        low_cutoff_UMI = 1000,high_cutoff_UMI = 20000)
p2 <- QC_Plot_GenevsFeature(seurat_object = HCC_18, feature1 = "percent_mito", 
                            low_cutoff_gene = 800,high_cutoff_gene = 12000, 
                            high_cutoff_feature = 20)
p1+p2
dev.off()
rm(p1,p2)
#### 表达量标准化 ####
HCC_18 <- NormalizeData(HCC_18, 
                        normalization.method = "LogNormalize",
                        scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
HCC_18 <- FindVariableFeatures(HCC_18, 
                               selection.method = "vst",
                               nfeatures = 2000)


#均一化#
HCC_18 <- ScaleData(
  object = HCC_18,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent_mito"))

#PCA计算
HCC_18 <- RunPCA(object = HCC_18, 
                 features = VariableFeatures(HCC_18),
                 verbose = F,npcs = 50)


ElbowPlot(HCC_18,ndims = 40)

HCC_18 <- FindNeighbors(HCC_18, dims = 1:20)
HCC_18 <- FindClusters(HCC_18, resolution = 1)



#### UMAP算法 ####
HCC_18<- RunUMAP(HCC_18, dims = 1:20)

pdf(paste0("./",sam.name,"/CellCluster_UMAP_split.pdf"),width = 11,height = 8)
DimPlot(HCC_18, pt.size = 0.75, reduction = "umap", label = T, group.by = "cell_type")
dev.off()


HCC_18 <- SetIdent(object = HCC_18, value = "cell_type")
levels(HCC_18@active.ident)
#[1] "C0_Tcell"   "C1_Tcell"   "C10_Tumor"  "C11_Mye."   "C12_Tumor"  "C13_Tumor"  "C14_Tumor" 
#[8] "C15_Mey."   "C16_Tumor"  "C17_Endo."  "C18_Epi."   "C19_Tcell"  "C2_Mye."    "C20_Mye"   
#[15] "C21_pDC"    "C22_Plasma" "C23_HSC"    "C3_Tcell"   "C4_NK"      "C5_Tcell"   "C6_Bcell"  
#[22] "C7_NK"      "C8_Mye."    "C9_Tumor"  
HCC_18 <- RenameIdents(object = HCC_18,
                       "C0_Tcell" = "T Cells",
                       "C1_Tcell" = "T Cells",
                       "C10_Tumor" = "Tumor",
                       "C11_Mye." = "Mye",
                       "C12_Tumor" = "Tumor",
                       "C13_Tumor" = "Tumor",
                       "C14_Tumor" = "Tumor",
                       "C15_Mey." = "Mye",
                       "C16_Tumor"= "Tumor",
                       "C17_Endo."= "Parenchyma Cells",
                       "C18_Epi." = "Parenchyma Cells",
                       "C19_Tcell"= "T Cells",
                       "C2_Mye." = "Mye",
                       "C20_Mye" = "Mye",
                       "C21_pDC" = "Mye",  
                       "C22_Plasma" = "B Cells" ,
                       "C23_HSC" = "Parenchyma Cells", 
                       "C3_Tcell"= "T Cells",
                       "C4_NK" = "NK",
                       "C5_Tcell" = "T Cells",
                       "C6_Bcell" = "B Cells", 
                       "C7_NK"  = "NK",
                       "C8_Mye."= "Mye",
                       "C9_Tumor"= "Tumor")
levels(HCC_18@active.ident)
#[1] "T Cells"          "Tumor"            "Mye"              "Parenchyma Cells"
#[5] "B Cells"          "NK"
HCC_18@active.ident <- factor(HCC_18@active.ident, levels = c("T Cells","NK","B Cells","Mye",
                                                              "Parenchyma Cells", "Tumor"))


pdf(paste0("./",sam.name,"/CellCluster_UMAP_celltype1.pdf"),width = 10,height = 8)
DimPlot(HCC_18, pt.size = 0.75, reduction = "umap", label = T)+
  scale_color_tableau(palette = "Classic 20")
dev.off()
#看看基因表达
pal <- viridis(n = 10, option = "D", direction = 1)

#FeaturePlot

pal_red <- c("#999999","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d")

####CK5 Analyze####
HCC_18<- readRDS("./220413HCC_18.rds")
HCC_18$tissue_source <- as.factor(HCC_18$tissue_source)
HCC_18 <- RunTSNE(HCC_18, dims = 1:20, do.fast = T)

Idents(HCC_18) <- 'cell_type'

pdf(paste0("./",sam.name,"/CellCluster_UMAP_celltype1.pdf"),width = 10,height = 8)
DimPlot(HCC_18, pt.size = 1, reduction = "tsne", label = F)+
  scale_color_igv()
dev.off()

m_featureplot <-FeaturePlot_scCustom(seurat_object = HCC_18, features = "SAMSN1", order = F,
                                     reduction = "tsne", 
                                     colors_use = pal_red,split.by ="orig.ident")+labs(fill=" ")
m_featureplot <- m_featureplot & annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)
m_featureplot <- m_featureplot & annotate(geom = 'segment', x = Inf, xend = Inf, color = 'black', y = -Inf, yend = Inf, size = 1)
m_featureplot <- m_featureplot & theme(axis.text = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
m_featureplot

table(HCC_18$tissue_source,HCC_18$Patients)

#                P01  P02  P03  P04  P05  P07  P08  P09  P10  P11  P12  P13  P14  P15  P16  P17
#Adjacent liver    0   58  125  365  100    0    0  124  171    0  332    0    0  289  173  248
#Tumor           119  898 1200  561  629  281 1436 1464 1446 1203  430  476  418 1037 1313  135

#                P18  P19
#Adjacent liver  238   39
#Tumor          1082  108

cell.use <- subset(HCC_18@meta.data,Patients %in% c("P01","P04","P05","P07","P08","P10",
                                                      "P11","P14","P16","P17","P19"))
HCC_sub <- subset(HCC_18, cells = row.names(cell.use))
rm(cell.use)
saveRDS(HCC_18, file = "./220413HCC_18.rds")

HCC <- readRDS("./220413HCC_18.rds")
saveRDS(HCC,"./220621HCC_18.rds") #后续再分析用这个
saveRDS(HCC_nk,"./220621HCC_ILC.rds")
rm(HCC)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(patchwork)
library(scales)
library(ggthemes)
library(viridis)
library(openxlsx)
library(scCustomize)

FeaturePlot(HCC,features = "CLEC12B",reduction = 'tsne')
pdf(file = "./2206/tsne_CLEC12B_Density.pdf",height = 8,width = 9)
Plot_Density_Custom(seurat_object = HCC, features = "CLEC12B",pt.size = 0.75,
                    reduction = 'tsne',viridis_palette="magma")+
  theme_par()
dev.off()

pdf(file = "./2206/tsne_minimal_cluster.pdf",height = 8,width = 10)
DimPlot(HCC, pt.size = 0.75, reduction = "tsne", label = F, group.by = 'cell_type1')+
  scale_color_igv()+
  theme_par()
dev.off()

####230210 for CHD ####
library(Seurat)
library(SeuratObject)
library(dplyr)
library(Matrix)
library(ggplot2)
library(ggsci)
library(patchwork)
library(scales)#调色板包
library(ggthemes)
library(viridis)
library(openxlsx)
library(scCustomize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggunchull)
library(ComplexHeatmap)
library(circlize)
library(scRNAtoolVis)

HCC <- readRDS("./220621HCC_18.rds")
rm(HCC)
pdf(file = "./2206/tsne_celltype.pdf",height = 6,width = 8)
DimPlot(HCC, pt.size = 0.75, reduction = "tsne", label = F, group.by = 'cell_type1')+
  scale_color_igv()+
  theme_par()+
  theme(plot.title = element_blank(),
      axis.title.y = element_text(color = "Black", size = 15, face = "plain"),
      axis.title.x  = element_text(color = "Black", size = 15, face = "plain"),
      legend.title = element_blank(), 
      legend.text = element_text(colour="Black", size=15 ,face="plain"),
      axis.text = element_text(colour="Black", size=15 ,face="plain"))+
  guides(color=guide_legend(override.aes = list(size=5,alpha=1)))
dev.off()

pdf(file = "./2206/tsne_origin.pdf",height = 6,width = 9)
DimPlot(HCC, pt.size = 0.75, reduction = "tsne", label = F)+
  scale_color_igv()+
  theme_par()+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 15, face = "plain"),
        axis.title.x  = element_text(color = "Black", size = 15, face = "plain"),
        legend.title = element_blank(), 
        legend.text = element_text(colour="Black", size=15 ,face="plain"),
        axis.text = element_text(colour="Black", size=15 ,face="plain"))+
  guides(color=guide_legend(override.aes = list(size=5,alpha=1)))
dev.off()

pdf(file = "./2206/tsne_where_is_SAMSN1.pdf",height = 6,width = 6.5)
FeaturePlot(HCC,pt.size = 0.75, reduction = "tsne", label = F,features = "SAMSN1",cols = c("lightgrey","red"))+
  theme_par()+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 15, face = "plain"),
        axis.title.x  = element_text(color = "Black", size = 15, face = "plain"),
       # legend.title = element_blank(), 
        legend.text = element_text(colour="Black", size=15 ,face="plain"),
        axis.text = element_text(colour="Black", size=15 ,face="plain"))
dev.off()

FeaturePlot(HCC,pt.size = 0.75, reduction = "tsne", label = F,split.by = "tissue_source",
            features = "SAMSN1",cols = c("lightgrey","red"))+
  theme_par()+
  theme(plot.title = element_blank(),
        axis.title.y = element_text(color = "Black", size = 15, face = "plain"),
        axis.title.x  = element_text(color = "Black", size = 15, face = "plain"),
        # legend.title = element_blank(), 
        legend.text = element_text(colour="Black", size=15 ,face="plain"),
        axis.text = element_text(colour="Black", size=15 ,face="plain"))



pdf("./2206/CK5PTvsIT.pdf",height = 5,width = 11)
FeatureCornerAxes(object = HCC,reduction = 'tsne',
                  groupFacet = 'tissue_source',
                  relLength = 0.5,
                  relDist = 0.1,
                  features = c("SAMSN1"),
                  show.legend = T,
                  cornerTextSize = 5,
                  base_size = 20
)
dev.off()


color1 <- colorRampPalette(c("#4575b4",'#ffffbf',"#d73027"))(100)
jjDotPlot(object = HCC,
          gene = c("SAMSN1"),
          id = 'cell_type1',
          split.by = 'tissue_source',
          split.by.aesGroup = T,
          ytree = F)+
  theme(axis.title.y = element_text(color = "Black", size = 12, face = "plain"),
        axis.title.x.bottom  = element_blank(),
        legend.title = element_text(color = "Black", size = 12, face = "plain"),
        legend.text = element_text(colour="Black", size=12 ,face="plain"),
        axis.text = element_text(colour="Black", size=12 ,face="plain"),
        axis.text.x = element_text(vjust = 1)
  )+
  RotatedAxis()+
  scale_size(range = c(1,9))+
  scale_fill_gradient2(high = color1[100],low = color1[1]) 
#breaks = c(-2,-1,0,1,2),labels = c(-2,-1,0,1,2),limits = c(-2,2)


a <- AverageExpression(HCC, assays = "RNA",slot = "data",group.by =c("tissue_source","cell_type1"))
a$RNA["SAMSN1",]
# Adjacent liver_T Cells               Adjacent liver_NK          Adjacent liver_B Cells 
#              3.8078647                       2.9788686                       1.9731207 
#     Adjacent liver_Mye Adjacent liver_Parenchyma Cells            Adjacent liver_Tumor 
#              3.0190963                       0.9124489                       0.4297335 
#          Tumor_T Cells                        Tumor_NK                   Tumor_B Cells 
#              4.5858959                       3.7964177                       2.7326824 
#              Tumor_Mye          Tumor_Parenchyma Cells                     Tumor_Tumor 
#              2.8797813                       0.8912194                       0.2741002 

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
