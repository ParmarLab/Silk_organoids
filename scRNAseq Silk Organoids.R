################################################################################
# Silk scaffolding drives self-assembly of functional and mature human brain organoids
################################################################################

################################################################################
# Data preparation 
################################################################################

# Packages
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library("escape")
library(cowplot)
library("Nebulosa")
library(voxhunt)

Organoids <- readRDS("Single-cell RNA seq/Organoids.rds")

################################################################################
# Figure 4 
################################################################################


# Figure 4 A and B
################################################################################
set.seed(2)
Organoids@meta.data$cell_id <- rownames(Organoids@meta.data)
new_df <- Organoids@meta.data %>% group_by(dataset) %>% sample_n(4462)
Organoids.randomsubset <- subset(Organoids,cells=new_df$cell_id)
UMAPcols <- c("#EB8C44","#80C3C6","#4391BA","#122E7B","#CF6EAB","#60A762","#255838")
figures4AandB <- DimPlot(Organoids.randomsubset, reduction = "umap", group.by="OrderedIdents",
                         split.by = "dataset", label = FALSE, 
                         pt.size = .1, cols=UMAPcols)
figures4AandB

# Figure 4 C
################################################################################
n_cells <- FetchData(Organoids, vars = c("dataset", "OrderedIdents")) %>%
  dplyr::count(dataset,OrderedIdents) %>%
  tidyr::spread(OrderedIdents, n)
n_cells_melt <- reshape2::melt(n_cells)
n_cells_melt$cellType <- n_cells_melt$variable
order <- c("Neural stem cell","Newborn neuron","Inhibitory neuron",
           "Excitatory neuron","Astroglia","OPC","Mural")
n_cells_melt$cellType <- factor(n_cells_melt$cellType, levels = order)
figure4C<-ggplot(n_cells_melt, aes(x=dataset, y=value, fill=cellType)) + 
  geom_bar(stat="identity", position = "fill") +
  scale_fill_manual(values=c("#EB8C44","#80C3C6","#4391BA","#122E7B",
                             "#CF6EAB","#60A762","#255838")) +
  coord_flip()
figure4C

# Figure 4 D
################################################################################
figure4D <- DotPlot(Organoids,
           features = c("TOP2A","CDK1","CENPF","PTTG1","ASCL1","PAX6","SOX2",
                        "GAP43","MAP2","NNAT","DCX","STMN2","SYT1","NSG2","SNAP25",
                        "GAD1","GABRA2","NEFL","SLC32A1","SST","ISL1","GRIN2B",
                        "BCL11B","SLC1A2","SIX3","NEUROD6","GFAP","AQP4","SOX9",
                        "SLC1A3","OLIG1","SOX10","PLP1","MPZ","COL1A1","PDGFRA",
                        "PDGFRB","IFITM1"), cols = c("white","#245D70"),group.by="OrderedIdents",
           dot.scale = 7)
figure4D

# Figure 4 E
################################################################################
NeuronDifferentiation.genes <-getGeneSets(library = "C5",
                                         gene.sets ="GOBP_NEURON_DIFFERENTIATION")
ES.seurat <- enrichIt(obj = Organoids, 
                      gene.sets = NeuronDifferentiation.genes, 
                      groups = 1000, cores = 2)
Organoids <- Seurat::AddMetaData(Organoids,ES.seurat)
figure4E <- ggplot(Organoids@meta.data,aes(x=dataset,
                               y=GOBP_NEURON_DIFFERENTIATION,
                               fill=dataset)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  scale_fill_manual(values = c("#7AAABA","#EECF82")) + 
  coord_flip()+
  NoLegend()+xlab("")
figure4E
summary(lm(GOBP_NEURON_DIFFERENTIATION~dataset,Organoids@meta.data))

# Figure 4 F
################################################################################
Organoids.matureneurons<-subset(Organoids, subset=seurat_clusters == c(2,3) )
NeuronMaturationGOterms <- c("GOBP_NEURON_MATURATION",
                             "GOBP_NEURON_MIGRATION",
                             "GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION",
                             "GOBP_POSTSYNAPSE_ASSEMBLY",
                             "GOBP_NEUROTRANSMITTER_SECRETION")
NeuronMaturation.genes <-getGeneSets(library = "C5",
                                    gene.sets = NeuronMaturationGOterms)
ES.seurat <- enrichIt(obj = Organoids.matureneurons, 
                      gene.sets = NeuronMaturation.genes, 
                      groups = 1000, cores = 2)
Organoids.matureneurons <- Seurat::AddMetaData(Organoids.matureneurons,ES.seurat)
figure4F <- ggplot(Organoids.matureneurons@meta.data,aes(x=dataset,
                                           y=GOBP_NEUROTRANSMITTER_SECRETION,
                                           fill=dataset)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  scale_fill_manual(values = c("#7AAABA","#EECF82")) + 
  coord_flip()+
  NoLegend()+xlab("")
figure4F
summary(lm(GOBP_NEUROTRANSMITTER_SECRETION~dataset,Organoids.matureneurons@meta.data))

# Figure 4 G
################################################################################
Organoids.matureneurons<-subset(Organoids, subset=seurat_clusters == c(2,3) )
figure4G <- VlnPlot(object=Organoids.matureneurons, 
        features = c("SYT1","NNAT","VAMP2","GABARAPL2",
                     "NSG2","ELAVL3","TUBB2B","PEBP1"), 
        group.by="dataset", pt.size=0,col=c("#7AAABA","#EECF82"),ncol=2)
figure4G

################################################################################
# Figure 5 
################################################################################

# Figure 5 A
################################################################################
StressResponse4A.genes <-getGeneSets(library = "C5",
                       gene.sets = c("GOBP_NECROPTOTIC_SIGNALING_PATHWAY",
                                     "GOBP_POSITIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS"))
ES.seurat <- enrichIt(obj = Organoids, 
                      gene.sets = StressResponse4A.genes, 
                      groups = 1000, cores = 2)
Organoids <- Seurat::AddMetaData(Organoids,ES.seurat)

StressResponse4B.genes <-getGeneSets(library = "H",
                                   gene.sets = c("HALLMARK_APOPTOSIS",
                                                 "HALLMARK_HYPOXIA"))
ES.seurat <- enrichIt(obj = Organoids, 
                      gene.sets = StressResponse4B.genes, 
                      groups = 1000, cores = 2)
Organoids <- Seurat::AddMetaData(Organoids,ES.seurat)

figure5A <- ggplot(Organoids@meta.data,aes(x=dataset,
                                           y=GOBP_POSITIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS,
                                           fill=dataset)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  scale_fill_manual(values = c("#7AAABA","#EECF82")) + 
  NoLegend()+xlab("")
figure5A
summary(lm(GOBP_POSITIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS~dataset,Organoids@meta.data))

# Figure 5 B
################################################################################
figure5B <- DotPlot(Organoids, features = c("TP53","CTSL","CASP7","GPX8","GSS",
                                            "CTSB","CAT","CAPN6",	"XIAP", "NAIP"), 
                    group.by="dataset",
                    cols = c("#D9E8ED","#245D70"),
                    dot.scale = 10) + coord_flip()
figure5B

################################################################################
# Supplementary Figure 3 
################################################################################

# Supplementary Figure 3 A                    
################################################################################
KantonOrg.seu <- readRDS("Single-cell RNA seq/KantonOrganoids.rds")

#Select only cells from Organoids at 4 months
KantonOrg.seu.4M<-subset(KantonOrg.seu, subset = Stage == c("Organoid-4M") )

#Select only clusters with 6 or more cells
Idents(KantonOrg.seu.4M)<-KantonOrg.seu.4M$PredCellType
KantonOrg.seu.4M.clean<-subset(KantonOrg.seu.4M, idents =  c("Astrocyte",
                                                             "Choroid",
                                                             "EN",
                                                             "Glyc",
                                                             "IN",
                                                             "IPC",
                                                             "Mural",
                                                             "OPC",
                                                             "RG"))

anchors <- FindTransferAnchors(reference = KantonOrg.seu.4M.clean,
                               query = Organoids, dims = 1:30)
predictions <- TransferData(anchorset = anchors, 
                            refdata = KantonOrg.seu.4M.clean$PredCellType, 
                            dims = 1:30) 
Organoids <- AddMetaData(Organoids, metadata = predictions)
Organoids$predicted.id  <- factor(Organoids$predicted.id, 
                                  levels = c("IPC","IN","Glyc","EN","RG",
                                             "Astrocyte","OPC","Mural", "Choroid"))
suppfigure3A<-DimPlot(Organoids,
                      group.by = "predicted.id",
                      label=FALSE, 
                      cols=c("#EB8C44","#80C3C6","#4391BA","#122E7B","#CF6EAB",
                             "#E8312D","#60A762","#255838","#9923A3")) 
  
suppfigure3A

# Supplementary Figure 3 B
################################################################################
colfunc <- colorRampPalette(c("white", "#122E7B"))
gradient<-colfunc(16)
suppfigure3B <- FeaturePlot(Organoids,features = c("prediction.score.max"), 
                            cols = gradient)
suppfigure3B

# Supplementary Figure 3 C
################################################################################
colfunc <- colorRampPalette(c("white", "#122E7B"))
gradient<-colfunc(16)

table1<-table(Organoids$predicted.id,Organoids$OrderedIdents)
correlationtable <- table1[, c(6, 1, 3, 4, 2, 7, 5)]
suppfigure3C <- heatmap(table1, Rowv = NA, Colv=NA, col = gradient )
suppfigure3C

# Supplementary Figure 3 D
################################################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Organoids <- CellCycleScoring(Organoids, s.features = s.genes, 
                              g2m.features = g2m.genes, set.ident = TRUE)
suppfigure3D <- DimPlot(Organoids)
suppfigure3D

# Supplementary Figure 3 E
################################################################################
densityplot.genes <- c("ASCL1","MKI67","CLSPN","STMN2","GRIA2","DLX5","SYN1","DLG4",
                       "NCAM1","GAD1","GABRA1","SST","SLC17A7","GNB1","SIX3","AQP4",
                       "GFAP","RFX4","OLIG2","S100B","MBP","COL1A2","FBLN1","S100A11")
suppfigure3E <- plot_density(Organoids,densityplot.genes) & 
  scale_colour_gradient(low="#A7D8D8",high="#10285E")
suppfigure3E

# Supplementary Figure 3 J                                                  
################################################################################
FiorenzanoOrg <- readRDS("Single-cell RNA seq/FiorenzanoOrganoids.rds")
FiorenzanoOrg.day120 <- subset(FiorenzanoOrg, subset = groups==c("day120"))
FiorenzanoOrg.day120.VLMC <- subset(FiorenzanoOrg.day120, 
                                    subset = celltypes2==c("VLMC"))
KantonOrg.seu.mural <-subset(KantonOrg.seu, subset = PredCellType == c("Mural"))
Organoids.mural <- subset(Organoids, subset = seurat_clusters ==c("4"))

ComparisonMuralVLMC <- merge(Organoids.mural, y = c(KantonOrg.seu.mural, FiorenzanoOrg.day120.VLMC), 
                             add.cell.ids = c("Sozzi", "Kanton", "Fiorenzano_VLMC"),
                             project = "ComparisonMuralVLMC")
ComparisonMuralVLMC <- NormalizeData(ComparisonMuralVLMC)
ComparisonMuralVLMC <- FindVariableFeatures(ComparisonMuralVLMC, 
                                            selection.method = "vst",
                                            nfeatures = 2000)
ComparisonMuralVLMC <- ScaleData(ComparisonMuralVLMC)
ComparisonMuralVLMC <- RunPCA(ComparisonMuralVLMC, 
                              features = VariableFeatures(object = ComparisonMuralVLMC))
ComparisonMuralVLMC$dataset <- ifelse(ComparisonMuralVLMC$orig.ident=="Silk1","Mural (Sozzi et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="Silk2","Mural (Sozzi et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="Non-Silk","Mural (Sozzi et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="barbara1","Mural (Kanton et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="1standardorgday120","VLMC (Fiorenzano et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="2standardorgday120","VLMC (Fiorenzano et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="3standardorgday120","VLMC (Fiorenzano et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="4standardorgday120","VLMC (Fiorenzano et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="5standardorgday120","VLMC (Fiorenzano et al.)",
                               ifelse(ComparisonMuralVLMC$orig.ident=="6standardorgday120","VLMC (Fiorenzano et al.)",NA))))))))))
ComparisonMuralVLMC <- RunHarmony(ComparisonMuralVLMC, "dataset", 
                                  plot_convergence = FALSE)
ComparisonMuralVLMC <- ComparisonMuralVLMC %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 0.01) %>% 
  identity()

suppfigure3J <- DimPlot(ComparisonMuralVLMC, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .1,
        group.by="dataset",cols=c("#EB8C44","#CF6EAB","#4391BA","red"))
suppfigure3J

# Supplementary Figure 3 K
################################################################################
suppfigure3K <- VlnPlot(ComparisonMuralVLMC, features = c("FOXP1","FOXC1","CD248",
                                                          "CREB3L1","FABP7","SOX2",
                                                          "SOX6","TUBB2B"),
                        group.by="dataset", 
                        pt.size=0,
                        col=c("#EB8C44","#CF6EAB","#4391BA"), ncol=2)
suppfigure3K
  
# Supplementary Figure 3 L
################################################################################
load_aba_data("Single-cell RNA seq/voxhunt_rds/")
genes_use <- variable_genes('E13', 300)$gene

SilkOrg <- subset(Organoids, subset=dataset == "Silk" )
NonSilkOrg <- subset(Organoids, subset=dataset == "Non-Silk" )

vox_mapSilk <- voxel_map(SilkOrg, genes_use=genes_use)
vox_mapNonSilk <- voxel_map(NonSilkOrg, genes_use=genes_use)

suppfigure3L <- plot_map(vox_mapSilk, scale=TRUE, guide = "none", groups = "Silk") +
               plot_map(vox_mapNonSilk, scale=TRUE, guide = "none", groups = "Non-Silk")
suppfigure3L

# Supplementary Figure 3 M
################################################################################
suppfigure3M <- VlnPlot(Organoids, features = c("MAP2","DCX","STMN2","SYT1","NNAT",
                                                "SNAP25","NCAM1","PAK3","COL1A1",
                                                "DCN","COL1A2","FBLN1"),
                        group.by="dataset", 
                        pt.size=0,
                        col=c("#7AAABA","#EECF82"))
suppfigure3M

# Supplementary Figure 3 N
################################################################################
Organoids.matureneurons<-subset(Organoids, subset=seurat_clusters == c(2,3))
suppfigure3N <- DotPlot(Organoids.matureneurons,
                        features = c("SCN3B","SCN3A","SCN2A","KCNK3","KCNB1",
                                     "KCNQ2","CACNG2","CACNG3","CACNG4","CACNA1B"), 
                        group.by="dataset",
                        cols = c("#D9E8ED","#245D70"),
                        dot.scale = 10) + coord_flip()
suppfigure3N

################################################################################
# Supplementary Figure 4 
################################################################################

# Supplementary Figure 4 A
################################################################################
StressResponseS4A.genes <-getGeneSets(library = "C5",
                                  gene.sets =c("GOBP_APOPTOTIC_CELL_CLEARANCE",
                                               "GOBP_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS"))
ES.seurat <- enrichIt(obj = Organoids, 
                      gene.sets = StressResponseS4A.genes, 
                      groups = 1000, cores = 2)
Organoids <- Seurat::AddMetaData(Organoids,ES.seurat)
suppfigure4A <- ggplot(Organoids@meta.data,aes(x=dataset,
                                           y=GOBP_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS,
                                           fill=dataset)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  scale_fill_manual(values = c("#7AAABA","#EECF82")) + 
  coord_flip()+
  NoLegend()+xlab("")
suppfigure4A
summary(lm(GOBP_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS~dataset,Organoids@meta.data))

sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] voxhunt_1.0.1      Nebulosa_1.2.0     patchwork_1.1.1    cowplot_1.1.1      escape_1.2.0       ggplot2_3.3.5     
#  [7] harmony_0.1.0      Rcpp_1.0.7         dplyr_1.0.7        SeuratObject_4.0.4 Seurat_4.0.6      

rm(list = ls())
