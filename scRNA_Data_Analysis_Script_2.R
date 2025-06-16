# Loading packages and objects --------------------------------------------
library(tidyverse)
library(dplyr)
library(tidyr)
library(UCell)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggrepel)
library(UpSetR)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(Seurat)
library(babelgene)
library(clusterProfiler)
library(dendextend)
library(clustree)
library(DESeq2)
library("org.Mm.eg.db")
library(rstatix)
library(ggpubr)

appendlist <- function(x,y,z){
  result <- list()
  length(result) <- 7
  names(result) <- names(x)
  for (i in names(result)){
    test1 <- unique(c(x[[i]],y[[i]],z[[i]]))
    result[[i]] <- test1
  }
  return(result)
}
intersectlist <- function(x,y){
  result <- list()
  length(result) <- 7
  names(result) <- names(x)
  for (i in names(result)){
    test1 <- intersect(x[[i]],y[[i]])
    result[[i]] <- test1
  }
  return(result)
}
setdifflist <- function(x,y){
  result <- list()
  length(result) <- 7
  names(result) <- names(x)
  for (i in  names(result)){
    test1 <- setdiff(x[[i]],y[[i]])
    result[[i]] <- test1
  }
  return(result)
}

"Data from GEO GSE282887"
VH1H2_late_pos <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/New_VHKO_VEKO_VHEKO_Late_Data.rds")

"Data from GEO GSE253168"
Vhl_Pax8_pos_batch_sex<- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Old_VKO_ConKO_Early_Late_Data.rds")

col_32 <- unique(c(brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"), brewer.pal(8, "Accent"),brewer.pal(8, "Set2") ))
celltypes_colour_2 <- c("PT S1" = "#08519c" ,"PT S2" = "#6baed6","PT S3" = "#bdd7e7","LoH" = "#bae4b3", "DCT" = "#74c476","CDPC" ="#31a354" ,"CDIC" = "#006d2c", "Podocyte" = "#dadaeb","PEC" = "#bcbddc","Fibroblast" = "#9e9ac8", "Endothelial" = "#807dba","VSMC" = "#6a51a3","Pericyte" = "#4a1486","Macrophage" = "#fcbba1","Neutrophil" = "#fc9272","NK Cell" = "#fb6a4a","Granulocyte" = "#ef3b2c","Monocyte" = "#cb181d","B Cell" = "#a50f15","T Cell"="#67000d")
sample_colors <- c("#9ecae1","#6baed6","#3182bd","#08519c","#fc9272","#fb6a4a","#de2d26","#a50f15","#bf8548","#7c3e04","#58300d","#41280f","#a1d99b","#74c476","#31a354","#006d2c","#bcbddc","#9e9ac8","#756bb1","#54278f")
Isoform_spec_colors <- c("Dependent on HIF1A alone" = "#bc8e54","Dependent on HIF2A alone" = "#208d4c","Dependent on both isoforms"= "#a50f15","Dependent on either isoform" = "#9467cf","Dependent on neither isoform" = "skyblue2","Ambiguous HIF dependence" = "gray")
condition_colour_2 <-  c("ConKO"="#595FC7","VKO"="#DA5E24", "VHKO" = "#bc8e54", "VEKO" = "#208d4c","VHEKO" = "#9467cf")
condition_colour <-  c("ConKO Early"="#93A0FF","VKO Early"="#F29163","ConKO Late"="#595FC7","VKO Late"="#DA5E24")
cellidentities_colour_3 <- c("PT S1" = "#08519c" ,"PT S2" = "#6baed6","PT S3" = "#bdd7e7","PT Class A" = "#444BF5", "PT Class B" = "#80f3FC","PT All" = "#3bb5ff","PT" = "#1b85b8", "CDIC" = "#006d2c","Non PT" = "#559e83", "PT like" = "#000000", "Mixed" = c("#AAAAAA"))
ptdentities_colour_4 <- c("S1 A"="#bf8548","S1 B"="#7c3e04","S2 A"="#74c476","S2 B"="#31a354","S3 A"="#9e9ac8","S3 B"="#756bb1","Non PT"="darkgrey")
sample_colors <- c("#9ecae1","#6baed6","#3182bd","#08519c","#fc9272","#fb6a4a","#de2d26","#a50f15","#bf8548","#7c3e04","#58300d","#41280f","#a1d99b","#74c476","#31a354","#006d2c","#bcbddc","#9e9ac8","#756bb1","#54278f")
Isoform_spec_colors <- c("HIF1A alone" = "#bc8e54","HIF2A alone" = "#208d4c","HIF1A or HIF2A"= "#a50f15","HIF1A + HIF2A" = "#9467cf","Dependent on neither isoform" = "skyblue2","Ambiguous HIF dependence" = "gray")


cluster_colour <- c(col_32, brewer.pal(n = 12,name = "Paired"))
names(cluster_colour) <- c(as.character(c(0:43)))

scale_fill_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)))
scale_fill_gradientn(colors = colorRampPalette(rev(c("#67001F","red","#F7F7F7","#41B6C4", "#081D58")))(100),values = scales::rescale(c(-1,-0.5,0,0.5,1)),limits = c(-1,1))


# Merging old and new data ------------------------------------------------
saveRDS(VH1H2_late_pos,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/New_VHKO_VEKO_VHEKO_Late_Data.rds")
saveRDS(Vhl_Pax8_pos_batch_sex,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Old_VKO_ConKO_Early_Late_Data.rds")

V <- subset(Vhl_Pax8_pos_batch_sex,cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$sample %in% c("OY622_pos", "OY615_pos", "OY1023_pos", "OY552_pos", "OY611_pos", "OY545_pos", "OY606_pos", "OY626_pos")])

VHE_late_pos <- merge(VH1H2_late_pos,V)

# Metadata Annotation  -----------------------------------------------------
#Adding sample as factor

VHE_late_pos[["Sample"]] <- factor(VHE_late_pos$sample, levels = c("OY552_pos", "OY611_pos","OY545_pos","OY622_pos",  "OY606_pos","OY626_pos", "OY615_pos", "OY1023_pos","JH186_pos", "JH187_pos","JH437_pos","JH444_pos", "JP120_pos", "JP122_pos", "JP123_pos", "JP165_pos","OP415_pos","OP404_pos","OP408_pos","OP437_pos"))

#Adding sex as factor
Idents(VHE_late_pos) <- "Sample"
Female <- WhichCells(VHE_late_pos, idents = c("OY552_pos", "OY611_pos", "OY545_pos", "OY606_pos", "OY626_pos","JH437_pos","JH444_pos","JP123_pos", "JP165_pos","OP415_pos"))
Male <- WhichCells(VHE_late_pos, idents = c("OY622_pos", "OY615_pos", "OY1023_pos", "JH186_pos", "JH187_pos", "JP120_pos", "JP122_pos","OP404_pos","OP408_pos","OP437_pos"))
Idents(VHE_late_pos, cells = Female) <- "Female"
Idents(VHE_late_pos, cells = Male) <- "Male"
VHE_late_pos[["Sex"]] <- factor(Idents(VHE_late_pos), levels = c("Male", "Female"))

#Adding Condition as factor
Idents(VHE_late_pos) <- "Sample"

Con_P <- WhichCells(VHE_late_pos, idents = c("OY552_pos", "OY611_pos", "OY545_pos", "OY622_pos"))
Vhl_P <- WhichCells(VHE_late_pos, idents = c("OY606_pos", "OY626_pos", "OY615_pos", "OY1023_pos"))
Hif1_P <- WhichCells(VHE_late_pos, idents = c("JH186_pos", "JH187_pos","JH437_pos","JH444_pos"))
Hif2_P <- WhichCells(VHE_late_pos, idents = c("JP120_pos", "JP122_pos", "JP123_pos", "JP165_pos"))
Hif12_P <- WhichCells(VHE_late_pos, idents = c("OP415_pos", "OP404_pos", "OP408_pos", "OP437_pos"))

Idents(VHE_late_pos, cells = Con_P) <- "ConKO"
Idents(VHE_late_pos, cells = Vhl_P) <- "VKO"
Idents(VHE_late_pos, cells = Hif1_P) <- "VHKO"
Idents(VHE_late_pos, cells = Hif2_P) <- "VEKO"
Idents(VHE_late_pos, cells = Hif12_P) <- "VHEKO"
VHE_late_pos[["Condition"]] <- factor(Idents(VHE_late_pos), levels = c("ConKO","VKO","VHKO","VEKO","VHEKO"))

# Assigning Cell Type  -----------------------------------------------------
"Deriving consensus lists of cell type specific marker genes"

"Getting cell type marker modules - Supplementary Table 3"
Modules <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_Type_Marker_list_male_female.csv",col_names = TRUE)

"Removing blank cell values"
remove_na <- function(x){x <- x[!is.na(x)]}
Modules <- as.list(Modules)
Modules <- lapply(Modules, remove_na)
remove_dup <- function(x){x <- x[!duplicated(x)]}
Modules <- lapply(Modules, remove_dup)
"Getting list of detected genes in object - all cells together"
All_genes <- VHE_late_pos@assays$RNA@data
Sum_genes <- rowSums(All_genes,dims = 1)
Expressed_genes <- Sum_genes[Sum_genes > 0]
Reference_genes_all <- names(Expressed_genes)
"Removing undetected genes from modules list before scoring"
removenoexpall <- function(x){x <- x[x %in% Reference_genes_all]}

"ORIGINAL MODULES"
Modules_all <- lapply(Modules, removenoexpall)


"Getting cell type marker modules"
Metamarkers <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Mouse_Kidney_Atlas_Meta_Markers.csv",col_names = TRUE)

Metamarkers_list <- list()
for(i in unique(Metamarkers$`Cell Population`)){
  test1 <- filter(Metamarkers, AUROC >= 0.8, `Fold Change Detection Rate` > 2, `Cell Population` == i)
  genes1 <- as.character(test1$`Gene symbol`)
  test2 <- filter(Metamarkers, AUROC >= 0.6, `Fold Change Detection Rate` > 3, `Cell Population` == i)
  genes2 <- as.character(test2$`Gene symbol`)
  test3 <- filter(Metamarkers, AUROC >= 0.5, `Fold Change Detection Rate` > 6, `Cell Population` == i)
  genes3 <- as.character(test3$`Gene symbol`)
  test4 <- unique(c(genes1,genes2,genes3))
  test5 <- arrange(filter(Metamarkers, `Cell Population` == i, `Gene symbol` %in% test4),Rank)$`Gene symbol`
  Metamarkers_list[[i]] <- unique(c(test5))}


"Mouse Kidney Atlas Metamarkers"
Metamarkers_all <- lapply(Metamarkers_list, removenoexpall)

"SYNERGISING METAMARKERS AND PREVOUS MODULES"

"Will not use metamarkers for S2 and S3 as they do not consider sex differences"
Metamarkers_previous_noS2S3 <- Metamarkers_all[c(3,5,8,13,14,18,19,20,23,24,25,26,27,28,29,33)]
Modules_S2S3 <- Modules_all[c(3,4,5,6)]

"Fusing ICA and ICB cell type marker genes so taking genes common to both subtytpes of CDIC"
Metamarkers_CDIC <- list("IC" = c(Reduce(intersect,x = Metamarkers_all[c(16,17)])))
Metamarkers_combine <- c(Metamarkers_previous_noS2S3,Metamarkers_CDIC, Modules_S2S3)

"Granuclocyte and VSMC markers are not present in the metamarkers - fine"

"MODULES VS METAMARKERS"
"Removing granulocyte and VSMC markers"
Modules_combine <- Modules_all[c(1:13,15:18,20:22)]
Metamarkers_rearrange <- Metamarkers_combine[c(12,15,18:21,6,3,11,17,14,13,4,5,7,10,9,8,1,16)]

"Retaining only those module genes that have passed metamarker filtered. Genes arranged by metamarker rank (not applicable to S2 and S3)"
Meta_modules <- Map(f = intersect,Metamarkers_rearrange,Modules_combine)

"Want to eliminate genes that are markers for multiple cell types. Note that S2/S3 male markers will never be compared to S2/S3 female markers so can retain common markers between the sexes"

Meta_modules_unique <- list()
for (i in seq_along(Meta_modules)){
  if(i == 3){ #S2 male
    test1 <- c(Meta_modules[i],Meta_modules[-c(i,5)]) #exclude S2 female and so on
    test2 <- Reduce(setdiff,test1)
    test3 <- names(Meta_modules)[i]
    Meta_modules_unique[[test3]] <- test2
  }
  else if(i == 4){
    test1 <- c(Meta_modules[i],Meta_modules[-c(i,6)])
    test2 <- Reduce(setdiff,test1)
    test3 <- names(Meta_modules)[i]
    Meta_modules_unique[[test3]] <- test2
  }
  else if(i == 5){
    test1 <- c(Meta_modules[i],Meta_modules[-c(i,3)])
    test2 <- Reduce(setdiff,test1)
    test3 <- names(Meta_modules)[i]
    Meta_modules_unique[[test3]] <- test2
  }
  else if(i == 6){
    test1 <- c(Meta_modules[i],Meta_modules[-c(i,4)])
    test2 <- Reduce(setdiff,test1)
    test3 <- names(Meta_modules)[i]
    Meta_modules_unique[[test3]] <- test2
  }
  else {
    test1 <- c(Meta_modules[i],Meta_modules[-i])
    test2 <- Reduce(setdiff,test1)
    test3 <- names(Meta_modules)[i]
    Meta_modules_unique[[test3]] <- test2
  }
}


max_length <- max(sapply(Meta_modules_unique, length))

# Pad each entry in the list to make them the same size
padded_list <- lapply(Meta_modules_unique, function(entry) {
  length_diff <- max_length - length(entry)
  c(entry, rep(NA, length_diff))
})

# Convert the list to a data frame if needed
Meta_modules_unique_df <- as.data.frame(padded_list)

write.table(x = Meta_modules_unique_df,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_type_modules_filtered_metamarkers.csv",sep = ",",row.names = FALSE,col.names = TRUE)

Meta_modules_unique_df <- read_csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_type_modules_filtered_metamarkers.csv")
"`````````````````````````````````````````````````````````````````````````````````"
"-----------------SCORING CELL TYPES---------------------"

cell_types <- names(Meta_modules_unique) #metamarkers that are also modules and filtered for unique genes
cell_types_meta <- paste(cell_types,"MetaUCell",sep = "_")

VHE_late_pos<- AddModuleScore_UCell(VHE_late_pos, features = Meta_modules_unique, name = "_MetaUCell",assay = "RNA",maxRank = 1500,ncores = 10,slot = "data",force.gc = TRUE)

# FINDING MAXIMUM SCORES
Cell_type_scores <- FetchData(VHE_late_pos, c("Sex", cell_types_meta))
cell_types_male <- cell_types_meta[c(1:4, 7:20)]
cell_types_female <- cell_types_meta[c(1:2, 5:20)]

# Use vectorized operations to find max
get_max <- function(y) {
  if (y["Sex"] == "Male") {
    max_indices <- which.max(y[cell_types_male])
    max_cell_type <- cell_types_male[max_indices]
  } else {
    max_indices <- which.max(y[cell_types_female])
    max_cell_type <- cell_types_female[max_indices]
  }
  return(max_cell_type)
}

# Filling in empty dataframe with max score cell type for each cell
Cells_celltype <- as.data.frame(matrix(nrow = nrow(Cell_type_scores), ncol = 1))
rownames(Cells_celltype) <- rownames(Cell_type_scores)
colnames(Cells_celltype) <- "Cellscore"
Cells_celltype[, 1] <- apply(Cell_type_scores, 1, get_max)

# Converting from module names to cell type names.
celltypenames <- data.frame(
  Cellscore = cell_types_meta,
  Cell_type = c("PEC", "PT S1", "PT S2", "PT S3", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte", "Pericyte", "Endothelial","Fibroblast", "Macrophage", "NK Cell", "Neutrophil", "Monocyte", "B Cell", "T Cell")
)

Cells_celltype_2 <- rownames_to_column(Cells_celltype,var = "Cell")
Cells_celltype_2 <- left_join(Cells_celltype_2, celltypenames, by = "Cellscore")

Cells_celltype_2$Cell_type <- factor(
  Cells_celltype_2$Cell_type,
  levels = c("PT S1", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte", "PEC", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "Neutrophil", "NK Cell", "Monocyte", "B Cell", "T Cell")
)

Cell_type_MetaUCell_idents <- Cells_celltype_2
rownames(Cell_type_MetaUCell_idents) <- Cell_type_MetaUCell_idents$Cell 

Idents(VHE_late_pos) <- Cell_type_MetaUCell_idents[,3]
VHE_late_pos[["Cell_types_MetaUCell"]] <- factor(Idents(VHE_late_pos), levels = c("PT S1", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte","PEC", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "Neutrophil", "NK Cell", "Monocyte", "B Cell", "T Cell"))


# Assigning PT Class ------------------------------------------------
"Data from Supplementary Table 1"

PT_modules_list <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Control_PT_Class_Module_A_B_correlations.rds")

Module_A_ordered <- PT_modules_list[[1]]
Module_B_ordered <- PT_modules_list[[2]]

PT_Modules <- list(Module_A_ordered, Module_B_ordered)

DefaultAssay(VHE_late_pos) <- "RNA"
VHE_late_pos <- AddModuleScore(VHE_late_pos, features = PT_Modules, name = "PT_Module_", nbin = 100, ctrl = 50)

colnames(VHE_late_pos@meta.data)[colnames(VHE_late_pos@meta.data) %in% c("PT_Module_1","PT_Module_2")] <- c("PT_Module_A","PT_Module_B")

Tub_A <- WhichCells(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3") & VHE_late_pos$PT_Module_A > 0.125])
Tub_B <- WhichCells(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3") & VHE_late_pos$PT_Module_A < 0.125])
Non_tub <- WhichCells(VHE_late_pos, cells = colnames(VHE_late_pos)[!VHE_late_pos$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3")])

Idents(VHE_late_pos, cells = Tub_A) <- "PT Class A"
Idents(VHE_late_pos, cells = Tub_B) <- "PT Class B"
Idents(VHE_late_pos, cells = Non_tub) <- "Non PT"

VHE_late_pos[["PT_Class_UCell"]] <- factor(Idents(VHE_late_pos), levels = c("PT Class A", "PT Class B", "Non PT"))

ptidents <- FetchData(VHE_late_pos,c("Cell_types_MetaUCell","PT_Class_UCell"))
ptidents[["PT_identity"]] <- apply(ptidents,1,function(x){
  if(x[["Cell_types_MetaUCell"]] %in% c("PT S1","PT S2","PT S3")){paste(x[["Cell_types_MetaUCell"]],x[["PT_Class_UCell"]])}
  else{"Non PT"}
})
ptidents$PT_identity <- factor(ptidents$PT_identity,levels = c("PT S1 PT Class A","PT S1 PT Class B","PT S2 PT Class A", "PT S2 PT Class B","PT S3 PT Class A", "PT S3 PT Class B","Non PT"),labels = c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B", "Non PT"))
VHE_late_pos[["PT_identity"]] <- ptidents[,3]

saveRDS(VHE_late_pos,file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_VHKO_VEKO_VHEKO_Late_pos.rds")

# VHKO, VEKO, VHEKO - UMAP based analysis-------------------------------------------------

table(VHE_late_pos$Sample, VHE_late_pos$Sex)

"Downsample to 4500 cells per sample for UMAP plotting"
set.seed(1234)
celllist_HE <- list()
for(i in levels(VHE_late_pos$Sample)){
  if(length(grep(pattern = "OY",x = i)) == 0){
     test <- sample(colnames(VHE_late_pos)[VHE_late_pos$Sample == i],4500,replace = FALSE)
  }
  else{test <- c()}
  celllist_HE[[i]] <- test
}

cells2 <- unlist(celllist_HE)
HE_late_pos_equal <- subset(VHE_late_pos,cells = cells2)
#seurat batch correction

batch_list=SplitObject(HE_late_pos_equal, split.by = "Sex")

batch_list <- lapply(X = batch_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
gc()
features <- SelectIntegrationFeatures(object.list = batch_list)
batch_list <- lapply(X = batch_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE,vars.to.regress=c("nFeature_RNA", "nCount_RNA", "percent.mito"))
  x <- RunPCA(x, features = features, verbose = FALSE)
})
gc()
anchors <- FindIntegrationAnchors(object.list = batch_list, reduction = "rpca",dims = 1:30)
rm(batch_list)
gc()
HE_late_pos_equal <- IntegrateData(anchorset = anchors, dims = 1:30)
rm(anchors)
#sex.rpca.male_ref.integrated <- ScaleData(sex.rpca.male_ref.integrated)
HE_late_pos_equal <- ScaleData(HE_late_pos_equal, vars.to.regress=c("nFeature_RNA", "nCount_RNA", "percent.mito"),verbose = FALSE)
HE_late_pos_equal <- RunPCA(HE_late_pos_equal)
HE_late_pos_equal <- RunUMAP(HE_late_pos_equal, dims = 1:30)
gc()
num.pcs=30
resolution <- c(1)
HE_late_pos_equal <- FindNeighbors(object=HE_late_pos_equal, dims=1:num.pcs)
HE_late_pos_equal <- FindClusters(object=HE_late_pos_equal, resolution=resolution)

saveRDS(HE_late_pos_equal,file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/New_VHKO_VEKO_VHEKO_equal_samples_integrated.rds")

gc()



# Supplementary Fig. 4a ---------------------------------------------------

HE_late_pos_equal<-readRDS(file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/New_VHKO_VEKO_VHEKO_equal_samples_integrated.rds")

umap <- FetchData(HE_late_pos_equal,vars = c("umap_1","umap_2","Condition","Sex","Sample"))
sample2 <- c()
for(i in unique(umap$Condition)){
  for (j in unique(umap$Sex)){
    test <- unique(filter(umap,Condition == i, Sex == j)$Sample)
    test2 <- c(1:length(test))
    names(test2) <- test
    sample2 <- c(sample2,test2)
  }
}
sample2 <- as.data.frame(sample2)%>%rownames_to_column(var = "Sample")
umap <- full_join(umap,sample2,by = "Sample")
umap[["Sample2"]] <- paste(umap$Condition,umap$Sex,umap$sample2)

sample_colors <- c("#9ecae1","#6baed6","#3182bd","#08519c","#fc9272","#fb6a4a","#de2d26","#a50f15","#bf8548","#7c3e04","#58300d","#41280f","#a1d99b","#74c476","#31a354","#006d2c","#bcbddc","#9e9ac8","#756bb1","#54278f")
names(sample_colors) <- c("ConKO F1", "ConKO F2", "ConKO F3", "ConKO M1", "VKO F1", "VKO F2", "VKO M1", "VKO M2","VHKO F1", "VHKO F2", "VHKO M1", "VHKO M2","VEKO F1", "VEKO F2", "VEKO M1", "VEKO M2", "VHEKO F1", "VHEKO M1", "VHEKO M2", "VHEKO M3")
umap$Sample2 <- factor(umap$Sample2,levels = c("ConKO Female 1", "ConKO Female 2", "ConKO Female 3", "ConKO Male 1", "VKO Female 1", "VKO Female 2", "VKO Male 1", "VKO Male 2","VHKO Female 1", "VHKO Female 2", "VHKO Male 1", "VHKO Male 2","VEKO Female 1", "VEKO Female 2", "VEKO Male 1", "VEKO Male 2", "VHEKO Female 1", "VHEKO Male 1", "VHEKO Male 2", "VHEKO Male 3"))
levels(umap$Sample2) <- c("ConKO F1", "ConKO F2", "ConKO F3", "ConKO M1", "VKO F1", "VKO F2", "VKO M1", "VKO M2","VHKO F1", "VHKO F2", "VHKO M1", "VHKO M2","VEKO F1", "VEKO F2", "VEKO M1", "VEKO M2", "VHEKO F1", "VHEKO M1", "VHEKO M2", "VHEKO M3")


p1 <- ggplot(arrange(umap,umap_2),aes(x = umap_1,y = umap_2))+
  geom_point(size = 0.1, shape = 21, color = "black")+
  geom_point(size = 0.01,shape = 16, aes(color = Sample2))+
  scale_color_manual(values = sample_colors)+
  theme_classic()+
  xlab(label = "UMAP 1")+
  ylab(label = "UMAP 2")+
  theme(legend.title = element_blank(),axis.title = element_text(size = 6,hjust = 0.5,margin = margin(0,0,0,0)), axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.size = unit(0.1,"cm"),legend.direction = "horizontal",legend.position = "right",legend.text = element_text(size = 7),legend.margin = margin(0, 0, 0, 0,"cm"),legend.key.spacing = unit(0.1, "cm"))+
  guides(color = guide_legend(override.aes = list(size =1.5),ncol = 4,byrow = TRUE))
p2 <- get_legend(p1)
p3 <- plot_grid(p1+NoLegend(),p2,nrow = 2,rel_heights = c(2,0.5))
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_4/Supplementary_Fig_4a/Supplementary_Fig_4a.pdf",width = 3.5,height = 2.5)
p3
dev.off()


# Supplementary Fig. 4b ---------------------------------------------------

umap <- FetchData(HE_late_pos_equal,vars = c("umap_1","umap_2","Cell_types_MetaUCell","PT_Class_UCell"))

umap[["PT_type"]] <- apply(X = umap,MARGIN = 1,FUN = function(x){
  if(x[["Cell_types_MetaUCell"]] %in% c("PT S1","PT S2","PT S3")){x[["Cell_types_MetaUCell"]]}
  else{"Non PT"}
})
umap$PT_type <- factor(umap$PT_type,levels = c("PT S1","PT S2","PT S3", "Non PT"))

umap$Cell_types_MetaUCell <- factor(umap$Cell_types_MetaUCell,levels = c("PT S1", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte", "PEC", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "Neutrophil", "NK Cell", "Monocyte", "B Cell", "T Cell"))
umap$PT_Class_UCell <- factor(umap$PT_Class_UCell,levels = c("PT Class A","PT Class B", "Non PT"))

PT_type_colors <- c("PT S1" = "#bf8548","PT S2" = "#74c476", "PT S3" = "#9e9ac8","Non PT" = "darkgrey" )


p0 <- ggplot(arrange(umap,umap_1),aes(x = umap_1,y = umap_2))+
  geom_point(mapping = aes(x = umap_1,y = umap_2, color = PT_type), size = 0.01,shape = 16)+
  scale_color_manual(values = PT_type_colors)+
  theme_classic()+
  xlab(label = "UMAP 1")+
  ylab(label = "UMAP 2")+
  theme(legend.title = element_blank(),axis.title = element_text(size = 6,hjust = 0.5,margin = margin(0,0,0,0)), axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.size = unit(0.1,"cm"),legend.direction = "horizontal",legend.position = "right",legend.text = element_text(size = 7),legend.margin = margin(0, 0, 0, 0,"cm"),legend.key.spacing = unit(0.1, "cm"))+
  guides(color = guide_legend(override.aes = list(size =1.5),ncol = 6,byrow = TRUE))
p1 <- ggplot(arrange(umap,umap_1),aes(x = umap_1,y = umap_2))+
  geom_point(size = 0.1, shape = 21, color = "black")+
  geom_point(data = filter(umap, PT_type == c("Non PT")),mapping = aes(x = umap_1,y = umap_2, color = PT_type), size = 0.01,shape = 16)+
  geom_point(data = filter(umap, PT_type != c("Non PT")),mapping = aes(x = umap_1,y = umap_2, color = PT_type), size = 0.01,shape = 16)+
  #geom_point(data = filter(umap, PT_type == "PT S3"),mapping = aes(x = umap_1,y = umap_2, color = PT_type), size = 0.01,shape = 16)+
  scale_color_manual(values = PT_type_colors)+
  theme_classic()+
  xlab(label = "UMAP 1")+
  ylab(label = "UMAP 2")+
  theme(legend.title = element_blank(),axis.title = element_text(size = 6,hjust = 0.5,margin = margin(0,0,0,0)), axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.size = unit(0.1,"cm"),legend.direction = "horizontal",legend.position = "right",legend.text = element_text(size = 6),legend.margin = margin(0, 0, 0, 0,"cm"),legend.key.spacing = unit(0.1, "cm"))+
  guides(color = guide_legend(override.aes = list(size =1),ncol = 6,byrow = TRUE))
p2 <- get_legend(p0)
p3 <- plot_grid(p1+NoLegend(),p2,nrow = 2,rel_heights = c(2,0.5))
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_4/Supplementary_Fig_4b/Supplementary_Fig_4b.pdf",width = 3.5,height = 2.5)
p3
dev.off()


# Supplementary Fig. 4c ---------------------------------------------------

umap <- FetchData(HE_late_pos_equal,vars = c("umap_1","umap_2","PT_Module_A"))

p1 <- ggplot(arrange(umap,umap_1),aes(x = umap_1,y = umap_2))+
  geom_point(size = 0.1, shape = 21,color = "black")+
  geom_point(size = 0.1,aes(colour = PT_Module_A),shape = 16)+
  scale_color_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)))+
  theme_classic()+
  labs(color = "PT Module A expression score")+
  xlab(label = "UMAP 1")+
  ylab(label = "UMAP 2")+
  theme(legend.title = element_text(size = 7,hjust = 1,vjust = 1),legend.text = element_text(size = 6), legend.position = "right",legend.key.height = unit(0.2, "cm"),legend.title.position = "left",legend.direction = "horizontal", axis.title = element_text(size = 6,hjust = 0.5), axis.text = element_blank(), strip.background = element_blank(),axis.ticks = element_blank())
p2 <- get_legend(p1)
p3 <- plot_grid(p1+NoLegend(),p2,nrow = 2,rel_heights = c(2,0.5))
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_4/Supplementary_Fig_4c/Supplementary_Fig_4c.pdf",width = 3.5,height = 2.5)
p3
dev.off()
# Supplementary Fig. 4d ---------------------------------------------------
umap <- FetchData(HE_late_pos_equal,vars = c("umap_1","umap_2","Cell_types_MetaUCell","PT_Class_UCell"))
p1 <- ggplot(arrange(umap,umap_1),aes(x = umap_1,y = umap_2))+
  geom_point(size = 0.1, shape = 21, color = "black")+
  geom_point(size = 0.01,shape = 16, aes(color = PT_Class_UCell))+
  scale_color_manual(values = cellidentities_colour_3)+
  theme_classic()+
  xlab(label = "UMAP 1")+
  ylab(label = "UMAP 2")+
  theme(legend.title = element_blank(),axis.title = element_text(size = 6,hjust = 0.5,margin = margin(0,0,0,0)), axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.size = unit(0.1,"cm"),legend.direction = "horizontal",legend.position = "right",legend.text = element_text(size = 7),legend.margin = margin(0, 0, 0, 0,"cm"),legend.key.spacing = unit(0.1, "cm"))+
  guides(color = guide_legend(override.aes = list(size =1.5),ncol = 3,byrow = TRUE))
p2 <- get_legend(p1)
p3 <- plot_grid(p1+NoLegend(),p2,nrow = 2,rel_heights = c(2,0.5))
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_4/Supplementary_Fig_4d/Supplementary_Fig_4d.pdf",width = 3.5,height = 2.5)
p3
dev.off()

# Supplementary Fig. 4f ---------------------------------------------------
celltypes <- FetchData(VHE_late_pos,vars = c("Sample","Condition","PT_identity"))

celltypes <- celltypes%>%
  group_by(Sample,Condition,PT_identity)%>%
  summarise("Count" = n())%>%
  complete(fill = list(Count = 0),PT_identity)%>%
  mutate("Total" = sum(Count))%>%
  mutate("Percent" = 100*Count/Total)

celltypes[["Dataset"]] <- ifelse(test = celltypes$Condition %in% c("ConKO","VKO"),yes = "ConKO + VKO",no = "VHKO + VEKO + VHEKO")

celltypes2 <- celltypes %>%
  group_by(Dataset,PT_identity)%>%
  summarise("Median" = quantile(Percent,0.5),"Q1" = quantile(Percent, 0.25),"Q3" = quantile(Percent,0.75))

p1 <- ggplot(celltypes2)+
  facet_wrap(facets = vars(Dataset),nrow = 1)+
  geom_col(aes(y = Median,x = PT_identity,fill = PT_identity),show.legend = FALSE)+
  geom_errorbar(aes(y = Median, x = PT_identity,ymin = Q1, ymax = Q3),color = "black",width = 0.3,linewidth = 0.2)+
  scale_fill_manual(values = ptdentities_colour_4)+
  scale_y_continuous(name = "Cells (%)",breaks = c(0,20,40,60))+
  theme_classic()+
  theme(axis.title = element_text(size = 9,colour = "black"),strip.text = element_text(size = 9, color = "black"),axis.text = element_text(size = 8, colour = "black"),legend.text = element_text(size = 6, colour = "black"),legend.title = element_blank(),legend.position = "bottom",legend.direction = "horizontal",axis.title.x = element_blank(),strip.background = element_blank())
pdf(file ="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_4/Supplementary_Fig_4f/Supplementary_Fig_4f.pdf",width = 7,height = 1)
p1
dev.off()





# PT Identities - Dimension Reduction and Clustering ----------------------

PT_ident_seurat <- SplitObject(object = VHE_late_pos,split.by = "PT_identity")
PT_ident_seurat <- PT_ident_seurat[-6]
integrateseurats <- function(a){  
  batch_list=SplitObject(a, split.by = "Sex")
  batch_list <- lapply(X = batch_list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  gc()
  features <- SelectIntegrationFeatures(object.list = batch_list)
  batch_list <- lapply(X = batch_list, FUN = function(y) {
    y <- ScaleData(y, features = features, verbose = FALSE,vars.to.regress=c("nFeature_RNA", "nCount_RNA", "percent.mito"))
    y <- RunPCA(y, features = features, verbose = FALSE)
  })
  gc()
  anchors <- FindIntegrationAnchors(object.list = batch_list, reduction = "rpca",dims = 1:30)
  rm(batch_list)
  gc()
  integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  rm(anchors)
  integrated <- ScaleData(integrated, vars.to.regress=c("nFeature_RNA", "nCount_RNA", "percent.mito"),verbose = FALSE)
  integrated <- RunPCA(integrated)
  integrated <- RunUMAP(integrated, dims = 1:30)
  gc()
  num.pcs=30
  resolution <- c(0.1,0.2,0.4,0.6,0.8,1)
  integrated <- FindNeighbors(object=integrated, dims=1:num.pcs)
  integrated <- FindClusters(object=integrated, resolution=resolution)
  return(integrated)}
PT_ident_seurat <- lapply(PT_ident_seurat, integrateseurats)

"`````````````````````````````````````````````````````````````````````````````"

saveRDS(PT_ident_seurat,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/PT_identities_separated_integrated.rds")

# Fig. 3a -----------------------------------------------------------------
PT_ident_seurat <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/PT_identities_separated_integrated.rds")

plotconditionumap <- function(x){
  umap <- FetchData(x,vars = c("umap_1","umap_2","Condition","Sex","Sample","integrated_snn_res.0.2","PT_identity"))
  umap$integrated_snn_res.0.2 <- factor(as.character(umap$integrated_snn_res.0.2), levels = as.character(c(0:length(unique(umap$integrated_snn_res.0.2))-1)))
  umap$Condition <- factor(umap$Condition,levels = c("ConKO","VKO","VHKO","VEKO","VHEKO"))
  sizes <- ifelse(test = unique(umap$PT_identity) %in% c("S3 A","S3 B"),yes = 0.5,no = 0.1)
  p1 <- ggplot(arrange(umap,umap_1),aes(x = umap_1,y = umap_2))+
    geom_point(size = sizes, shape = 21, color = "black")+
    geom_point(size = sizes,shape = 16, aes(color = Condition))+
    scale_color_manual(values = condition_colour_2)+
    theme_classic()+
    xlab(label = "UMAP 1")+
    ylab(label = "UMAP 2")+
    labs(title = unique(umap$PT_identity))+ 
    theme(legend.title = element_blank(),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0)), axis.text = element_blank(),axis.ticks = element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),legend.key.size = unit(0.1,"cm"),legend.direction = "horizontal",legend.position = "bottom",legend.text = element_text(size = 6),legend.margin = margin(0, 0, 0, 0),legend.spacing = unit(0, "cm"),plot.title = element_text(size = 7,hjust = 0.5,vjust = 0.5),legend.key.spacing = unit(0.1,"cm"))+
    guides(color = guide_legend(override.aes = list(size =1.5),nrow = 1,byrow= TRUE))
  p2 <- p1+NoLegend()
  return(p2)}
condition_plots <- lapply(PT_ident_seurat,plotconditionumap)
names(condition_plots) <- paste("umap",gsub(pattern = " ",replacement = "_",x = names(condition_plots)),sep = "_")
list2env(condition_plots,envir = globalenv())

cluster_condition_color_composition <- function(x){
  clustercomp <- FetchData(x, vars = c("Sample","Condition","PT_identity", "integrated_snn_res.0.2","integrated_snn_res.0.4") )
  nos2 <- length(unique(clustercomp$integrated_snn_res.0.2))-1
  nos4 <- length(unique(clustercomp$integrated_snn_res.0.4))-1
  clustercomp$integrated_snn_res.0.2 <- factor(as.character(clustercomp$integrated_snn_res.0.2),levels = c(0:nos2))
  clustercomp$integrated_snn_res.0.4 <- factor(as.character(clustercomp$integrated_snn_res.0.4),levels = c(0:nos4))
  clusters <- ifelse(test = unique(clustercomp$PT_identity) %in% c("S3 A","S3 B"),yes = "integrated_snn_res.0.4",no = "integrated_snn_res.0.2")
  clustercomp <- clustercomp[,c("Sample","Condition","PT_identity",clusters)]
  colnames(clustercomp) <- c("Sample","Condition","PT_identity", "Clusters")
  clustercomp2 <- clustercomp%>%
    group_by(Sample,Condition, Clusters)%>%
    summarise(Count = n())%>%
    complete(fill = list(Count = 0),Clusters)
  clustercomp2 <- group_by(clustercomp2, Sample,Condition)%>%
    mutate(Total = sum(Count))
  clustercomp2[["Percent"]] <- 100*clustercomp2$Count/clustercomp2$Total
  clustercomp2 <- group_by(clustercomp2,Condition,Clusters)%>%
    summarise("Median" = median(Percent),"Q25" = quantile(Percent,0.25),"Q75" = quantile(Percent,0.75))
  cluster_colour <- c(col_32, brewer.pal(n = 12,name = "Paired"))
  names(cluster_colour) <- c(as.character(c(0:43)))
  clustercomp2$Condition <- factor(clustercomp2$Condition,levels = c("ConKO","VKO","VHKO","VEKO","VHEKO"))
  p1 <- ggplot(clustercomp2)+
    geom_col(aes(x = Clusters,y = Median,fill =Condition),show.legend = FALSE)+
    geom_errorbar(aes(x = Clusters, y = Median,ymin = Q25, ymax = Q75),color = "black",width = 0.3,linewidth = 0.2)+
    scale_fill_manual(values = condition_colour_2)+
    facet_wrap(facets = vars(Condition),nrow = 5,scales = "fixed",strip.position = "right",axis.labels = "margins",axes = "margins",dir = "v")+
    theme_classic()+
    theme(strip.text.y.right = element_text(size = 7,margin = margin(0.1,0.1,0.1,0.1,unit = "cm"),angle = 0),strip.background = element_blank(),legend.title = element_blank(),legend.text = element_blank(),axis.title = element_text(size = 7,hjust = 0.5), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 6,vjust = 0.5,angle = 0,hjust = 1),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(size = 6,hjust = 0.5,vjust = 0.5),axis.line = element_line(linewidth = 0.3),axis.ticks.length = unit(0.01,"cm"))+
    xlab(label = "Clusters")+
    ylab(label = "Proportion in cluster (%)")+
    labs(title = unique(clustercomp$PT_identity))+
    scale_y_continuous(n.breaks = 2)
  p2 <- p1+NoLegend()
  return(p2)}
cluster_composition_condition_plots <- lapply(PT_ident_seurat, cluster_condition_color_composition)
names(cluster_composition_condition_plots) <- paste("cluster",gsub(pattern = " ",replacement = "_",x = names(cluster_composition_condition_plots)),sep = "_")
list2env(cluster_composition_condition_plots,envir = globalenv())

p1 <- plot_grid(umap_S1_A,cluster_S1_A,umap_S1_B,cluster_S1_B,umap_S2_A,cluster_S2_A,umap_S2_B,cluster_S2_B,umap_S3_A,cluster_S3_A,umap_S3_B,cluster_S3_B,nrow = 3)
pdf(file ="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_3/Fig_3a/Fig_3a.pdf",width = 18/2.54,height = 13.5/2.54)
p1
dev.off()


# HIF-dependent Gene Expression -------------------------------------------

sexcondition <- FetchData(VHE_late_pos,c("Sample","Sex","Condition"))%>%group_by(Sample,Sex,Condition)%>%summarise("Count" = n())
rownames(sexcondition) <- gsub(pattern = "_",replacement = "-",x = sexcondition$Sample)

"```````````````````````````````````````````````"
"VKO vs ConKO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(VHE_late_pos$PT_identity)){
  test <- subset(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i & VHE_late_pos$Condition %in% c("VKO", "ConKO")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_V <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_V <- sexcondition[colnames(cts_V$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_V,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","VKO","ConKO"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

V_DE_list <- lapply(cts_V, deseqlist)
saveRDS(V_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VKOvConKO_DE_pseudobulk.rds")

"````````````````````````````````````````````````````````````````````"
"VHEKO vs VKO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(VHE_late_pos$PT_identity)){
  test <- subset(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i & VHE_late_pos$Condition %in% c("VKO", "VHEKO")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_HE <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_HE <- sexcondition[colnames(cts_HE$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_HE,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","VHEKO","VKO"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

HE_DE_list <- lapply(cts_HE, deseqlist)
saveRDS(HE_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHEKOvVKO_DE_pseudobulk.rds")
"````````````````````````````````````````````````````````````````````"
"VHKO vs VKO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(VHE_late_pos$PT_identity)){
  test <- subset(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i & VHE_late_pos$Condition %in% c("VKO", "VHKO")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_H <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_H <- sexcondition[colnames(cts_H$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_H,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","VHKO","VKO"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

H_DE_list <- lapply(cts_H, deseqlist)
saveRDS(H_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHKOvVKO_DE_pseudobulk.rds")
"``````````````````````````````````````````````````````````````````````````````"
"VEKO vs VKO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(VHE_late_pos$PT_identity)){
  test <- subset(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i & VHE_late_pos$Condition %in% c("VKO", "VEKO")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_E <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_E <- sexcondition[colnames(cts_E$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_E,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","VEKO","VKO"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

E_DE_list <- lapply(cts_E, deseqlist)
saveRDS(E_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VEKOvVKO_DE_pseudobulk.rds")

"```````````````````````````````````````````````"
"VHEKO vs ConKO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(VHE_late_pos$PT_identity)){
  test <- subset(VHE_late_pos, cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i & VHE_late_pos$Condition %in% c("VHEKO", "ConKO")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_VHE <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_VHE <- sexcondition[colnames(cts_VHE$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_VHE,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","VHEKO","ConKO"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

VHE_DE_list <- lapply(cts_VHE, deseqlist)
saveRDS(VHE_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHEKOvConKO_DE_pseudobulk.rds")


# HIF-dependence of Vhl-dependent genes -----------------------------------

"Getiing list of differentially regulated genes"

V_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VKOvConKO_DE_pseudobulk.rds")
H_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHKOvVKO_DE_pseudobulk.rds")
E_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VEKOvVKO_DE_pseudobulk.rds")
HE_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHEKOvVKO_DE_pseudobulk.rds")

"Defining HIF dependent genes - in each PT identity, genes that are up/down in VKO vs ConKO, down/up after ANY HIF KO vs VKO, and whose extent of regulation following any HIF KO is atleast half of the regulation following VKO"

"```````````````````````````````````````````````````````````````````````````````"
"Applying significance thresholds"

V_DE_up <- lapply(V_DE_list, function(x){filter(x,log2FoldChange > 1,padj < 0.05)$gene})
V_DE_down <- lapply(V_DE_list, function(x){filter(x,log2FoldChange < -1,padj < 0.05)$gene})

HE_DE_up <- lapply(HE_DE_list, function(x){filter(x,log2FoldChange > 1,padj < 0.05)$gene})
HE_DE_down <- lapply(HE_DE_list, function(x){filter(x,log2FoldChange < -1,padj < 0.05)$gene})
H_DE_up <- lapply(H_DE_list, function(x){filter(x,log2FoldChange > 1,padj < 0.05)$gene})
H_DE_down <- lapply(H_DE_list, function(x){filter(x,log2FoldChange < -1,padj < 0.05)$gene})
E_DE_up <- lapply(E_DE_list, function(x){filter(x,log2FoldChange > 1,padj < 0.05)$gene})
E_DE_down <- lapply(E_DE_list, function(x){filter(x,log2FoldChange < -1,padj < 0.05)$gene})

"`````````````````````````````````````````````````````````````````````````````````"
"HIF dependent genes by significance thresholds. Sig in VKO vs ConKO + significant in reverse direction in ANY HIF KO vs VKO"

V_down_Hif_dep_sig <- intersectlist(V_DE_down,appendlist(HE_DE_up,H_DE_up,E_DE_up))
V_up_Hif_dep_sig <- intersectlist(V_DE_up,appendlist(HE_DE_down,H_DE_down,E_DE_down))
"`````````````````````````````````````````````````````````````````````````````````"
"HIF dependent genes by L2FC. If L2FC VKO vs ConKO is x, then magnitude of L2FC in any HIF KO vs VKO has to be more than x/2 in opposite direction. e.g. if gene A has L2FC 4.2 in VKO vs ConKO, it has to have atleast L2FC -2.1 in VHKO vs VKO, VEKO vs VKO, or VHEKO vs VKO"

V_HE_H_E_FC <- list()
length(V_HE_H_E_FC) <- 7
names(V_HE_H_E_FC) <- names(V_DE_list)
for(i in names(V_HE_H_E_FC)){
  test1 <- V_DE_list[[i]][,c(1,3,7)]
  names(test1) <- c("Gene","VKO vs ConKO - L2FC","VKO vs ConKO - p_adj")
  test2 <- HE_DE_list[[i]][,c(1,3,7)]
  names(test2) <- c("Gene","VHEKO vs VKO - L2FC","VHEKO vs VKO - p_adj")
  test3 <- H_DE_list[[i]][,c(1,3,7)]
  names(test3) <- c("Gene","VHKO vs VKO - L2FC","VHKO vs VKO - p_adj")
  test4 <- E_DE_list[[i]][,c(1,3,7)]
  names(test4) <- c("Gene","VEKO vs VKO - L2FC","VEKO vs VKO - p_adj")
  test5 <- VHE_DE_list[[i]][,c(1,3,7)]
  names(test5) <- c("Gene","VHEKO vs ConKO - L2FC","VHEKO vs ConKO - p_adj")
  test <- full_join(test1,test2,by = "Gene")
  test <- full_join(test,test3,by = "Gene")
  test <- full_join(test,test4,by = "Gene")
  test <- full_join(test,test5,by = "Gene")
  test[is.na(test)] <- 0
  test[["PT_identity"]] <- i
  test$PT_identity <- factor(test$PT_identity,levels = c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B","Non PT"))
  test <- test[,c(1,12,2:11)]
  V_HE_H_E_FC[[i]] <- test
}

L2FC_Hif_dep <- function(x){
  result <- list()
  length(result) <- 7
  for (i in 1:7){
    genes <- appendlist(V_DE_up,V_DE_down,NULL)[[i]]
    genes_1 <- c()
    for(j in genes){
      test11 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VHEKO vs VKO - L2FC`
      test12 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VHKO vs VKO - L2FC`
      test13 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VEKO vs VKO - L2FC`
      test2 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VKO vs ConKO - L2FC`
      diff <- c()
      if(test2 > 0){
        test1 <- min(test11,test12,test13)
      }
      else{test1 <- max(test11,test12,test13)}
      if(test1*test2 >= 0){diff <- NA}
      else if(test2 > 0){diff <- ifelse(test = test1 < -test2/2,yes = j,no = NA)}
      else{diff <- ifelse(test = test1 > -test2/2,yes = j,no = NA)}
      genes_1 <-  c(genes_1,diff)
    }
    if(length(genes_1) == 0){genes_1 <- NA}
    result[[i]] <- genes_1[!is.na(genes_1)]
  }
  names(result) <- names(V_DE_list)
  return(result)
}
V_Hif_dep_l2fc <- L2FC_Hif_dep(HE_DE_list)

"```````````````````````````````````````````````````````````````````````````````"
"Combining significance and L2FC thresholding"

V_down_Hif_dep <- intersectlist(V_down_Hif_dep_sig,V_Hif_dep_l2fc)
V_up_Hif_dep <- intersectlist(V_up_Hif_dep_sig,V_Hif_dep_l2fc)

"``````````````````````````````````````````````````````````````````````````````"
"Identifying genes that are unaltered with single or combined HIF deletion. These are significantly up/down in VKO vs ConKO AND in VHEKO vs ConKO. Additionally, they are NOT significantly down/up in any HIF KO vs VKO. Furthermore, the magnitude of L2FC in VHEKO vs ConKO should be more than half the L2FC VKO vs ConKO in the same direction. Finally, the magnitude of L2FC in ANY HIF KO vs VKO should be less than half the L2FC VKO vs ConKO in the opposite direction."

"Significance thresholds"
V_down_Hif_ind_sig <- setdifflist(intersectlist(V_DE_down, VHE_DE_down),appendlist(HE_DE_up,H_DE_up,E_DE_up))
V_up_Hif_ind_sig <- setdifflist(intersectlist(V_DE_up, VHE_DE_up),appendlist(HE_DE_down,H_DE_down,E_DE_down))

"L2FC thresholds"

L2FC_Hif_ind <- function(x){
  result <- list()
  length(result) <- 7
  for (i in 1:7){
    genes <- appendlist(V_DE_up,V_DE_down,NULL)[[i]]
    genes_1 <- c()
    for(j in genes){
      test11 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VHEKO vs VKO - L2FC`
      test12 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VHKO vs VKO - L2FC`
      test13 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VEKO vs VKO - L2FC`
      test14 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VHEKO vs ConKO - L2FC`
      test2 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VKO vs ConKO - L2FC`
      diff <- c()
      if(test2 >= 0){
        if(test14 < test2/2){diff <- NA}
        else{
          test1 <- min(test11,test12,test13)
          if(test1*test2 >= 0){diff <- j}
          else{diff <- ifelse(test = test1 > -test2/2,yes = j,no = NA)}
        }
      }
      else{
        if(test14 > test2/2){diff <- NA}
        else{
          test1 <- max(test11,test12,test13)
          if(test1*test2 >= 0){diff <- j}
          else{diff <- ifelse(test = test1 < -test2/2,yes = j,no = NA)}
        }
      }
      genes_1 <-  c(genes_1,diff)
      }
    if(length(genes_1) == 0){genes_1 <- NA}
    result[[i]] <- genes_1[!is.na(genes_1)]
    }
  names(result) <- names(V_DE_list)
  return(result)
}
V_Hif_ind_l2fc <- L2FC_Hif_ind(HE_DE_list)

"Applying significance and L2FC thresholds"
V_down_Hif_ind <- intersectlist(V_down_Hif_ind_sig,V_Hif_ind_l2fc)
V_up_Hif_ind <- intersectlist(V_up_Hif_ind_sig,V_Hif_ind_l2fc)

"````````````````````````````````````````````````````````````````````````````````"


# HIF isoform specificity of HIF dependent genes --------------------------
"Need HIF1A Alone = HIF-dependent, HIF1A-dependent, but not HIF2A-dependent"
"Need HIF2A Alone = HIF-dependent, HIF2A-dependent, but not HIF1A-dependent"
"Need both isoforms = HIF-dependent, HIF1A-dependent, and HIF2A-dependent"
"Need either isoforms = HIF-dependent, but not HIF1A-dependent or HIF2A-dependent"

"Significance thresholds"
V_up_H_dep_sig <- intersectlist(V_up_Hif_dep,H_DE_down)
V_down_H_dep_sig <- intersectlist(V_down_Hif_dep,H_DE_up)
V_up_E_dep_sig <- intersectlist(V_up_Hif_dep,E_DE_down)
V_down_E_dep_sig <- intersectlist(V_down_Hif_dep,E_DE_up)

"L2FC thresholds"
L2FC_dep <- function(x){
  result <- x
  result <- list()
  length(result) <- 7
  names(result) <- names(V_DE_list)
  for (i in names(result)){
    genes <- appendlist(V_up_Hif_dep,V_down_Hif_dep,NULL)[[i]]
    genes_1 <- c()
    for(j in genes){
      test1 <- filter(V_HE_H_E_FC[[i]], Gene == j)[[x]]
      test2 <- filter(V_HE_H_E_FC[[i]], Gene == j)$`VKO vs ConKO - L2FC`
      diff <- c()
      if(test1*test2 >= 0){diff <- NA}
      if(test2 > 0){diff <- ifelse(test = test1 < -test2/2,yes = j,no = NA)}
      else{diff <- ifelse(test = test1 > -test2/2,yes = j,no = NA)}
      genes_1 <-  c(genes_1,diff)
    }
    if(length(genes_1) == 0){genes_1 <- NA}
    result[[i]] <- genes_1[!is.na(genes_1)]
  }
  return(result)
}

V_up_H_dep <- intersectlist(L2FC_dep("VHKO vs VKO - L2FC"),V_up_H_dep_sig)
V_down_H_dep <- intersectlist(L2FC_dep("VHKO vs VKO - L2FC"),V_down_H_dep_sig)
V_up_E_dep <- intersectlist(L2FC_dep("VEKO vs VKO - L2FC"),V_up_E_dep_sig)
V_down_E_dep <- intersectlist(L2FC_dep("VEKO vs VKO - L2FC"),V_down_E_dep_sig)

Need_H_up <- setdifflist(V_up_H_dep,V_up_E_dep)
Need_E_up <- setdifflist(V_up_E_dep,V_up_H_dep)
Need_both_up <- intersectlist(V_up_H_dep,V_up_E_dep)
Need_either_up <- setdifflist(V_up_Hif_dep,appendlist(V_up_H_dep,V_up_E_dep,NULL))
Need_H_down <- setdifflist(V_down_H_dep,V_down_E_dep)
Need_E_down <- setdifflist(V_down_E_dep,V_down_H_dep)
Need_both_down <- intersectlist(V_down_H_dep,V_down_E_dep)
Need_either_down <- setdifflist(V_down_Hif_dep,appendlist(V_down_H_dep,V_down_E_dep,NULL))

saveRDS(list("Need_H_up" = Need_H_up,"Need_E_up"= Need_E_up,"Need_both_up"=Need_both_up,"Need_either_up"=Need_either_up,"Need_H_down" = Need_H_down,"Need_E_down"= Need_E_down,"Need_both_down"=Need_both_down,"Need_either_down"=Need_either_down,"Need_neither_up"=V_up_Hif_ind, "Need_neither_down"=V_down_Hif_ind,"All_up" = V_DE_up,"All_down" = V_DE_down,"Hif_dep_up" = V_up_Hif_dep,"Hif_dep_down" = V_down_Hif_dep),file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")

"`````````````````````````````````````````````````````````````````````````````"
"Isoform specificity of HIF dependent genes across PT identities"
addhifspec <- function(x){
  x[["HIF_dependence"]] <- apply(x,1,function(y){
    j <- y[["Gene"]]
    i <- y[["PT_identity"]]
    if(j %in% c(Need_H_up[[i]],Need_H_down[[i]])){"Dependent on HIF1A alone"}
    else if(j %in% c(Need_E_up[[i]],Need_E_down[[i]])){"Dependent on HIF2A alone"}
    else if(j %in% c(Need_both_up[[i]],Need_both_down[[i]])){"Dependent on both isoforms"}
    else if(j %in% c(Need_either_up[[i]],Need_either_down[[i]])){"Dependent on either isoform"}
    #else if(j %in% c(Need_neither_up[[i]])){PT_identity_HIF_dependence_up[j,i] <- "Dependent on neither isoform"}
    else if(j %in% c(V_DE_up[[i]],V_DE_down[[i]])){"Ambiguous HIF dependence"}
    else{NA}
  })
  return(x)}

V_HE_H_E_FC_hif <- lapply(V_HE_H_E_FC, addhifspec)
saveRDS(V_HE_H_E_FC_hif,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/V_H_E_HE_DE.rds")


# Renaming isoform-specificity of genes -----------------------------------
V_HE_H_E_FC_hif <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/V_H_E_HE_DE.rds")

V_HE_H_E_FC_hif <- lapply(V_HE_H_E_FC_hif, function(x){
  test <- x
  test$HIF_dependence <- factor(test$HIF_dependence,levels = c("Ambiguous HIF dependence", "Dependent on either isoform","Dependent on HIF1A alone", "Dependent on HIF2A alone","Dependent on both isoforms"),labels = c("Ambiguous HIF dependence","HIF1A + HIF2A","HIF1A alone","HIF2A alone","HIF1A or HIF2A"))
  return(test)
})

saveRDS(V_HE_H_E_FC_hif,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/V_H_E_HE_DE_rename.rds")
# Supplementary Table 2 ---------------------------------------------------
test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())
V_HE_H_E_FC_hif <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/V_H_E_HE_DE_rename.rds")

V_HE_H_E_FC_hif_df <- do.call(rbind,V_HE_H_E_FC_hif) 
write.table(x =V_HE_H_E_FC_hif_df,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Papers/Vhl_HIF_Paper/Nature_Communications_Revisions/Supp_Table_2.csv",sep = ",",row.names = FALSE,col.names = TRUE)
write.table(x =V_HE_H_E_FC_hif_df,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Table_2/Supplemetary_Table_2.csv",sep = ",",row.names = FALSE,col.names = TRUE)


# Fig. 3b -----------------------------------------------------------------
V_HE_H_E_FC_hif <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/V_H_E_HE_DE_rename.rds")

Hif_dep_up <- Reduce(f = rbind,lapply(V_HE_H_E_FC_hif, function(x){filter(x,!is.na(HIF_dependence),Gene %in% Hif_dep_up[[unique(x[["PT_identity"]])]])%>%
    group_by(HIF_dependence)%>%
    summarise("Count" = n())%>%
    complete(fill = list(Count = 0),HIF_dependence)%>%
    mutate("Total" = sum(Count))%>%
    mutate("Percent" = 100*Count/Total,"PT_identity" = unique(x[["PT_identity"]]))
}))
Hif_dep_down <- Reduce(f = rbind,lapply(V_HE_H_E_FC_hif, function(x){filter(x,!is.na(HIF_dependence),Gene %in% Hif_dep_down[[unique(x[["PT_identity"]])]])%>%
    group_by(HIF_dependence)%>%
    summarise("Count" = n())%>%
    complete(fill = list(Count = 0),HIF_dependence)%>%
    mutate("Total" = sum(Count))%>%
    mutate("Percent" = 100*Count/Total,"PT_identity" = unique(x[["PT_identity"]]))
}))
   
Hif_dep_up$HIF_dependence <- factor(Hif_dep_up$HIF_dependence,levels = c("Dependent on either isoform","Dependent on HIF1A alone", "Dependent on HIF2A alone","Dependent on both isoforms"))
Hif_dep_up$PT_identity <- factor(Hif_dep_up$PT_identity,levels = c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B","Non PT"))
Hif_dep_down$HIF_dependence <- factor(Hif_dep_down$HIF_dependence,levels = c("Dependent on either isoform","Dependent on HIF1A alone", "Dependent on HIF2A alone","Dependent on both isoforms"))
Hif_dep_down$PT_identity <- factor(Hif_dep_down$PT_identity,levels = c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B", "Non PT"))

Hif_dep_up[["Category"]] <- "Upregulated"
Hif_dep_down[["Category"]] <- "Downregulated"
Hif_dep_count <- rbind(Hif_dep_up,Hif_dep_down)

Hif_dep_count$Category <- factor(Hif_dep_count$Category,levels = c("Upregulated","Downregulated"))

p1 <- ggplot(filter(Hif_dep_count, PT_identity %in% c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B"), HIF_dependence != "Ambiguous HIF dependence"),aes(x = HIF_dependence,y = Percent,fill = HIF_dependence))+
  facet_grid(rows = vars(Category),cols = vars(PT_identity),scales = "free_y",space = "fixed",axes = "all",axis.labels = "margins")+
  geom_col()+
  scale_fill_manual(values = Isoform_spec_colors)+
  theme_classic()+
  ylab(label = expression(paste("HIF", alpha, "-dependent genes (%)",sep = "")))+
  xlab(label = NULL)+
  theme_classic()+
  theme(legend.title = element_blank(),legend.text = element_blank(),axis.title = element_text(size = 9,hjust = 0.5,margin = margin(0,0,0,0)),axis.text = element_text(size = 6),axis.text.x = element_text(size = 7,angle = 90,hjust = 1,vjust = 0.5),plot.margin = unit(c(0,0.1,0.1,0.1),"cm"),strip.background = element_blank(),strip.text.y.right = element_text(size = 8,hjust = 0.5,angle = 270),plot.title = element_text(size = 8,hjust = 0.5,vjust = 0.5,margin = margin(0,0,0,0)))+
  NoLegend()
pdf(file = "//Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_3/Fig_3b/Fig_3b.pdf",width = 18/2.54,height = 3.3)
p1
dev.off()


# Supplementary Fig 5a -----------------------------------------------------
V_HE_H_E_FC_hif_df <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Table_2/Supplemetary_Table_2.csv")
V_HE_H_E_FC_hif_df$HIF_dependence <- factor(V_HE_H_E_FC_hif_df$HIF_dependence,levels = c("Dependent on both isoforms","Dependent on HIF1A alone","Dependent on HIF2A alone", "Dependent on either isoform", "Ambiguous HIF dependence"))
test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())

Isoform_spec_colors <- c(`Dependent on HIF1A alone` = "#bc8e54", `Dependent on HIF2A alone` = "#208d4c",`Dependent on both isoforms` = "#a50f15", `Dependent on either isoform` = "#9467cf")
Isoform_spec_colors_2 <- c("Dependent" = "skyblue2", "Not dependent" ="grey")

test <- arrange(V_HE_H_E_FC_hif_df,desc(HIF_dependence))
test[["hif_Dep"]] <- apply(test,1,FUN = function(x){
  if(is.na(x[["HIF_dependence"]])){NA}
  else if(x[["HIF_dependence"]] == "Ambiguous HIF dependence"){"Not dependent"}
  else {"Dependent"}
})

test <- filter(test, PT_identity != "Non PT")
test <- arrange(test,Gene)
p1 <- ggplot()+
  scale_color_manual(values = Isoform_spec_colors_2)+
  scale_x_continuous(name = "L2FC VKO vs ConKO",limits = c(-6.5,11.5),breaks = c(-5,-1,3,7,11),expand = c(0,0))+
  scale_y_continuous(name = "L2FC VHEKO vs VKO",limits = c(-7.6,5.5),breaks = c(-7,-4,-1,2,5),expand = c(0,0))+
  labs(color = expression(paste("HIF", alpha, "-dependence:",sep = "")))+
  geom_hline(yintercept = 0,linetype = "solid")+
  geom_vline(xintercept = 0,linetype = "solid")+
  geom_point(data = filter(test,!is.na(HIF_dependence)), aes(x = `VKO vs ConKO - L2FC`,y = `VHEKO vs VKO - L2FC`,color = hif_Dep),size = 0.5)+
  facet_wrap(facets = vars(PT_identity),nrow = 2,axes = "all")+
  theme_classic()+
  theme(legend.title = element_text(size = 9,hjust = 0,vjust = 0.5),legend.text = element_text(size = 8),plot.title = element_text(size = 9,color = "black",hjust = 0.5), axis.title = element_text(size = 8,hjust = 0.5,margin = margin(0,0,0,0)),axis.text = element_text(size = 6),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),strip.background = element_blank(),strip.text.y.right = element_text(size = 7,hjust = 0.5,angle = 0),legend.margin = margin(0, 0, 0, 0.1),legend.spacing = unit(0, "cm"),legend.key.size = unit(0.2,"cm"),legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size =1),nrow = 1,title.position = "left"))
p2 <- get_legend(p1)
p5 <- plot_grid(p1+NoLegend(),p2,nrow = 2,rel_heights = c(11,0.5))
pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_5/Supplementary_Fig_5a/Supplementary_Fig_5a.pdf",width = 18/2.54,height = 11.5/2.54)
p5
dev.off()

correlations <- matrix(ncol = 2, nrow = 6)
rownames(correlations) <- unique(test$PT_identity)
colnames(correlations) <- c("Not dependent","Dependent")
for(i in unique(test$PT_identity)){
  for(j in c("Not dependent","Dependent")){
    test1 <- filter(test, PT_identity == i,hif_Dep == j)
    correlation <- cor.test(x = test1$`VKO vs ConKO - L2FC`,y = test1$`VHEKO vs VKO - L2FC`,method = "spearman",alternative = "t",exact = FALSE)
    if(correlation$p.value < 0.05){correlations[i,j] <- correlation$estimate}
    else{correlations[i,j] <- 0}
  }}
correlations_df <- correlations%>%
  as.data.frame()%>%
  rownames_to_column(var = "PT_identity")

write.table(correlations_df,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_5/Supplementary_Fig_5a/Spearman_correlations.csv",sep = ",",row.names = FALSE,col.names = TRUE)


# Supplementary Fig 5b ----------------------------------------------------
V_HE_H_E_FC_hif_df <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Table_2/Supplemetary_Table_2.csv")

test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())

test <- arrange(filter(V_HE_H_E_FC_hif_df, PT_identity != "Non PT"),desc(HIF_dependence))
Isoform_spec_colors <- c(`HIF1A alone` = "#bc8e54", `HIF2A alone` = "#208d4c",`HIF1A or HIF2A` = "#a50f15", `HIF1A + HIF2A` = "#9467cf")
test$HIF_dependence <- factor(test$HIF_dependence,levels = names(Isoform_spec_colors))



p3 <- ggplot()+
  scale_color_manual(values = Isoform_spec_colors)+
  scale_y_continuous(name = "L2FC VEKO vs VKO",limits = c(-7.6,5.5),breaks = c(-7,-4,-1,2,5),expand = c(0,0))+
  scale_x_continuous(name = "L2FC VHKO vs VKO",limits = c(-7.6,5.5),breaks = c(-7,-4,-1,2,5),expand = c(0,0))+
  #labs(color = expression(paste("HIF", alpha, "-isoform specificity: Dependence on",sep = "")))+
  labs(color = "Regulation requires deletion of:")+
  geom_hline(yintercept = 0,linetype = "solid")+
  geom_vline(xintercept = 0,linetype = "solid")+
  geom_point(data = filter(test,!is.na(HIF_dependence)), aes(y = `VEKO vs VKO - L2FC`,x = `VHKO vs VKO - L2FC`),color = "grey",size = 0.5)+
  geom_point(data = filter(test,!is.na(HIF_dependence),HIF_dependence %in% c("HIF1A alone","HIF2A alone")), aes(y = `VEKO vs VKO - L2FC`,x = `VHKO vs VKO - L2FC`,color = HIF_dependence),size = 0.7)+
  facet_wrap(facets = vars(PT_identity),nrow = 2,axes = "all")+
  theme_classic()+
  theme(legend.title = element_text(size = 8,hjust = 0.5,vjust = 0.5),legend.text = element_text(size = 8),plot.title = element_text(size = 9,color = "black",hjust = 0.5), axis.title = element_text(size = 8,hjust = 0.5,margin = margin(0,0,0,0)),axis.text = element_text(size = 6),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),strip.background = element_blank(),strip.text.y.right = element_text(size = 7,hjust = 0.5,angle = 0),legend.margin = margin(0, 0, 0, 0),legend.spacing = unit(0, "cm"),legend.key.size = unit(0.2,"cm"),legend.position = "right")+
  guides(color = guide_legend(override.aes = list(size =1),nrow = 1,title.position = "left"))
p4 <- get_legend(p3)

p6 <- plot_grid(p3+NoLegend(),p4,nrow = 2,rel_heights = c(11,0.5))
p6
pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_5/Supplementary_Fig_5b/Supplementary_Fig_5b.pdf",width = 18/2.54,height = 11.5/2.54)
p6
dev.off()

correlations <- matrix(ncol = 3, nrow = 6)
rownames(correlations) <- unique(test$PT_identity)
colnames(correlations) <- c("HIF1A alone","HIF2A alone","All")

for(i in unique(test$PT_identity)){
  for(j in c("HIF1A alone","HIF2A alone")){
    test1 <- filter(test, PT_identity == i, HIF_dependence %in% j)
    if(nrow(test1) > 4){
      correlation <- cor.test(x = test1$`VHKO vs VKO - L2FC`,y = test1$`VEKO vs VKO - L2FC`,method = "spearman",alternative = "t",exact = FALSE)
      #if(correlation$p.value < 0.01){correlations[i,j] <- correlation$estimate}
      #else{correlations[i,j] <- 0}
      correlations[i,j] <- correlation$estimate
    }
    else{correlations[i,j] <- NA}}}
for(i in unique(test$PT_identity)){
    test1 <- filter(test, PT_identity == i, !is.na(HIF_dependence))
    if(nrow(test1) > 4){
      correlation <- cor.test(x = test1$`VHKO vs VKO - L2FC`,y = test1$`VEKO vs VKO - L2FC`,method = "spearman",alternative = "t",exact = FALSE)
      #if(correlation$p.value < 0.01){correlations[i,j] <- correlation$estimate}
      #else{correlations[i,j] <- 0}
      correlations[i,"All"] <- correlation$estimate
    }
    else{correlations[i,"All"] <- NA}}

correlations_df <- correlations%>%
  as.data.frame()%>%
  rownames_to_column(var = "PT_identity")

write.table(correlations_df,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_5/Supplementary_Fig_5b/Spearman_correlations_alltogether.csv",sep = ",",row.names = FALSE,col.names = TRUE)

# Fig 4a ------------------------------------------------------------------
test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())

all_need_hif1_up <- unique(unlist(Need_H_up))
all_need_hif2_up <- unique(unlist(Need_E_up))

all_up <- list("Dependent on HIF1A alone" = all_need_hif1_up,"Dependent on HIF2A alone" = all_need_hif2_up)
names(all_up) <- c("HIF1A","HIF2A")
p1 <- ggvenn::ggvenn(data = all_up,fill_color = c("white","white","white"),show_percentage = FALSE,set_name_color = "black",stroke_color = "black",set_name_size = 2,text_size = 3,stroke_size = 0.4,auto_scale = FALSE)
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4a/Fig_4a_up.pdf",width = 1.75,height = 1.75)
p1
dev.off()

all_need_hif1_down <- unique(unlist(Need_H_down))
all_need_hif2_down <- unique(unlist(Need_E_down))

all_down <- list("Dependent on HIF1A alone" = all_need_hif1_down,"Dependent on HIF2A alone" = all_need_hif2_down)
names(all_down) <- c("HIF1A","HIF2A")
p1 <- ggvenn::ggvenn(data = all_down,fill_color = c("white","white","white"),show_percentage = FALSE,set_name_color = "black",stroke_color = "black",set_name_size = 2,text_size = 3,stroke_size = 0.4,auto_scale = FALSE)
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4a/Fig_4a_down.pdf",width = 1.75,height = 1.75)
p1
dev.off()



# Fig 4c ------------------------------------------------------------------
VHE_late_pos <- readRDS(file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_VHKO_VEKO_VHEKO_Late_pos.rds")

"Want to have < 30,000 cells for ggplot to handle the heatmap accurately. WIll take 5999 cells per condition"

table(VHE_late_pos$Sex,VHE_late_pos$Condition)
"Downsample to 1500 cells per sample"
set.seed(1234)
celllist_VHE <- list()
for(i in levels(VHE_late_pos$Condition)){
  cells <- colnames(VHE_late_pos)[VHE_late_pos$Condition == i & VHE_late_pos$Sex == "Female"]
  test <- sample(cells,5999,replace = FALSE)
  celllist_VHE[[i]] <- test
}

cells2 <- unlist(celllist_VHE)

test <- subset(VHE_late_pos,cells = cells2,features = unique(c(unlist(Need_H_up),unlist(Need_E_up),unlist(Need_E_down))))

Hif1_glyco_genes <- c("Pfkl",'Gpi1',"Pgam1","Eno1","Pgk1","Pdk1","Car9","Bnip3")

quantile(test@assays[["RNA"]]@data[c(Hif1_glyco_genes),])

scale_fill_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)))

p1 <- DoHeatmap(object = test,features = c(Hif1_glyco_genes),group.by = c("Condition"),slot = "data",assay = "RNA",draw.lines = TRUE,group.colors = condition_colour_2,size = 4,hjust = 0.5,vjust = 0.5,angle = 0,lines.width = 500,raster = FALSE,disp.min = 0,disp.max = 2)+
  scale_fill_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)),na.value = "white")
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4c/Fig_4c.pdf",width = 5,height = 3.5)
p1
dev.off()


# Fig. 4d -----------------------------------------------------------------

VHE_late_pos <- readRDS(file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_VHKO_VEKO_VHEKO_Late_pos.rds")

"Want to have < 30,000 cells for ggplot to handle the heatmap accurately. WIll take 5999 cells per condition"

table(VHE_late_pos$Sex,VHE_late_pos$Condition)
"Downsample to 1500 cells per sample"
set.seed(1234)
celllist_VHE <- list()
for(i in levels(VHE_late_pos$Condition)){
  cells <- colnames(VHE_late_pos)[VHE_late_pos$Condition == i & VHE_late_pos$Sex == "Female"]
  test <- sample(cells,5999,replace = FALSE)
  celllist_VHE[[i]] <- test
}

cells2 <- unlist(celllist_VHE)

test <- subset(VHE_late_pos,cells = cells2,features = unique(c(unlist(Need_H_up),unlist(Need_E_up),unlist(Need_E_down))))

Hif2_up_genes <- c("Angptl3","Igfbp5", "Apela","Col4a1","P4ha1","Npnt","Fabp5","Slc3a1")

quantiles <- quantile(test@assays[["RNA"]]@data[c(Hif2_up_genes),])

p1 <- DoHeatmap(object = test,features = c(Hif2_up_genes),group.by = c("Condition"),slot = "data",assay = "RNA",draw.lines = TRUE,group.colors = condition_colour_2,size = 4,hjust = 0.5,vjust = 0.5,angle = 0,lines.width = 500,raster = FALSE,disp.min = 0,disp.max = 2)+
  scale_fill_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)),na.value = "white")
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4d/Fig_4d.pdf",width = 5,height = 3.5)
p1
dev.off()


# Fig. 4f -----------------------------------------------------------------

VHE_late_pos <- readRDS(file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_VHKO_VEKO_VHEKO_Late_pos.rds")

"Want to have < 30,000 cells for ggplot to handle the heatmap accurately. WIll take 5999 cells per condition"

table(VHE_late_pos$Sex,VHE_late_pos$Condition)
"Downsample to 1500 cells per sample"
set.seed(1234)
celllist_VHE <- list()
for(i in levels(VHE_late_pos$Condition)){
  cells <- colnames(VHE_late_pos)[VHE_late_pos$Condition == i & VHE_late_pos$Sex == "Female"]
  test <- sample(cells,5999,replace = FALSE)
  celllist_VHE[[i]] <- test
}

cells2 <- unlist(celllist_VHE)

test <- subset(VHE_late_pos,cells = cells2,features = unique(c(unlist(Need_H_up),unlist(Need_E_up),unlist(Need_E_down))))

Hif2_down_genes <- down_genes <- c("Slc5a12","Slc6a19","Inmt","Cyp4a14","Cyp4a10","Slc22a8","Gsta3")

quantile(test@assays[["RNA"]]@data[c(Hif2_down_genes),])

p1 <- DoHeatmap(object = test,features = c(Hif2_down_genes),group.by = c("Condition"),slot = "data",assay = "RNA",draw.lines = TRUE,group.colors = condition_colour_2,size = 4,hjust = 0.5,vjust = 0.5,angle = 0,lines.width = 500,raster = FALSE,disp.min = 0,disp.max = 2)+
  scale_fill_gradientn(colours = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),values = scales::rescale(c(0,50, 100)),na.value = "white")
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4f/Fig_4f.pdf",width = 5,height = 3.5)
p1
dev.off()

# GO Term Analysis - Isoform SpecificGenes --------------------------
V_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VKOvConKO_DE_pseudobulk.rds")
H_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHKOvVKO_DE_pseudobulk.rds")
E_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VEKOvVKO_DE_pseudobulk.rds")
HE_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHEKOvVKO_DE_pseudobulk.rds")
test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())
V_HE_H_E_FC_hif_df <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Table_2/Supplemetary_Table_2.csv")

Reference_genes_all <- Reduce(intersect,list(unlist(lapply(V_DE_list,function(x){x[["gene"]]})),unlist(lapply(HE_DE_list,function(x){x[["gene"]]})),unlist(lapply(H_DE_list,function(x){x[["gene"]]})),unlist(lapply(E_DE_list,function(x){x[["gene"]]}))))

"Some GO terms have comments advising against use:"
avoidGO <- c("GO:0009628","GO:0001101","GO:0010243","GO:1901698","GO:0071417","GO:1901699")

"Get enriched GO terms in list of genes that H_only, E_only. Filter terms such that adjusted p value is < 0.01 and there are atleast 10 target genes in each term"

Need_H_up_enrichGO <- enrichGO(gene = unique(unlist(Need_H_up)),OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "fdr",universe = Reference_genes_all,qvalueCutoff = 0.05,minGSSize = 50,maxGSSize = 1000,readable = TRUE)
orsize <- ifelse(as.numeric(unlist(lapply(str_split(Need_H_up_enrichGO@result$GeneRatio,pattern = "/",n = 2),function(x)x[[1]]))) > 7,yes = TRUE,no = FALSE)
test2 <- Need_H_up_enrichGO@result[orsize,c("ID","Description","p.adjust","geneID")]
colnames(test2) <- c("ID","Description", "Need_H","Need_H_genes")
terms2 <- filter(test2, Need_H < 0.01)$ID

Need_E_up_enrichGO <- enrichGO(gene = unique(unlist(Need_E_up)),OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "fdr",universe = Reference_genes_all,qvalueCutoff = 0.05,minGSSize = 50,maxGSSize = 1000,readable = TRUE)
orsize <- ifelse(as.numeric(unlist(lapply(str_split(Need_E_up_enrichGO@result$GeneRatio,pattern = "/",n = 2),function(x)x[[1]]))) > 7,yes = TRUE,no = FALSE)
test3 <- Need_E_up_enrichGO@result[orsize,c("ID","Description","p.adjust","geneID")]
colnames(test3) <- c("ID","Description", "Need_E","Need_E_genes")
terms3 <- filter(test3, Need_E < 0.01)$ID


test <- full_join(test2,test3,by = c("ID","Description"))
isoform_up_enrichgo <- filter(test, ID %in% unique(c(terms2,terms3)))
isoform_up_enrichgo <- filter(isoform_up_enrichgo, ! ID %in% avoidGO)
rownames(isoform_up_enrichgo) <- isoform_up_enrichgo$Description

write.table(isoform_up_enrichgo,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Overrepresentated_GO_terms_upregulated.csv",sep = ",",row.names = FALSE,col.names = TRUE)

Need_H_down_enrichGO <- enrichGO(gene = unique(unlist(Need_H_down)),OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "fdr",universe = Reference_genes_all,qvalueCutoff = 0.05,minGSSize = 50,maxGSSize = 1000,readable = TRUE)
orsize <- ifelse(as.numeric(unlist(lapply(str_split(Need_H_down_enrichGO@result$GeneRatio,pattern = "/",n = 2),function(x)x[[1]]))) > 7,yes = TRUE,no = FALSE)
test2 <- Need_H_down_enrichGO@result[orsize,c("ID","Description","p.adjust","geneID")]
colnames(test2) <- c("ID","Description", "Need_H","Need_H_genes")
terms2 <- filter(test2, Need_H < 0.01)$ID

Need_E_down_enrichGO <- enrichGO(gene = unique(unlist(Need_E_down)),OrgDb = "org.Mm.eg.db",keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "fdr",universe = Reference_genes_all,qvalueCutoff = 0.05,minGSSize = 50,maxGSSize = 1000,readable = TRUE)
orsize <- ifelse(as.numeric(unlist(lapply(str_split(Need_E_down_enrichGO@result$GeneRatio,pattern = "/",n = 2),function(x)x[[1]]))) > 7,yes = TRUE,no = FALSE)
test3 <- Need_E_down_enrichGO@result[orsize,c("ID","Description","p.adjust","geneID")]
colnames(test3) <- c("ID","Description", "Need_E","Need_E_genes")
terms3 <- filter(test3, Need_E < 0.01)$ID


test <- full_join(test2,test3,by = c("ID","Description"))
isoform_down_enrichgo <- filter(test, ID %in% unique(c(terms2,terms3)))
isoform_down_enrichgo <- filter(isoform_down_enrichgo, ! ID %in% avoidGO)
rownames(isoform_down_enrichgo) <- isoform_down_enrichgo$Description

write.table(isoform_down_enrichgo,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Overrepresentated_GO_terms_downregulated.csv",sep = ",",row.names = FALSE,col.names = TRUE)

saveRDS(list("Need_H_up_enrichGO"=Need_H_up_enrichGO,"Need_E_up_enrichGO"=Need_E_up_enrichGO,"Need_H_down_enrichGO"=Need_H_down_enrichGO,"Need_E_down_enrichGO"=Need_E_down_enrichGO),file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/EnrichGo_Hif1a_Epas1_specific_genes.rds")


# Fig. 4b -----------------------------------------------------------------

"Getting list of HIF1A specific and HIF2A specific genes that are members of each significantly enriched GO term"

Hif2_enriched_genes_up <- as.list(arrange(filter(isoform_up_enrichgo,Need_E < 0.01),Need_E)[["Need_E_genes"]])
Hif2_enriched_genes_up <- lapply(Hif2_enriched_genes_up, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif2_enriched_genes_up <- lapply(Hif2_enriched_genes_up,as.character)
names(Hif2_enriched_genes_up) <- as.character(arrange(filter(isoform_up_enrichgo,Need_E < 0.01),Need_E)[["Description"]])

Hif1_enriched_genes_up <- as.list(arrange(filter(isoform_up_enrichgo,Need_H < 0.01),Need_H)[["Need_H_genes"]])
Hif1_enriched_genes_up <- lapply(Hif1_enriched_genes_up, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif1_enriched_genes_up <- lapply(Hif1_enriched_genes_up,as.character)
names(Hif1_enriched_genes_up) <- as.character(arrange(filter(isoform_up_enrichgo,Need_H < 0.01),Need_H)[["Description"]])

intersect(names(Hif2_enriched_genes_up),names(Hif1_enriched_genes_up))

H_E_GO_terms_l2fc <- matrix(nrow = nrow(isoform_up_enrichgo),ncol = 2)
rownames(H_E_GO_terms_l2fc) <- rownames(isoform_up_enrichgo)
colnames(H_E_GO_terms_l2fc) <- c("VHKO\nvs VKO","VEKO\nvs VKO")
for(i in rownames(H_E_GO_terms_l2fc)){
  genes <- c(Hif1_enriched_genes_up[[i]],Hif2_enriched_genes_up[[i]])
  l2fc_H <- c()
  l2fc_E <- c()
  for(j in genes){
    H <- 2^filter(V_HE_H_E_FC_hif_df, !is.na(HIF_dependence), Gene == j)[["VHKO vs VKO - L2FC"]]
    l2fc_H <- c(l2fc_H,H)
  }
  for(j in genes){
    E <- 2^filter(V_HE_H_E_FC_hif_df, !is.na(HIF_dependence), Gene == j)[["VEKO vs VKO - L2FC"]]
    l2fc_E <- c(l2fc_E,E)
  }
  l2fc_H <- log(mean(l2fc_H),base = 2)
  l2fc_E <- log(mean(l2fc_E),base = 2)
  H_E_GO_terms_l2fc[i,"VHKO\nvs VKO"] <- l2fc_H
  H_E_GO_terms_l2fc[i,"VEKO\nvs VKO"] <- l2fc_E
}

test6 <- H_E_GO_terms_l2fc
"Some GO term descriptions are too long"
rownames(test6)[nchar(rownames(test6)) > 30] <- c("ribonucleotide catabolism", "response to decreased oxygen", "nucleoside phosphate catabolism", "organophosphate catabolism","carboxylate metabolism", "monocarboxylate metabolism","carbohydrate derivative catabolism", "monosaccharide metabolism","generating precursor metabolites", "ribonucleotide metabolism","ribose phosphate metabolism", "purine metabolism","nucleoside metabolism", "nucleobase metabolism","purine metabolic process", "nucleobase catabolism","organophosphate metabolism", "purine-containing metabolism","carbohydrate metabolic process", "nitrogen catabolism","anatomical morphogenesis", "striated muscle differentiation", "organic acid transport", "carboxylate transmembrane transport", "organic acid biosynthesis", "negative regulation of transport","ECM organization", "ECM structure organization","external capsule organization")

scale_matrix <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 2, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 2, max_abs_values, "/")
  return(scaled_mat)
}

"Final plot"
clustrow_l2fc <- hclust(dist(scale_matrix(test6)),method = "ward.D2")
pheatmap::pheatmap(mat = as.matrix(test6),color = colorRampPalette(rev(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58")))(100),cluster_cols = FALSE,cluster_rows = clustrow_l2fc,fontsize = 6,height = 6.5,width = 3.5,treeheight_row = 10,show_colnames = TRUE,angle_col = 0,legend = TRUE,fontsize_col = 6.8,fontsize_row = 6.5,breaks = seq(-2,2,length.out = 101),legend_breaks = c(-2,-1,0,0.5),filename = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4b/Fig_4b.pdf")




# Fig. 4e -----------------------------------------------------------------

"Getting list of HIF1A specific and HIF2A specific genes that are members of each significantly enriched GO term"

Hif2_enriched_genes_down <- as.list(arrange(filter(isoform_down_enrichgo,Need_E < 0.01),Need_E)[["Need_E_genes"]])
Hif2_enriched_genes_down <- lapply(Hif2_enriched_genes_down, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif2_enriched_genes_down <- lapply(Hif2_enriched_genes_down,as.character)
names(Hif2_enriched_genes_down) <- as.character(arrange(filter(isoform_down_enrichgo,Need_E < 0.01),Need_E)[["Description"]])

Hif1_enriched_genes_down <- as.list(arrange(filter(isoform_down_enrichgo,Need_H < 0.01),Need_H)[["Need_H_genes"]])
Hif1_enriched_genes_down <- lapply(Hif1_enriched_genes_down, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif1_enriched_genes_down <- lapply(Hif1_enriched_genes_down,as.character)
names(Hif1_enriched_genes_down) <- as.character(arrange(filter(isoform_down_enrichgo,Need_H < 0.01),Need_H)[["Description"]])

intersect(names(Hif2_enriched_genes_down),names(Hif1_enriched_genes_down))

H_E_GO_terms_l2fc <- matrix(nrow = nrow(isoform_down_enrichgo),ncol = 2)
rownames(H_E_GO_terms_l2fc) <- rownames(isoform_down_enrichgo)
colnames(H_E_GO_terms_l2fc) <- c("VHKO\nvs VKO","VEKO\nvs VKO")
for(i in rownames(H_E_GO_terms_l2fc)){
  genes <- c(Hif1_enriched_genes_down[[i]],Hif2_enriched_genes_down[[i]])
  l2fc_H <- c()
  l2fc_E <- c()
  for(j in genes){
    H <- 2^filter(V_HE_H_E_FC_hif_df, !is.na(HIF_dependence), Gene == j)[["VHKO vs VKO - L2FC"]]
    l2fc_H <- c(l2fc_H,H)
  }
  for(j in genes){
    E <- 2^filter(V_HE_H_E_FC_hif_df, !is.na(HIF_dependence), Gene == j)[["VEKO vs VKO - L2FC"]]
    l2fc_E <- c(l2fc_E,E)
  }
  l2fc_H <- log(mean(l2fc_H),base = 2)
  l2fc_E <- log(mean(l2fc_E),base = 2)
  H_E_GO_terms_l2fc[i,"VHKO\nvs VKO"] <- l2fc_H
  H_E_GO_terms_l2fc[i,"VEKO\nvs VKO"] <- l2fc_E
}

test6 <- H_E_GO_terms_l2fc
"Some GO term descriptions are too long"
rownames(test6)[nchar(rownames(test6)) > 30] <- c("carboxylate metabolism", "cellular xenobiotic response","xenobiotic stimulus response", "sulfur compound metabolism","olefinic compound metabolism", "monocarboxylate metabolism","cellular lipid metabolism", "unsaturated fat metabolism","vascular process", "long-chain fat metabolism"
)

scale_matrix <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 2, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 2, max_abs_values, "/")
  return(scaled_mat)
}

"Final plot"
clustrow_l2fc <- hclust(dist(scale_matrix(test6)),method = "ward.D2")
pheatmap::pheatmap(mat = as.matrix(test6),color = colorRampPalette(c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858',"#081D58"))(100),cluster_cols = FALSE,cluster_rows = clustrow_l2fc,fontsize = 6,height = 6/2.54,width = 3.5,treeheight_row = 10,show_colnames = TRUE,angle_col = 0,legend = TRUE,fontsize_col = 6.8,fontsize_row = 6.5,breaks = seq(-0.1,2,length.out = 101),legend_breaks = c(0,1,2),filename = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_4/Fig_4e/Fig_4e.pdf")


# Supplementary Fig 7 -----------------------------------------------------
isoform_up_enrichgo <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Overrepresentated_GO_terms_upregulated.csv")
Hif1_enriched_genes_up <- as.list(arrange(filter(isoform_up_enrichgo,Need_H < 0.01),Need_H)[["Need_H_genes"]])
Hif1_enriched_genes_up <- lapply(Hif1_enriched_genes_up, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif1_enriched_genes_up <- lapply(Hif1_enriched_genes_up,as.character)
names(Hif1_enriched_genes_up) <- as.character(arrange(filter(isoform_up_enrichgo,Need_H < 0.01),Need_H)[["Description"]])

makeupset <- function(x){
  genes <- unique(unlist(x))
  gene_matrix <- matrix(nrow = length(genes),ncol = length(x))
  rownames(gene_matrix) <- genes
  colnames(gene_matrix) <- names(x)
  for(i in rownames(gene_matrix)){for(j in colnames(gene_matrix)){
    test <- ifelse(test = length(x[[j]][x[[j]] == i]) == 1,yes = 1,no = 0)
    gene_matrix[i,j] <- test
  }}
  return(as.data.frame(gene_matrix))
}
Hif1_enriched_genes_mat <- makeupset(Hif1_enriched_genes_up)

scale_matrix_row <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 2, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 2, max_abs_values, "/")
  return(scaled_mat)
}
scale_matrix_col <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 1, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 1, max_abs_values, "/")
  return(scaled_mat)
}

colnames(Hif1_enriched_genes_mat)[nchar(colnames(Hif1_enriched_genes_mat)) > 30] <- c("ribonucleotide catabolism", "response to decreased oxygen", "nucleoside phosphate catabolism", "organophosphate catabolism","carboxylate metabolism", "monocarboxylate metabolism","carbohydrate derivative catabolism", "monosaccharide metabolism","generating precursor metabolites", "ribonucleotide metabolism","ribose phosphate metabolism", "purine metabolism","nucleoside metabolism", "nucleobase metabolism","purine metabolic process", "nucleobase catabolism","organophosphate metabolism", "purine-containing metabolism","carbohydrate metabolic process", "nitrogen catabolism")


clustrow <- hclust(dist(scale_matrix_row(as.matrix(Hif1_enriched_genes_mat[,]))),method = "ward.D2")
clustcol <- hclust(dist(t(scale_matrix_col(as.matrix(Hif1_enriched_genes_mat[,])))),method = "ward.D2")
pheatmap::pheatmap(mat = as.matrix(Hif1_enriched_genes_mat[,]),fontsize = 4,cluster_cols = clustcol,cluster_rows = clustrow,scale = "none",color = colorRampPalette(c('#fff7f7',"#081D58"))(2),legend = FALSE,                 filename = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_7/Supplementary_Fig_7.pdf",width = 7,height = 9,treeheight_col = 10,treeheight_row = 10,border_color = NA)


# Supplementary Fig. 8---------------------------
isoform_up_enrichgo <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Overrepresentated_GO_terms_upregulated.csv")
Hif2_enriched_genes_up <- as.list(arrange(filter(isoform_up_enrichgo,Need_E < 0.01),Need_E)[["Need_E_genes"]])
Hif2_enriched_genes_up <- lapply(Hif2_enriched_genes_up, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif2_enriched_genes_up <- lapply(Hif2_enriched_genes_up,as.character)
names(Hif2_enriched_genes_up) <- as.character(arrange(filter(isoform_up_enrichgo,Need_E < 0.01),Need_E)[["Description"]])

makeupset <- function(x){
  genes <- unique(unlist(x))
  gene_matrix <- matrix(nrow = length(genes),ncol = length(x))
  rownames(gene_matrix) <- genes
  colnames(gene_matrix) <- names(x)
  for(i in rownames(gene_matrix)){for(j in colnames(gene_matrix)){
    test <- ifelse(test = length(x[[j]][x[[j]] == i]) == 1,yes = 1,no = 0)
    gene_matrix[i,j] <- test
  }}
  return(as.data.frame(gene_matrix))
}
Hif2_enriched_genes_mat <- makeupset(Hif2_enriched_genes_up)

scale_matrix_row <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 2, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 2, max_abs_values, "/")
  return(scaled_mat)
}
scale_matrix_col <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 1, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 1, max_abs_values, "/")
  return(scaled_mat)
}

colnames(Hif2_enriched_genes_mat)[nchar(colnames(Hif2_enriched_genes_mat)) > 30] <- c("striated muscle differentiation","anatomical morphogenesis", "organic acid transport", "carboxylate transmembrane transport","carboxylate metabolism", "organic acid biosynthesis", "negative regulation of transport","ECM organization", "ECM structure organization","external capsule organization")


clustrow <- hclust(dist(scale_matrix_row(as.matrix(Hif2_enriched_genes_mat[,]))),method = "ward.D2")
clustcol <- hclust(dist(t(scale_matrix_col(as.matrix(Hif2_enriched_genes_mat[,])))),method = "ward.D2")
pheatmap::pheatmap(mat = as.matrix(Hif2_enriched_genes_mat[,]),fontsize = 4,cluster_cols = clustcol,cluster_rows = clustrow,scale = "none",color = colorRampPalette(c('#fff7f7',"#081D58"))(2),legend = FALSE,                 filename = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_8/Supplementary_Fig_8.pdf",width = 7,height = 9,treeheight_col = 10,treeheight_row = 10,border_color = NA)



# Supplementary Fig 9 -----------------------------------------------------

isoform_down_enrichgo <- read_csv(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Overrepresentated_GO_terms_downregulated.csv")
Hif2_enriched_genes_down <- as.list(arrange(filter(isoform_down_enrichgo,Need_E < 0.01),Need_E)[["Need_E_genes"]])
Hif2_enriched_genes_down <- lapply(Hif2_enriched_genes_down, function(x) x <- str_split(string = x,pattern = "/",simplify = TRUE))
Hif2_enriched_genes_down <- lapply(Hif2_enriched_genes_down,as.character)
names(Hif2_enriched_genes_down) <- as.character(arrange(filter(isoform_down_enrichgo,Need_E < 0.01),Need_E)[["Description"]])

makedownset <- function(x){
  genes <- unique(unlist(x))
  gene_matrix <- matrix(nrow = length(genes),ncol = length(x))
  rownames(gene_matrix) <- genes
  colnames(gene_matrix) <- names(x)
  for(i in rownames(gene_matrix)){for(j in colnames(gene_matrix)){
    test <- ifelse(test = length(x[[j]][x[[j]] == i]) == 1,yes = 1,no = 0)
    gene_matrix[i,j] <- test
  }}
  return(as.data.frame(gene_matrix))
}
Hif2_enriched_genes_mat <- makedownset(Hif2_enriched_genes_down)

scale_matrix_row <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 2, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 2, max_abs_values, "/")
  return(scaled_mat)
}
scale_matrix_col <- function(mat) {
  centered_mat <- scale(mat, center = TRUE, scale = FALSE)
  max_abs_values <- apply(centered_mat, 1, function(x) max(abs(x)))
  scaled_mat <- sweep(centered_mat, 1, max_abs_values, "/")
  return(scaled_mat)
}

colnames(Hif2_enriched_genes_mat)[nchar(colnames(Hif2_enriched_genes_mat)) > 30] <- c("carboxylate metabolism", "cellular xenobiotic response","xenobiotic stimulus response", "sulfur compound metabolism","olefinic compound metabolism", "monocarboxylate metabolism","cellular lipid metabolism", "unsaturated fat metabolism","vascular process", "long-chain fat metabolism")


clustrow <- hclust(dist(scale_matrix_row(as.matrix(Hif2_enriched_genes_mat[,]))),method = "ward.D2")
clustcol <- hclust(dist(t(scale_matrix_col(as.matrix(Hif2_enriched_genes_mat[,])))),method = "ward.D2")
pheatmap::pheatmap(mat = as.matrix(Hif2_enriched_genes_mat[,]),fontsize = 4,cluster_cols = clustcol,cluster_rows = clustrow,scale = "none",color = colorRampPalette(c('#fff7f7',"#081D58"))(2),legend = FALSE,                 filename = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_9/Supplementary_Fig_9.pdf",width = 7,height = 9,treeheight_col = 10,treeheight_row = 10,border_color = NA)

# Hypergeomteric test for PT marker genes--------------------------------------------------------

"Getting PT Markers"

Meta_modules_unique_df <- read_csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_type_modules_filtered_metamarkers.csv")
Meta_modules_unique <- as.list(Meta_modules_unique_df)
Meta_modules_unique <- lapply(Meta_modules_unique,function(x){x <-x[!is.na(x)]})
ptmarkers <- unique(unlist(lapply(Meta_modules_unique[2:6],function(x){x <-x[!is.na(x)]})))


"Getiing list of differentially regulated genes"

V_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VKOvConKO_DE_pseudobulk.rds")
H_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHKOvVKO_DE_pseudobulk.rds")
E_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VEKOvVKO_DE_pseudobulk.rds")
HE_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VHEKOvVKO_DE_pseudobulk.rds")
test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())
"`````````````````````````````````````````````````````````````````````````````````"
"Performing a hypergeomtric test for enrichment - https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/"

universe <- Reduce(intersect,list(unlist(lapply(V_DE_list,function(x){x[["gene"]]})),unlist(lapply(HE_DE_list,function(x){x[["gene"]]})),unlist(lapply(H_DE_list,function(x){x[["gene"]]})),unlist(lapply(E_DE_list,function(x){x[["gene"]]}))))

E_down_genes <- unique(unlist(Need_E_down))

PT_E_p <- phyper(q = length(intersect(E_down_genes,ptmarkers))-1,m = length(ptmarkers),n = length(universe) - length(ptmarkers),k = length(E_down_genes),lower.tail = FALSE)

# Fig 5a ------------------------------------------------------------------
Meta_modules_unique_df <- read_csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_type_modules_filtered_metamarkers.csv")
Meta_modules_unique <- as.list(Meta_modules_unique_df)
Meta_modules_unique <- lapply(Meta_modules_unique,function(x){x <-x[!is.na(x)]})

PT_score_plots <- list()

for(i in levels(VHE_late_pos$Cell_types_MetaUCell)[1:3]){
  if(i == "PT S1"){
    test <- AddModuleScore(subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell == i]),features = Meta_modules_unique[c(2)],nbin = 100,ctrl = 50,name = "Module_")
    test <- FetchData(test,c("Module_1","Condition","PT_identity"))
    colnames(test) <- c("Score","Condition","PT_identity")
  }
  else if(i == "PT S2"){
    test1 <- AddModuleScore(subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell == i & VHE_late_pos$Sex == "Male"]),features = Meta_modules_unique[c(3)],nbin = 100,ctrl = 50,name = "Module_")
    test1 <- FetchData(test1,c("Module_1","Condition","PT_identity"))
    test2 <- AddModuleScore(subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell == i & VHE_late_pos$Sex == "Female"]),features = Meta_modules_unique[c(5)],nbin = 100,ctrl = 50,name = "Module_")
    test2 <- FetchData(test2,c("Module_1","Condition","PT_identity"))
    colnames(test1) <- c("Score","Condition","PT_identity")
    colnames(test2) <- c("Score","Condition","PT_identity")
    test <- rbind(test1,test2)
    colnames(test) <- c("Score","Condition","PT_identity")
  }
  else if(i == "PT S3"){
    test1 <- AddModuleScore(subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell == i & VHE_late_pos$Sex == "Male"]),features = Meta_modules_unique[c(4)],nbin = 100,ctrl = 50,name = "Module_")
    test1 <- FetchData(test1,c("Module_1","Condition","PT_identity"))
    test2 <- AddModuleScore(subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$Cell_types_MetaUCell == i & VHE_late_pos$Sex == "Female"]),features = Meta_modules_unique[c(6)],nbin = 100,ctrl = 50,name = "Module_")
    test2 <- FetchData(test2,c("Module_1","Condition","PT_identity"))
    colnames(test1) <- c("Score","Condition","PT_identity")
    colnames(test2) <- c("Score","Condition","PT_identity")
    test <- rbind(test1,test2)
    colnames(test) <- c("Score","Condition","PT_identity")
  }
  summary_data <- test %>%
    group_by(Condition) %>%
    summarise(lower = quantile(Score, 0.05),q1 = quantile(Score, 0.25),median = median(Score),q3 = quantile(Score, 0.75),upper = quantile(Score, 0.95),.groups = "drop")
  p1 <- ggplot(summary_data, aes(x = Condition)) +
    geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
    geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
    scale_color_manual(values = condition_colour_2)+
    labs(title = paste(i,"cells"),y = paste(i,"markers expression score"),x=NULL)+
    theme_classic()+
    theme(strip.text = element_text(size = 7),strip.background = element_blank(),legend.title = element_blank(),legend.text = element_blank(),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 0.5,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"),plot.title = element_text(size = 7,hjust = 0.5,vjust = 0.5))
  p2 <- p1+NoLegend()
  PT_score_plots[[i]] <- p2
}
p1 <- plot_grid(PT_score_plots$`PT S1`,PT_score_plots$`PT S2`,PT_score_plots$`PT S3`,nrow = 1)
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5a/Fig_5a.pdf",width = 18/2.54,height = 2)
p1
dev.off()

# Fig 5c ------------------------------------------------------------------

"Lists of genes that are specifically regulated by HIF1A or HIF2A were used an input in LISA and CheA3 analysis online. The results were exported as .tsv files that were then saved as .xlsx files"

Lisa_H <- readxl::read_xlsx("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5c/Need_H_Lisa_results.xlsx")
colnames(Lisa_H) <- c("TF", "Rank", "1st Sample p-value", "2nd Sample p-value", "3rd Sample p-value", "4th Sample p-value", "5th Sample p-value")
Chea3_H <- readxl::read_xlsx("//Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5c/Need_H_ChEA3_results.xlsx")
Lisa_Chea3_H <- full_join(Lisa_H[,c("TF","Rank")],Chea3_H[,c("TF","Rank")],by = "TF")
colnames(Lisa_Chea3_H) <- c("TF","Lisa","ChEA3")
Lisa_Chea3_H$Lisa[is.na(Lisa_Chea3_H$Lisa)] <- nrow(Lisa_H)+1
Lisa_Chea3_H$ChEA3[is.na(Lisa_Chea3_H$ChEA3)] <- nrow(Chea3_H)+1

TF_H <- list("LISA" = filter(Lisa_Chea3_H, Lisa < 51)[["TF"]],"CHEA3" = filter(Lisa_Chea3_H, ChEA3 < 51)[["TF"]])
Consensus_TF_H <- Reduce(intersect, TF_H)

Lisa_E <- readxl::read_xlsx("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5c/Need_E_Lisa_results.xlsx")
colnames(Lisa_E) <- c("TF", "Rank", "1st Sample p-value", "2nd Sample p-value", "3rd Sample p-value", "4th Sample p-value", "5th Sample p-value")
Chea3_E <- readxl::read_xlsx("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5c/Need_E_ChEA3_results.xlsx")
Lisa_Chea3_E <- full_join(Lisa_E[,c("TF","Rank")],Chea3_E[,c("TF","Rank")],by = "TF")
colnames(Lisa_Chea3_E) <- c("TF","Lisa","ChEA3")
Lisa_Chea3_E$Lisa[is.na(Lisa_Chea3_E$Lisa)] <- nrow(Lisa_E)+1
Lisa_Chea3_E$ChEA3[is.na(Lisa_Chea3_E$ChEA3)] <- nrow(Chea3_E)+1

TF_E <- list("LISA" = filter(Lisa_Chea3_E, Lisa < 51)[["TF"]],"CHEA3" = filter(Lisa_Chea3_E, ChEA3 < 51)[["TF"]])
Consensus_TF_E <- Reduce(intersect, TF_E)

Lisa <- full_join(Lisa_H[,c(1,3)],Lisa_E[,c(1,3)],by = "TF",suffix = c("H","E"))
colnames(Lisa) <- c("TF","H_LISA","E_LISA")
Lisa$H_LISA <- as.numeric(Lisa$H_LISA)
Lisa$E_LISA <- as.numeric(Lisa$E_LISA)

CheA3 <- full_join(Chea3_H[,c(3,4)],Chea3_E[,c(3,4)],by = "TF",suffix = c("H","E"))
colnames(CheA3) <- c("TF","H_CheA3","E_CheA3")


CheA3_Lisa <- full_join(CheA3,Lisa,by = "TF")
CheA3_Lisa$H_CheA3[is.na(CheA3_Lisa$H_CheA3)] <- 1583
CheA3_Lisa$E_CheA3[is.na(CheA3_Lisa$E_CheA3)] <- 1583
CheA3_Lisa$H_LISA[is.na(CheA3_Lisa$H_LISA)] <- 1
CheA3_Lisa$E_LISA[is.na(CheA3_Lisa$E_LISA)] <- 1

p1 <- ggplot()+
  geom_point(data = CheA3_Lisa,aes(y = -1*log10(H_LISA),x = -1*log10(H_CheA3)),color = "grey",size = 0.8)+
  geom_point(data = filter(CheA3_Lisa,TF %in% Consensus_TF_H), aes(y = -1*log10(H_LISA),x = -1*log10(H_CheA3)),color = "#bc8e54",size = 1.2)+
  geom_text_repel(data = filter(CheA3_Lisa,TF %in% Consensus_TF_H), aes(y = -1*log10(H_LISA),x = -1*log10(H_CheA3),label = TF), colour = "#bc8e54",max.overlaps = 10,force = 2,force_pull = 1, size = 2.5,ylim = c(10,30),xlim = c(-2,-1))+
  labs(y = "-log(p) - LISA",x = "-log(Score) - CheA3",title = "TFs for HIF1A-specific genes")+
  theme_classic()+
  theme(plot.title = element_text(size = 11,hjust = 0.5),axis,title = element_text(size = 9),axis.text = element_text(size = 7))

p2 <- ggplot()+
  geom_point(data = CheA3_Lisa,aes(y = -1*log10(E_LISA),x = -1*log10(E_CheA3)),color = "grey",size = 0.8)+
  geom_point(data = filter(CheA3_Lisa,TF %in% Consensus_TF_E), aes(y = -1*log10(E_LISA),x = -1*log10(E_CheA3)),color = "#208d4c",size = 1.2)+
  geom_text_repel(data = filter(CheA3_Lisa,TF %in% Consensus_TF_E), aes(y = -1*log10(E_LISA),x = -1*log10(E_CheA3),label = TF), colour = "#208d4c",max.overlaps = 10,force = 10,force_pull = 1, size = 2.5,ylim = c(10,30),xlim = c(-2,0))+
  labs(y = "-log(p) - LISA",x = "-log(Score) - CheA3",title = "TFs for HIF2A-specific genes")+
  theme_classic()+
  theme(plot.title = element_text(size = 11,hjust = 0.5),axis,title = element_text(size = 9),axis.text = element_text(size = 7))

p3 <- plot_grid(p1+NoLegend(),p2+NoLegend(),ncol = 2)
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_5/Fig_5c/Fig_5c.pdf",width = 18/2.54,height = 2)
p3
dev.off()

# Supplementary Fig 10 ----------------------------------------------------
"Getting lists of HIF1 and HIF2 dependent genes from Schonenberger et al 2016 Cancer Res. This is microarray data from primary RPTECs with Vhl, Hif1a, and Epas1 recombined by adenoviral cre and analysed 4 days after recombination. Only upregulated genes have been reported"

Scho_genes <- readxl::read_xlsx(path = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_10/Schonenberger_RPTEC_hif1_hif2.xlsx",sheet = 2,col_names = TRUE)

Scho_genes_list <- lapply(as.list(Scho_genes),function(x)x <- x[!is.na(x)])

"GSEA on our data"

V_DE_list <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/VHE_late_pos_VKOvConKO_DE_pseudobulk.rds")

V_DE_rank <- lapply(V_DE_list, function(x){
  x <- arrange(x,desc(log2FoldChange))
  values <- as.numeric(x[["log2FoldChange"]])
  names(values) <- as.character(x[["gene"]])
  return(values)
})

V_DE_scho_gsea <- lapply(V_DE_rank, function(x){
  fgsea::fgsea(pathways = Scho_genes_list[1:2],stats = x)
})

scho_gsea_hif1_plots <- Map(x = V_DE_rank,z = names(V_DE_rank),f = function(x,z){fgsea::plotEnrichment(pathway = Scho_genes_list[[1]],stats = x)+labs(x=NULL,y="Enrichement Score",title = z)+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(size = 8,margin = unit(c(0,0,0,0),"cm")))})
scho_gsea_hif2_plots <- Map(x = V_DE_rank,z = names(V_DE_rank),f = function(x,z){fgsea::plotEnrichment(pathway = Scho_genes_list[[2]],stats = x)+labs(x=NULL,y="Enrichement Score",title = z)+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(size = 8,margin = unit(c(0,0,0,0),"cm")))})

V_DE_scho_gsea_NES_p <- lapply(V_DE_scho_gsea,function(x){
  x <- x[,c("pathway","padj","NES")]
})
V_DE_scho_gsea_NES_p <- bind_rows(V_DE_scho_gsea_NES_p,.id = "PT_identity")

colnames(V_DE_scho_gsea_NES_p) <- c("PT Identity","Geneset","p","NES")

write.table(arrange(filter(V_DE_scho_gsea_NES_p,`PT Identity` != "Non PT"),Geneset,`PT Identity`),file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_10/Schonenberger_HIF_genes_NES_p_values.csv",sep = ",",row.names = FALSE,col.names = TRUE)

p1 <- plot_grid(scho_gsea_hif1_plots$`S1 A`,scho_gsea_hif1_plots$`S1 B`,scho_gsea_hif1_plots$`S2 A`,scho_gsea_hif1_plots$`S2 B`,scho_gsea_hif1_plots$`S3 A`,scho_gsea_hif1_plots$`S3 B`,ncol = 3)
p2 <- plot_grid(scho_gsea_hif2_plots$`S1 A`,scho_gsea_hif2_plots$`S1 B`,scho_gsea_hif2_plots$`S2 A`,scho_gsea_hif2_plots$`S2 B`,scho_gsea_hif2_plots$`S3 A`,scho_gsea_hif2_plots$`S3 B`,ncol = 3)
p3 <- plot_grid(p1,p2,nrow = 2)
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_10/Supplementary_Fig_10_HIF1A.pdf",width = 17/2.54,height = 8/2.54)
p1
dev.off()
pdf(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_10/Supplementary_Fig_10_HIF2A.pdf",width = 17/2.54,height = 8/2.54)
p2
dev.off()

# Early and Adaptive Gene Expression --------------------------------------

Vhl_Pax8_pos_batch_sex<- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Old_VKO_ConKO_Early_Late_Data.rds")

Vhl_Pax8_pos_batch_sex[["Sample"]] <- factor(Vhl_Pax8_pos_batch_sex$sample, levels = c("OY725_pos", "OY954_pos", "OY955_pos", "OY1004_pos", "OY733_pos", "OY978_pos", "OY1166_pos", "OY1021_pos", "OY552_pos", "OY611_pos","OY545_pos","OY622_pos",  "OY606_pos","OY626_pos", "OY615_pos", "OY1023_pos"))

#Adding sex as factor
Idents(Vhl_Pax8_pos_batch_sex) <- "Sample"
Female <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY725_pos", "OY954_pos", "OY955_pos","OY733_pos", "OY978_pos","OY552_pos", "OY611_pos", "OY545_pos", "OY606_pos", "OY626_pos"))
Male <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY1004_pos","OY1021_pos","OY1166_pos", "OY622_pos", "OY615_pos", "OY1023_pos"))
Idents(Vhl_Pax8_pos_batch_sex, cells = Female) <- "Female"
Idents(Vhl_Pax8_pos_batch_sex, cells = Male) <- "Male"
Vhl_Pax8_pos_batch_sex[["Sex"]] <- factor(Idents(Vhl_Pax8_pos_batch_sex), levels = c("Male", "Female"))
Vhl_Pax8_pos_batch_sex$Sex <- factor(Vhl_Pax8_pos_batch_sex$Sex, levels = c("Male", "Female"))

#Adding Condition as factor
Idents(Vhl_Pax8_pos_batch_sex) <- "Sample"
ConEp <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY725_pos", "OY954_pos", "OY955_pos", "OY1004_pos"))
KOEp <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY733_pos", "OY978_pos", "OY1021_pos", "OY1166_pos"))
ConLp <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY552_pos", "OY611_pos", "OY545_pos", "OY622_pos"))
KOLp <- WhichCells(Vhl_Pax8_pos_batch_sex, idents = c("OY606_pos", "OY626_pos", "OY615_pos", "OY1023_pos"))

Idents(Vhl_Pax8_pos_batch_sex, cells = ConEp) <- "Control Early Pos"
Idents(Vhl_Pax8_pos_batch_sex, cells = KOEp) <- "KO Early Pos"
Idents(Vhl_Pax8_pos_batch_sex, cells = ConLp) <- "Control Late Pos"
Idents(Vhl_Pax8_pos_batch_sex, cells = KOLp) <- "KO Late Pos"

Vhl_Pax8_pos_batch_sex[["Condition"]] <- factor(Idents(Vhl_Pax8_pos_batch_sex), levels = c("Control Early Neg","Control Late Neg","KO Early Neg","KO Late Neg","Control Early Pos","Control Late Pos", "KO Early Pos", "KO Late Pos"))

Idents(Vhl_Pax8_pos_batch_sex) <- "Condition"
DimPlot(Vhl_Pax8_pos_batch_sex)

"```````````Assigning Cell Type``````````````````````````````````````````````````"

Meta_modules_unique_df <- read_csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Cell_type_modules_filtered_metamarkers.csv")

cell_types <- names(Meta_modules_unique) #metamarkers that are also modules and filtered for unique genes
cell_types_meta <- paste(cell_types,"MetaUCell",sep = "_")

Vhl_Pax8_pos_batch_sex<- AddModuleScore_UCell(Vhl_Pax8_pos_batch_sex, features = Meta_modules_unique, name = "_MetaUCell",assay = "RNA",maxRank = 1500,ncores = 10,slot = "data",force.gc = TRUE)

# FINDING MAXIMUM SCORES
Cell_type_scores <- FetchData(Vhl_Pax8_pos_batch_sex, c("Sex", cell_types_meta))
cell_types_male <- cell_types_meta[c(1:4, 7:20)]
cell_types_female <- cell_types_meta[c(1:2, 5:20)]

# Use vectorized operations to find max
get_max <- function(y) {
  if (y["Sex"] == "Male") {
    max_indices <- which.max(y[cell_types_male])
    max_cell_type <- cell_types_male[max_indices]
  } else {
    max_indices <- which.max(y[cell_types_female])
    max_cell_type <- cell_types_female[max_indices]
  }
  return(max_cell_type)
}

# Filling in empty dataframe with max score cell type for each cell
Cells_celltype <- as.data.frame(matrix(nrow = nrow(Cell_type_scores), ncol = 1))
rownames(Cells_celltype) <- rownames(Cell_type_scores)
colnames(Cells_celltype) <- "Cellscore"
Cells_celltype[, 1] <- apply(Cell_type_scores, 1, get_max)

# Converting from module names to cell type names.
celltypenames <- data.frame(
  Cellscore = cell_types_meta,
  Cell_type = c("PEC", "PT S1", "PT S2", "PT S3", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte", "Pericyte", "Endothelial","Fibroblast", "Macrophage", "NK Cell", "Neutrophil", "Monocyte", "B Cell", "T Cell")
)

Cells_celltype_2 <- rownames_to_column(Cells_celltype,var = "Cell")
Cells_celltype_2 <- left_join(Cells_celltype_2, celltypenames, by = "Cellscore")

Cells_celltype_2$Cell_type <- factor(
  Cells_celltype_2$Cell_type,
  levels = c("PT S1", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte", "PEC", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "Neutrophil", "NK Cell", "Monocyte", "B Cell", "T Cell")
)

Cell_type_MetaUCell_idents <- Cells_celltype_2
rownames(Cell_type_MetaUCell_idents) <- Cell_type_MetaUCell_idents$Cell 

Idents(Vhl_Pax8_pos_batch_sex) <- Cell_type_MetaUCell_idents[,3]
Vhl_Pax8_pos_batch_sex[["Cell_types_MetaUCell"]] <- factor(Idents(Vhl_Pax8_pos_batch_sex), levels = c("PT S1", "PT S2", "PT S3", "LoH", "DCT", "CDPC", "CDIC", "Podocyte","PEC", "Fibroblast", "Endothelial", "Pericyte", "Macrophage", "Neutrophil", "NK Cell", "Monocyte", "B Cell", "T Cell"))

"``````````````````Assigning PT Class````````````````````````````````````````````"

PT_modules_list <- readRDS(file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Control_PT_Class_Module_A_B_correlations.rds")

Module_A_ordered <- PT_modules_list[[1]]
Module_B_ordered <- PT_modules_list[[2]]

PT_Modules <- list(Module_A_ordered, Module_B_ordered)

DefaultAssay(Vhl_Pax8_pos_batch_sex) <- "RNA"
Vhl_Pax8_pos_batch_sex <- AddModuleScore(Vhl_Pax8_pos_batch_sex, features = PT_Modules, name = "PT_Module_", nbin = 100, ctrl = 50)

colnames(Vhl_Pax8_pos_batch_sex@meta.data)[colnames(Vhl_Pax8_pos_batch_sex@meta.data) %in% c("PT_Module_1","PT_Module_2")] <- c("PT_Module_A","PT_Module_B")

Tub_A <- WhichCells(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3") & Vhl_Pax8_pos_batch_sex$PT_Module_A > 0.125])
Tub_B <- WhichCells(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3") & Vhl_Pax8_pos_batch_sex$PT_Module_A < 0.125])
Non_tub <- WhichCells(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[!Vhl_Pax8_pos_batch_sex$Cell_types_MetaUCell %in% c("PT S1","PT S2","PT S3")])

Idents(Vhl_Pax8_pos_batch_sex, cells = Tub_A) <- "PT Class A"
Idents(Vhl_Pax8_pos_batch_sex, cells = Tub_B) <- "PT Class B"
Idents(Vhl_Pax8_pos_batch_sex, cells = Non_tub) <- "Non PT"

Vhl_Pax8_pos_batch_sex[["PT_Class_UCell"]] <- factor(Idents(Vhl_Pax8_pos_batch_sex), levels = c("PT Class A", "PT Class B", "Non PT"))

ptidents <- FetchData(Vhl_Pax8_pos_batch_sex,c("Cell_types_MetaUCell","PT_Class_UCell"))
ptidents[["PT_identity"]] <- apply(ptidents,1,function(x){
  if(x[["Cell_types_MetaUCell"]] %in% c("PT S1","PT S2","PT S3")){paste(x[["Cell_types_MetaUCell"]],x[["PT_Class_UCell"]])}
  else{"Non PT"}
})
ptidents$PT_identity <- factor(ptidents$PT_identity,levels = c("PT S1 PT Class A","PT S1 PT Class B","PT S2 PT Class A", "PT S2 PT Class B","PT S3 PT Class A", "PT S3 PT Class B","Non PT"),labels = c("S1 A","S1 B","S2 A","S2 B","S3 A","S3 B", "Non PT"))
Vhl_Pax8_pos_batch_sex[["PT_identity"]] <- ptidents[,3]

saveRDS(Vhl_Pax8_pos_batch_sex,file="/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_ConKO_Early_Late_Data.rds")

"`````````````````````DIFFERENTIAL GENE EXPRESSION````````````````````````````````"

sexcondition <- FetchData(Vhl_Pax8_pos_batch_sex,c("Sample","Sex","Condition"))%>%group_by(Sample,Sex,Condition)%>%summarise("Count" = n())
rownames(sexcondition) <- gsub(pattern = "_",replacement = "-",x = sexcondition$Sample)

"```````````````````````````````````````````````"
"VKO vs ConKO - Early"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(Vhl_Pax8_pos_batch_sex$PT_identity)){
  test <- subset(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$PT_identity == i & Vhl_Pax8_pos_batch_sex$Condition %in% c("KO Early Pos", "Control Early Pos")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_Early <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_Early <- sexcondition[colnames(cts_Early$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_Early,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","KO Early Pos","Control Early Pos"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

Early_DE_list <- lapply(cts_Early, deseqlist)
saveRDS(Early_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Vhl_Pax8_pos_batch_sex_VKOvConKO_Early_DE_pseudobulk.rds")

"``````````````````````````````````````````````````````````````````````````````"
"Late vs Early - KO"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(Vhl_Pax8_pos_batch_sex$PT_identity)){
  test <- subset(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$PT_identity == i & Vhl_Pax8_pos_batch_sex$Condition %in% c("KO Late Pos", "KO Early Pos")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_Adaptive_KO <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_Adaptive_KO <- sexcondition[colnames(cts_Adaptive_KO$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_Adaptive_KO,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","KO Late Pos","KO Early Pos"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

Adaptive_KO_DE_list <- lapply(cts_Adaptive_KO, deseqlist)
saveRDS(Adaptive_KO_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Vhl_Pax8_pos_batch_sex_VKO_LatevEarly_DE_pseudobulk.rds")

"``````````````````````````````````````````````````````````````````````````````"
"Late vs Early - Control"

"Making list of aggregated counts data for different PT types and classes by sample"
cts <- list()
for(i in levels(Vhl_Pax8_pos_batch_sex$PT_identity)){
  test <- subset(Vhl_Pax8_pos_batch_sex, cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$PT_identity == i & Vhl_Pax8_pos_batch_sex$Condition %in% c("Control Late Pos", "Control Early Pos")])
  test2 <- AggregateExpression(test,assays = "RNA",return.seurat = FALSE,group.by = "Sample")
  cts[[i]] <- test2
}
"Converting to matrices"
cts_Adaptive_Control <- lapply(cts, function(x){x <- as.matrix(x$RNA)})

"Creating condition and sex metadata columns for each sample"
colData_Adaptive_Control <- sexcondition[colnames(cts_Adaptive_Control$`S1 A`),c("Sex","Condition")]

"DESeq for each PT type and PT Class"
deseqlist <- function(x){
  dds <- DESeqDataSetFromMatrix(countData = x,colData = colData_Adaptive_Control,design = ~ Sex+ Condition)
  keep <- rowSums(counts(dds)) > 10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds,contrast = c("Condition","Control Late Pos","Control Early Pos"),independentFiltering = FALSE, alpha = 0.01,lfcThreshold= 0, altHypothesis="greaterAbs")
  res_tbl <- res %>%data.frame() %>%rownames_to_column(var = "gene") %>%as_tibble() %>%arrange(desc(log2FoldChange))
}

Adaptive_Control_DE_list <- lapply(cts_Adaptive_Control, deseqlist)
saveRDS(Adaptive_Control_DE_list,"/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Vhl_Pax8_pos_batch_sex_ConKO_LatevEarly_DE_pseudobulk.rds")

"`````````````````````````````````````````````````````````````````````````````"

Early_DE_up <- lapply(Early_DE_list, function(x){filter(x,log2FoldChange > 2,padj < 0.05)$gene})
Early_DE_down <- lapply(Early_DE_list, function(x){filter(x,log2FoldChange < -2,padj < 0.05)$gene})
Adaptive_KO_DE_up <- lapply(Adaptive_KO_DE_list, function(x){filter(x,log2FoldChange > 2,padj < 0.05)$gene})
Adaptive_KO_DE_down <- lapply(Adaptive_KO_DE_list, function(x){filter(x,log2FoldChange < -2,padj < 0.05)$gene})
Adaptive_Control_DE_up <- lapply(Adaptive_Control_DE_list, function(x){filter(x,log2FoldChange > 2,padj < 0.05)$gene})
Adaptive_Control_DE_down <- lapply(Adaptive_Control_DE_list, function(x){filter(x,log2FoldChange < -2,padj < 0.05)$gene})

Adaptive_DE_up <- setdifflist(Adaptive_KO_DE_up,Adaptive_Control_DE_up)
Adaptive_DE_down <- setdifflist(Adaptive_KO_DE_down,Adaptive_Control_DE_down)

Only_Early_up <- setdifflist(Early_DE_up,Adaptive_DE_up)
Only_Early_down <- setdifflist(Early_DE_down,Adaptive_DE_down)
Only_Adaptive_up <- setdifflist(Adaptive_DE_up,Early_DE_up)
Only_Adaptive_down <- setdifflist(Adaptive_DE_down,Early_DE_down)

saveRDS(list("Only_Early_up" = Only_Early_up,"Only_Early_down" = Only_Early_down,"Only_Adaptive_up" = Only_Adaptive_up,"Only_Adaptive_down" = Only_Adaptive_down),file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Early_Adaptive_Gene_lists.rds")

# Fig 6a ------------------------------------------------------------------
test <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Early_Adaptive_Gene_lists.rds")
list2env(test,envir = globalenv())

for(i in levels(VHE_late_pos$PT_identity)){
  modules <- list(Only_Early_up[[i]],Only_Early_down[[i]],Only_Adaptive_up[[i]],Only_Adaptive_down[[i]])
  names(modules) <- c("Only_Early_up","Only_Early_down","Only_Adaptive_up","Only_Adaptive_down")
  modules <- modules[sapply(modules,FUN = length) != 0]
  test <- subset(VHE_late_pos,cells = colnames(VHE_late_pos)[VHE_late_pos$PT_identity == i])
  name <- gsub(pattern = " ",replacement = "_",x = i)
  test <- AddModuleScore(test,features = modules,nbin = 100,ctrl = 50,name = "Module_")
  colnames(test@meta.data)[grep(pattern = "^Module_",x = colnames(test@meta.data))] <- names(modules)
  assign(x = name,value = test,envir = globalenv())
}

modulescores <- lapply(list(S1_A,S1_B,S2_A,S2_B,S3_A,S3_B),function(x){x <- FetchData(x,c("PT_identity","Condition","Only_Early_up","Only_Early_down","Only_Adaptive_up","Only_Adaptive_down"))%>%pivot_longer(cols = -c(1:2),names_to = "Modules",values_to = "Score")})
modulescores <- Reduce(rbind,modulescores)
modulescores$Modules <- factor(modulescores$Modules,levels = c("Only_Early_up","Only_Early_down","Only_Adaptive_up","Only_Adaptive_down"),labels = c("Only Early Up","Only Early Down","Only Adaptive Up","Only Adaptive Down"))

saveRDS(modulescores,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6a/VHE_Early_Adaptive_Module_scores.rds")

"````````````````````````````````````````````````````````````````````````````"
"Scale scores so that median in ConKO is 0 and median in VKO is 1. Then plot only VHKO, VEKO and VHEKO"

modulescores_VHE <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6a/VHE_Early_Adaptive_Module_scores.rds")
mod_sumary <- group_by(filter(modulescores_VHE,PT_identity != "Non PT",Modules %in% c("Only Early Up","Only Adaptive Up","Only Adaptive Down")),Condition,Modules)%>% summarise("Median" = quantile(Score,0.5),"Q1" = quantile(Score,0.25),"Q3" = quantile(Score,0.75))

scaled_VHE <- filter(modulescores_VHE,PT_identity != "Non PT",Modules %in% c("Only Early Up","Only Adaptive Up","Only Adaptive Down"))

scalescores <- function(x){
  module <- x[["Modules"]]
  condition <- x[["Condition"]]
  score <- as.numeric(x[["Score"]])
  conkoscore <- as.numeric(mod_sumary[mod_sumary$Condition == "ConKO" & mod_sumary$Modules == module,"Median"])
  vkoscore <- as.numeric(mod_sumary[mod_sumary$Condition == "VKO" & mod_sumary$Modules == module,"Median"])
  scale <- (score - conkoscore)/(vkoscore - conkoscore)
  if(module == "Only Adaptive Down"){scale <- -1*scale}
  return(scale)
}
scaled_VHE[["Scaled_score"]] <- apply(scaled_VHE,1,scalescores)

scaled_VHE$Modules <- factor(scaled_VHE$Modules,levels = c("Only Early Up","Only Adaptive Up","Only Adaptive Down"),labels = c("Early Up","Adaptive Up","Adaptive Down"))

condition_colour_2 <-  c("ConKO"="#595FC7","VKO"="#DA5E24", "VHKO" = "#bc8e54", "VEKO" = "#208d4c","VHEKO" = "#9467cf")
summary_data <- filter(scaled_VHE, Modules == "Early Up") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")

p1 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour_2)+
  labs(title = "Early Up",y = "Scaled expression score",x=NULL)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
summary_data <- filter(scaled_VHE, Modules == "Adaptive Up") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")
p2 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour_2)+
  labs(title = "Adaptive Up",y = "Scaled expression score",x=NULL)+
  scale_fill_manual(values = condition_colour_2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
summary_data <- filter(scaled_VHE, Modules == "Adaptive Down") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")
p3 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour_2)+
  labs(title = "Adaptive Down",y = "Scaled expression score",x=NULL)+
  scale_fill_manual(values = condition_colour_2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
p4 <- plot_grid(p1,p2,p3,nrow = 1)
pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6a/Fig_6a.pdf",width = 18/2.54,height = 2)
p4
dev.off()


# Fig 6b ------------------------------------------------------------------

Vhl_Pax8_pos_batch_sex<- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/ConKO_VKO_ConKO_Early_Late_Data.rds")

hif_modules <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")[c(1,2,6)]

for(i in levels(Vhl_Pax8_pos_batch_sex$PT_identity)){
  modules <-lapply(hif_modules, function(x){x <- x[[i]]})
  modules <- modules[sapply(modules,FUN = length) != 0]
  test <- subset(Vhl_Pax8_pos_batch_sex,cells = colnames(Vhl_Pax8_pos_batch_sex)[Vhl_Pax8_pos_batch_sex$PT_identity == i])
  name <- gsub(pattern = " ",replacement = "_",x = i)
  test <- AddModuleScore(test,features = modules,nbin = 100,ctrl = 50,name = "Module_")
  colnames(test@meta.data)[grep(pattern = "^Module_",x = colnames(test@meta.data))] <- names(modules)
  assign(x = name,value = test,envir = globalenv())
}

modulescores <- lapply(list(S1_A,S1_B,S2_A,S2_B,S3_A,S3_B),function(x){x <- FetchData(x,c("PT_identity","Condition","Need_H_up","Need_E_up","Need_E_down"))%>%pivot_longer(cols = -c(1:2),names_to = "Modules",values_to = "Score")})
modulescores <- Reduce(rbind,modulescores)
modulescores$Modules <- factor(modulescores$Modules,levels = c("Need_H_up","Need_E_up","Need_E_down"),labels = c("HIF1A Up","HIF2A Up","HIF2A Down"))
modulescores$Condition <- factor(modulescores$Condition,levels = c("Control Early Pos", "Control Late Pos", "KO Early Pos", "KO Late Pos"),labels = c("ConKO E","ConKO L","VKO E","VKO L"))

saveRDS(modulescores,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6b/Vhl_Pax8_pos_batch_sex_HIF1A_HIF2A_Module_scores.rds")

"`````````````````````````````````````````````````````````````````````````````"
"Scale so that ConKO E is 0 and VKO L is 1"
modulescores_Vhl <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6b/Vhl_Pax8_pos_batch_sex_HIF1A_HIF2A_Module_scores.rds")

mod_sumary <- group_by(filter(modulescores_Vhl,PT_identity != "Non PT",Modules %in% c("HIF1A Up","HIF2A Up","HIF2A Down")),Condition,Modules)%>% summarise("Median" = quantile(Score,0.5),"Q1" = quantile(Score,0.25),"Q3" = quantile(Score,0.75))

scaled_Vhl <- filter(modulescores_Vhl,PT_identity != "Non PT",Modules %in% c("HIF1A Up","HIF2A Up","HIF2A Down"))

scalescores <- function(x){
  module <- x[["Modules"]]
  condition <- x[["Condition"]]
  score <- as.numeric(x[["Score"]])
  conkoscore <- as.numeric(mod_sumary[mod_sumary$Condition == "ConKO E" & mod_sumary$Modules == module,"Median"])
  vkoscore <- as.numeric(mod_sumary[mod_sumary$Condition == "VKO L" & mod_sumary$Modules == module,"Median"])
  scale <- (score - conkoscore)/(vkoscore - conkoscore)
  if(module == "HIF2A Down"){scale <- -1*scale}
  return(scale)
}
scaled_Vhl[["Scaled_score"]] <- apply(scaled_Vhl,1,scalescores)

scaled_Vhl$Modules <- factor(scaled_Vhl$Modules,levels = c("HIF1A Up","HIF2A Up","HIF2A Down"))

levels(scaled_Vhl$Condition) <- c("ConKO\nEarly","ConKO\nLate","VKO\nEarly","VKO\nLate")

condition_colour <-  c("ConKO\nEarly"="#93A0FF","VKO\nEarly"="#F29163","ConKO\nLate"="#595FC7","VKO\nLate"="#DA5E24")
summary_data <- filter(scaled_Vhl, Modules == "HIF1A Up") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")
p1 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour)+
  labs(title = "HIF1A Up",y = "Scaled expression score",x=NULL)+
  scale_fill_manual(values = condition_colour)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
summary_data <- filter(scaled_Vhl, Modules == "HIF2A Up") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")
p2 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour)+
  labs(title = "HIF2A Up",y = "Scaled expression score",x=NULL)+
  scale_fill_manual(values = condition_colour)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
summary_data <- filter(scaled_Vhl, Modules == "HIF2A Down") %>%
  group_by(Condition) %>%
  summarise(lower = quantile(Scaled_score, 0.05),q1 = quantile(Scaled_score, 0.25),median = median(Scaled_score),q3 = quantile(Scaled_score, 0.75),upper = quantile(Scaled_score, 0.95),.groups = "drop")
p3 <- ggplot(summary_data, aes(x = Condition)) +
  geom_errorbar(aes(x = Condition, ymin = median, ymax = median, color = Condition),width = 0.5, linewidth = 1 ) +
  geom_errorbar(aes(x = Condition, ymin = q1,ymax = q3,y = median,color = Condition),width = 0.2,linewidth = 0.5)+
  scale_color_manual(values = condition_colour)+
  labs(title = "HIF2A Down",y = "Scaled expression score",x=NULL)+
  scale_fill_manual(values = condition_colour)+
  theme_classic()+
  theme(plot.title = element_text(size = 8,hjust = 0.5),axis.title = element_text(size = 7,hjust = 0.5,margin = margin(0,0,0,0,"cm")), axis.text = element_text(size = 6,hjust = 0.5),axis.text.x = element_text(size = 7,vjust = 1,hjust = 0.5,angle = 0,colour = "black"),plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))+NoLegend()
p4 <- plot_grid(p1,p2,p3,nrow = 1)
pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6b/Fig_6b.pdf",width = 18/2.54,height = 2)
p4
dev.off()


# Supplementary Fig 11 ------------------------------------------------------------------

"Accessed differential expression data from Supplementary Information provided in
1. Harlander et al., 2017; Nat Med
2. Hoefflin et al., 2022; Nat Commun
3. Nargund et al., 2017; Cell Rep.
Saved the data as .xlsx files"

Harlander <- readxl::read_xlsx(path = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_11/Harlander_ccrcc_rnaseq.xlsx",sheet = 2)
Hoefflin <- readxl::read_xlsx(path = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_11/Hoefflin_ccrcc_rnaseq.xlsx",sheet = 2)
Nargund <- readxl::read_xlsx(path = "//Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_11/Nargund_ccrcc_rnaseq.xlsx",sheet = 2)

test <- readRDS(file =  "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Hif_isoform_dependent_genes.rds")
list2env(x = test,envir = globalenv())
test <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/Early_Adaptive_Gene_lists.rds")
list2env(x = test,envir = globalenv())

H_early_up <- intersectlist(Need_H_up,Only_Early_up)
H_adaptive_up <- intersectlist(Need_H_up,Only_Adaptive_up)
E_early_up <- intersectlist(Need_E_up,Only_Early_up)
E_adaptive_up <- intersectlist(Need_E_up,Only_Adaptive_up)
E_adaptive_down <- intersectlist(Need_E_down,Only_Adaptive_down)


pathways <- list("H_early_up" = unique(unlist(H_early_up)),"E_early_up" =unique(unlist(E_early_up)),"E_adaptive_up" = unique(unlist(E_adaptive_up)),"E_adaptive_down"= unique(unlist(E_adaptive_down)))

ccRCC_gsea <- list()

"Harlander"
Harlander_rank <- as.numeric(arrange(Harlander, desc(L2FC))$L2FC)
names(Harlander_rank) <- as.character(arrange(Harlander, desc(L2FC))$Gene)
harlanderfgsea <- fgsea::fgsea(pathways = pathways,stats = Harlander_rank)
ccRCC_gsea[["Harlander"]] <- harlanderfgsea
p1 <- fgsea::plotEnrichment(pathway = pathways[[1]],stats = Harlander_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p2 <- fgsea::plotEnrichment(pathway = pathways[[2]],stats = Harlander_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p3 <- fgsea::plotEnrichment(pathway = pathways[[3]],stats = Harlander_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p4 <- fgsea::plotEnrichment(pathway = pathways[[4]],stats = Harlander_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p5 <- plot_grid(p1,p2,p3,p4,nrow = 4)

"Hoefflin"
Hoefflin_rank <- as.numeric(arrange(Hoefflin, desc(L2FC_VPR))$L2FC_VPR)
names(Hoefflin_rank) <- as.character(arrange(Hoefflin, desc(L2FC_VPR))$Gene)
Hoefflinfgsea <- fgsea::fgsea(pathways = pathways,stats = Hoefflin_rank)
ccRCC_gsea[["Hoefflin_ccRCC"]] <- Hoefflinfgsea
p1 <- fgsea::plotEnrichment(pathway = pathways[[1]],stats = Hoefflin_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p2 <- fgsea::plotEnrichment(pathway = pathways[[2]],stats = Hoefflin_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p3 <- fgsea::plotEnrichment(pathway = pathways[[3]],stats = Hoefflin_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p4 <- fgsea::plotEnrichment(pathway = pathways[[4]],stats = Hoefflin_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p6 <- plot_grid(p1,p2,p3,p4,nrow = 4)

"Nargund"
annots <- orthologs(genes = Nargund$EnsemblID,species = "mouse",human = FALSE)
Nargund <- left_join(Nargund,annots[,c("ensembl","symbol")],by = join_by("EnsemblID" == "ensembl"))
Nargund <- Nargund[!is.na(Nargund$L2FC),c("symbol","L2FC")]
Nargund_rank <- as.numeric(arrange(Nargund, desc(L2FC))$L2FC)
names(Nargund_rank) <- as.character(arrange(Nargund, desc(L2FC))$symbol)
Nargundfgsea <- fgsea::fgsea(pathways = pathways,stats = Nargund_rank)
ccRCC_gsea[["Nargund"]] <- Nargundfgsea
p1 <- fgsea::plotEnrichment(pathway = pathways[[1]],stats = Nargund_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p2 <- fgsea::plotEnrichment(pathway = pathways[[2]],stats = Nargund_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p3 <- fgsea::plotEnrichment(pathway = pathways[[3]],stats = Nargund_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p4 <- fgsea::plotEnrichment(pathway = pathways[[4]],stats = Nargund_rank)+labs(x=NULL,y="Enrichement Score")+theme_classic()+theme(title = element_blank(),text = element_text(size = 6,colour = "black"),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),plot.title = element_text(margin = unit(c(0,0,0,0),"cm")))
p7 <- plot_grid(p1,p2,p3,p4,nrow = 4)
p8 <- plot_grid(p5,p6,p7,nrow = 1)

pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_11/Supplementary_Fig_11_notitle.pdf",width = 17/2.54,height = 14/2.54)
p8
dev.off()

ccRCC_gsea_NES_p <- lapply(ccRCC_gsea, function(x){
  x <- x[,c("pathway","padj","NES")]
})
ccRCC_gsea_NES_p <- bind_rows(ccRCC_gsea_NES_p,.id = "Dataset")

colnames(ccRCC_gsea_NES_p) <- c("Dataset","Geneset","p","NES")

ccRCC_gsea_NES_p$Geneset <- factor(ccRCC_gsea_NES_p$Geneset,levels = c("H_early_up","E_early_up","E_adaptive_up","E_adaptive_down"),labels = c("HIF1A Early Up","HIF2A Early Up","HIF2A Adaptive Up","HIF2A Adaptive Down"))
ccRCC_gsea_NES_p$Dataset <- factor(ccRCC_gsea_NES_p$Dataset,levels = c("Harlander", "Hoefflin_ccRCC", "Nargund"),labels = c("Harlander_2017", "Hoefflin_2020", "Nargund_2017"))


write.table(arrange(ccRCC_gsea_NES_p,Geneset,Dataset),file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Supplementary_Fig_11/NES_p_values.csv",sep = ",",row.names = FALSE,col.names = TRUE)




# Scoring TCGA KIRC samples for early and adaptive genes - David Mole ------------------------------------------------------------------

library(survival)
library(survminer)
library(RTCGA.clinical)



# ensemble gene converter
# https://www.biotools.fr/mouse/ensembl_symbol_converter

"`````````````````````````NOTE`````````````````````````````````````````````````"
"David has downloaded TCGA RNA expression data for all cancer types and matched normal samples. From this data, he has then filtered the RNA seq to include only KIRC tumour and normal samples and only HIF1A Early Up, HIF2A Early Up, HIF2A Adaptive Up, and HIF2A Adaptive Down genes. 
For this he converted the gene lists produced in our data to human orthologs. The genes were filtered twice - once to remove antisense transcripts (AST) and then to remove genes with a DT annotation"


head -1 /Volumes/Extreme_Pro/cBioportal_data/kirc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt > /Users/drmole/Documents/Group_results/Samvid_Kurlekar/data_mrna_seq_v2_rsem_zscores_ref_all_samples_HIF1A_early_up.txt
grep -wf  /Users/drmole/Documents/Group_results/Samvid_Kurlekar/HIF1A_early_up.txt  /Volumes/Extreme_Pro/cBioportal_data/kirc_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt >> /Users/drmole/Documents/Group_results/Samvid_Kurlekar/data_mrna_seq_v2_rsem_zscores_ref_all_samples_HIF1A_early_up.txt

"HIF1A Early Up - Tumour"
# grep to extract data for cancer samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF1A_early_up.txt /arc/ratcliffe/data/TCGApanCancer/files/TCGApanCancer_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF1A_early_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF1A_early_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF1A_early_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF1A_early_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF1A_early_up_genes3.tsv

"HIF1A Early Up - Normal"
# grep to extract data for normal samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF1A_early_up.txt /arc/ratcliffe/data/TCGApanNormal/files/TCGApanNormal_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF1A_early_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF1A_early_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF1A_early_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF1A_early_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF1A_early_up_genes3.tsv

"HIF2A Early Up - Tumour"
# grep to extract data for cancer samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_early_up.txt /arc/ratcliffe/data/TCGApanCancer/files/TCGApanCancer_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_early_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_early_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_early_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_early_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_early_up_genes3.tsv

"HIF2A Early Up - Normal"
# grep to extract data for normal samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_early_up.txt /arc/ratcliffe/data/TCGApanNormal/files/TCGApanNormal_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_early_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_early_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_early_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_early_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_early_up_genes3.tsv

"HIF2A Adaptive Up - Tumour"
# grep to extract data for cancer samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_adaptive_up.txt /arc/ratcliffe/data/TCGApanCancer/files/TCGApanCancer_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_up_genes3.tsv

"HIF2A Adaptive Up - Normal"
# grep to extract data for normal samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_adaptive_up.txt /arc/ratcliffe/data/TCGApanNormal/files/TCGApanNormal_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_up_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_up_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_up_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_up_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_up_genes3.tsv

"HIF2A Adaptive Down - Tumour"
# grep to extract data for cancer samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_adaptive_down.txt /arc/ratcliffe/data/TCGApanCancer/files/TCGApanCancer_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_down_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_down_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_down_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_down_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanCancer_HIF2A_adaptive_down_genes3.tsv

"HIF2A Adaptive Down - Normal"
# grep to extract data for normal samples 
grep -w -f /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/HIF2A_adaptive_down.txt /arc/ratcliffe/data/TCGApanNormal/files/TCGApanNormal_combined_names.txt > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_down_genes.tsv
# remove antisense(-AS1) transcripts
grep -v AS1 /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_down_genes.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_down_genes2.tsv
# remove (-DT) transcripts
grep -v DT /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_down_genes2.tsv > /arc/ratcliffe/data/TCGApanCancer/files/adaptive_genes/TCGApanNormal_HIF2A_adaptive_down_genes3.tsv

"``````````````````````````````````````````````````````````````````````````````"
"Reference for sample ID and cancer types"
# import data for the tumour type of each sample
tumour_types <- read.csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/files.2021-06-07-tumour.csv", as.is=T, header = TRUE)
tumour_types <- as_tibble(tumour_types)
tumour_types <- tumour_types %>%
  separate(file_name, into = c("sample", "b"), sep=".FPKM")
tumour_types <- transmute(tumour_types, sample=sample, tumour_type=cases.0.project.project_id)

normal_types <- read.csv("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/files.2021-06-07-normal.csv", as.is=T, header = TRUE)
normal_types <- as_tibble(normal_types)
normal_types <- normal_types %>%
  separate(file_name, into = c("sample", "b"), sep=".FPKM")
normal_types <- transmute(normal_types, sample=sample, tumour_type=cases.0.project.project_id)

"``````````````````````READING TUMOR DATA`````````````````````````````````````"

"HIF1A EARLY UP"

TCGA_tumor_HIF1A_Early_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanCancer_HIF1A_early_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_tumor_HIF1A_Early_Up <- as_tibble(TCGA_tumor_HIF1A_Early_Up)
TCGA_tumor_HIF1A_Early_Up <- TCGA_tumor_HIF1A_Early_Up %>%
  add_column(sample_type = "Tumour")
TCGA_tumor_HIF1A_Early_Up <- transmute(TCGA_tumor_HIF1A_Early_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_tumor_HIF1A_Early_Up <- tumour_types %>%
  inner_join(TCGA_tumor_HIF1A_Early_Up, by="sample")

"HIF2A EARLY UP"

TCGA_tumor_HIF2A_Early_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanCancer_HIF2A_early_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_tumor_HIF2A_Early_Up <- as_tibble(TCGA_tumor_HIF2A_Early_Up)
TCGA_tumor_HIF2A_Early_Up <- TCGA_tumor_HIF2A_Early_Up %>%
  add_column(sample_type = "Tumour")
TCGA_tumor_HIF2A_Early_Up <- transmute(TCGA_tumor_HIF2A_Early_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_tumor_HIF2A_Early_Up <- tumour_types %>%
  inner_join(TCGA_tumor_HIF2A_Early_Up, by="sample")

"HIF2A Adaptive UP"

TCGA_tumor_HIF2A_Adaptive_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanCancer_HIF2A_adaptive_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_tumor_HIF2A_Adaptive_Up <- as_tibble(TCGA_tumor_HIF2A_Adaptive_Up)
TCGA_tumor_HIF2A_Adaptive_Up <- TCGA_tumor_HIF2A_Adaptive_Up %>%
  add_column(sample_type = "Tumour")
TCGA_tumor_HIF2A_Adaptive_Up <- transmute(TCGA_tumor_HIF2A_Adaptive_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_tumor_HIF2A_Adaptive_Up <- tumour_types %>%
  inner_join(TCGA_tumor_HIF2A_Adaptive_Up, by="sample")

"HIF2A Adaptive Down"

TCGA_tumor_HIF2A_Adaptive_Down <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanCancer_HIF2A_adaptive_down_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_tumor_HIF2A_Adaptive_Down <- as_tibble(TCGA_tumor_HIF2A_Adaptive_Down)
TCGA_tumor_HIF2A_Adaptive_Down <- TCGA_tumor_HIF2A_Adaptive_Down %>%
  add_column(sample_type = "Tumour")
TCGA_tumor_HIF2A_Adaptive_Down <- transmute(TCGA_tumor_HIF2A_Adaptive_Down, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_tumor_HIF2A_Adaptive_Down <- tumour_types %>%
  inner_join(TCGA_tumor_HIF2A_Adaptive_Down, by="sample")


"Combining all"

TCGA_tumor <- rbind(TCGA_tumor_HIF1A_Early_Up,TCGA_tumor_HIF2A_Early_Up,TCGA_tumor_HIF2A_Adaptive_Up,TCGA_tumor_HIF2A_Adaptive_Down)

"```````````````````````````````````READING Normal DATA``````````````````````````"

"HIF1A EARLY UP"

TCGA_Normal_HIF1A_Early_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanNormal_HIF1A_early_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_Normal_HIF1A_Early_Up <- as_tibble(TCGA_Normal_HIF1A_Early_Up)
TCGA_Normal_HIF1A_Early_Up <- TCGA_Normal_HIF1A_Early_Up %>%
  add_column(sample_type = "Normal")
TCGA_Normal_HIF1A_Early_Up <- transmute(TCGA_Normal_HIF1A_Early_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_Normal_HIF1A_Early_Up <- normal_types %>%
  inner_join(TCGA_Normal_HIF1A_Early_Up, by="sample")


"HIF2A EARLY UP"

TCGA_Normal_HIF2A_Early_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanNormal_HIF2A_early_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_Normal_HIF2A_Early_Up <- as_tibble(TCGA_Normal_HIF2A_Early_Up)
TCGA_Normal_HIF2A_Early_Up <- TCGA_Normal_HIF2A_Early_Up %>%
  add_column(sample_type = "Normal")
TCGA_Normal_HIF2A_Early_Up <- transmute(TCGA_Normal_HIF2A_Early_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_Normal_HIF2A_Early_Up <- normal_types %>%
  inner_join(TCGA_Normal_HIF2A_Early_Up, by="sample")


"HIF2A Adaptive UP"

TCGA_Normal_HIF2A_Adaptive_Up <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanNormal_HIF2A_adaptive_up_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_Normal_HIF2A_Adaptive_Up <- as_tibble(TCGA_Normal_HIF2A_Adaptive_Up)
TCGA_Normal_HIF2A_Adaptive_Up <- TCGA_Normal_HIF2A_Adaptive_Up %>%
  add_column(sample_type = "Normal")
TCGA_Normal_HIF2A_Adaptive_Up <- transmute(TCGA_Normal_HIF2A_Adaptive_Up, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_Normal_HIF2A_Adaptive_Up <- normal_types %>%
  inner_join(TCGA_Normal_HIF2A_Adaptive_Up, by="sample")


"HIF2A Adaptive Down"

TCGA_Normal_HIF2A_Adaptive_Down <- read.table("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGApanNormal_HIF2A_adaptive_down_genes3.tsv", sep=" ", as.is=T, header=F)
TCGA_Normal_HIF2A_Adaptive_Down <- as_tibble(TCGA_Normal_HIF2A_Adaptive_Down)
TCGA_Normal_HIF2A_Adaptive_Down <- TCGA_Normal_HIF2A_Adaptive_Down %>%
  add_column(sample_type = "Normal")
TCGA_Normal_HIF2A_Adaptive_Down <- transmute(TCGA_Normal_HIF2A_Adaptive_Down, sample=V1, sample_type=sample_type, gene=V2, expression=V3)
TCGA_Normal_HIF2A_Adaptive_Down <- normal_types %>%
  inner_join(TCGA_Normal_HIF2A_Adaptive_Down, by="sample")


"Combining all"

TCGA_Normal <- rbind(TCGA_Normal_HIF1A_Early_Up,TCGA_Normal_HIF2A_Early_Up,TCGA_Normal_HIF2A_Adaptive_Up,TCGA_Normal_HIF2A_Adaptive_Down)

"``````````````````````COMBINING TUMOR AND NORMAL````````````````````````````````"
TCGA <- rbind(TCGA_Normal,TCGA_tumor)%>%as.data.frame()%>%filter(!duplicated(paste(sample,gene)))
TCGA_matrix <- pivot_wider(data = TCGA,names_from = gene,values_from = expression)

"Scaling so that for each gene across all samples, mean is 0 and SD is 1"
TCGA_norm <- scale(TCGA_matrix[,(4:148)])
TCGA_matrix <- cbind(TCGA_matrix[,(1:3)], TCGA_norm[,(1:145)])

"Score for a geneset is the sum of the scaled gene expression for all genes in the geneset"

TCGA_sums <- TCGA_matrix[,(1:3)]
TCGA_sums[["HIF1A Early Up"]] <- rowSums(TCGA_norm[,unique(TCGA_tumor_HIF1A_Early_Up$gene)], na.rm = TRUE)
TCGA_sums[["HIF2A Early Up"]] <- rowSums(TCGA_norm[,unique(TCGA_tumor_HIF2A_Early_Up$gene)], na.rm = TRUE)
TCGA_sums[["HIF2A Adaptive Up"]] <- rowSums(TCGA_norm[,unique(TCGA_tumor_HIF2A_Adaptive_Up$gene)], na.rm = TRUE)
TCGA_sums[["HIF2A Adaptive Down"]] <- rowSums(TCGA_norm[,unique(TCGA_tumor_HIF2A_Adaptive_Down$gene)], na.rm = TRUE)

TCGA_sums <- pivot_longer(TCGA_sums,cols = 4:7,names_to = "Metagene",values_to = "Score")
colnames(TCGA_sums) <- c("ID", "Cancer", "Sample", "Metagene")

saveRDS(TCGA_sums,file = "/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGA_HIF_Module_Scores.rds")

# Fig 6c ------------------------------------------------------------------

"```````````````````````FILTERING CCRCC SAMPLES`````````````````````````````````"

TCGA_sums <- readRDS("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/scRNA_seq_Data_Analysis/TCGA_HIF_Module_Scores.rds")
KIRC <- filter(TCGA_sums,tumour_type == "TCGA-KIRC")[,-2]
colnames(KIRC) <- c("ID","Type","Metagene","Score")
n_numbers <- group_by(KIRC,Type,Metagene)%>%summarise("Count" = n())
KIRC$Type <- factor(KIRC$Type,levels = c("Normal","Tumour"),labels = c("Normal\n(n = 72)","ccRCC\n(n = 538)"))
KIRC$Metagene <- factor(KIRC$Metagene,levels = c("HIF1A Early Up", "HIF2A Early Up", "HIF2A Adaptive Up", "HIF2A Adaptive Down"))

KIRC_summary <- group_by(KIRC,Type,Metagene)%>%
  summarise("Median" = quantile(Score,0.5),"Q1" = quantile(Score,0.25),"Q3" = quantile(Score,0.75))%>%ungroup()

KIRC_kruskal <- KIRC %>%
  group_by(Metagene) %>%
  kruskal_test(Score ~ Type)%>%
  adjust_pvalue(method = "bonferroni")

pvals <-KIRC_kruskal%>%
  mutate(group1 = "Normal\n(n = 72)", group2 = "ccRCC\n(n = 538)",y.position = case_when(Metagene == "HIF1A Early Up" ~ 25,Metagene == "HIF2A Early Up" ~ 17, Metagene == "HIF2A Adaptive Up" ~ 18, Metagene == "HIF2A Adaptive Down" ~ 11),padj = ifelse(p.adj < 0.001,yes = "p < 0.001",no = paste("p = ",p.adj,sep = "")))

p1 <- ggplot(KIRC_summary,aes(x = Type,color = Type))+
  geom_errorbar(aes(ymin = Median, ymax = Median),width = 0.5, linewidth = 1) +
  geom_errorbar(aes(ymin = Q1,ymax = Q3),width = 0.2,linewidth = 0.5)+
  facet_wrap(facets = vars(Metagene),scales = "free_y",nrow = 1)+
  labs(x=NULL,y = "Expression Score",title = NULL)+
  scale_color_manual(values = c("Normal\n(n = 72)"="#595FC7","ccRCC\n(n = 538)"="#DA5E24"))+
  stat_pvalue_manual(pvals,label = "padj", x = "group1",y.position = "y.position",size = 2.5)+
  theme_classic()+
  theme(strip.text = element_text(size = 11,hjust = 0.5,color = "black"),title = element_text(size = 8),axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 7),strip.background = element_blank())+
  NoLegend()
pdf("/Users/samvid_k/Library/CloudStorage/OneDrive-Nexus365/DPhil_Clinical_Medicine/Lima_et_al_2025/Raw_Data/Fig_6/Fig_6c/Fig_6c.pdf",width = 18/2.54,height = 7.5/2.54)
p1
dev.off()
