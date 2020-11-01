###### Loading necessary packages   ########
library(SingleR)
library(dplyr)
library(Seurat)
library(reshape)
library(data.table)
library(readr)
library(hash)
library(ggplot2)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)
library(patchwork)
library(gtable)
library(plotly)

samples <- hash()

####################### Definig paths and variables ################################
samples[["sample1"]] <- "/path/to/sample1"
samples[["sample2"]] <- "/path/to/sample2"


merged_project <- "project_name"
R_object_folder <- "/path/to/output/dir"
cellcounts_csv <- "/path/to/output/dir/all_cellcounts.csv"
mergedcounts_csv <- "/path/to/output/dir/all_mergedcounts.csv"

myColors = c("#E31A1C", "#CAB2D6",  "#6A3D9A", "#FFFF00", "#1F78B4", "#A6CEE3", 
             "#33A02C", "#FB9A99", "#B2DF8A", "#FF7F00", "#FDBF6F", "#B15928")
myColors1 = c("#E31A1C", "#CAB2D6",  "#6A3D9A", "#1F78B4", "#A6CEE3", "#33A02C", 
              "#FB9A99", "#B2DF8A", "#FF7F00", "#FDBF6F", "#B15928")

###### the file that contains the rules for cell type merging (also found in GitHub repo)
ct <- read_delim("/path/to/cell_type_merging.txt", "\t", escape_double = FALSE, trim_ws = TRUE)



### Load reference set, create and merge SingleR objects ###
BP_ENCODE <- BlueprintEncodeData()
for (s in ls(samples)) {
  X10 <- Read10X(data.dir = samples[[s]])
  seurat <- CreateSeuratObject(counts = X10, project = s, min.cells = 2, min.features = 100)
  seurat$orig.ident <- s
  singleR <- SingleR(test = X10, ref = BP_ENCODE, labels = BP_ENCODE$label.fine)
  
  celltypes <- data.frame("barcode" = rownames(singleR), "cell_type" = singleR$labels)
  rownames(celltypes) <- celltypes$barcode
  seurat@meta.data <- merge(x = seurat@meta.data, y = celltypes, by = 0, all.x = TRUE)
  rownames(seurat@meta.data) <- seurat@meta.data$barcode
  seurat@meta.data <- subset(seurat@meta.data, select=-c(Row.names))
  assign(s, seurat)
}

################################# All SAMPLES ####################################
seurat_merged <- merge(eval(parse(text = keys(samples)[1])), 
                       y = c(eval(parse(text = keys(samples)[2]))), # add as many samples as you want to load
                       add.cell.ids = c(keys(samples)),
                       project = merged_project)

seurat_merged@meta.data$barcode <- rownames(seurat_merged@meta.data)

seurat_merged@meta.data <- merge(seurat_merged@meta.data, ct, by.x = "cell_type", by.y = "cell_type", all.x=TRUE)
rownames(seurat_merged@meta.data) <- seurat_merged@meta.data$barcode

#################################################################################

############################ Labelling samples ##################################

seurat_merged@meta.data$days <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$days[seurat_merged@meta.data$days == "sample1"] <- 0
seurat_merged@meta.data$days[seurat_merged@meta.data$days == "sample2"] <- 1

seurat_merged@meta.data$patient <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$patient[seurat_merged@meta.data$patient == "sample1"] <- "P1"
seurat_merged@meta.data$patient[seurat_merged@meta.data$patient == "sample2"] <- "P2"

seurat_merged@meta.data$treatment <- seurat_merged@meta.data$orig.ident
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "sample1"] <- "injected"
seurat_merged@meta.data$treatment[seurat_merged@meta.data$treatment == "sample2"] <- "non-injected"

seurat_merged@meta.data$disease <- seurat_merged@meta.data$patient
seurat_merged@meta.data$disease[seurat_merged@meta.data$patient == "P1"] <- "P1 - pCFCL"
seurat_merged@meta.data$disease[seurat_merged@meta.data$patient == "P2"] <- "P2 - pCDLBCL-LT"

###################### Creating cell-type output files ############################
x <- as.data.frame(table(seurat_merged$orig.ident, seurat_merged$cell_type))
y <- as.data.frame(table(seurat_merged$orig.ident, seurat_merged$merged_type))
write.csv(x, cellcounts_csv)
write.csv(y, mergedcounts_csv)


###################### Loading cell cycle markers ############################
cc.genes <- readLines("/home/ubuntu/scRNA/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

###################### Calculating mitochondrial fraction ############################
seurat_merged@meta.data[order(seurat_merged@meta.data$barcode),]
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")

seurat_merged <- subset(seurat_merged, subset = nFeature_RNA > 99, )
seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(seurat_merged)
seurat_merged <- ScaleData(seurat_merged, features = all.genes)
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged))
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.5)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:30)
seurat_merged <- CellCycleScoring(seurat_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

############################## UMAP of clusters #########################################
DimPlot(seurat_merged, reduction = "umap", label = T, label.size = 12)

########################### Identifying pDCs ##################################### 
FeaturePlot(seurat_merged, c("CLEC4C", "LILRA4"))
seurat_merged@meta.data$merged_type[seurat_merged@meta.data$seurat_clusters == 23] <- "pDC" # or whichever cluster expresses those markers the most

################### Identifying the malignant populations ############################
markers <- c("MS4A1", "LY75", "SIGLEC1", "CD14", "KLRF1", "CLEC4C", "CD38", "CD4", "CD8A", "FOXP3")
c <- VlnPlot(seurat_merged, markers, cols = myColors1, pt.size = 0, ncol = 5)
seurat_merged@meta.data$merged_type[seurat_merged@meta.data$seurat_clusters == 12] <- "Malignant" # change to appropriate cluster
seurat_merged@meta.data$merged_type[seurat_merged@meta.data$seurat_clusters == 14] <- "Malignant" # change to appropriate cluster
seurat_merged@meta.data$merged_type[seurat_merged@meta.data$seurat_clusters == 16] <- "Malignant" # change to appropriate cluster

########################### Adding VDJ sequencing data #########################

sample1_BCR <- "/path/to/sample1_BCR/filtered_contig_annotations.csv"
sample2_BCR <- "/path/to/sample2_BCR/filtered_contig_annotations.csv"

csv_sample1_BCR <- read.csv(sample1_BCR, header=TRUE, sep=",")
csv_sample2_BCR <- read.csv(sample2_BCR, header=TRUE, sep=",")


sample1_TCR <- "/path/to/sample1_TCR/filtered_contig_annotations.csv"
sample2_TCR <- "/path/to/sample2_TCR/filtered_contig_annotations.csv"

csv_sample1_TCR <- read.csv(sample1_TCR, header=TRUE, sep=",")
csv_sample2_TCR <- read.csv(sample2_TCR, header=TRUE, sep=",")


############# Define immune receptor functions ##############

igh <- function(x, csv){
  bc <- x["barcode"]
  c(nrow(csv[csv$barcode == bc & csv$chain == "IGH", ]))
}

igk <- function(x, csv){
  bc <- x["barcode"]
  c(nrow(csv[csv$barcode == bc & csv$chain == "IGK", ]))
}

igl <- function(x, csv){
  bc <- x["barcode"]
  c(nrow(csv[csv$barcode == bc & csv$chain == "IGL", ]))
}

light <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & (csv$chain == "IGL" | csv$chain == "IGK"), ]$v_gene,
        csv[csv$barcode == bc & (csv$chain == "IGL" | csv$chain == "IGK"), ]$j_gene,
        csv[csv$barcode == bc & (csv$chain == "IGL" | csv$chain == "IGK"), ]$c_gene)
}

heavy <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "IGH",]$v_gene,
        csv[csv$barcode == bc & csv$chain == "IGH",]$d_gene,
        csv[csv$barcode == bc & csv$chain == "IGH",]$j_gene,
        csv[csv$barcode == bc & csv$chain == "IGH",]$c_gene)
}

lcdr3 <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & (csv$chain == "IGL" | csv$chain == "IGK"), ]$cdr3)
}

hcdr3 <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "IGH",]$cdr3)
}

bcr <- function(csv){
  csv <- subset(csv, csv$raw_consensus_id != "None")
  current <- data.frame(unique(csv$barcode))
  names(current) <- "barcode"
  current <- distinct(merge(current, select(csv, c("barcode", "raw_clonotype_id")),
                            by.x = "barcode", by.y = "barcode"))
  current$IGH <- apply(current, 1, igh, csv=csv)
  current$IGK <- apply(current, 1, igk, csv=csv)
  current$IGL <- apply(current, 1, igl, csv=csv)
  current$Light_chain <- apply(current, 1, light, csv=csv)
  current$Heavy_chain <- apply(current, 1, heavy, csv=csv)
  current$Lcdr3 <- apply(current, 1, lcdr3, csv=csv)
  current$Hcdr3 <- apply(current, 1, hcdr3, csv=csv)
  current <- merge(current, as.data.frame(table(current$raw_clonotype_id)), by.x = "raw_clonotype_id", by.y = "Var1")
  names(current)[names(current) == "Freq"] <- "BCR_clonotype_freq"
  names(current)[names(current) == "raw_clonotype_id"] <- "BCR_clonotype_id"
  return(current)
}

tra <- function(x, csv){
  bc <- x["barcode"]
  c(nrow(csv[csv$barcode == bc & csv$chain == "TRA", ]))
}

trb <- function(x, csv){
  bc <- x["barcode"]
  c(nrow(csv[csv$barcode == bc & csv$chain == "TRB", ]))
}

alpha_chain <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "TRA",]$v_gene,
        csv[csv$barcode == bc & csv$chain == "TRA",]$j_gene,
        csv[csv$barcode == bc & csv$chain == "TRA",]$c_gene)
}


beta_chain <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "TRB",]$v_gene,
        csv[csv$barcode == bc & csv$chain == "TRB",]$d_gene,
        csv[csv$barcode == bc & csv$chain == "TRB",]$j_gene,
        csv[csv$barcode == bc & csv$chain == "TRB",]$c_gene)
}

acdr3 <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "TRA",]$cdr3)
}

bcdr3 <- function(x, csv){
  bc <- x["barcode"]
  paste(csv[csv$barcode == bc & csv$chain == "TRB",]$cdr3)
}

tcr <- function(csv){
  csv <- subset(csv, csv$raw_consensus_id != "None")
  current <- data.frame(unique(csv$barcode))
  names(current) <- "barcode"
  current <- distinct(merge(current, select(csv, c("barcode", "raw_clonotype_id")),
                            by.x = "barcode", by.y = "barcode"))
  current$TRA <- apply(current, 1, tra, csv=csv)
  current$TRB <- apply(current, 1, trb, csv=csv)
  current$Alpha_chain <- apply(current, 1, alpha_chain, csv=csv)
  current$Beta_chain <- apply(current, 1, beta_chain, csv=csv)
  current$Acdr3 <- apply(current, 1, acdr3, csv=csv)
  current$Bcdr3 <- apply(current, 1, bcdr3, csv=csv)
  current <- merge(current, as.data.frame(table(current$raw_clonotype_id)), by.x = "raw_clonotype_id", by.y = "Var1")
  names(current)[names(current) == "Freq"] <- "TCR_clonotype_freq"
  names(current)[names(current) == "raw_clonotype_id"] <- "TCR_clonotype_id"
  return(current)
}

##############################################################################################################
####################################                BCRs               #######################################
##############################################################################################################

### RUN FOR sample1 ###
per_cell_csv_sample1_BCR <- bcr(csv_sample1_BCR)
per_cell_csv_sample1_BCR$barcode <- paste0("sample1_" ,per_cell_csv_sample1_BCR$barcode)


### RUN FOR sample2 ###
per_cell_csv_sample2_BCR <- bcr(csv_sample2_BCR)
per_cell_csv_sample2_BCR$barcode <- paste0("sample2_" ,per_cell_csv_sample2_BCR$barcode)


### Combine and merge ###
all_bcr <- rbind(per_cell_csv_sample1_BCR, per_cell_csv_sample2_BCR)
rownames(all_bcr) <- all_bcr$barcode
all_bcr <- subset(all_bcr, select=-barcode)
#rownames(seurat_merged@meta.data) <- seurat_merged$Row.names
seurat_merged@meta.data <- merge(x = seurat_merged@meta.data, y = all_bcr, by = 0, all.x = TRUE)
rownames(seurat_merged@meta.data) <- seurat_merged$Row.names
seurat_merged@meta.data <- subset(seurat_merged@meta.data, select=-Row.names)

##############################################################################################################
####################################                TCRs               #######################################
##############################################################################################################

### RUN FOR sample1 ###
per_cell_csv_sample1_TCR <- tcr(csv_sample1_TCR)
per_cell_csv_sample1_TCR$barcode <- paste0("sample1_" ,per_cell_csv_sample1_TCR$barcode)

### RUN FOR sample2 ###
per_cell_csv_sample2_TCR <- tcr(csv_sample2_TCR)
per_cell_csv_sample2_TCR$barcode <- paste0("sample2_" ,per_cell_csv_sample2_TCR$barcode)

### Combine and merge ###
all_tcr <- rbind(per_cell_csv_sample1_TCR, per_cell_csv_sample2_TCR)
rownames(all_tcr) <- all_tcr$barcode
all_tcr <- subset(all_tcr, select=-barcode)
#rownames(seurat_merged@meta.data) <- seurat_merged$Row.names
seurat_merged@meta.data <- merge(x = seurat_merged@meta.data, y = all_tcr, by = 0, all.x = TRUE)
rownames(seurat_merged@meta.data) <- seurat_merged$Row.names
seurat_merged@meta.data <- subset(seurat_merged@meta.data, select=-Row.names)


#################### Identifying potentially malignant clonotypes ####################

seurat_merged@meta.data$aplusb <- paste(seurat_merged@meta.data$Alpha_chain, seurat_merged@meta.data$Beta_chain, sep=", ")
seurat_merged@meta.data$lplush <- paste(seurat_merged@meta.data$Light_chain, seurat_merged@meta.data$Heavy_chain, sep=", ")
sort(table(seurat_merged@meta.data$lplush[seurat_merged@meta.data$BCR_clonotype_freq > 3 & seurat_merged@meta.data$patient == "P1"]),decreasing=T)
sort(table(seurat_merged@meta.data$lplush[seurat_merged@meta.data$BCR_clonotype_freq > 3 & seurat_merged@meta.data$patient == "P2"]),decreasing=T)

################ Picking malignant clonotypes based on previous statistics ##############
seurat_merged@meta.data$malignancy[seurat_merged@meta.data$patient == "P1" & grepl("IGKV1-6 IGKJ1", seurat_merged@meta.data$lplush, fixed = TRUE)] <- "P15 - MZL"
seurat_merged@meta.data$malignancy[seurat_merged@meta.data$patient == "P1" & grepl("IGHV1-2 None IGHJ4", seurat_merged@meta.data$lplush, fixed = TRUE)] <- "P15 - MZL"

seurat_merged@meta.data$malignancy[seurat_merged@meta.data$patient == "P2" & grepl("IGKV3-11 IGKJ4", seurat_merged@meta.data$lplush, fixed = TRUE)] <- "P12 - FCL"
seurat_merged@meta.data$malignancy[seurat_merged@meta.data$patient == "P2" & grepl("IGHV1-46 IGHD2-21 IGHJ5", seurat_merged@meta.data$lplush, fixed = TRUE)] <- "P12 - FCL"


############ We check how the clonotype and gene expression data aligns for the malignant population ##############
DimPlot(seurat_merged, reduction = "umap", cols = myColors, group.by = "merged_type")
DimPlot(seurat_merged, reduction = "umap", group.by = "malignancy")

############### Creating barplot of cell_type data #################
reduced_data <- seurat_merged@meta.data[seurat_merged@meta.data$patient %in% c("P1", "P2"),]
c <- ggplot() + 
  geom_bar(data = reduced_data, aes(x = days, fill = first_type),
           position = "fill") +  
  theme_classic() + 
  scale_fill_manual(name = "cell_type" ,values = myColors1) + 
  scale_y_continuous(labels=percent) +
  xlab("Days p.i.") + 
  ylab("Percentage of cells") +
  labs(title = "Cell type composition changes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(.~patient, scales = "free_x") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
c