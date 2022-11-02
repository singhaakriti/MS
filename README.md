# MS

#### This README.md contains full code for all analyses performed for Aakriti Singh's MS thesis. 

### Versions:
FASTQ files downloaded from Globus Connect Personal (version 3.1.6)
FASTQ files processed through 10X Genomics Cell Ranger Cloud (version 7.0.0) : Single Cell Multiome ATAC and 3' Gene Expression
filtered.feature.bc.matrix files utilized and analyzed:
  - through "hdf5r" package (version 1.3.5)
  - Seurat (version 4.1.1)
  - R (version 4.2.1)
  - ScType (2022)
  - Azimuth (version 0.4.6) using packages "SeuratData" (0.0.2) and "lungrefSeuratData" (2.0.0)
  - ggplot (version 3.36)
  - patchwork (version 1.1.2)
 atac_fragments.tsv.gz files utilized through "Signac" (version 1.7.0)
 Top differentiall expressed markers manually analyzed using Human Protein Atlas (version 21.1)
 
 
 ## Four supervised clusters: analysis with BG18, BG19, BG20, BG22, BG26
 
insert libraries:
```{r}
library(dplyr)
library(Seurat)
library(Matrix)
library(patchwork)
```

load in data sets and create seurat objects
```{r}
NBG18k <- read.csv("/Users/asingh/Downloads/BG18_normalized_excelOpened.csv", header = TRUE, row.names = 1)  
NBG18k_seurat <-CreateSeuratObject(counts = NBG18k, project = "NBG18k")  

NBG19k <- read.csv("/Users/asingh/PycharmProjects/research/practiceResearch/12.21/BG19_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG19k_seurat <- CreateSeuratObject(counts = NBG19k, project = "NBG19k")

NBG20k <- read.csv("/Users/asingh/Downloads/BG20_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG20k_seurat <- CreateSeuratObject(counts = NBG20k, project = "NBG20k")

NBG22k <- read.csv("/Users/asingh/Downloads/BG22_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG22k_seurat <- CreateSeuratObject(counts = NBG22k, project = "NBG22k")

NBG26k <- read.csv("/Users/asingh/Desktop/jan 22/BG26_normalized_excelOpened (1).csv", header = TRUE, row.names = 1)
NBG26k_seurat <- CreateSeuratObject(counts = NBG26k, project = "NBG26k")
```

remove original files to save memory:
```{r}
rm(NBG19k)
rm(NBG20k)
rm(NBG18k)
rm(NBG22k)
rm(NBG26k)
```

merge the seurat objects:
```{r}
mergedNBG <- merge(NBG18k_seurat, y = c(NBG19k_seurat, NBG20k_seurat, NBG22k_seurat, NBG26k_seurat), add.cell.ids = c("NBG_18", "NBG_19", "NBG_20", "NBG_22", "NBG_26"), project = "NBG_CF")
```

remove original seurat objects:
```{r}
rm(NBG18k_seurat)
rm(NBG19k_seurat)
rm(NBG20k_seurat)
rm(NBG22k_seurat)
rm(NBG26k_seurat)
```

```{r}
mergedNBG[["percent.mt"]] <- PercentageFeatureSet(mergedNBG, pattern = "^MT-")
GROVlnPlot(mergedNBG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
scale the data:
```{r}
all.genes <- rownames(mergedNBG)
mergedNBG <- ScaleData(mergedNBG, features = all.genes)
# old: mergedNBG <- ScaleData(object = mergedNBG)
```

Normalizing (must run before running PCA):
```{r}
mergedNBG <- FindVariableFeatures(object = mergedNBG, selection.method = "vst", nfeatures = 2000)
```


run PCA (must perform this before clustering)
```{r}
mergedNBG <- RunPCA(mergedNBG, features = VariableFeatures(object = mergedNBG))
# old: mergedNBG <- RunPCA(mergedNBG, features = mergedNBG@assays$RNA@var.features,  ndims.print = 1:5, nfeatures.print = 5)
```
^^^CAN do multiple PCA analysis on this i.e. dim plot, heat map .. 


cluster the cells:
```{r}
mergedNBG <- FindNeighbors(mergedNBG, dims = 1:10)
mergedNBG <- FindClusters(mergedNBG, resolution = 0.15)
# change resolution -- change number of clusters (want less clusters)
# 0.5 = 12
# ADJUST RESOLUTION AND LOOK AT UMAP CLUSTER -- we like the 5 cluster groups umap 
```

run UMAP:
```{r}
#reticulate::py_install(packages ='umap-learn')
mergedNBG <- RunUMAP(mergedNBG, dims = 1:10)
DimPlot(mergedNBG, reduction = "umap")
#pip install umap-learn
#below is full UMAP with all clusters
```
create a new seurat object selecting the clusters you want to keep (0, 1, 3, and 4)
remove 2, 4, 5 and 7 (want to keep the basal cells which are the clusters on the right)

```{r}
#mergedNBG.basal <- subset(mergedNBG, idents = 0) # just basal cells 
#mergedNBG.basal

mergedNBG.basal <- subset(mergedNBG, idents = c(0, 1, 3, 4))
mergedNBG.basal <- RunUMAP(mergedNBG.basal, dims = 1:10)
DimPlot(mergedNBG.basal, reduction = "umap")
```
https://satijalab.org/seurat/archive/v3.0/interaction_vignette.html
now i have just the basal cells
rerun clusters 
now it is only looking at differences between basal cells (before was looking between cell types)
look at cf vs nonCF
```{r}

```

looking at where specific cell markers are as a double check; we know that BPIFA1 is in ciliated cells 
```{r}
FeaturePlot(object = mergedNBG, features = c("BPIFA1"))
```

TOP2A is in fibroblasts.
```{r}
FeaturePlot(object = mergedNBG, features = c("TOP2A"))
```

KRT5 is a basal cell indicator 
```{r}
FeaturePlot(object = mergedNBG, features = c("KRT5"))
```


cluster 0 
```{r}
cluster0.markers <- FindMarkers(mergedNBG, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)
```
UMAP of just cluster 0
```{r}
mergedNBG.basal <- subset(mergedNBG, idents = c(0))
mergedNBG.basal <- RunUMAP(mergedNBG.basal, dims = 1:10)
DimPlot(mergedNBG.basal, reduction = "umap")
```


cluster1
```{r}
cluster1.markers <- FindMarkers(mergedNBG, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)
```

cluster 3
```{r}
cluster3.markers <- FindMarkers(mergedNBG, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)
```

cluster 4
```{r}
cluster4.markers <- FindMarkers(mergedNBG, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
```

change cluster identities
```{r}
new.cluster.ids <- c("Ciliated/Club", "Basal 1", "a", "Basal 2", "Club", "b", "c", "d")
names(new.cluster.ids) <- levels(mergedNBG)
mergedNBG <- RenameIdents(mergedNBG, new.cluster.ids)

```



named UMAP
```{r}
DimPlot(mergedNBG, reduction = "umap") 
```
labeled on map 
```{r}
DimPlot(mergedNBG, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```
labeling only the subtype:
```{r}
new.cluster.ids <- c("Ciliated/Club", "Basal 1", "Basal 2", "Club")
names(new.cluster.ids) <- levels(mergedNBG.basal)
mergedNBG.basal <- RenameIdents(mergedNBG.basal, new.cluster.ids)
```

```{r}
DimPlot(mergedNBG.basal, reduction = "umap", label = TRUE, pt.size = 0.5) 
```



```{r}
VlnPlot(object = mergedNBG.basal, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"))
```

```{r}
VlnPlot(object = mergedNBG.basal, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"))
```

```{r}
FeaturePlot(object = mergedNBG.basal, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"))

```

CLUSTER VIOLIN/UMAP PLOTS 
0
```{r}
VlnPlot(object = mergedNBG.basal, features = c("TSPAN1", "KRT17", "F3", "MTRNR2L8", "TXNIP", "MTRNR2L2", "NUCB2", "ALDH3A1", "ANXA1", "ATPIF1"))
```
```{r}
FeaturePlot(object = mergedNBG.basal, features = c("TSPAN1", "KRT17", "F3", "MTRNR2L8", "TXNIP", "MTRNR2L2", "NUCB2", "ALDH3A1", "ANXA1", "ATPIF1"))
```


1
```{r}
VlnPlot(object = mergedNBG.basal, features = c("S100A9", "S100A8", "TXN", "MKI67", "TOP2A", "KRT6A", "STMN1", "LY6D", "S100A2", "KRT19"))
```

```{r}
FeaturePlot(object = mergedNBG.basal, features = c("S100A9", "S100A8", "TXN", "MKI67", "TOP2A", "KRT6A", "STMN1", "LY6D", "S100A2", "KRT19"))

```

3
```{r}
VlnPlot(object = mergedNBG.basal, features = c("KRT15", "MMP10", "RPL13A", "S100A2", "RPL27A", "KRT5", "EEF1A1", "MALAT1", "RPS28", "ALDH3A1"))

```
```{r}
FeaturePlot(object = mergedNBG.basal, features = c("KRT15", "MMP10", "RPL13A", "S100A2", "RPL27A", "KRT5", "EEF1A1", "MALAT1", "RPS28", "ALDH3A1"))
```


4
```{r}
VlnPlot(object = mergedNBG.basal, features = c("GABRP", "SERPINB3", "NTS", "IGFBP3", "HES1", "AQP3", "KRT7", "SERPINB4", "HSPB1", "SLPI"))

```

```{r}
FeaturePlot(object = mergedNBG.basal, features = c("GABRP", "SERPINB3", "NTS", "IGFBP3", "HES1", "AQP3", "KRT7", "SERPINB4", "HSPB1", "SLPI"))
```


FINDING DISTINGUISHING MARKERS:

```{r}
cluster0.markers <- FindMarkers(mergedNBG, ident.1 = 0, ident.2 = c(1, 3, 4), min.pct = 0.25)
head(cluster0.markers, n = 10)
```

checks:
```{r}
FeaturePlot(object = mergedNBG, features = c("TSPAN1"))
```

trying to seperate clusters into 5 from 4:

```{r}
mergedNBG <- FindNeighbors(mergedNBG, dims = 1:10)
mergedNBG <- FindClusters(mergedNBG, resolution = 0.25)
```


```{r}
mergedNBG <- RunUMAP(mergedNBG, dims = 1:10)
DimPlot(mergedNBG, reduction = "umap")
```
new UMAP with the 5:
```{r}
mergedNBG.basal <- subset(mergedNBG, idents = c(0, 2, 3, 4, 5))
mergedNBG.basal <- RunUMAP(mergedNBG.basal, dims = 1:10)
DimPlot(mergedNBG.basal, reduction = "umap")
```

new clusters with 5:

cluster 0
```{r}
clusterb0.markers <- FindMarkers(mergedNBG, ident.1 = 0, min.pct = 0.25)
head(clusterb0.markers, n = 10)
```

cluster 2
```{r}
clusterb2.markers <- FindMarkers(mergedNBG, ident.1 = 2, min.pct = 0.25)
head(clusterb2.markers, n = 10)
```

cluster 3
```{r}
clusterb3.markers <- FindMarkers(mergedNBG, ident.1 = 3, min.pct = 0.25)
head(clusterb3.markers, n = 10)
```


cluster 4: CLUB
```{r}
clusterb4.markers <- FindMarkers(mergedNBG, ident.1 = 4, min.pct = 0.25)
head(clusterb4.markers, n = 10)
```


cluster 5
```{r}
clusterb5.markers <- FindMarkers(mergedNBG, ident.1 = 5, min.pct = 0.25)
head(clusterb5.markers, n = 10)
```

### Six supervised clusters: analysis with BG18, BG19, BG20, BG22, BG26
insert libraries
```{r}
library(dplyr)
library(Seurat)
library(Matrix)
library(patchwork)
```

import data sets and create seurat objects
```{r}
NBG18k <- read.csv("/Users/asingh/Downloads/BG18_normalized_excelOpened.csv", header = TRUE, row.names = 1)  
NBG18k_seurat <-CreateSeuratObject(counts = NBG18k, project = "NBG18k")  

NBG19k <- read.csv("/Users/asingh/PycharmProjects/research/practiceResearch/12.21/BG19_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG19k_seurat <- CreateSeuratObject(counts = NBG19k, project = "NBG19k")

NBG20k <- read.csv("/Users/asingh/Downloads/BG20_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG20k_seurat <- CreateSeuratObject(counts = NBG20k, project = "NBG20k")

NBG22k <- read.csv("/Users/asingh/Downloads/BG22_normalized_excelOpened(2).csv", header = TRUE, row.names = 1)
NBG22k_seurat <- CreateSeuratObject(counts = NBG22k, project = "NBG22k")

#NBG26k <- read.csv("/Users/asingh/Desktop/jan 22/BG26_normalized_excelOpened (1).csv", header = TRUE, row.names = 1)
#NBG26k_seurat <- CreateSeuratObject(counts = NBG26k, project = "NBG26k")
```

```{r}
NBG26k <- read.csv("/Users/asingh/Desktop/June22/jan 22/BG26_normalized_excelOpened (1).csv", header = TRUE, row.names = 1)
NBG26k_seurat <- CreateSeuratObject(counts = NBG26k, project = "NBG26k")
```


remove og files 
```{r}
rm(NBG19k)
rm(NBG20k)
rm(NBG18k)
rm(NBG22k)
rm(NBG26k)
```

merge seurat objects 
```{r}
mergedNBG <- merge(NBG18k_seurat, y = c(NBG19k_seurat, NBG20k_seurat, NBG22k_seurat, NBG26k_seurat), add.cell.ids = c("NBG_18", "NBG_19", "NBG_20", "NBG_22", "NBG_26"), project = "NBG_CF")

```

remove og seurat objects 
```{r}
rm(NBG18k_seurat)
rm(NBG19k_seurat)
rm(NBG20k_seurat)
rm(NBG22k_seurat)
rm(NBG26k_seurat)
```

```{r}
NBG26k_seurat[[]]
```


```{r}
mergedNBG[["percent.mt"]] <- PercentageFeatureSet(mergedNBG, pattern = "^MT-")

```

```{r}
VlnPlot(mergedNBG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "orig.ident")
```
```{r}
grep ("^MT-", rownames(mergedNBG[["RNA"]]),value = T)
```

```{r}
subsetmt <- subset(x = mergedNBG, subset = MT-CYB > 0)
subsetmt[[]]
```

scale data
```{r}
all.genes <- rownames(mergedNBG)
mergedNBG <- ScaleData(mergedNBG, features = all.genes)
```

normalize
```{r}
mergedNBG <- FindVariableFeatures(object = mergedNBG, selection.method = "vst", nfeatures = 2000)

```

run PCA
```{r}
mergedNBG <- RunPCA(mergedNBG, features = VariableFeatures(object = mergedNBG))
```

```{r}
VlnPlot(object = mergedNBG, features = c("KRT5", "KRT14", "NGFR", "ITGA6", "P63", "FOXJ1", "SCGB1A1", "MUC5B", "SPDEF", "CD56"), group.by = "orig.ident")

```

```{r}
mergedNBG <- FindNeighbors(mergedNBG, dims = 1:10)
mergedNBG <- FindClusters(mergedNBG, resolution = 0.33)
```

look at UMAP clusters 
```{r}
mergedNBG <- RunUMAP(mergedNBG, dims = 1:10)
DimPlot(mergedNBG, reduction = "umap")
```


removed ciliated cells

look at 6 supervised clusters 
```{r}
mergedNBG.basal <- subset(mergedNBG, idents = c(0, 2, 3, 4, 5, 6))
mergedNBG.basal <- RunUMAP(mergedNBG.basal, dims = 1:10)
DimPlot(mergedNBG.basal, reduction = "umap")
```

```{r}
VlnPlot(object = mergedNBG.basal, features = c("KRT5", "KRT14", "NGFR", "ITGA6", "P63", "FOXJ1", "SCGB1A1", "MUC5B", "SPDEF", "CD56"))
```
```{r}
mergedNBG <- JackStraw(mergedNBG, num.replicate = 100)
mergedNBG <- ScoreJackStraw(mergedNBG, dims = 1:20)
JackStrawPlot(mergedNBG, dims = 1:15)
```
```{r}
ElbowPlot(mergedNBG)
```


cluster markers, repeat per cluster
```{r}
cluster0.markers <- FindMarkers(mergedNBG, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)
```

```{r}
cluster2.markers <- FindMarkers(mergedNBG, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)
```

```{r}
cluster3.markers <- FindMarkers(mergedNBG, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)
```

```{r}
cluster4.markers <- FindMarkers(mergedNBG, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
```

```{r}
cluster5.markers <- FindMarkers(mergedNBG, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 10)
```

```{r}
cluster6.markers <- FindMarkers(mergedNBG, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 10)
```

```{r}
mergedNBG.basal <- subset(mergedNBG, idents = c(0, 2, 3, 4, 5, 6))
mergedNBG.basal <- RunUMAP(mergedNBG.basal, dims = 1:10)
DimPlot(mergedNBG.basal, reduction = "umap", label = "TRUE", pt.size = 0.5)
```

```{r}
new.cluster.ids <- c("Transitional Basal", "Basal (TOP2A+) Proliferating", "Club", "Basal (KRT5+) Differentiated 1", "Club (SCGB3A1+)", "Basal (KRT5+) Differentiated 2")
names(new.cluster.ids) <- levels(mergedNBG.basal)
mergedNBG.basal <- RenameIdents(mergedNBG.basal, new.cluster.ids)
DimPlot(mergedNBG.basal, reduction = "umap", label = "TRUE", pt.size = 0.5)
```

```{r}
DimPlot(mergedNBG.basal, reduction = "umap", pt.size = 0.5, group.by = "orig.ident")
```

## Examining total 24 samples from published data (Carerro et. al., 2021)
#### Attempt 1
insert libraries
```{r}
library(dplyr)
library(Seurat)
library(Matrix)
library(patchwork)
```


incorrect 
```{r}
#data_dir <- "/Users/asingh/Downloads/GSE150674_Seurat_Object.rds"
#mergedcompleteBG <- Read10X(data.dir = data_dir)
```

pull in the already merged RDS file from GEO
```{r}
mergedcompleteBG <- readRDS(file = "/Users/asingh/Downloads/GSE150674_Seurat_Object.rds.gz")
```

i do see percent.mt here (finally)
```{r}
mergedcompleteBG[["percent.mt"]] <- PercentageFeatureSet(mergedcompleteBG, pattern = "^MT-")
VlnPlot(mergedcompleteBG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

scale the data:
```{r}
all.genes <- row.names(mergedcompleteBG)
mergedcompleteBG <- ScaleData(mergedcompleteBG, features = all.genes)
```

Normalizing
```{r}
mergedcompleteBG <- FindVariableFeatures(object = mergedcompleteBG, selection.method = "vst", nfeatures = 2000)
```

run PCA
```{r}
mergedcompleteBG <- RunPCA(mergedcompleteBG, features = VariableFeatures(object = mergedcompleteBG))
```

cluster cells
```{r}
mergedcompleteBG <- FindNeighbors(mergedcompleteBG, dims = 1:10)
mergedcompleteBG <- FindClusters(mergedcompleteBG, resolution = 0.2)
```

runUMAP
```{r}
mergedcompleteBG <- RunUMAP(mergedcompleteBG, dims = 1:10)
DimPlot(mergedcompleteBG, reduction = "umap")
```
#### Attempt 2
insert libraries:
```{r echo=TRUE}
library(Seurat)
library(SeuratObject)
```

import already created seurat object››
```{r import}
mergedBGcomplete <- readRDS(file = "/Users/asingh/Downloads/CSMC.CFF.UCLA.CO.CF.cleaned.reclust.slim.rds")
```

scale data
```{r scale data, message=FALSE}
all.genes <- rownames(mergedBGcomplete)
mergedBGcomplete <- ScaleData(mergedBGcomplete, features = all.genes)
```

normalize the data 
```{r normalize}
mergedBGcomplete <- FindVariableFeatures(object = mergedBGcomplete, selection.method = "vst", nfeatures = 2000)
```

run PCA
```{r}
mergedBGcomplete <- RunPCA(mergedBGcomplete, features = VariableFeatures(object = mergedBGcomplete))
```

find # clusters 
```{r}
mergedBGcomplete <- FindNeighbors(mergedBGcomplete, dims = 1:10)
mergedBGcomplete <- FindClusters(mergedBGcomplete, resolution = 0.2)
```

look at umap clusters
```{r}
mergedBGcomplete <- RunUMAP(mergedBGcomplete, dims = 1:10)
DimPlot(mergedBGcomplete, reduction = "umap")
```


inspect metadata 
```{r}
mergedBGcomplete[[]]
```

```{r}
jedumap <- readRDS(file = "/Users/asingh/Downloads/CO_CF_ALI_integrated.umap.done.rds")
```

```{r}
jedumap <- RunUMAP(jedumap, dims = 1:10)
DimPlot(jedumap, reduction = "umap")
```
#### Attempt 3

```{r}
library(dplyr)
library(Seurat)
library(Matrix)
library(patchwork)
```


```{r}
cff_seurat <- readRDS("/Volumes/aakdisc/GSE150674_Seurat_Object.rds")

```

```{r}

```


```{r}
cells1 <- CellsByIdentities(cff_seurat)
```

```{r}
cells1
```

```{r}
cff_seurat[[]]
```

```{r}
DimPlot(cff_seurat, reduction = "umap")
```

```{r}
cff_seurat <- FindVariableFeatures(cff_seurat, selection.method = "vst", nfeatures = 2000)
```

```{r}
cff_seurat <- FindNeighbors(cff_seurat, dims = 1:20)

```

```{r}
cff_seurat <- FindClusters(cff_seurat, resolution = 0.5)
```

```{r}
DimPlot(cff_seurat, reduction = "umap", group.by = "orig.ident")
```

## Batch 1 samples grouped by P0/P1

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

```{r}
CFdatafiltered <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CFfiltered")
```


```{r}
all.genes <- rownames(CFdatafiltered)
CFdatafiltered <- ScaleData(CFdatafiltered, features = all.genes)
```

```{r}
CFdatafiltered <- FindVariableFeatures(CFdatafiltered, selection.method = "vst", nfeatures = 2000)
```

```{r}
CFdatafiltered <- RunPCA(CFdatafiltered, features = VariableFeatures(object = CFdatafiltered))
```

```{r}
CFdatafiltered <- FindNeighbors(CFdatafiltered, dims = 1:10)
CFdatafiltered <- FindClusters(CFdatafiltered, resolution = 0.5)
```


```{r}
CFdatafiltered <- RunUMAP(CFdatafiltered, dims = 1:10)
DimPlot(CFdatafiltered)
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap", group.by = "orig.ident")
```

Batch 1 samples grouped by CF/CO
---
```{r}
library(Seurat)
```



read file 
```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "CF")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

```{r}
CFdatafiltered <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CFfiltered")
```



```{r}
all.genes <- rownames(CFdatafiltered)
CFdatafiltered <- ScaleData(CFdatafiltered, features = all.genes)
```

```{r}
CFdatafiltered <- FindVariableFeatures(CFdatafiltered, selection.method = "vst", nfeatures = 2000)
```

```{r}
CFdatafiltered <- RunPCA(CFdatafiltered, features = VariableFeatures(object = CFdatafiltered))
```

```{r}
CFdatafiltered <- FindNeighbors(CFdatafiltered, dims = 1:10)
CFdatafiltered <- FindClusters(CFdatafiltered, resolution = 0.5)
```


```{r}
CFdatafiltered <- RunUMAP(CFdatafiltered, dims = 1:10)
DimPlot(CFdatafiltered)

```


```{r}
DimPlot(CFdatafiltered, reduction = "umap", label = TRUE)
```




```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

```{r}
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

```{r}
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 
```

```{r}
gs_list = gene_sets_prepare(db_, tissue)
```

```{r}
es.max = sctype_score(scRNAseqData = CFdatafiltered[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
```

```{r}
cL_resutls = do.call("rbind", lapply(unique(CFdatafiltered@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(CFdatafiltered@meta.data[CFdatafiltered@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(CFdatafiltered@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 
```

```{r}
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
```

```{r}
CFdatafiltered@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  CFdatafiltered@meta.data$customclassif[CFdatafiltered@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```
```{r}
Idents(object = CFdatafiltered) <- CFdatafiltered@meta.data$'customclassif' 
blah2 <- FindMarkers(object = batch2, ident.1 = "Pulmonary alveolar type I cells")
```
```{r}
CFdatafiltered[[]]
```
```{r}
Idents(object = CFdatafiltered) <- CFdatafiltered@meta.data$'seurat_clusters' 
blah2 <- FindMarkers(object = batch2, ident.1 = 1)

cluster2.markers <- FindMarkers(CFdatafiltered, ident.1 = 5, min.pct = 0.25)
head(cluster2.markers, n = 5)
```
```{r}
head(cluster2.markers, n = 10)
```


```{r}
cluster0.markers <- FindMarkers(CFdatafiltered, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
```

```{r}
cluster1.markers <- FindMarkers(CFdatafiltered, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

```{r}
cluster2.markers <- FindMarkers(CFdatafiltered, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```

```{r}
cluster3.markers <- FindMarkers(CFdatafiltered, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
```

```{r}
cluster4.markers <- FindMarkers(CFdatafiltered, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
```

```{r}
cluster5.markers <- FindMarkers(CFdatafiltered, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```{r}
cluster6.markers <- FindMarkers(CFdatafiltered, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)
```

```{r}
cluster7.markers <- FindMarkers(CFdatafiltered, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 5)
```

```{r}
cluster8.markers <- FindMarkers(CFdatafiltered, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 5)
```

```{r}
cluster9.markers <- FindMarkers(CFdatafiltered, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)
```

```{r}
cluster10.markers <- FindMarkers(CFdatafiltered, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)
```

```{r}
cluster11.markers <- FindMarkers(CFdatafiltered, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 5)
```

```{r}
cluster12.markers <- FindMarkers(CFdatafiltered, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 5)
```

```{r}
cluster13.markers <- FindMarkers(CFdatafiltered, ident.1 = 13, min.pct = 0.25)
head(cluster13.markers, n = 5)
```

```{r}
cluster14.markers <- FindMarkers(CFdatafiltered, ident.1 = 14, min.pct = 0.25)
head(cluster14.markers, n = 5)
```

```{r}
FeaturePlot(CFdatafiltered, features = c("KRT6B", "KRT16", "SPRR1", "IVL"))
```
```{r}
FeaturePlot(CFdatafiltered, features = c("TOP2A", "MK167", "KRT15", "KRT5", "SCGB3A1"))
```

```{r}
FeaturePlot(CFdatafiltered, features = c("KRT6B", "KRT16", "SPRR1", "IVL"), max.cutoff = 5)
```
```{r}
FeaturePlot(object = CFdatafiltered, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"), max.cutoff = 5)
```

## Batch 1 and batch 2 with ScType

```{r}
library(Seurat)
```

import all batch 1:

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "bothbatches")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

import all batch 2: 


```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "bothbatches")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```

merge both 
```{r}
bothbatches <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0, so2.2, so2.4, so2.6, so2.7, so2.8, so2.1), project = "CFfiltered")
```

run through seurat pipeline

```{r}
all.genes <- rownames(bothbatches)
bothbatches <- ScaleData(bothbatches, features = all.genes)
```

```{r}
bothbatches <- FindVariableFeatures(bothbatches, selection.method = "vst", nfeatures = 2000)
```

```{r}
bothbatches <- RunPCA(bothbatches, features = VariableFeatures(object = bothbatches))
```

```{r}
bothbatches <- FindNeighbors(bothbatches, dims = 1:10)
bothbatches <- FindClusters(bothbatches, resolution = 0.5)
```


```{r}
bothbatches <- RunUMAP(bothbatches, dims = 1:10)
DimPlot(bothbatches)
```

```{r}
DimPlot(bothbatches, reduction = "umap", group.by = "orig.ident")
```


```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = bothbatches[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(bothbatches@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(bothbatches@meta.data[bothbatches@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(bothbatches@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

bothbatches@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  bothbatches@meta.data$customclassif[bothbatches@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}

DimPlot(bothbatches, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```

```{r}
bothbatches[[]]
```

```{r}
Idents(bothbatches) <- bothbatches$customclassif
meso_markers <- FindMarkers(bothbatches, ident.1="Mesothelial cells")
head(meso_markers, n = 10)
```

```{r}
Idents(bothbatches) <- bothbatches$customclassif
pul1_markers <- FindMarkers(bothbatches, ident.1="Pulmonary alveolar type I cells")
head(pul1_markers, n = 10)
```

```{r}
Idents(bothbatches) <- bothbatches$customclassif
pul2_markers <- FindMarkers(bothbatches, ident.1="Pulmonary alveolar type II cells")
head(pul2 _markers, n = 10)
```

## Batch 1 and Batch 2 grouped by passage
---
combining batch 1 & 2
seperate by P0/P1/P
---

```{r}
library(Seurat)
```

import all batch 1:

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

import all batch 2: 


```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "P1")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "P5")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "P1")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "P5")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "P5")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "P1")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```

merge both 
```{r}
bothbatches <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0, so2.2, so2.4, so2.6, so2.7, so2.8, so2.1), project = "CFfiltered")
```

```{r}
VlnPlot(bothbatches, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


run through seurat pipeline

```{r}
all.genes <- rownames(bothbatches)
bothbatches <- ScaleData(bothbatches, features = all.genes)
```

```{r}
bothbatches <- FindVariableFeatures(bothbatches, selection.method = "vst", nfeatures = 2000)
```

```{r}
bothbatches <- RunPCA(bothbatches, features = VariableFeatures(object = bothbatches))
```

```{r}
bothbatches <- FindNeighbors(bothbatches, dims = 1:10)
bothbatches <- FindClusters(bothbatches, resolution = 0.5)
```


```{r}
bothbatches <- RunUMAP(bothbatches, dims = 1:10)
DimPlot(bothbatches)
```

```{r}
DimPlot(bothbatches, reduction = "umap", group.by = "orig.ident")
```



```{r eval=FALSE, include=FALSE}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = bothbatches[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(bothbatches@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(bothbatches@meta.data[bothbatches@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(bothbatches@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

bothbatches@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  bothbatches@meta.data$customclassif[bothbatches@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(bothbatches, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```

```{r}
#basal stem cell genes
FeaturePlot(object = bothbatches, features = c("KRT5", "KRT14", "KRT15", "TP63", "ITGA6"), max.cutoff = 5)

```

```{r}
#basal squamous cells
FeaturePlot(object = bothbatches, features = c("SPRR1B", "SPRR3", "KRT13", "IVL"), max.cutoff = 5)

```


```{r}
# basal -- inflamatory cells
FeaturePlot(object = bothbatches, features = c("S100A9", "S110A8", "SERPINEB3", "SLPI"), max.cutoff = 5)

```

```{r}
# basal -- cycling genes
FeaturePlot(object = bothbatches, features = c("TOP2A", "MKI67"), max.cutoff = 5)

```

```{r}
# basal phenotypes of interest to them -- is this up in p5? 
FeaturePlot(object = bothbatches, features = c("KRT6B", "KRT16", "KRT17", "KRT23"), max.cutoff = 5)

```

```{r}
# ciliated 
FeaturePlot(object = bothbatches, features = c("FOXJ1", "DNAH5", "MCIDAS", "CCNO", "DNAH11", "SPAG1", "LRRC6"), max.cutoff = 5)

```


```{r}
# secretory cells 
FeaturePlot(object = bothbatches, features = c("SCGB1A1", "SCGB3A2", "MUC5AC", "MUC5B", "SPDEF", "SPDEF"), max.cutoff = 5)

```

```{r}
# ciliated progenitors 
FeaturePlot(object = bothbatches, features = c("FOXN4"), max.cutoff = 5)

```


```{r}
# ionocytes
FeaturePlot(object = bothbatches, features = c("FOXI1"), max.cutoff = 5)

```

```{r}
#neuroendocrine
FeaturePlot(object = bothbatches, features = c("ASCL1"), max.cutoff = 5)

```

```{r}
# cf gene 
FeaturePlot(object = bothbatches, features = c("CFTR"), max.cutoff = 5)

```

## Batch 1 and batch 2 grouped by CO/CF
---
combining batch 1 & batch 2
sep CF vs CO
---

```{r}
library(Seurat)
```

import all batch 1:

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "CF")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

import all batch 2: 


```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "Healthy")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "Healthy")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "CF")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "CF")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "Healthy")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220801T183018Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "Healthy")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```

merge both 
```{r}
bothbatches <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0, so2.2, so2.4, so2.6, so2.7, so2.8, so2.1), project = "CFfiltered")
```

run through seurat pipeline

```{r}
all.genes <- rownames(bothbatches)
bothbatches <- ScaleData(bothbatches, features = all.genes)
```

```{r}
bothbatches <- FindVariableFeatures(bothbatches, selection.method = "vst", nfeatures = 2000)
```

```{r}
bothbatches <- RunPCA(bothbatches, features = VariableFeatures(object = bothbatches))
```

```{r}
bothbatches <- FindNeighbors(bothbatches, dims = 1:10)
bothbatches <- FindClusters(bothbatches, resolution = 0.5)
```


```{r}
bothbatches <- RunUMAP(bothbatches, dims = 1:10)
DimPlot(bothbatches)
```

```{r}
DimPlot(bothbatches, reduction = "umap", group.by = "orig.ident")
```

## Batch 2 alone, filtered
---
batch2 filtered
---

```{r}
library(Seurat)
```



read file 

```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "1")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "2")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "4")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "6")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "7")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "8")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```



```{r}
batch2 <- merge(so2.1, y = c(so2.2, so2.4, so2.6, so2.7, so2.8), project = "batch2")
```

```{r}
batch2[[]]
```


```{r}

rm(batch2_1)
rm(batch2_2)
rm(batch2_4)
rm(batch2_6)
rm(batch2_7)
rm(batch2_8)
```


```{r}
batch2 <- NormalizeData(batch2)
```

```{r}
batch2 <- FindVariableFeatures(batch2, selection.method = "vst", nfeatures = 2000)
```


```{r}
all.genes <- rownames(batch2)
batch2 <- ScaleData(batch2, features = all.genes)
```


```{r}
batch2 <- RunPCA(batch2, features = VariableFeatures(object = batch2))
```

```{r}
ElbowPlot(batch2)
```



```{r}
batch2 <- FindNeighbors(batch2, dims = 1:10)
batch2 <- FindClusters(batch2, resolution = 0.5)
```

```{r}
batch2 <- FindNeighbors(batch2, dims = 1:15)
batch2 <- FindClusters(batch2, resolution = 0.5)
```


```{r}
batch2 <- RunUMAP(batch2, dims = 1:10)
DimPlot(batch2, reduction = "umap")
```

```{r}
DimPlot(batch2, reduction = "umap", group.by = "orig.ident")
```


```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

```{r}
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

```{r}
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 
```

```{r}
gs_list = gene_sets_prepare(db_, tissue)
```

```{r}
es.max = sctype_score(scRNAseqData = batch2[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
```

```{r}
cL_resutls = do.call("rbind", lapply(unique(batch2@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(batch2@meta.data[batch2@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(batch2@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 
```

```{r}
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
```

```{r}
batch2@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  batch2@meta.data$customclassif[batch2@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
```

```{r}
DimPlot(batch2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```

```{r}
FeaturePlot(object = batch2, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"), max.cutoff = 5)
```

```{r}
FeaturePlot(batch2, features = c("KRT6B", "KRT16", "SPRR1", "IVL"), max.cutoff = 5)
```

```{r}
batch2[[]
```

```{r}
Idents(object = batch2) <- batch2@meta.data$'customclassif'
blah <- FindMarkers(object = batch2, ident.1 = "Ionocytes")
```

```{r}
head(blah, n = 10)
```
```{r}
Idents(object = batch2) <- batch2@meta.data$'customclassif'
blah3 <- FindMarkers(object = batch2, ident.1 = "Pulmonary alveolar type I cells")
```

```{r}
head(blah2, n = 10)
```
```{r}
cluster2.markers <- FindMarkers(batch2, ident.1 = "2", min.pct = 0.25)
head(cluster2.markers, n = 5)
```

### Batch 2 unfiltered

```{r}
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(rhdf5)

```

```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "1")
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "2")
```

```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "4")
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "6")
```

```{r}
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "7")
```

```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "8")
```


```{r}
batch2 <- merge(so2.1, y = c(so2.2, so2.4, so2.6, so2.7, so2.8), project = "batch2")
```

```{r}
batch2
```

rm old datasets to save space
```{r}
rm(batch2_1)
rm(batch2_2)
rm(batch2_4)
rm(batch2_6)
rm(batch2_7)
rm(batch2_8)
```

```{r}

```

```{r}
batch2[[]]
```


```{r}
batch2[["percent.mt"]] <- PercentageFeatureSet(batch2, pattern = "^MT")
```

```{r}
VlnPlot(batch2, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
```

```{r}
batch2 <- NormalizeData(batch2, normalization.method = "LogNormalize", scale.factor = 10000)
```


```{r}
batch2 <- FindVariableFeatures(batch2, selection.method = "vst", nfeatures = 2000)
```

```{r}
all.genes <- rownames(batch2)
batch2 <- ScaleData(batch2, features = all.genes)
```

```{r}
batch2 <- RunPCA(batch2, features = VariableFeatures(object = batch2))
```


```{r}
batch2 <- FindNeighbors(batch2, dims = 1:10)
batch2 <- FindClusters(batch2, resolution = 0.5)
```

```{r}
batch2 <- RunUMAP(batch2, dims = 1:10)
DimPlot(batch2, reduction = "umap")
```
```{r}
DimPlot(batch2, reduction = "umap", group.by = "orig.ident")
```

```{r}
FeaturePlot(object = batch2, features = c("TP63", "KRT5", "NGFR", "KRT8", "ITGA6", "KRT14", "ITGB4", "SCGB1A1", "FOXJ1", "DNAH5", "DNAH11", "PCNT", "MUC5B", "PDPN", "AQP5", "GPRC5A"))
```

```{r}
FeaturePlot(object = batch2, features = "ASCL3")
```

## Batch 1 and batch 2 merged, filtered, grouped by CF/CO

checking the file 
```{r}
h5ls("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5")
```
read file 
```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```

```{r}

```

```{r}

```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 25)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 25)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "CF")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 25)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 25)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 25)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 25)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 25)
SO_H18_P0[[]]
```

```{r}
CFdatafiltered <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CFfiltered")
```


```{r}
all.genes <- rownames(CFdatafiltered)
CFdatafiltered <- ScaleData(CFdatafiltered, features = all.genes)
```
```{r}
CFdatafiltered <- FindVariableFeatures(CFdatafiltered, selection.method = "vst", nfeatures = 2000)
```

```{r}
CFdatafiltered <- RunPCA(CFdatafiltered, features = VariableFeatures(object = CFdatafiltered))
```

```{r}
CFdatafiltered <- FindNeighbors(CFdatafiltered, dims = 1:10)
CFdatafiltered <- FindClusters(CFdatafiltered, resolution = 0.5)
```


```{r}
CFdatafiltered <- RunUMAP(CFdatafiltered, dims = 1:10)
DimPlot(CFdatafiltered)
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap", group.by = "orig.ident")
```

## Batch 1 unfiltered 
---
title: "MERGED FILES"
output: html_notebook
---

```{r}
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(rhdf5)
```

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "CF")
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "CF")
```

```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "CF")
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "CF")
```

```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "CF")
```

```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "CF")
```

```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "CF")
```

```{r}
CFdata <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CF")
```

```{r}
CFdata
```

rm old datasets to save space
```{r}
rm(CF13_P1)
rm(CF10_P1)
rm(CF16_P0)
rm(H14_P1)
rm(H15_P0)
rm(H15_P1)
rm(H18_P0)
```



```{r}
head(colnames(CFdata))
```

```{r}
table(CFdata$orig.ident)
```

```{r}
CFdata[["percent.mt"]] <- PercentageFeatureSet(CFdata, pattern = "^MT")
```

```{r}
VlnPlot(CFdata, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
```

```{r}
CFdata <- NormalizeData(CFdata, normalization.method = "LogNormalize", scale.factor = 10000)
```

```{r}
CFdata <- FindVariableFeatures(CFdata, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(CFdata), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(CFdata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

```{r}
all.genes <- rownames(CFdata)
CFdata <- ScaleData(CFdata, features = all.genes)
```

```{r}
CFdata <- RunPCA(CFdata, features = VariableFeatures(object = CFdata))
```

```{r}
# Examine and visualize PCA results a few different ways
print(CFdata[["pca"]], dims = 1:5, nfeatures = 5)
```

```{r}
VizDimLoadings(CFdata, dims = 1:2, reduction = "pca")
```

```{r}
DimPlot(CFdata, reduction = "pca")
```

```{r}
DimHeatmap(CFdata, dims = 1:15, cells = 500, balanced = TRUE)
```

```{r}
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
CFdata <- JackStraw(CFdata, num.replicate = 100)
CFdata <- ScoreJackStraw(CFdata, dims = 1:20)
```

```{r}
JackStrawPlot(CFdata, dims = 1:15)
```

```{r}
ElbowPlot(CFdata)
```

```{r}
CFdata <- FindNeighbors(CFdata, dims = 1:10)
CFdata <- FindClusters(CFdata, resolution = 0.5)
```

```{r}
CFdata <- RunUMAP(CFdata, dims = 1:10)
```

```{r}
DimPlot(CFdata, reduction = "umap")
```

```{r}
cluster0.markers <- FindMarkers(CFdata, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
```

```{r}
cluster1.markers <- FindMarkers(CFdata, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

```{r}
cluster2.markers <- FindMarkers(CFdata, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```

```{r}
cluster2.markers <- FindMarkers(CFdata, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```

```{r}
cluster3.markers <- FindMarkers(CFdata, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
```

```{r}
cluster4.markers <- FindMarkers(CFdata, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
```

```{r}
cluster5.markers <- FindMarkers(CFdata, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```{r}
cluster6.markers <- FindMarkers(CFdata, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)
```

```{r}
cluster7.markers <- FindMarkers(CFdata, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 5)
```

```{r}
cluster8.markers <- FindMarkers(CFdata, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 5)
```

```{r}
cluster9.markers <- FindMarkers(CFdata, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)
```

```{r}
cluster10.markers <- FindMarkers(CFdata, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)
```

```{r}
cluster11.markers <- FindMarkers(CFdata, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 5)
```

```{r}
cluster12.markers <- FindMarkers(CFdata, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 5)
```

```{r}
cluster13.markers <- FindMarkers(CFdata, ident.1 = 13, min.pct = 0.25)
head(cluster13.markers, n = 5)
```

```{r}
cluster14.markers <- FindMarkers(CFdata, ident.1 = 14, min.pct = 0.25)
head(cluster14.markers, n = 5)
```

```{r}
cluster15.markers <- FindMarkers(CFdata, ident.1 = 15, min.pct = 0.25)
head(cluster15.markers, n = 5)
```


```{r}
CFdata[[]]
```

```{r}
CFdata <- FindNeighbors(CFdata, dims = 1:10)
CFdata <- FindClusters(CFdata, resolution = 0.1)

```

```{r}
CFdata <- RunUMAP(CFdata, dims = 1:10)
DimPlot(CFdata, reduction = "umap")
```

```{r}
CFdata <- FindNeighbors(CFdata, dims = 1:10)
CFdata <- FindClusters(CFdata, resolution = 0.3)
```

```{r}
CFdata <- RunUMAP(CFdata, dims = 1:10)
DimPlot(CFdata, reduction = "umap")
```
## Batch 1 filtered
---
title: "LOOK AT DATASETS INDIVIDUALLY"
output: html_notebook
---

don't need 
```{r}
# install.packages("BiocManager")
# BiocManager::install("rhdf5")
```

deleted and reinstalled seurat 
```{r}
# install.packages("Seurat", repo = 'https://mac.R-project.org')
```

load libraries
```{r}
library(Seurat)
library(hdf5r)
library(dplyr)
library(patchwork)
library(rhdf5)
```

checking the file 
```{r}
h5ls("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5")
```
read file 
```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "SO_CF13_P1")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```

```{r}

```

```{r}

```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```
```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 25)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "SO_CF10_P1")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 25)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "SO_CF16_P0")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 25)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "SO_H14_P1")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 25)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "SO_H15_P0")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 25)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "SO_H15_P1")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 25)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "SO_H18_P0")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 25)
SO_H18_P0[[]]
```

```{r}
CFdatafiltered <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CFfiltered")
```

```{r}
CFdatafiltered[[]]
```

```{r}
CFdatafiltered[[]]
```


```{r}
VlnPlot(CFdatafiltered, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)
```

```{r}
all.genes <- rownames(CFdatafiltered)
CFdatafiltered <- ScaleData(CFdatafiltered, features = all.genes)
```

```{r}
CFdatafiltered <- FindVariableFeatures(CFdatafiltered, selection.method = "vst", nfeatures = 2000)
```


```{r}
CFdatafiltered <- RunPCA(CFdatafiltered, features = VariableFeatures(object = CFdatafiltered))
```

```{r}
DimPlot(CFdatafiltered, reduction = "pca", group.by = "orig.ident")
```

```{r}
ElbowPlot(CFdatafiltered)
```

```{r}
CFdatafiltered <- FindNeighbors(CFdatafiltered, dims = 1:10)
CFdatafiltered <- FindClusters(CFdatafiltered, resolution = 0.5)

```

```{r}
CFdatafiltered <- RunUMAP(CFdatafiltered, dims = 1:10)
DimPlot(CFdatafiltered)
```
```{r}
cluster0.markers <- FindMarkers(CFdatafiltered, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
```

```{r}
cluster1.markers <- FindMarkers(CFdatafiltered, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```


```{r}
cluster2.markers <- FindMarkers(CFdatafiltered, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```


```{r}
cluster2.markers <- FindMarkers(CFdatafiltered, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```

```{r}
cluster3.markers <- FindMarkers(CFdatafiltered, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
```

```{r}
cluster4.markers <- FindMarkers(CFdatafiltered, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 5)
```

```{r}
cluster5.markers <- FindMarkers(CFdatafiltered, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```{r}
cluster6.markers <- FindMarkers(CFdatafiltered, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)
```

```{r}
cluster7.markers <- FindMarkers(CFdatafiltered, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 5)
```

```{r}
cluster8.markers <- FindMarkers(CFdatafiltered, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 5)
```

```{r}
cluster9.markers <- FindMarkers(CFdatafiltered, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 5)
```

```{r}
cluster10.markers <- FindMarkers(CFdatafiltered, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 5)
```

```{r}
cluster11.markers <- FindMarkers(CFdatafiltered, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 5)
```

```{r}
cluster12.markers <- FindMarkers(CFdatafiltered, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 5)
```

```{r}
cluster13.markers <- FindMarkers(CFdatafiltered, ident.1 = 13, min.pct = 0.25)
head(cluster13.markers, n = 5)
```

```{r}
cluster14.markers <- FindMarkers(CFdatafiltered, ident.1 = 14, min.pct = 0.25)
head(cluster14.markers, n = 5)
```

```{r}
DimPlot(CFdatafiltered, reduction = "pca", group.by = "orig.ident")
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap")
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap", label = "TRUE", pt.size = 0.2)
```


```{r}
DimPlot(CFdatafiltered, reduction = "umap", group.by = "orig.ident")
```

```{r}

```

```{r}
install.packages("BiocManager")
BiocManager::install("monocle")
```

```{r}
library(monocle)
```

```{r}
as.CellDataSet(CFdatafiltered)
```

```{r}
CFdatafiltered
```

```{r}

```



## Both batches merged, grouped by P0.P1
---
group by P0 vs P1
---

read file 
```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```

```{r}

```

```{r}

```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 25)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 25)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 25)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 25)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 25)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 25)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 25)
SO_H18_P0[[]]
```

```{r}
CFdatafiltered <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0), project = "CFfiltered")
```


```{r}
all.genes <- rownames(CFdatafiltered)
CFdatafiltered <- ScaleData(CFdatafiltered, features = all.genes)
```

```{r}
CFdatafiltered <- FindVariableFeatures(CFdatafiltered, selection.method = "vst", nfeatures = 2000)
```

```{r}
CFdatafiltered <- RunPCA(CFdatafiltered, features = VariableFeatures(object = CFdatafiltered))
```

```{r}
CFdatafiltered <- FindNeighbors(CFdatafiltered, dims = 1:10)
CFdatafiltered <- FindClusters(CFdatafiltered, resolution = 0.5)
```


```{r}
CFdatafiltered <- RunUMAP(CFdatafiltered, dims = 1:10)
DimPlot(CFdatafiltered)
```

```{r}
DimPlot(CFdatafiltered, reduction = "umap", group.by = "orig.ident")
```

## Using Azimuth for ATAC data
---
azimuth pipeline
---


Load libraries.
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```



```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

CF16_P0 <- CreateSeuratObject(counts = rna_counts)
CF16_P0[["percent.mt"]] <- PercentageFeatureSet(CF16_P0, pattern = "^MT-")

CF16_P0 <- RunAzimuth(CF16_P0, reference = "lungref")


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
CF16_P0[["ATAC"]] <- chrom_assay

VlnPlot(CF16_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
CF16_P0 <- subset(
  x = CF16_P0,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(CF16_P0) <- "RNA"
CF16_P0 <- SCTransform(CF16_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(CF16_P0) <- "ATAC"
CF16_P0 <- RunTFIDF(CF16_P0)
CF16_P0 <- FindTopFeatures(CF16_P0, min.cutoff = 'q0')
CF16_P0 <- RunSVD(CF16_P0)
CF16_P0 <- RunUMAP(CF16_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

CF16_P0 <- FindMultiModalNeighbors(CF16_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CF16_P0 <- RunUMAP(CF16_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CF16_P0 <- FindClusters(CF16_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(CF16_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF16_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF16_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(CF16_P0, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p4 

```
```{r}
p1 <- DimPlot(CF16_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF16_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF16_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```


```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

CF13_P1 <- CreateSeuratObject(counts = rna_counts)
CF13_P1[["percent.mt"]] <- PercentageFeatureSet(CF13_P1, pattern = "^MT-")

CF13_P1 <- RunAzimuth(CF13_P1, reference = "lungref")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
CF13_P1[["ATAC"]] <- chrom_assay

VlnPlot(CF13_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

CF13_P1 <- subset(
  x = CF13_P1,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(CF13_P1) <- "RNA"
CF13_P1 <- SCTransform(CF13_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(CF13_P1) <- "ATAC"
CF13_P1 <- RunTFIDF(CF13_P1)
CF13_P1 <- FindTopFeatures(CF13_P1, min.cutoff = 'q0')
CF13_P1 <- RunSVD(CF13_P1)
CF13_P1 <- RunUMAP(CF13_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

CF13_P1 <- FindMultiModalNeighbors(CF13_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CF13_P1 <- RunUMAP(CF13_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CF13_P1 <- FindClusters(CF13_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(CF13_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF13_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF13_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(CF16_P0, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)
 
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p4
```

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H14_P1 <- CreateSeuratObject(counts = rna_counts)
H14_P1[["percent.mt"]] <- PercentageFeatureSet(H14_P1, pattern = "^MT-")

H14_P1 <- RunAzimuth(H14_P1, reference = "lungref")


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H14_P1[["ATAC"]] <- chrom_assay

VlnPlot(H14_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

H14_P1 <- subset(
  x = H14_P1,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(H14_P1) <- "RNA"
H14_P1 <- SCTransform(H14_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H14_P1) <- "ATAC"
H14_P1 <- RunTFIDF(H14_P1)
H14_P1 <- FindTopFeatures(H14_P1, min.cutoff = 'q0')
H14_P1 <- RunSVD(H14_P1)
H14_P1 <- RunUMAP(H14_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H14_P1 <- FindMultiModalNeighbors(H14_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H14_P1 <- RunUMAP(H14_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H14_P1 <- FindClusters(H14_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H14_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H14_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H14_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(H14_P1, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p4
```

--- NEXT SAMPLE --- 
H15_P0
error: H15_P0 <- FindMultiModalNeighbors(H15_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H15_P0 <- CreateSeuratObject(counts = rna_counts)
H15_P0[["percent.mt"]] <- PercentageFeatureSet(H15_P0, pattern = "^MT-")

H15_P0 <- RunAzimuth(H15_P0, reference = "lungref")



grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H15_P0[["ATAC"]] <- chrom_assay

VlnPlot(H15_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

H15_P0 <- subset(
  x = H15_P0,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(H15_P0) <- "RNA"
H15_P0 <- SCTransform(H15_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H15_P0) <- "ATAC"
H15_P0 <- RunTFIDF(H15_P0)
H15_P0 <- FindTopFeatures(H15_P0, min.cutoff = 'q0')
H15_P0 <- RunSVD(H15_P0)
H15_P0 <- RunUMAP(H15_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H15_P0 <- FindMultiModalNeighbors(H15_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H15_P0 <- RunUMAP(H15_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H15_P0 <- FindClusters(H15_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H15_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H15_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H15_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(H15_P0, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p4
```

--- NEXT SAMPLE --- 
H15_P1

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H15_P1 <- CreateSeuratObject(counts = rna_counts)
H15_P1[["percent.mt"]] <- PercentageFeatureSet(H15_P1, pattern = "^MT-")

H15_P1 <- RunAzimuth(H15_P1, reference = "lungref")



grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H15_P1[["ATAC"]] <- chrom_assay

VlnPlot(H15_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

H15_P1 <- subset(
  x = H15_P1,
    percent.mt < 20
)
# RNA analysis
DefaultAssay(H15_P1) <- "RNA"
H15_P1 <- SCTransform(H15_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H15_P1) <- "ATAC"
H15_P1 <- RunTFIDF(H15_P1)
H15_P1 <- FindTopFeatures(H15_P1, min.cutoff = 'q0')
H15_P1 <- RunSVD(H15_P1)
H15_P1 <- RunUMAP(H15_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H15_P1 <- FindMultiModalNeighbors(H15_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H15_P1 <- RunUMAP(H15_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H15_P1 <- FindClusters(H15_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H15_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H15_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H15_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(H15_P1, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p4
```

--- NEXT SAMPLE --- 
H18_P0

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H18_P0 <- CreateSeuratObject(counts = rna_counts)
H18_P0[["percent.mt"]] <- PercentageFeatureSet(H18_P0, pattern = "^MT-")
H18_P0 <- RunAzimuth(H18_P0, reference = "lungref")




grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H18_P0[["ATAC"]] <- chrom_assay

VlnPlot(H18_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
H18_P0 <- subset(
  x = H18_P0,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(H18_P0) <- "RNA"
H18_P0 <- SCTransform(H18_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H18_P0) <- "ATAC"
H18_P0 <- RunTFIDF(H18_P0)
H18_P0 <- FindTopFeatures(H18_P0, min.cutoff = 'q0')
H18_P0 <- RunSVD(H18_P0)
H18_P0 <- RunUMAP(H18_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H18_P0 <- FindMultiModalNeighbors(H18_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H18_P0 <- RunUMAP(H18_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H18_P0 <- FindClusters(H18_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H18_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H18_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H18_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(H18_P0, reduction = "wnn.umap", group.by = "predicted.ann_level_3", label = TRUE, label.size = 3)

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p4
```

## ATAC WNN analysis
---
merging objects together -- atac wnn analysis with azimuth
---

Load libraries
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```

Load in samples -- NOTE: missing CF10_P1 tsv.gz file. Not included in this merge. 
```{r}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
inputdata.CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H14_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")
#inputdata.H15_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5files/HEALTHY-18-015_P0_filtered_feature_bc_matrix (1).h5")
inputdata.H15_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H18_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")
```

Read RNA and ATAC data; create seurat objects
```{r}
rna_counts_CF16_P0 <- inputdata.CF16_P0$`Gene Expression`
atac_counts_CF16_P0 <- inputdata.CF16_P0$Peaks
CF16_P0 <- CreateSeuratObject(counts = rna_counts_CF16_P0)

rna_counts_CF13_P1 <- inputdata.CF13_P1$`Gene Expression`
atac_count_CF13_P1 <- inputdata.CF13_P1$Peaks
CF13_P1 <- CreateSeuratObject(counts = rna_counts_CF13_P1)

rna_counts_H14_P1 <- inputdata.H14_P1$`Gene Expression`
atac_counts_H14_P1 <- inputdata.H14_P1$Peaks
H14_P1 <- CreateSeuratObject(counts = rna_counts_H14_P1)

#rna_counts_H15_P0 <- inputdata.H15_P0$`Gene Expression`
#atac_counts_H15_P0 <- inputdata.H15_P0$Peaks
#H15_P0 <- CreateSeuratObject(counts = rna_counts_H15_P0)

rna_counts_H15_P1 <- inputdata.H15_P1$`Gene Expression`
atac_counts_H15_P1 <- inputdata.H15_P1$Peaks
H15_P1 <- CreateSeuratObject(counts = rna_counts_H15_P1)

rna_counts_H18_P0 <- inputdata.H18_P0$`Gene Expression`
atac_counts_H18_P0 <- inputdata.H18_P0$Peaks
H18_P0 <- CreateSeuratObject(counts = rna_counts_H18_P0)
```

Merge RNA objects together.
```{r}
wnnmerged <- merge(CF16_P0, y = c(CF13_P1, H14_P1, H15_P1, H18_P0), project = "batch1wnn")
```
```{r}
ncol(wnnmerged)
```


Load atac fragment files 
grange.counts <- StringToGRanges(rownames(rna_counts_H15_P0), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
rna_counts_H15_P0 <- rna_counts_H15_P0[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P0_atac_fragments.tsv.gz"
chrom_assay_H15_P0 <- CreateChromatinAssay(
   counts = rna_counts_H15_P0,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
```{r}
grange.counts <- StringToGRanges(rownames(atac_counts_CF16_P0), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts_CF16_P0 <- atac_counts_CF16_P0[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file_CF16_P0 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
chrom_assay_CF16_P0 <- CreateChromatinAssay(
   counts = atac_counts_CF16_P0,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file_CF16_P0,
   min.cells = 10,
   annotation = annotations
 )

grange.counts <- StringToGRanges(rownames(atac_count_CF13_P1), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_count_CF13_P1 <- atac_count_CF13_P1[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file_CF13_P1 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
chrom_assay_CF13_P1 <- CreateChromatinAssay(
   counts = atac_count_CF13_P1,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file_CF13_P1,
   min.cells = 10,
   annotation = annotations
 )

grange.counts <- StringToGRanges(rownames(atac_counts_H14_P1), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts_H14_P1 <- atac_counts_H14_P1[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file_H14_P1 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
chrom_assay_H14_P1 <- CreateChromatinAssay(
   counts = atac_counts_H14_P1,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file_H14_P1,
   min.cells = 10,
   annotation = annotations
 )



grange.counts <- StringToGRanges(rownames(atac_counts_H15_P1), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts_H15_P1 <- atac_counts_H15_P1[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file_H15_P1 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
chrom_assay_H15_P1 <- CreateChromatinAssay(
   counts = atac_counts_H15_P1,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file_H15_P1,
   min.cells = 10,
   annotation = annotations
 )

grange.counts <- StringToGRanges(rownames(atac_counts_H18_P0), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts_H18_P0 <- atac_counts_H18_P0[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
frag.file_H18_P0 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"
chrom_assay_H18_P0 <- CreateChromatinAssay(
   counts = atac_counts_H18_P0,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file_H18_P0,
   min.cells = 10,
   annotation = annotations
 )
```

Merge atac objects together
```{r}
chrom_assay <- merge(
  x = chrom_assay_CF16_P0, 
  y = list(chrom_assay_CF13_P1, chrom_assay_H14_P1, chrom_assay_H14_P1, chrom_assay_H15_P1, chrom_assay_H18_P0)
)
```

```{r}
ncol(chrom_assay)
```


```{r}
atacmerged <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
```


```{r}
wnnmerged[["ATAC"]] <- chrom_assay

VlnPlot(wnnmerged, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
```


```{r}
atacmerged[["RNA"]] <- wnnmerged
```

```{r}
atacmerged[["RNA"]] <- CreateAssayObject(counts = counts$"Gene Expression"[,colnames(atacmerged)])
```

## ATAC merged and annotated batch 2
---
ATAC MERGED ANNOTATED FILTERED FULL BATCH1
---
load libraries
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```

```{r eval=FALSE, include=FALSE}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
frag.file160 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)

object160_2 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`)
object160_2 <- RunAzimuth(object160_2, reference = "lungref")

object160_2[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)

object160_2 <- subset(
  x = object160_2,
    percent.mt < 20
)
```


```{r}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
inputdata.CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H14_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")
inputdata.H15_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H18_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")


frag.file160 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
frag.file131 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
frag.file141 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
frag.file151 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
frag.file180 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"

chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)
chromassay131 <- CreateChromatinAssay(counts = inputdata.CF13_P1$Peaks, sep = c(":", "-"), fragments = frag.file131, min.cells = 1, min.features = -1)
chromassay141 <- CreateChromatinAssay(counts = inputdata.H14_P1$Peaks, sep = c(":", "-"), fragments = frag.file141, min.cells = 1, min.features = -1)
chromassay151 <- CreateChromatinAssay(counts = inputdata.H15_P1$Peaks, sep = c(":", "-"), fragments = frag.file151, min.cells = 1, min.features = -1)
chromassay180 <- CreateChromatinAssay(counts = inputdata.H18_P0$Peaks, sep = c(":", "-"), fragments = frag.file180, min.cells = 1, min.features = -1)

#CHECK MERGE
#chromassay <- merge(
 # x = chromassay160, 
  #y = list(chromassay131, chromassay141, chromassay151, chromassay180)
#)

object160 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`)
object160 <- RunAzimuth(object160, reference = "lungref")
object160[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)
object160[["percent.mt"]] <- PercentageFeatureSet(object160, pattern = "MT")
object160 <- subset(x = object160, percent.mt < 20 )

object131 <- CreateSeuratObject(counts = inputdata.CF13_P1$`Gene Expression`)
object131 <- RunAzimuth(object131, reference = "lungref")
object131[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF13_P1$`Peaks`)
object131[["percent.mt"]] <- PercentageFeatureSet(object131, pattern = "MT")
object131 <- subset(x = object131, percent.mt < 20 )

object141 <- CreateSeuratObject(counts = inputdata.H14_P1$`Gene Expression`)
object141 <- RunAzimuth(object141, reference = "lungref")
object141[["ATAC"]] <- CreateAssayObject(counts = inputdata.H14_P1$`Peaks`)
object141[["percent.mt"]] <- PercentageFeatureSet(object141, pattern = "MT")
object141 <- subset(x = object141, percent.mt < 20 )

object151 <- CreateSeuratObject(counts = inputdata.H15_P1$`Gene Expression`)
object151 <- RunAzimuth(object151, reference = "lungref")
object151[["ATAC"]] <- CreateAssayObject(counts = inputdata.H15_P1$`Peaks`)
object151[["percent.mt"]] <- PercentageFeatureSet(object151, pattern = "MT")
object151 <- subset(x = object151, percent.mt < 20 )

#object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
#object5[["RNA"]] < CreateAssayObject(counts = inputdata.H18_P0$`Gene Expression`)

#mergedobject <- merge(object, y = c(object2, object141, object151, object5), project = "atacmerge_b1")
```

error: Cannot find RNA in this Seurat object 
```{r eval=FALSE, include=FALSE}
object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
inputdata.H18_P0$`Gene Expression`
```

```{r}
mergedobject <- merge(object160, y = c(object131, object141, object151), project = "atacmerge_b1")
```

```{r}
# RNA analysis
DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)
```

```{r}
mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

```{r}
p4 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p5 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p6 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)

p4 + p5 + p6
```
```{r}
mergedobject[[]]
```

## ATAC merged and annotated batch 1
---
ATAC MERGED ANNOTATED FILTERED FULL BATCH1
---
load libraries
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```

```{r eval=FALSE, include=FALSE}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
frag.file160 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)

object160_2 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`)
object160_2 <- RunAzimuth(object160_2, reference = "lungref")

object160_2[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)

object160_2 <- subset(
  x = object160_2,
    percent.mt < 20
)
```


```{r}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
inputdata.CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H14_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")
inputdata.H15_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H18_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")


frag.file160 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
frag.file131 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
frag.file141 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
frag.file151 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
frag.file180 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"

chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)
chromassay131 <- CreateChromatinAssay(counts = inputdata.CF13_P1$Peaks, sep = c(":", "-"), fragments = frag.file131, min.cells = 1, min.features = -1)
chromassay141 <- CreateChromatinAssay(counts = inputdata.H14_P1$Peaks, sep = c(":", "-"), fragments = frag.file141, min.cells = 1, min.features = -1)
chromassay151 <- CreateChromatinAssay(counts = inputdata.H15_P1$Peaks, sep = c(":", "-"), fragments = frag.file151, min.cells = 1, min.features = -1)
chromassay180 <- CreateChromatinAssay(counts = inputdata.H18_P0$Peaks, sep = c(":", "-"), fragments = frag.file180, min.cells = 1, min.features = -1)

#CHECK MERGE
#chromassay <- merge(
 # x = chromassay160, 
  #y = list(chromassay131, chromassay141, chromassay151, chromassay180)
#)

object160 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`)
object160 <- RunAzimuth(object160, reference = "lungref")
object160[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)
object160[["percent.mt"]] <- PercentageFeatureSet(object160, pattern = "MT")
object160 <- subset(x = object160, percent.mt < 20 )

object131 <- CreateSeuratObject(counts = inputdata.CF13_P1$`Gene Expression`)
object131 <- RunAzimuth(object131, reference = "lungref")
object131[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF13_P1$`Peaks`)
object131[["percent.mt"]] <- PercentageFeatureSet(object131, pattern = "MT")
object131 <- subset(x = object131, percent.mt < 20 )

object141 <- CreateSeuratObject(counts = inputdata.H14_P1$`Gene Expression`)
object141 <- RunAzimuth(object141, reference = "lungref")
object141[["ATAC"]] <- CreateAssayObject(counts = inputdata.H14_P1$`Peaks`)
object141[["percent.mt"]] <- PercentageFeatureSet(object141, pattern = "MT")
object141 <- subset(x = object141, percent.mt < 20 )

object151 <- CreateSeuratObject(counts = inputdata.H15_P1$`Gene Expression`)
object151 <- RunAzimuth(object151, reference = "lungref")
object151[["ATAC"]] <- CreateAssayObject(counts = inputdata.H15_P1$`Peaks`)
object151[["percent.mt"]] <- PercentageFeatureSet(object151, pattern = "MT")
object151 <- subset(x = object151, percent.mt < 20 )

#object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
#object5[["RNA"]] < CreateAssayObject(counts = inputdata.H18_P0$`Gene Expression`)

#mergedobject <- merge(object, y = c(object2, object141, object151, object5), project = "atacmerge_b1")
```

error: Cannot find RNA in this Seurat object 
```{r eval=FALSE, include=FALSE}
object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
inputdata.H18_P0$`Gene Expression`
```

```{r}
mergedobject <- merge(object160, y = c(object131, object141, object151), project = "atacmerge_b1")
```

```{r}
# RNA analysis
DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)
```

```{r}
mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

```{r}
p4 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p5 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p6 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)

p4 + p5 + p6
```
```{r}
mergedobject[[]]
```

## WNN analysis 
---
WNN analysis -- all samples from batch2
---

Load packages
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
```

Load h5 file and extract RNA and ATAC data; create Seurat object with ATAC data
```{r}
inputdata.10x <- Read10X_h5("")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample1b <- CreateSeuratObject(counts = rna_counts)
sample1b[["percent.mt"]] <- PercentageFeatureSet(sample1b, pattern = "^MT-")

```

Add ATAC data
```{r}
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- ""
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample1b[["ATAC"]] <- chrom_assay
```

Basic QC
```{r}
VlnPlot(sample1b, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
```
Can perform QC according to above, i didn;t yet 

```{r}
##
```

Pre-processing
```{r}
# RNA analysis
DefaultAssay(sample1b) <- "RNA"
sample1b <- SCTransform(sample1b, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample1b) <- "ATAC"
sample1b <- RunTFIDF(sample1b)
sample1b <- FindTopFeatures(sample1b, min.cutoff = 'q0')
sample1b <- RunSVD(sample1b)
sample1b <- RunUMAP(sample1b, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
```

Calculate WNN graph
```{r}
sample1b <- FindMultiModalNeighbors(sample1b, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample1b <- RunUMAP(sample1b, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample1b <- FindClusters(sample1b, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
```

-- CAN LABEL CELL IDENTS HERE --

```{r}
p1 <- DimPlot(sample1b, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample1b, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample1b, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
Sample 2 / Batch 2

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample2_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample2 <- CreateSeuratObject(counts = rna_counts)
sample2[["percent.mt"]] <- PercentageFeatureSet(sample1b, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/ATAC fragments/Batch2/bacth2_atac_fragments/Batch2_Sample2_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample2[["ATAC"]] <- chrom_assay

VlnPlot(sample2, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample2) <- "RNA"
sample2 <- SCTransform(sample2, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample2) <- "ATAC"
sample2 <- RunTFIDF(sample2)
sample2 <- FindTopFeatures(sample2, min.cutoff = 'q0')
sample2 <- RunSVD(sample2)
sample2 <- RunUMAP(sample2, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample2 <- FindMultiModalNeighbors(sample2, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample2 <- RunUMAP(sample2, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample2 <- FindClusters(sample2, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample2, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample2, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample2, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 

Sample 4 // Batch 2
```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample4_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample4 <- CreateSeuratObject(counts = rna_counts)
sample4[["percent.mt"]] <- PercentageFeatureSet(sample4, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/ATAC fragments/Batch2/bacth2_atac_fragments/Batch2_Sample4_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample4[["ATAC"]] <- chrom_assay

VlnPlot(sample4, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample4) <- "RNA"
sample4 <- SCTransform(sample4, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample4) <- "ATAC"
sample4 <- RunTFIDF(sample4)
sample4 <- FindTopFeatures(sample4, min.cutoff = 'q0')
sample4 <- RunSVD(sample4)
sample4 <- RunUMAP(sample4, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample4 <- FindMultiModalNeighbors(sample4, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample4 <- RunUMAP(sample4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample4 <- FindClusters(sample4, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample4, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample4, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample4, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 

Sample 6 // Batch 2
```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample6_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample6 <- CreateSeuratObject(counts = rna_counts)
sample6[["percent.mt"]] <- PercentageFeatureSet(sample6, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/ATAC fragments/Batch2/bacth2_atac_fragments/Batch2_Sample6_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample6[["ATAC"]] <- chrom_assay

VlnPlot(sample6, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample6) <- "RNA"
sample6 <- SCTransform(sample6, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample6) <- "ATAC"
sample6 <- RunTFIDF(sample6)
sample6 <- FindTopFeatures(sample6, min.cutoff = 'q0')
sample6 <- RunSVD(sample6)
sample6 <- RunUMAP(sample6, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample6 <- FindMultiModalNeighbors(sample6, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample6 <- RunUMAP(sample6, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample6 <- FindClusters(sample6, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample6, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample6, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample6, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
Sample 7 // Batch 2

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample7_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample7 <- CreateSeuratObject(counts = rna_counts)
sample7[["percent.mt"]] <- PercentageFeatureSet(sample7, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/ATAC fragments/Batch2/bacth2_atac_fragments/Batch2_Sample7_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample7[["ATAC"]] <- chrom_assay

VlnPlot(sample7, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample7) <- "RNA"
sample7 <- SCTransform(sample7, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample7) <- "ATAC"
sample7 <- RunTFIDF(sample7)
sample7 <- FindTopFeatures(sample7, min.cutoff = 'q0')
sample7 <- RunSVD(sample7)
sample7 <- RunUMAP(sample7, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample7 <- FindMultiModalNeighbors(sample7, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample7 <- RunUMAP(sample7, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample7 <- FindClusters(sample7, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample7, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample7, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample7, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
Sample 8 // Batch 2

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample8_filtered_feature_bc_matrix.h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample8 <- CreateSeuratObject(counts = rna_counts)
sample8[["percent.mt"]] <- PercentageFeatureSet(sample8, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/ATAC fragments/Batch2/bacth2_atac_fragments/Batch2_Sample8_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample8[["ATAC"]] <- chrom_assay

VlnPlot(sample8, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample8) <- "RNA"
sample8 <- SCTransform(sample8, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample8) <- "ATAC"
sample8 <- RunTFIDF(sample8)
sample8 <- FindTopFeatures(sample8, min.cutoff = 'q0')
sample8 <- RunSVD(sample8)
sample8 <- RunUMAP(sample8, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample8 <- FindMultiModalNeighbors(sample8, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample8 <- RunUMAP(sample8, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample8 <- FindClusters(sample8, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample8, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample8, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample8, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
BATCH #1 STARTING HERE 
CF10_P1

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-010_P1_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

CF10_P1 <- CreateSeuratObject(counts = rna_counts)
CF10_P1[["percent.mt"]] <- PercentageFeatureSet(CF10_P1, pattern = "^MT-")
CF10_P1 <- subset(CF10_P1, subset = percent.mt <= 20)

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-010_P1_atac_fragments.tsv"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
CF10_P1[["ATAC"]] <- chrom_assay

VlnPlot(CF10_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(CF10_P1) <- "RNA"
CF10_P1 <- SCTransform(CF10_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(CF10_P1) <- "ATAC"
CF10_P1 <- RunTFIDF(CF10_P1)
CF10_P1 <- FindTopFeatures(CF10_P1, min.cutoff = 'q0')
CF10_P1 <- RunSVD(CF10_P1)
CF10_P1 <- RunUMAP(CF10_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

CF10_P1 <- FindMultiModalNeighbors(CF10_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CF10_P1 <- RunUMAP(CF10_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CF10_P1 <- FindClusters(CF10_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(CF10_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF10_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF10_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```
```{r}
file.exists("/Users/asingh/Downloads/Batch1/atac files /CF-18-010_P1_atac_fragments.tsv")
file.exists(paste0("/Users/asingh/Downloads/Batch1/atac files /CF-18-010_P1_atac_fragments.tsv", '.tbi'))
```


--- NEXT SAMPLE --- 
CF13_P1 // this sample has complete pipeline of filtered as well as ScType annotation + cluster markers 

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

CF13_P1 <- CreateSeuratObject(counts = rna_counts)
CF13_P1[["percent.mt"]] <- PercentageFeatureSet(CF13_P1, pattern = "^MT-")

CF13_P1[[]]

all.genes <- rownames(CF13_P1)
CF13_P1 <- ScaleData(CF13_P1, features = all.genes)

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
CF13_P1[["ATAC"]] <- chrom_assay

VlnPlot(CF13_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

CF13_P1 <- subset(
  x = CF13_P1,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(CF13_P1) <- "RNA"
CF13_P1 <- SCTransform(CF13_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(CF13_P1) <- "ATAC"
CF13_P1 <- RunTFIDF(CF13_P1)
CF13_P1 <- FindTopFeatures(CF13_P1, min.cutoff = 'q0')
CF13_P1 <- RunSVD(CF13_P1)
CF13_P1 <- RunUMAP(CF13_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

CF13_P1 <- FindMultiModalNeighbors(CF13_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CF13_P1 <- RunUMAP(CF13_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CF13_P1 <- FindClusters(CF13_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)




lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = CF13_P1[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(CF13_P1@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(CF13_P1@meta.data[CF13_P1@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(CF13_P1@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

CF13_P1@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  CF13_P1@meta.data$customclassif[CF13_P1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}




p1 <- DimPlot(CF13_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF13_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF13_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p4 <- DimPlot(CF13_P1, reduction = "wnn.umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
p1 + p2 + p3 + p4 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

cluster0.markers <- FindMarkers(CF13_P1, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(CF13_P1, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(CF13_P1, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
cluster3.markers <- FindMarkers(CF13_P1, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 5)
```
Annotating using Azimuth
```{r}

#devtools::install_github("satijalab/seurat-data")
#devtools::install_github("satijalab/azimuth")
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
```
```{r}
install.packages("remotes")
remotes::install_github("satijalab/azimuth")1
```

```{r}
library(Seurat)
library(Azimuth)
library(SeuratData)

options(timeout = 1000)
InstallData("lungref.SeuratData")

CF13_P1 <- RunAzimuth(CF13_P1, reference = "lungref")

p1 <- DimPlot(CF13_P1, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p1
```

```{r}
CF13_P1 <- RunAzimuth(CF13_P1, reference = "lungref")

p1 <- DimPlot(CF13_P1, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p1
```

```{r}

```


```{r}
available_data <- AvailableData()
available_data[grep("Azimuth", available_data[, 3]), 1:3]


```

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
rna_counts <- inputdata.10x$`Gene Expression`
CF16_P0 <- CreateSeuratObject(counts = rna_counts)

CF16_P0 <- RunAzimuth(CF16_P0, reference = "lungref")

p1 <- DimPlot(CF16_P0, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p1
```

```{r}

p1 <- DimPlot(CF16_P0, group.by = "predicted.ann_level_4", label = TRUE, label.size = 3) + NoLegend()
p1
```
```{r}
CF16_P0[[]]
```

--- NEXT SAMPLE --- 
CF16_P0

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

CF16_P0 <- CreateSeuratObject(counts = rna_counts)
CF16_P0[["percent.mt"]] <- PercentageFeatureSet(CF16_P0, pattern = "^MT-")


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
CF16_P0[["ATAC"]] <- chrom_assay

VlnPlot(CF16_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
CF16_P0 <- subset(
  x = CF16_P0,
    percent.mt < 20
)

# RNA analysis
DefaultAssay(CF16_P0) <- "RNA"
CF16_P0 <- SCTransform(CF16_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(CF16_P0) <- "ATAC"
CF16_P0 <- RunTFIDF(CF16_P0)
CF16_P0 <- FindTopFeatures(CF16_P0, min.cutoff = 'q0')
CF16_P0 <- RunSVD(CF16_P0)
CF16_P0 <- RunUMAP(CF16_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

CF16_P0 <- FindMultiModalNeighbors(CF16_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
CF16_P0 <- RunUMAP(CF16_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CF16_P0 <- FindClusters(CF16_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(CF16_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(CF16_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(CF16_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
H14_P1

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H14_P1 <- CreateSeuratObject(counts = rna_counts)
H14_P1[["percent.mt"]] <- PercentageFeatureSet(H14_P1, pattern = "^MT-")



grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H14_P1[["ATAC"]] <- chrom_assay

VlnPlot(H14_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(H14_P1) <- "RNA"
H14_P1 <- SCTransform(H14_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H14_P1) <- "ATAC"
H14_P1 <- RunTFIDF(H14_P1)
H14_P1 <- FindTopFeatures(H14_P1, min.cutoff = 'q0')
H14_P1 <- RunSVD(H14_P1)
H14_P1 <- RunUMAP(H14_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H14_P1 <- FindMultiModalNeighbors(H14_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H14_P1 <- RunUMAP(H14_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H14_P1 <- FindClusters(H14_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H14_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H14_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H14_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
H15_P0

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H15_P0 <- CreateSeuratObject(counts = rna_counts)
H15_P0[["percent.mt"]] <- PercentageFeatureSet(H15_P0, pattern = "^MT-")
H15_P0 <- subset(H15_P0, subset = percent.mt <= 20)


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H15_P0[["ATAC"]] <- chrom_assay

VlnPlot(H15_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(H15_P0) <- "RNA"
H15_P0 <- SCTransform(H15_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H15_P0) <- "ATAC"
H15_P0 <- RunTFIDF(H15_P0)
H15_P0 <- FindTopFeatures(H15_P0, min.cutoff = 'q0')
H15_P0 <- RunSVD(H15_P0)
H15_P0 <- RunUMAP(H15_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H15_P0 <- FindMultiModalNeighbors(H15_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H15_P0 <- RunUMAP(H15_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H15_P0 <- FindClusters(H15_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H15_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H15_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H15_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
H15_P1

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H15_P1 <- CreateSeuratObject(counts = rna_counts)
H15_P1[["percent.mt"]] <- PercentageFeatureSet(H15_P1, pattern = "^MT-")
H15_P1 <- subset(H15_P1, subset = percent.mt <= 20)


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H15_P1[["ATAC"]] <- chrom_assay

VlnPlot(H15_P1, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(H15_P1) <- "RNA"
H15_P1 <- SCTransform(H15_P1, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H15_P1) <- "ATAC"
H15_P1 <- RunTFIDF(H15_P1)
H15_P1 <- FindTopFeatures(H15_P1, min.cutoff = 'q0')
H15_P1 <- RunSVD(H15_P1)
H15_P1 <- RunUMAP(H15_P1, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H15_P1 <- FindMultiModalNeighbors(H15_P1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H15_P1 <- RunUMAP(H15_P1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H15_P1 <- FindClusters(H15_P1, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H15_P1, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H15_P1, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H15_P1, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 
H18_P0

```{r}
inputdata.10x <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

H18_P0 <- CreateSeuratObject(counts = rna_counts)
H18_P0[["percent.mt"]] <- PercentageFeatureSet(H18_P0, pattern = "^MT-")
H18_P0 <- subset(H18_P0, subset = percent.mt <= 20)


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
H18_P0[["ATAC"]] <- chrom_assay

VlnPlot(H18_P0, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(H18_P0) <- "RNA"
H18_P0 <- SCTransform(H18_P0, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(H18_P0) <- "ATAC"
H18_P0 <- RunTFIDF(H18_P0)
H18_P0 <- FindTopFeatures(H18_P0, min.cutoff = 'q0')
H18_P0 <- RunSVD(H18_P0)
H18_P0 <- RunUMAP(H18_P0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

H18_P0 <- FindMultiModalNeighbors(H18_P0, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
H18_P0 <- RunUMAP(H18_P0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
H18_P0 <- FindClusters(H18_P0, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(H18_P0, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(H18_P0, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(H18_P0, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

--- NEXT SAMPLE --- 





```{r eval=FALSE, include=FALSE}
inputdata.10x <- Read10X_h5("")

rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

sample1b <- CreateSeuratObject(counts = rna_counts)
sample1b[["percent.mt"]] <- PercentageFeatureSet(sample1b, pattern = "^MT-")
H18_P0 <- subset(H18_P0, subset = percent.mt <= 20)


grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- ""
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )
sample1b[["ATAC"]] <- chrom_assay

VlnPlot(sample1b, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()

# RNA analysis
DefaultAssay(sample1b) <- "RNA"
sample1b <- SCTransform(sample1b, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(sample1b) <- "ATAC"
sample1b <- RunTFIDF(sample1b)
sample1b <- FindTopFeatures(sample1b, min.cutoff = 'q0')
sample1b <- RunSVD(sample1b)
sample1b <- RunUMAP(sample1b, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

sample1b <- FindMultiModalNeighbors(sample1b, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
sample1b <- RunUMAP(sample1b, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sample1b <- FindClusters(sample1b, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(sample1b, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(sample1b, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(sample1b, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

## To examine genes in a feature plot
<FeaturePlot(object = bothbatchesrnramy, features = c("KRT5", "KRT14", "KRT15", "TP63", "ITGA6"), max.cutoff = 5)>

## Getting differential markers comparing P5 to P1 
---
combining batch 1 & 2
seperate by P0/P1/P
---

```{r}
library(Seurat)
```

import all batch 1:

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "P1")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "P0")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

import all batch 2: 


```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "P1")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "P5")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "P1")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "P5")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "P5")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "P1")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```

merge both 
```{r}
bothbatches <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0, so2.2, so2.4, so2.6, so2.7, so2.8, so2.1), project = "CFfiltered")
```

```{r}
VlnPlot(bothbatches, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


run through seurat pipeline

```{r}
all.genes <- rownames(bothbatches)
bothbatches <- ScaleData(bothbatches, features = all.genes)
```

```{r}
bothbatches <- FindVariableFeatures(bothbatches, selection.method = "vst", nfeatures = 2000)
```

```{r}
bothbatches <- RunPCA(bothbatches, features = VariableFeatures(object = bothbatches))
```

```{r}
bothbatches <- FindNeighbors(bothbatches, dims = 1:10)
bothbatches <- FindClusters(bothbatches, resolution = 0.5)
```


```{r}
bothbatches <- RunUMAP(bothbatches, dims = 1:10)
DimPlot(bothbatches)
```

```{r}
DimPlot(bothbatches, reduction = "umap", group.by = "orig.ident")
```

```{r}
FeaturePlot(object = bothbatches, features = c("NR3C1", "FOXK1"), max.cutoff = 20)
```

```{r}
Idents(bothbatches) <- bothbatches@meta.data$orig.ident
diffmarkerstxt <- FindMarkers(bothbatches, ident.1 = "P5", ident.2 = c("P1"), min.pct = 0.25)
head(diffmarkerstxt)
```

```{r}
Idents(bothbatches) <- bothbatches@meta.data$orig.ident
diffmarkerstxt1 <- FindMarkers(bothbatches, ident.1 = "P5", ident.2 = "P1", min.pct = 0.25)
head(diffmarkerstxt1, n = 10)
```

```{r}
Idents(bothbatches) <- bothbatches@meta.data$orig.ident
diffmarkerstxt2 <- FindMarkers(bothbatches, ident.1 = "P5", ident.2 = "P1", min.pct = 0.25)
write.csv(diffmarkerstxt2, "diffmarkerstxt2")
```

## Patient sample sets (CF16)
---
title: "R Notebook"
output: html_notebook
---


```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Azimuth)
library(SeuratData)
library(patchwork)
```

```{r}
CFP0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")
CF16P1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
soCF16P0 <- CreateSeuratObject(counts = CFP0$`Gene Expression`, project = "P0")
soCF16P1 <- CreateSeuratObject(counts = CF16P1$`Gene Expression`, project = "P1")
merged_samples <- merge(soCF16P0, y = c(soCF16P1), project = "timepointCF16")
```

```{r}
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
all.genes <- rownames(merged_samples)
merged_samples <- ScaleData(merged_samples, features = all.genes)
merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)
merged_samples <- RunPCA(merged_samples, features = VariableFeatures(object = merged_samples))
ElbowPlot(merged_samples)
merged_samples <- FindNeighbors(merged_samples, dims = 1:10)
merged_samples <- FindClusters(merged_samples, resolution = 0.5)
merged_samples <- RunUMAP(merged_samples, dims = 1:10)
DimPlot(merged_samples)
DimPlot(merged_samples, reduction = "umap", group.by = "orig.ident")
```

```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = merged_samples[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(merged_samples@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(merged_samples@meta.data[merged_samples@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_samples@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

merged_samples@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  merged_samples@meta.data$customclassif[merged_samples@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}

DimPlot(merged_samples, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```

## Patient sample set (H18)
Time points by sample 

```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Azimuth)
library(SeuratData)
library(patchwork)
```

```{r}
H18P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")
H18P5 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
H18P1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
```

```{r}
soH18P0 <- CreateSeuratObject(counts = H18P0$`Gene Expression`, project = "P0")
soH18P5 <- CreateSeuratObject(counts = H18P5$`Gene Expression`, project = "P5")
soH18P1 <- CreateSeuratObject(counts = H18P1$`Gene Expression`, project = "P1")

merged_samples <- merge(soH18P0, y = c(soH18P5, soH18P1), project = "timepointH18")

```


```{r}
seurat <- function(merged_samples) {
  merged_samples[["percent.mt"]] <- PercentageFeatureSet(merged_samples, pattern = "^MT-")
  merged_samples <- subset(merged_samples, subset = percent.mt <= 20)
  VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  all.genes <- rownames(merged_samples)
  merged_samples <- ScaleData(merged_samples, features = all.genes)
  merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)
  merged_samples <- RunPCA(merged_samples, features = VariableFeatures(object = merged_samples))
  ElbowPlot(merged_samples)
  merged_samples <- FindNeighbors(merged_samples, dims = 1:10)
  merged_samples <- FindClusters(merged_samples, resolution = 0.5)
  merged_samples <- RunUMAP(merged_samples, dims = 1:10)
  DimPlot(merged_samples)
  DimPlot(merged_samples, reduction = "umap", group.by = "orig.ident")
  return(merged_samples)
}

seurat(merged_samples)
```

```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = merged_samples[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(merged_samples@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(merged_samples@meta.data[merged_samples@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_samples@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

merged_samples@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  merged_samples@meta.data$customclassif[merged_samples@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])

DimPlot(merged_samples, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```

```{r}
merged_samples[["percent.mt"]] <- PercentageFeatureSet(merged_samples, pattern = "^MT-")
VlnPlot(merged_samples, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
all.genes <- rownames(merged_samples)
merged_samples <- ScaleData(merged_samples, features = all.genes)
merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)
merged_samples <- RunPCA(merged_samples, features = VariableFeatures(object = merged_samples))
ElbowPlot(merged_samples)
merged_samples <- FindNeighbors(merged_samples, dims = 1:10)
merged_samples <- FindClusters(merged_samples, resolution = 0.5)
merged_samples <- RunUMAP(merged_samples, dims = 1:10)
DimPlot(merged_samples)
DimPlot(merged_samples, reduction = "umap", group.by = "orig.ident")
```
```{r}
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Lung" 

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = merged_samples[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(merged_samples@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(merged_samples@meta.data[merged_samples@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(merged_samples@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

merged_samples@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  merged_samples@meta.data$customclassif[merged_samples@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])}

DimPlot(merged_samples, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')  
```
```{r}
FindMarkers()
```

## Filtered, merged ATAC and WNN
---
title: "R Notebook"
output: html_notebook
---

load libraries
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```


```{r}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
inputdata.CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H14_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")
inputdata.H15_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H18_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")


frag.file160 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-016_P0_atac_fragments.tsv.gz"
frag.file131 <- "/Users/asingh/Downloads/Batch1/atac files /CF-18-013_P1_atac_fragments.tsv.gz"
frag.file141 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-014-P1_atac_fragments.tsv.gz"
frag.file151 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-015_P1_atac_fragments.tsv.gz"
frag.file180 <- "/Users/asingh/Downloads/Batch1/atac files /HEALTHY-18-018_P0_atac_fragments.tsv.gz"

chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)
chromassay131 <- CreateChromatinAssay(counts = inputdata.CF13_P1$Peaks, sep = c(":", "-"), fragments = frag.file131, min.cells = 1, min.features = -1)
chromassay141 <- CreateChromatinAssay(counts = inputdata.H14_P1$Peaks, sep = c(":", "-"), fragments = frag.file141, min.cells = 1, min.features = -1)
chromassay151 <- CreateChromatinAssay(counts = inputdata.H15_P1$Peaks, sep = c(":", "-"), fragments = frag.file151, min.cells = 1, min.features = -1)
chromassay180 <- CreateChromatinAssay(counts = inputdata.H18_P0$Peaks, sep = c(":", "-"), fragments = frag.file180, min.cells = 1, min.features = -1)

#CHECK MERGE
#chromassay <- merge(
 # x = chromassay160,
  #y = list(chromassay131, chromassay141, chromassay151, chromassay180)
#)

object160 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`)
object160 <- RunAzimuth(object160, reference = "lungref")
object160[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)
object160[["percent.mt"]] <- PercentageFeatureSet(object160, pattern = "MT")
object160 <- subset(x = object160, percent.mt < 20 )

object131 <- CreateSeuratObject(counts = inputdata.CF13_P1$`Gene Expression`)
object131 <- RunAzimuth(object131, reference = "lungref")
object131[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF13_P1$`Peaks`)
object131[["percent.mt"]] <- PercentageFeatureSet(object131, pattern = "MT")
object131 <- subset(x = object131, percent.mt < 20 )

object141 <- CreateSeuratObject(counts = inputdata.H14_P1$`Gene Expression`)
object141 <- RunAzimuth(object141, reference = "lungref")
object141[["ATAC"]] <- CreateAssayObject(counts = inputdata.H14_P1$`Peaks`)
object141[["percent.mt"]] <- PercentageFeatureSet(object141, pattern = "MT")
object141 <- subset(x = object141, percent.mt < 20 )

object151 <- CreateSeuratObject(counts = inputdata.H15_P1$`Gene Expression`)
object151 <- RunAzimuth(object151, reference = "lungref")
object151[["ATAC"]] <- CreateAssayObject(counts = inputdata.H15_P1$`Peaks`)
object151[["percent.mt"]] <- PercentageFeatureSet(object151, pattern = "MT")
object151 <- subset(x = object151, percent.mt < 20 )

#object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
#object5[["RNA"]] < CreateAssayObject(counts = inputdata.H18_P0$`Gene Expression`)

#mergedobject <- merge(object, y = c(object2, object141, object151, object5), project = "atacmerge_b1")
```

error: Cannot find RNA in this Seurat object
```{r eval=FALSE, include=FALSE}
object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
inputdata.H18_P0$`Gene Expression`
```

```{r}
mergedobject <- merge(object160, y = c(object131, object141, object151), project = "atacmerge_b1")
```

```{r}
# RNA analysis
DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)
```

```{r}
mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

```{r}
p4 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p5 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p6 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)

p4 + p5 + p6
```
```{r}
mergedobject[[]]
```

---
ATAC MERGED ANNOTATED FILTERED FULL BATCH2
---
load libraries
```{r}
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

library(Azimuth)
library(SeuratData)
library(patchwork)
```

```{r}
rm(mergedNBG)
rm(mergedNBG.basal)
rm(pbmc)
rm(pbmc.data)
rm(plot2)
rm(plot1)
```

```{r}
rm(mergedNBG)
rm(mergedNBG.basal)
rm(pbmc)
rm(pbmc.data)
rm(plot2)
rm(plot1)
```



THIS:
```{r}
inputdata.CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-016_P0_filtered_feature_bc_matrix (1).h5")
inputdata.CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/CF-18-013_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H14_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-014-P1_filtered_feature_bc_matrix (2).h5")
inputdata.H15_P1 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-015_P1_filtered_feature_bc_matrix (1).h5")
inputdata.H18_P0 <- Read10X_h5("/Users/asingh/Downloads/Batch1/h5 files/HEALTHY-18-018_P0_filtered_feature_bc_matrix (1).h5")
inputdata.1 <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample1_filtered_feature_bc_matrix.h5")
#inputdata.2 <- Read10X_h5("/Users/asingh/Downloads/ATAC #fragments/Batch2/Batch2_h5_files/Sample2_filtered_feature_bc_matrix.h5")
inputdata.4 <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample4_filtered_feature_bc_matrix.h5")
inputdata.6 <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample6_filtered_feature_bc_matrix.h5")
inputdata.7 <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample7_filtered_feature_bc_matrix.h5")
inputdata.8 <- Read10X_h5("/Users/asingh/Downloads/ATAC fragments/Batch2/Batch2_h5_files/Sample8_filtered_feature_bc_matrix.h5")


#frag.file1 <- ""
#frag.file2 <- ""
#frag.file4 <- ""
#frag.file6 <- ""
#frag.file7 <- ""
#frag.file8 <- ""

#chromassay160 <- CreateChromatinAssay(counts = inputdata.CF16_P0$Peaks, sep = c(":", "-"), fragments = frag.file160, min.cells = 1, min.features = -1)
#chromassay131 <- CreateChromatinAssay(counts = inputdata.CF13_P1$Peaks, sep = c(":", "-"), fragments = frag.file131, min.cells = 1, min.features = -1)
#chromassay141 <- CreateChromatinAssay(counts = inputdata.H14_P1$Peaks, sep = c(":", "-"), fragments = frag.file141, min.cells = 1, min.features = -1)
#chromassay151 <- CreateChromatinAssay(counts = inputdata.H15_P1$Peaks, sep = c(":", "-"), fragments = frag.file151, min.cells = 1, min.features = -1)
#chromassay180 <- CreateChromatinAssay(counts = inputdata.H18_P0$Peaks, sep = c(":", "-"), fragments = frag.file180, min.cells = 1, min.features = -1)

#CHECK MERGE
#chromassay <- merge(
 # x = chromassay160,
  #y = list(chromassay131, chromassay141, chromassay151, chromassay180)
#)

object160 <- CreateSeuratObject(counts = inputdata.CF16_P0$`Gene Expression`, project = "P0")
object160 <- RunAzimuth(object160, reference = "lungref")
object160[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF16_P0$Peaks)
object160[["percent.mt"]] <- PercentageFeatureSet(object160, pattern = "MT")
object160 <- subset(x = object160, percent.mt < 20 )

object131 <- CreateSeuratObject(counts = inputdata.CF13_P1$`Gene Expression`, project = "P1")
object131 <- RunAzimuth(object131, reference = "lungref")
object131[["ATAC"]] <- CreateAssayObject(counts = inputdata.CF13_P1$`Peaks`)
object131[["percent.mt"]] <- PercentageFeatureSet(object131, pattern = "MT")
object131 <- subset(x = object131, percent.mt < 20 )

object141 <- CreateSeuratObject(counts = inputdata.H14_P1$`Gene Expression`, project = "P1")
object141 <- RunAzimuth(object141, reference = "lungref")
object141[["ATAC"]] <- CreateAssayObject(counts = inputdata.H14_P1$`Peaks`)
object141[["percent.mt"]] <- PercentageFeatureSet(object141, pattern = "MT")
object141 <- subset(x = object141, percent.mt < 20 )

object151 <- CreateSeuratObject(counts = inputdata.H15_P1$`Gene Expression`, project = "P1")
object151 <- RunAzimuth(object151, reference = "lungref")
object151[["ATAC"]] <- CreateAssayObject(counts = inputdata.H15_P1$`Peaks`)
object151[["percent.mt"]] <- PercentageFeatureSet(object151, pattern = "MT")
object151 <- subset(x = object151, percent.mt < 20 )

object1 <- CreateSeuratObject(counts = inputdata.1$`Gene Expression`, project = "P1")
object1 <- RunAzimuth(object1, reference = "lungref")
object1[["ATAC"]] <- CreateAssayObject(counts = inputdata.1$Peaks)
object1[["percent.mt"]] <- PercentageFeatureSet(object1, pattern = "MT")
object1 <- subset(x = object1, percent.mt < 20 )

#object2 <- CreateSeuratObject(counts = inputdata.2$`Gene Expression`)
#object2 <- RunAzimuth(object2, reference = "lungref")
#object2[["ATAC"]] <- CreateAssayObject(counts = inputdata.2$`Peaks`)
#object2[["percent.mt"]] <- PercentageFeatureSet(object2, pattern = "MT")
#object2 <- subset(x = object2, percent.mt < 20 )

object4 <- CreateSeuratObject(counts = inputdata.4$`Gene Expression`, project = "P1")
object4 <- RunAzimuth(object4, reference = "lungref")
object4[["ATAC"]] <- CreateAssayObject(counts = inputdata.4$`Peaks`)
object4[["percent.mt"]] <- PercentageFeatureSet(object4, pattern = "MT")
object4 <- subset(x = object4, percent.mt < 20 )

object6 <- CreateSeuratObject(counts = inputdata.6$`Gene Expression`, project = "P5")
object6 <- RunAzimuth(object6, reference = "lungref")
object6[["ATAC"]] <- CreateAssayObject(counts = inputdata.6$`Peaks`)
object6[["percent.mt"]] <- PercentageFeatureSet(object6, pattern = "MT")
object6 <- subset(x = object6, percent.mt < 20 )

object7 <- CreateSeuratObject(counts = inputdata.7$`Gene Expression`, project = "P5")
object7 <- RunAzimuth(object7, reference = "lungref")
object7[["ATAC"]] <- CreateAssayObject(counts = inputdata.7$`Peaks`)
object7[["percent.mt"]] <- PercentageFeatureSet(object7, pattern = "MT")
object7 <- subset(x = object7, percent.mt < 20 )

object8 <- CreateSeuratObject(counts = inputdata.8$`Gene Expression`, project = "P1")
object8 <- RunAzimuth(object8, reference = "lungref")
object8[["ATAC"]] <- CreateAssayObject(counts = inputdata.8$`Peaks`)
object8[["percent.mt"]] <- PercentageFeatureSet(object8, pattern = "MT")
object8 <- subset(x = object8, percent.mt < 20 )

rm(inputdata.CF13_P1)
rm(inputdata.8)
rm(inputdata.7)
rm(inputdata.6)
rm(inputdata.4)
rm(inputdata.1)
rm(inputdata.H18_P0)
rm(inputdata.H15_P1)
rm(inputdata.H14_P1)
rm(inputdata.CF16_P0)

mergedobject <- merge(object1, y = c(object4, object6, object7, object8, object160, object151, object141, object131), project = "atacmerge_oct31")

#object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
#object5[["RNA"]] < CreateAssayObject(counts = inputdata.H18_P0$`Gene Expression`)

#mergedobject <- merge(object, y = c(object2, object4, object6, object5), project = "atacmerge_b1")



DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)

mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```
```{r}
p1 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 2.5,  repel = TRUE) + ggtitle("WNN")

p1  
p2 
p3 
```
```{r}
p1 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p1
p2 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5,  repel = TRUE) + ggtitle("WNN")
p2
p3
```


error: Cannot find RNA in this Seurat object
```{r eval=FALSE, include=FALSE}
object5 <- CreateSeuratObject(chromassay180, assay = "ATAC")
inputdata.H18_P0$`Gene Expression`
```

```{r}
mergedobject <- merge(object1, y = c(object2, object4, object6, object7, object8), project = "atacmerge_b2")
```

```{r}
# RNA analysis
DefaultAssay(mergedobject) <- "RNA"
mergedobject <- SCTransform(mergedobject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mergedobject) <- "ATAC"
mergedobject <- RunTFIDF(mergedobject)
mergedobject <- FindTopFeatures(mergedobject, min.cutoff = 'q0')
mergedobject <- RunSVD(mergedobject)
```

```{r}
mergedobject <- RunUMAP(mergedobject, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mergedobject <- FindMultiModalNeighbors(mergedobject, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mergedobject <- RunUMAP(mergedobject, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mergedobject <- FindClusters(mergedobject, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mergedobject, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mergedobject, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mergedobject, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
```

```{r}
p4 <- DimPlot(mergedobject, reduction = "umap.rna", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p5 <- DimPlot(mergedobject, reduction = "umap.atac", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)
p6 <- DimPlot(mergedobject, reduction = "wnn.umap", group.by = "predicted.ann_level_4", label = TRUE, label.size = 3)

p4 + p5 + p6
```
## Getting differential markers comparing P0 to P1
---
title: "R Notebook"
output: html_notebook
---

```{r}
library(Seurat)
```

import all batch 1:

```{r}
CF13_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-013_P1_filtered_feature_bc_matrix_1.h5", use.names = TRUE, unique.features = TRUE)
```

create seurat object 
```{r}
SO_CF13_P1 <- CreateSeuratObject(counts = CF13_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF13_P1
```

```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF13_P1, pattern = "^MT-")
```

```{r}
SO_CF13_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF13_P1, 
                                    pattern = "^RP[SL]")

```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1[[]]
```


```{r}
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.ribo <= 60)
SO_CF13_P1 <- subset(SO_CF13_P1, subset = percent.mt <= 20)
```


```{r}
SO_CF13_P1[[]]
```

```{r}
VlnPlot(SO_CF13_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF13_P1[[]]
```

```{r}
CF10_P1 <- Read10X_h5("/Users/asingh/Downloads/CF-18-010_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF10_P1 <- CreateSeuratObject(counts = CF10_P1$`Gene Expression`, project = "CF")
```

```{r}
SO_CF10_P1[["percent.mt"]] <- PercentageFeatureSet(SO_CF10_P1, pattern = "^MT-")
VlnPlot(SO_CF10_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
SO_CF10_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_CF10_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF10_P1[[]]
```


```{r}
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.ribo <= 60)
SO_CF10_P1 <- subset(SO_CF10_P1, subset = percent.mt <= 20)
```

```{r}
SO_CF10_P1[[]]
```


```{r}
CF16_P0 <- Read10X_h5("/Users/asingh/Downloads/CF-18-016_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_CF16_P0 <- CreateSeuratObject(counts = CF16_P0$`Gene Expression`, project = "CF")
```

```{r}
SO_CF16_P0[["percent.mt"]] <- PercentageFeatureSet(SO_CF16_P0, pattern = "^MT-")
VlnPlot(SO_CF16_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_CF16_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_CF16_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_CF16_P0[[]]
```


```{r}
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.ribo <= 60)
SO_CF16_P0 <- subset(SO_CF16_P0, subset = percent.mt <= 20)
SO_CF16_P0[[]]
```

```{r}
H14_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-014-P1_filtered_feature_bc_matrix (1).h5")
```

```{r}
SO_H14_P1 <- CreateSeuratObject(counts = H14_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H14_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H14_P1, pattern = "^MT-")
VlnPlot(SO_H14_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H14_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H14_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H14_P1[[]]
```

```{r}
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.ribo <= 60)
SO_H14_P1 <- subset(SO_H14_P1, subset = percent.mt <= 20)
SO_H14_P1[[]]
```


```{r}
H15_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P0_filtered_feature_bc_matrix.h5")
```


```{r}
SO_H15_P0 <- CreateSeuratObject(counts = H15_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P0, pattern = "^MT-")
VlnPlot(SO_H15_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P0[[]]
```

```{r}
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.ribo <= 60)
SO_H15_P0 <- subset(SO_H15_P0, subset = percent.mt <= 20)
SO_H15_P0[[]]
```


```{r}
H15_P1 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-015_P1_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H15_P1 <- CreateSeuratObject(counts = H15_P1$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H15_P1[["percent.mt"]] <- PercentageFeatureSet(SO_H15_P1, pattern = "^MT-")
VlnPlot(SO_H15_P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H15_P1[["percent.ribo"]] <- PercentageFeatureSet(SO_H15_P1, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H15_P1[[]]
```

```{r}
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.ribo <= 60)
SO_H15_P1 <- subset(SO_H15_P1, subset = percent.mt <= 20)
SO_H15_P1[[]]
```


```{r}
H18_P0 <- Read10X_h5("/Users/asingh/Downloads/HEALTHY-18-018_P0_filtered_feature_bc_matrix.h5")
```

```{r}
SO_H18_P0 <- CreateSeuratObject(counts = H18_P0$`Gene Expression`, project = "Healthy")
```

```{r}
SO_H18_P0[["percent.mt"]] <- PercentageFeatureSet(SO_H18_P0, pattern = "^MT-")
VlnPlot(SO_H18_P0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
SO_H18_P0[["percent.ribo"]] <- PercentageFeatureSet(SO_H18_P0, 
                                    pattern = "^RP[SL]")
```

```{r}
SO_H18_P0[[]]
```


```{r}
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.ribo <= 60)
SO_H18_P0 <- subset(SO_H18_P0, subset = percent.mt <= 20)
SO_H18_P0[[]]
```

import all batch 2: 


```{r}
batch2_1 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
so2.1 <- CreateSeuratObject(counts = batch2_1$`Gene Expression`, project = "Healthy")
```

```{r}
so2.1
```

```{r}
so2.1[[]]
```


```{r}
so2.1[["percent.mt"]] <- PercentageFeatureSet(so2.1, pattern = "^MT-")
```

```{r}
so2.1[["percent.ribo"]] <- PercentageFeatureSet(so2.1, 
                                    pattern = "^RP[SL]")

```


```{r}
so2.1[[]]
```


```{r}
so2.1[[]]
```


```{r}
so2.1 <- subset(so2.1, subset = percent.ribo <= 60)
so2.1 <- subset(so2.1, subset = percent.mt <= 20)
```


```{r}
so2.1[[]]
```

```{r}
VlnPlot(so2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.1[[]]
```

```{r}
batch2_2 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample2_filtered_feature_bc_matrix.h5")
so2.2 <- CreateSeuratObject(counts = batch2_2$`Gene Expression`, project = "Healthy")
```

```{r}
so2.2[["percent.mt"]] <- PercentageFeatureSet(so2.2, pattern = "^MT-")
VlnPlot(so2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

```

```{r}
so2.2[["percent.ribo"]] <- PercentageFeatureSet(so2.2, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.2[[]]
```


```{r}
so2.2 <- subset(so2.2, subset = percent.ribo <= 60)
so2.2 <- subset(so2.2, subset = percent.mt <= 20)
```

```{r}
so2.2[[]]
```


```{r}
batch2_4 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample4_filtered_feature_bc_matrix.h5")
so2.4 <- CreateSeuratObject(counts = batch2_4$`Gene Expression`, project = "CF")
```

```{r}
so2.4[["percent.mt"]] <- PercentageFeatureSet(so2.4, pattern = "^MT-")
VlnPlot(so2.4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.4[["percent.ribo"]] <- PercentageFeatureSet(so2.4, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.4[[]]
```


```{r}
so2.4 <- subset(so2.4, subset = percent.ribo <= 60)
so2.4 <- subset(so2.4, subset = percent.mt <= 20)
so2.4[[]]
```

```{r}
batch2_6 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample6_filtered_feature_bc_matrix.h5")
so2.6 <- CreateSeuratObject(counts = batch2_6$`Gene Expression`, project = "CF")
```

```{r}
so2.6[["percent.mt"]] <- PercentageFeatureSet(so2.6, pattern = "^MT-")
VlnPlot(so2.6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.6[["percent.ribo"]] <- PercentageFeatureSet(so2.6, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.6[[]]
```

```{r}
so2.6 <- subset(so2.6, subset = percent.ribo <= 60)
so2.6 <- subset(so2.6, subset = percent.mt <= 20)
so2.6[[]]
```


```{r}
batch2_7 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample7_filtered_feature_bc_matrix.h5")
so2.7 <- CreateSeuratObject(counts = batch2_7$`Gene Expression`, project = "Healthy")
```

```{r}
so2.7[["percent.mt"]] <- PercentageFeatureSet(so2.7, pattern = "^MT-")
VlnPlot(so2.7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.7[["percent.ribo"]] <- PercentageFeatureSet(so2.7, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.7[[]]
```

```{r}
so2.7 <- subset(so2.7, subset = percent.ribo <= 60)
so2.7 <- subset(so2.7, subset = percent.mt <= 20)
so2.7[[]]
```


```{r}
batch2_8 <- Read10X_h5("/Users/asingh/Downloads/drive-download-20220920T213314Z-001/Sample8_filtered_feature_bc_matrix.h5")
so2.8 <- CreateSeuratObject(counts = batch2_8$`Gene Expression`, project = "Healthy")
```

```{r}
so2.8[["percent.mt"]] <- PercentageFeatureSet(so2.8, pattern = "^MT-")
VlnPlot(so2.8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
so2.8[["percent.ribo"]] <- PercentageFeatureSet(so2.8, 
                                    pattern = "^RP[SL]")
```

```{r}
so2.8[[]]
```

```{r}
so2.8 <- subset(so2.8, subset = percent.ribo <= 60)
so2.8 <- subset(so2.8, subset = percent.mt <= 20)
so2.8[[]]
```

merge both 
```{r}
bothbatches <- merge(SO_CF13_P1, y = c(SO_CF10_P1, SO_CF16_P0, SO_H14_P1, SO_H15_P0, SO_H15_P1, SO_H18_P0, so2.2, so2.4, so2.6, so2.7, so2.8, so2.1), project = "CFfiltered")
```

```{r}
VlnPlot(bothbatches, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```


run through seurat pipeline

```{r}
all.genes <- rownames(bothbatches)
bothbatches <- ScaleData(bothbatches, features = all.genes)
```

```{r}
bothbatches <- FindVariableFeatures(bothbatches, selection.method = "vst", nfeatures = 2000)
```

```{r}
bothbatches <- RunPCA(bothbatches, features = VariableFeatures(object = bothbatches))
```

```{r}
bothbatches <- FindNeighbors(bothbatches, dims = 1:10)
bothbatches <- FindClusters(bothbatches, resolution = 0.5)
```


```{r}
bothbatches <- RunUMAP(bothbatches, dims = 1:10)
DimPlot(bothbatches)
```

```{r}
DimPlot(bothbatches, reduction = "umap", group.by = "orig.ident")
```

```{r}
FeaturePlot(object = bothbatches, features = c("NR3C1", "FOXK1"), max.cutoff = 20)
```

```{r}
#Idents(bothbatches) <- bothbatches@meta.data$orig.ident
#diffmarkerstxt <- FindMarkers(bothbatches, ident.1 = "P5", ident.2 = c("P1"), min.pct = 0.25)
#head(diffmarkerstxt)
```

```{r}
#Idents(bothbatches) <- bothbatches@meta.data$orig.ident

#diffmarkerstxt1 <- FindMarkers(bothbatches, ident.1 = "P5", ident.2 = "P1", min.pct = 0.25)
#head(diffmarkerstxt1, n = 10)
```

```{r}
Idents(bothbatches) <- bothbatches@meta.data$orig.ident
diffmarkerstxt2 <- FindMarkers(bothbatches, ident.1 = "CF", ident.2 = "Healthy", min.pct = 0.25)
write.csv(diffmarkerstxt2, "diffmarkerstxt3")
```

```{r}
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
#library(loomR)
library(Seurat)
library(patchwork)
```

```{r}
BiocManager::install("RColorBrewer")
```

```{r}
bb.sce <- as.SingleCellExperiment(bothbatches)
```

```{r}
library(RColorBrewer)
```

```{r}
geneFilter <- apply(assays(bb.sce)$counts,1,function(x){
    sum(x >= 3) >= 10
})
bb.sce <- bb.sce[geneFilter, ]
```
```{r}
rm(batch2_1)
rm(batch2_7)
rm(batch2_8)
```





