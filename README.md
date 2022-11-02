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
 
 
 ## Practice Tutorial using Seurat's vignette
 # insert all the libraries
library(dplyr)
library(Seurat)
library(patchwork) 

# load in the data set 
pbmc.data <- Read10X(data.dir = "/Users/asingh/Downloads/filtered_gene_bc_matrices/hg19")

# initialize package with raw data 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc  

# visualize into violin plot 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# visualize into a feature scatter plot 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# normalize the data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# is also the same as: pbmc <- NormalizeData(pbmc)

# want to calculate features that show high variation 
# to do this, use the FindVariableFeatures() function

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# to identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot1 + plot2

# to scale the data, use ScaleData function 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# perform PCA on scaled data
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# there are many ways to visualize this: VizDimReduction(), DimPlot(), DimHeatMap()
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# determine dimensionality 
# this getss ride of technical noise 

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

# use JackStrawPlot for p values 
JackStrawPlot(pbmc, dims = 1:15)  # takes a long time!

# can use elbow plot 
ElbowPlot(pbmc)

# cluster the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#UMAP: look at individual clusters 
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# save plot using saveRDS

# cluster biomarkers 
# ex: find markers of cluste 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

#look at marker expression 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# look at raw data
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()  # did not work 

# assign cell type to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
