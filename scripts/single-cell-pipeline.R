library(Seurat)
library(patchwork)
library(monocle)

#loading dataset
e2.data <- Read10X(data.dir = "C:/Users/NIKSHEP/OneDrive/Documents/human/E2/")

#Creating Seurat object
e2<-CreateSeuratObject(e2.data,project ="GSM5293866",min.cells =3,min.features = 500)

#Adding Mitochondrial contamination percentage
e2[["percent.mt"]]<-PercentageFeatureSet(e2,pattern = "^MT-")

#plotting to see viable subsets
e2vln <- e2
VlnPlot(e2vln,features=c("nCount_RNA","nFeature_RNA","percent.mt"),ncol = 3)

plot1 <- FeatureScatter(e2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(e2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#subset
e2 <- subset(e2, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)

#normalizing
e2<-NormalizeData(e2)

#finding variable features
e2<-FindVariableFeatures(e2)
plot1<-VariableFeaturePlot(e2)
top10<-head(VariableFeatures(e2),10)
plot2<-LabelPoints(plot=plot1,points=top10,repel=TRUE)
plot1 + plot2

#scaling genes

all.genes<-rownames(e2)
e2<-ScaleData(e2, features= all.genes)

#PCA
e2<-RunPCA(e2,features=VariableFeatures(object=e2))

ElbowPlot(e2)

#UMAP
e2<-FindNeighbors(e2,dims = 1:15)
e2<-FindClusters(e2,resolution=0.4)
e2<-RunUMAP(e2,dims=1:15)

DimPlot(e2,reduction = "umap",label = TRUE)

#finding markers for each cluster
cluster2.markers <- FindMarkers(n2, ident.1 = 2, min.pct = 0.25, test.use ="bimod")
rownames(cluster2.markers)
length(cluster2.markers[ , 1])

#Feature plot for specific gene
FeaturePlot(e2,features=c("CCL4","HES4"))


#sample N2 ------------
n2.data <- Read10X(data.dir = "C:/Users/NIKSHEP/OneDrive/Documents/human/N2")
n2 <- CreateSeuratObject(n2.data, project = "sample-n2", min.cells = 3, min.features = 500)

#
n2[["percent.mt"]] <- PercentageFeatureSet(n2, pattern = "^MT-")

#
n2vln <- n2
VlnPlot(n2vln, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3)

#
n2 <- subset(n2, subset = nFeature_RNA<6000 & nFeature_RNA>500 & percent.mt<10)

#
n2 <- NormalizeData(n2)

#
n2 <- FindVariableFeatures(n2, nfeatures = 2000)
top10 <- head(VariableFeatures(n2), n = 10)
plot <- LabelPoints(plot = VariableFeaturePlot(n2), points = top10, repel = TRUE)
plot

#
all.genes <- rownames(n2)
n2 <- ScaleData(n2, features = all.genes)

#
n2 <- RunPCA(n2, features = VariableFeatures(n2))
DimPlot(n2, reduction = "pca")

#
ElbowPlot(n2)

#
n2 <- FindNeighbors(n2)
n2 <- FindClusters(n2, resolution = 0.2)

#
n2 <- RunUMAP(n2, dims = 1:10)
DimPlot(n2, reduction = "umap", label = TRUE)

#
a1.markers.0 <- FindMarkers(n2, ident.1 = 0, test.use = "bimod")
a1.markers.1 <- FindMarkers(n2, ident.1 = 1, test.use = "bimod")
a1.markers.2 <- FindMarkers(n2, ident.1 = 2, test.use = "bimod")
a1.markers.3 <- FindMarkers(n2, ident.1 = 3, test.use = "bimod")

#create seurat object
a1.data <- Read10X(data.dir = "C:/Users/NIKSHEP/OneDrive/Documents/human/A1/")
a1 <- CreateSeuratObject(counts = a1.data, project = "sample-a1", min.cells = 3, min.features = 500)

#separate mitochondrial genes
a1[["percent.mt"]] <- PercentageFeatureSet(a1, pattern = "^MT-")

#violin plot stuff (especially mt genes -> cutoff decided = 20)
a1vln <- a1
VlnPlot(a1vln, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
plot1 <- FeatureScatter(a1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(a1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(a1, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot1+plot2

#filter out cells
a1 <- subset(a1, subset = nFeature_RNA>500 & nFeature_RNA<6000 & percent.mt<40)

#Normalize
a1 <- NormalizeData(a1)

#feature selection
a1 <- FindVariableFeatures(a1, nfeatures = 3000)
top10 <- head(VariableFeatures(a1), n=10)
plot1 <- LabelPoints(plot = VariableFeaturePlot(a1), points = top10, repel = TRUE)
plot1

#scaling the data
all.genes <- rownames(a1)
a1 <- ScaleData(a1, features = all.genes)

#PCA reduction
a1 <- RunPCA(a1, features = VariableFeatures(a1))
print(a1[["pca"]], dims=1:10, nfeatures=10)
DimPlot(a1, reduction = "pca")

#elbow plots
ElbowPlot(a1)

#cluster after PCA
a1 <- FindNeighbors(a1, dims = 1:11)
a1 <- FindClusters(a1, resolution = 0.5)
head(Idents(a1))

#non linear dimensional reduction
a1 <- RunUMAP(a1, dims = 1:11)
DimPlot(a1, reduction = "umap", label = TRUE)

#biomarkers
a1.markers <- FindAllMarkers(a1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.markers <- a1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

n2.markers <- FindAllMarkers(n2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.markers.n2 <- n2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

e2.markers <- FindAllMarkers(e2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.markers.e2 <-e2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

markers <- c(top2.markers$gene, top2.markers.n2$gene)
markers <- c(markers, top2.markers.e2$gene)
markers

feature_genes <- top2.markers$gene
feature_genes

FeaturePlot(a1, features = c("C1QB", "C1QA"))

#
FeaturePlot(n2, features = c("CCL4","CCL3"))

umap <- function(a){
  DimPlot(a, reduction = "umap", label = TRUE)
}

feature_plot <- function(sample, gene){
  FeaturePlot(sample, features = gene)
}

singlecell_plot <- function(d, i, k){
  if(i == 1){
    scatter_plot(k)
  }
  if(i == 2){
    ElbowPlot(d)
  }
}

scatter_a <- FindVariableFeatures(a1, nfeatures = 2000)
top10 <- head(VariableFeatures(scatter_a), n = 10)

scatter_e <- FindVariableFeatures(placenta.E2, nfeatures = 2000)
top10 <- head(VariableFeatures(scatter_e), n = 10)

scatter_n <- FindVariableFeatures(n2, nfeatures = 2000)
top10 <- head(VariableFeatures(scatter_n), n = 10)

scatter_plot(scatter_e)

plot_a

scatter_plot <- function(k, i){
  plot <- LabelPoints(plot = VariableFeaturePlot(k), points = top10, repel = TRUE)
  plot
}

violin_plot <- function(a){
  VlnPlot(a,features=c("nCount_RNA","nFeature_RNA","percent.mt"),ncol = 3)
}

violin_text <- function(k){
  if(k == "Control"){"The cutoffs used are 6000, 500 and 10"}
  else if(k == "Early Pregnancy Loss"){"The cutoffs used are 6000, 500 and 40"}
  else if(k == "Endometriosis"){"The cutoffs used are 6000, 500 and 15"}
}

violin_final <- function(a, gene) {
  VlnPlot(a, features = gene)

}

cellData <- as.CellDataSet(a1)

cellData <- estimateSizeFactors(cellData)
cellData <- estimateDispersions(cellData)

cellData <- monocle::detectGenes(cellData, min_expr = 0.1)
print(head(fData(cellData)))

expressed_genes <- row.names(subset(fData(cellData),
                                    num_cells_expressed >= 10))
print(head(pData(cellData)))
