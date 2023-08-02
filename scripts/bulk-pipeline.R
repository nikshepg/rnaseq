# loading packages --------
library(here)
library(tidyverse)
library(limma)
library(edgeR)
library(Homo.sapiens)
library(MASS)
library(reshape2)
library(reshape)
library(RColorBrewer)
library(ggVennDiagram)
library(gplots)


# extracting data from files--------
sample1 <- read.csv(here("data", "readCounts-sample-1.txt"), sep = '\t', header = FALSE)
sample1 <- dplyr::select(sample1, c(V1, V7))
colnames(sample1) <- c("EntrezID", "Sample1")
summary(sample1)

sample2 <- read.csv(here("data", "readCounts-sample-2.txt"), sep = '\t', header = FALSE)
sample2 <- dplyr::select(sample2, c(V1, V7))
colnames(sample2) <- c("EntrezID", "Sample2")
summary(sample2)

sample3 <- read.csv(here("data", "readCounts-sample-3.txt"), sep = '\t', header = FALSE)
sample3 <- dplyr::select(sample3, c(V1, V7))
colnames(sample3) <- c("EntrezID", "Sample3")
summary(sample3)

sample4 <- read.csv(here("data", "readCounts-sample-4.txt"), sep = '\t', header = FALSE)
sample4 <- dplyr::select(sample4, c(V1, V7))
colnames(sample4) <- c("EntrezID", "Sample4")
summary(sample4)

sample5 <- read.csv(here("data", "readCounts-sample-5.txt"), sep = '\t', header = FALSE)
sample5 <- dplyr::select(sample5, c(V1, V7))
colnames(sample5) <- c("EntrezID", "Sample5")
summary(sample5)

sample6 <- read.csv(here("data", "readCounts-sample-6.txt"), sep = '\t', header = FALSE)
sample6 <- dplyr::select(sample6, c(V1, V7))
colnames(sample6) <- c("EntrezID", "Sample6")
summary(sample6)

sample7 <- read.csv(here("data", "readCounts-sample-7.txt"), sep = '\t', header = FALSE)
sample7 <- dplyr::select(sample7, c(V1, V7))
colnames(sample7) <- c("EntrezID", "Sample7")
summary(sample7)

sample8 <- read.csv(here("data", "readCounts-sample-8.txt"), sep = '\t', header = FALSE)
sample8 <- dplyr::select(sample8, c(V1, V7))
colnames(sample8) <- c("EntrezID", "Sample8")
summary(sample8)

sample9 <- read.csv(here("data", "readCounts-sample-9.txt"), sep = '\t', header = FALSE)
sample9 <- dplyr::select(sample9, c(V1, V7))
colnames(sample9) <- c("EntrezID", "Sample9")
summary(sample9)

sample10 <- read.csv(here("data", "readCounts-sample-10.txt"), sep = '\t', header = FALSE)
sample10 <- dplyr::select(sample10, c(V1, V7))
colnames(sample10) <- c("EntrezID", "Sample10")
summary(sample10)

sample11 <- read.csv(here("data", "readCounts-sample-11.txt"), sep = '\t', header = FALSE)
sample11 <- dplyr::select(sample11, c(V1, V7))
colnames(sample11) <- c("EntrezID", "Sample11")
summary(sample11)

sample12 <- read.csv(here("data", "readCounts-sample-12.txt"), sep = '\t', header = FALSE)
sample12 <- dplyr::select(sample12, c(V1, V7))
colnames(sample12) <- c("EntrezID", "Sample12")
summary(sample12)

sample13 <- read.csv(here("data", "readCounts-sample-13.txt"), sep = '\t', header = FALSE)
sample13 <- dplyr::select(sample13, c(V1, V7))
colnames(sample13) <- c("EntrezID", "Sample13")
summary(sample13)

sample14 <- read.csv(here("data", "readCounts-sample-14.txt"), sep = '\t', header = FALSE)
sample14 <- dplyr::select(sample14, c(V1, V7))
colnames(sample14) <- c("EntrezID", "Sample14")
summary(sample14)

sample13 == sample14



rownames(sample1) <- sample1[, 1]
sample1[, 1] <- NULL
data <- sample1
head(data)

data <- cbind(data, sample2)
data <- cbind(data, sample3)
data <- cbind(data, sample4)
data <- cbind(data, sample5)
data <- cbind(data, sample6)
data <- cbind(data, sample7)
data <- cbind(data, sample8)
data <- cbind(data, sample9)
data <- cbind(data, sample10)
data <- cbind(data, sample11)
data <- cbind(data, sample12)
data <- cbind(data, sample13)
data <- cbind(data, sample14)

data <- dplyr::select(data, c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8", "Sample9", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14"))
head(data)

# DGE list ----------
x <- DGEList(data)
head(x)
x

group <- as.factor(c("first", "first", "first", "first", "first", "first", "first", "first", "second", "second", "second", "second", "second", "second"))
x$samples$group <- group
x$samples

geneid <- rownames(x)
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "TXCHROM"),
                keytype="ENSEMBL")

genes <- genes[!duplicated(genes$ENSEMBL),]
dim(genes)
head(genes)
x$genes <- genes
x$genes



# preprocessing data --------
lcpm <- cpm(x, log=TRUE)
summary(lcpm)

table(rowSums(x$counts==0)==14)
keep.exprs <- filterByExpr(x)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
samplenames <- colnames(x)

# graphs -------
x1 <- x
lcpm1 <- lcpm

lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(x1)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm1[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm1[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")

# plots to visualize data -------
rownames <- rownames(lcpm)

sortmax <- order(lcpm[ , 1], decreasing = TRUE)

new <- as.data.frame(lcpm)
new <- arrange(new, desc(Sample1))

new["GeneID"] <- rownames(new)
new
new <- melt(new, id.vars="GeneID", value.name= "Counts", variable.name= "Sample")
head(new)
new[18000, ]
colnames(new) <- c("GeneID", "Sample", "Count")
nrow(new)


new %>%
  ggplot(aes(x = GeneID, y = Count, color = Sample)) +
  geom_bar(stat = "identity") +
  facet_grid(~Sample)

# MDS plotting ------
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

col.group
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
g <- ggbiplot(pc, group = genes)
pc <- prcomp(lcpm,
             center = TRUE,
             scale. = TRUE)
pc$scale
attributes(pc)
print(pc)
head(pc)
print(g)
g

# Design matrix -----
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  FirstVsSecond = first - second,
  levels = colnames(design))
contr.matrix

# voom plot -------
v <- voom(x, design, plot=TRUE)
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(tfit, p.value = 0.05))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit, p.value = 0.05)
summary(dt)
rownames(lcpm[which(dt[,1] == -1), ])

which(tfit$p.value < 0.0005)
tfit

de.common <- which(dt[,1]!=0)
length(de.common)

vennDiagram(dt[,1], circle.col=c("turquoise", "salmon"))
write.fit(tfit, dt, file="results.txt")

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))

FirstVsSecond1 <- topTreat(tfit, coef=1, n=Inf)
which(FirstVsSecond1$adj.P.Val < 0.03)
v$genes$ENSEMBL
head(FirstVsSecond)


FirstVsSecond <- FirstVsSecond1$ENSEMBL[1:100]
i <- which(v$genes$ENSEMBL %in% FirstVsSecond)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group,
          col=mycol, trace="none", density.info="none",
          margin=c(8,8), dendrogram="column")

which(v$genes$ENSEMBL %in% FirstVsSecond)
rownames[which(v$genes$ENSEMBL %in% FirstVsSecond)]

rownames(lcpm[which(v$genes$ENSEMBL %in% FirstVsSecond), ])
x$genes$SYMBOL[which(v$genes$ENSEMBL %in% FirstVsSecond)]

UpRegulated <- rownames(lcpm[which(dt[,1] == 1), ])
DownRegulated <- rownames(lcpm[which(dt[,1] == -1), ])

length(UpRegulated)

FirstVsSecond$ENSEMBL

dev.off()
par("mar")

GoPlotBio <- read.csv(here::here("data", "GO_Biological_Process_2023_table.txt"), sep = '\t')
GoPlotBio <- dplyr::select(GoPlotBio, c("Term", "Combined.Score"))
GoPlotBio <- GoPlotBio[order(GoPlotBio$Combined.Score, decreasing = TRUE), ]

GoPlotCell <- read.csv(here::here("data", "GO_Cellular_Component_2023_table.txt"), sep = '\t')
GoPlotCell <- dplyr::select(GoPlotCell, c("Term", "Combined.Score"))
GoPlotCell <- GoPlotCell[order(GoPlotCell$Combined.Score, decreasing = TRUE), ]

GoPlotMol <- read.csv(here::here("data", "GO_Molecular_Function_2023_table.txt"), sep = '\t')
GoPlotMol <- dplyr::select(GoPlotMol, c("Term", "Combined.Score"))
GoPlotMol <- GoPlotMol[order(GoPlotMol$Combined.Score, decreasing = TRUE), ]

kegg <- read.csv(here("data", "kegg.txt"), sep = '\t')
kegg <- dplyr::select(kegg, c("Term", "Combined.Score"))
kegg <- kegg[order(kegg$Combined.Score, decreasing = TRUE), ]

wiki <- read.csv(here("data", "wiki.txt"), sep = '\t')
wiki <- dplyr::select(wiki, c("Term", "Combined.Score"))
wiki <- wiki[order(wiki$Combined.Score, decreasing = TRUE), ]

GoPlotMol$ontology <- "Molecular Function"
GoPlotCell$ontology <- "Cellular Component"
GoPlotBio$ontology <- "Biological Process"

GoPlotBio

all <- rbind(GoPlotMol, GoPlotBio)
all <- rbind(all, GoPlotCell)
all

ggplot(all, aes(y = fct_inorder(Term), x = Combined.Score, fill = ontology)) +
  geom_bar(stat = 'identity')

ggplot(GoPlotBio[1:10, ], aes(y = fct_inorder(Term), x = Combined.Score, fill = Combined.Score)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low="black",high="blue")

ggplot(GoPlotCell[1:10, ], aes(y = fct_inorder(Term), x = Combined.Score, fill = Combined.Score)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low="black",high="red")

ggplot(GoPlotMol[1:10, ], aes(y = fct_inorder(Term), x = Combined.Score, fill = Combined.Score)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low="black",high="red")

trimester(3, 50)

trimester <- function(i, n) {
  #First Trimester

Ngenes<-n
firsttrim <- dplyr::select(as.data.frame(lcpm), c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7", "Sample8"))
scores1<-rowMeans(firsttrim)
scores1<-as.data.frame(scores1)

scores1desc<-dplyr::arrange(scores1,desc(scores1))
barplot1<-head(scores1desc,Ngenes)
head(barplot1)
rownames1<-rownames(barplot1)

#---
genes1<-cbind(x$genes,scores1)
genes1desc<-dplyr::arrange(genes1,desc(scores1))
genes1<-head(genes1desc,Ngenes)
rownames1<-genes1$SYMBOL
#barplot(barplot1$scores1, main = "First Trimester", names.arg =rownames1 , xlab = "GeneId", ylab = "log-CPM", las=2, cex.names=.25,col=rainbow(1))

Ngenes <- n
secondtrim <- dplyr::select(as.data.frame(lcpm),c( "Sample9", "Sample10", "Sample11", "Sample12", "Sample13", "Sample14"))
scores2<-rowMeans(secondtrim)
scores2<-as.data.frame(scores2)
scores2desc<-dplyr::arrange(scores2,desc(scores2))
barplot2<-head(scores2desc,Ngenes)
rownames2<-rownames(barplot2)
head(barplot2)

genes2<-cbind(x$genes,scores2)
genes2desc<-dplyr::arrange(genes2,desc(scores2))
genes2<-head(genes2desc,Ngenes)
rownames2<-genes2$SYMBOL

#barplot(barplot1$scores1, main = "First Trimester", names.arg =rownames1 , xlab = "GeneId", ylab = "log-CPM", las=2, cex.names=.25,col=rainbow(1))

if (i == 1) {
  ggplot(barplot1, aes(x=rownames1,y=barplot1$scores1)) +
  geom_segment(aes(x = rownames1, xend = rownames1, y = 0, yend = barplot1$scores1),color='grey')+
  geom_point(size = 3, color = "darkorange") +
  labs(title = "First Trimester",
         x = "GeneId",
         y = "Log-CPM")+
  theme(axis.text.x = element_text(angle = 90, size =10))
}
#barplot(barplot$scores, main = "second Trimester", names.arg =rownames2 , xlab = "GeneId", ylab = "log-CPM")

else if (i ==2) {
ggplot(barplot2, aes(x=rownames2,y=barplot2$scores2)) +
  geom_segment(aes(x = rownames2, xend = rownames2, y = 0, yend = barplot2$scores2),color='grey')+
  geom_point(size = 3, color = "darkred") +
  labs(title = "Second Trimester",
       x = "GeneId",
       y = "Log-CPM")+
  theme(axis.text.x = element_text(angle = 90, size =10))}

#Venn Diagram

else if (i == 3) {
Ngenes <- n
ggVennDiagram(
  x = list(rownames1,rownames2),
  category.names = c("First Trimester" , "Second Trimester "),


  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,

  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",

  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)}
}
library(extrafont)
font_import()
loadfonts(device = "win")
  ggplot(kegg[1:10, ], aes(y = fct_inorder(Term), x = Combined.Score, fill = Combined.Score)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low = "black", high = "orange") +
  scale_y_discrete(labels = label_wrap(40)) +
  labs(x = "Combined Score", y = "Process") +
  theme(text = element_text(size = 12, family = "Arial"))

library(FactoMineR)
library(factoextra)
res.pca = PCA(lcpm, scale.unit = TRUE, ncp = 5, quali.sup = c(1), graph = TRUE)
fviz_pca_ind(res.pca, geom.ind = "point", col.ind = group, axes = c(1, 2))

#single gene analysis ----------

geneplot <- function(gene, colour_first, colour_second) {
  gene_lcpm <- lcpm[match(gene, rownames2), ]
  gene_lcpm <- melt(gene_lcpm, value.name = "sample", id = "columns")
  gene_lcpm$sample <- rownames(gene_lcpm)
  ggplot(gene_lcpm, aes(x = fct_inorder(sample), y = value, fill = group)) +
    geom_bar(stat = "identity") +
    labs(x = "Sample", y = "log(CPM)") +
    scale_fill_manual(values = c(colour_first, colour_second)) +
    theme(text = element_text(size = 10)) +
    scale_x_discrete(guide = guide_axis(n.dodge=1.2))
}

geneplot("PEX10", "orange", "blue")