library(FactoMineR)
library(factoextra)
library(ballgown)
library(clusterProfiler)
library(pathview)

# Load experimental design
experimental.design <- read.csv("experimental_design.csv",as.is=T)
experimental.design$sample

## Load results from hisat2 + stringtie
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=experimental.design)

## Extract gene expression and name columns with sample.labels
gene.expression <- gexpr(bg.data)
colnames(gene.expression) <- c("NO3_4mM_1","NO3_4mM_2","NO3_12mM_1","NO3_12mM_2")
head(gene.expression)
nrow(gene.expression)

## Remove genes never expressed
gene.expression <- gene.expression[apply(X = gene.expression,MARGIN = 1,sum) != 0,]
head(gene.expression)
nrow(gene.expression)

## Remove genes with very different expression between replicates
gene.expression <- gene.expression[!(((gene.expression[,"NO3_4mM_1"] == 0 & gene.expression[,"NO3_4mM_2"] > 0) | (gene.expression[,"NO3_4mM_1"] > 0 & gene.expression[,"NO3_4mM_2"] == 0)) |
                                     ((gene.expression[,"NO3_12mM_1"] == 0 & gene.expression[,"NO3_12mM_2"] > 0) | (gene.expression[,"NO3_12mM_1"] > 0 & gene.expression[,"NO3_12mM_2"] == 0))),]
nrow(gene.expression)
head(gene.expression)

## Write and read gne expression matrix
write.table(x = gene.expression,file = "haematococcus_gene_expression.tsv",quote = F,sep = "\t")
gene.expression <- read.table(file = "haematococcus_gene_expression.tsv",header = T,sep = "\t")
head(gene.expression)

## PCA
pca.gene.expression <- data.frame(colnames(gene.expression),t(gene.expression))
colnames(pca.gene.expression)[1] <- "Sample"

res.pca <- PCA(pca.gene.expression, graph = FALSE,scale.unit = TRUE,quali.sup = 1 )
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70),main = "")

## Hierarchical clustering
res.hcpc <- HCPC(res.pca, graph=FALSE,nb.clust = 2)    

png(filename = "figures/hierarchical_tree.png",width = 400,height = 400)
fviz_dend(res.hcpc,k=2,
          cex = 1,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 4500      # Augment the room for labels
)
dev.off()

lowN <- rowMeans(gene.expression[,1:2])
highN <- rowMeans(gene.expression[,3:4])

png(filename = "figures/haematococcus_scatter.png")
plot(log2(highN+1),log2(lowN+1),pch=19,cex=0.7,xlab="12mM",ylab="4mM",xlim=c(0,15),ylim=c(0,15),col="darkgrey",cex.lab=1.5)
lines(x=c(-1,20),y=c(-1,20),col="red",lwd=3)
dev.off()

png(filename = "figures/replicates_scatter_4mM.png")
plot(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1),pch=19,cex=0.7,xlab="NO3_4mM_1",ylab="NO3_4mM_2",col="darkgrey",cex.lab=1.5)
lines(x=c(-1,20),y=c(-1,20),col="red",lwd=3)
text(x = 2,y=11,labels = paste(round(100*cor(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1)),digits=2),"%"),cex=1.5)
dev.off()

png(filename = "figures/replicates_scatter_12mM.png")
plot(log2(gene.expression[,3]+1),log2(gene.expression[,4]+1),pch=19,cex=0.7,xlab="NO3_12mM_1",ylab="NO3_12mM_2",col="darkgrey",cex.lab=1.5)
lines(x=c(-1,20),y=c(-1,20),col="red",lwd=3)
text(x = 2,y=11,labels = paste(round(100*cor(log2(gene.expression[,3]+1),log2(gene.expression[,4]+1)),digits=2),"%"),cex=1.5)
dev.off()

## Correspondence between genome notation and NCBI notation.
halan.gfh.correspondence <- read.table(file="halaN_GFH_correspondence/gfh_halan_correspondence.tsv",sep="\t",as.is=T)
head(halan.gfh.correspondence)

halan <- halan.gfh.correspondence$V2
gfh <- halan.gfh.correspondence$V1
names(halan) <- gfh
names(gfh) <- halan

## Some genes of interest
hala <- read.table(file="enzimas Haematococcus.csv",header=T,sep="\t",as.is=T)
for(i in 1:nrow(hala))
{
  current.gene.name <- hala[i,1]
  gene.ids <- strsplit(hala[i,2],split=",")[[1]]
  if(gene.ids != "no encontrada")
  {
    for(j in 1:length(gene.ids))
    {
      current.gene <- halan[gene.ids[j]]
      if(current.gene %in% rownames(gene.expression))
      {
        sds <- c(sd(gene.expression[current.gene,3:4]),
                 sd(gene.expression[current.gene,1:2]))
        means <- c(highN[current.gene],lowN[current.gene])
        png(filename = paste("barplots/",paste(gene.ids[j],".png",sep = ""),sep=""))
        xpos <- barplot(means,ylim = c(0,1.1*max(means+sds)),col=c("blue","red"),main=current.gene.name,names.arg = c("12mM","4mM"),ylab="FPKM")
        arrows(x0 = xpos,y0 =means-sds,x1 = xpos,y1 = means+sds,code=3,angle=90,length=0.1,lwd=2 )
        dev.off()
      }
    }
  }
}

## Log2 transformation
log.gene.expression <- log2(gene.expression+1)
head(log.gene.expression)
nrow(log.gene.expression)

## Differential expression analysis
library(limma)

## Specification of the experimental design
factor.experimental.design <- c(1,1,2,2)
limma.experimental.design <- model.matrix(~ -1+factor(factor.experimental.design))
colnames(limma.experimental.design) <- c("NO3_4mM", "NO3_12mM")

## Linear model fit
linear.fit <- lmFit(log.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(NO3_4mM-NO3_12mM,
                                 levels=c("NO3_4mM", "NO3_12mM"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract results
degs.results <- topTable(contrast.results, number=nrow(log.gene.expression),coef=1,sort.by="logFC")
head(degs.results)

fold.change <- degs.results$logFC
p.values <- degs.results$P.Value
q.values <- degs.results$adj.P.Val
genes.ids <- rownames(degs.results)

names(fold.change) <- genes.ids
names(q.values) <- genes.ids
names(p.values) <- genes.ids

fc.threshold <- 2
q.val.threshold <- 0.05

activated.genes <- genes.ids[fold.change > log2(fc.threshold) & q.values < q.val.threshold]
repressed.genes <- genes.ids[fold.change < - log2(fc.threshold) & q.values < q.val.threshold]

length(activated.genes)
length(repressed.genes)

log10.qval <- -log10(q.values)
log10.pval <- -log10(p.values)

png(filename = "figures/volcano_plot.png")
plot(fold.change,log10.qval,pch=19,cex=0.7,col="grey", xlab="Fold Change", ylab="-log10(p-value)",cex.lab=1.5)
points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)
dev.off()

#install.packages("org.Hlacustris.eg.db",repos = NULL)
library("org.Hlacustris.eg.db")

enrich.go <- enrichGO(gene          = gfh[activated.genes],
                      OrgDb         = org.Hlacustris.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE,
                      keyType = "GID")

write.table(x = as.data.frame(enrich.go),file = "go_enrichment/go_activated.tsv",quote = F,sep = "\t",row.names = F)

png(filename = "go_enrichment/barplot_go_activated.png")
barplot(enrich.go)
dev.off()

png(filename = "go_enrichment/emap_go_activated.png")
emapplot(enrich.go)
dev.off()

enrich.go <- enrichGO(gene          = gfh[repressed.genes],
                      OrgDb         = org.Hlacustris.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE,
                      keyType = "GID")

write.table(x = as.data.frame(enrich.go),file = "go_enrichment/go_repressed.tsv",quote = F,sep = "\t",row.names = F)

png(filename = "go_enrichment/barplot_go_repressed.png")
barplot(enrich.go)
dev.off()

png(filename = "go_enrichment/emap_go_repressed.png")
emapplot(enrich.go)
dev.off()

## KEEG pathway enrichment
hlac.ko <- select(org.Hlacustris.eg.db,columns = c("KO"),keys=keys(org.Hlacustris.eg.db,keytype = "GID"))
ko.universe <- hlac.ko$KO
ko.universe <- ko.universe[!is.na(ko.universe)]

activated.target.ko <- subset(hlac.ko,GID %in% gfh[activated.genes])$KO
activated.target.ko <- activated.target.ko[!is.na(activated.target.ko)]

pathway.enrichment <- as.data.frame(enrichKEGG(gene = activated.target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05))

write.table(x = pathway.enrichment,file = "pathway_enrichment/pathway_activated.tsv",quote = F,sep = "\t",row.names = F)

ko.activated <- pathway.enrichment$ID

repressed.target.ko <- subset(hlac.ko,GID %in% gfh[repressed.genes])$KO
repressed.target.ko <- repressed.target.ko[!is.na(repressed.target.ko)]

pathway.enrichment <- as.data.frame(enrichKEGG(gene = repressed.target.ko, organism = "ko", universe = ko.universe,qvalueCutoff = 0.05))

## Graphical representation of pathways
genes.pathway <- rep(0, length(ko.universe))
names(genes.pathway) <- ko.universe

genes.pathway[activated.target.ko] <- 1
genes.pathway[repressed.target.ko] <- -1

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00906",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00020",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00061",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00620",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00640",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00640",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00010",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko01212",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))


pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00480",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))



pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko04141",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00500",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))


## Some barplots
current.gene <- "HaLaN_05256"#"HaLaN_01176"#"HaLaN_03235"#"HaLaN_11280"
current.gene.name <- "pyruvate dehydrogenase"

sds <- c(sd(gene.expression[current.gene,3:4]),
         sd(gene.expression[current.gene,1:2]))
means <- c(highN[current.gene],lowN[current.gene])
png(filename = paste(current.gene,".png",sep = ""))
xpos <- barplot(means,ylim = c(0,1.1*max(means+sds)),col=c("blue","red"),main=current.gene.name,names.arg = c("12mM","4mM"),ylab="FPKM")
arrows(x0 = xpos,y0 =means-sds,x1 = xpos,y1 = means+sds,code=3,angle=90,length=0.1,lwd=2 )
dev.off()
