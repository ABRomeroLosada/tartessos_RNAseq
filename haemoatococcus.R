library(FactoMineR)
library(factoextra)
library(ballgown)
library(clusterProfiler)
library(pathview)

`if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ballgown")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("limma")


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
write.table(x = gene.expression,file = "haematococcus_complete_gene_expression.tsv",quote = F,sep = "\t")

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

length(activated.genes) #505
length(repressed.genes) #5603

log10.qval <- -log10(q.values)
log10.pval <- -log10(p.values)

png(filename = "figures/volcano_plot.png")
plot(fold.change,log10.qval,pch=19,cex=0.7,col="grey", xlab="Fold Change", ylab="-log10(p-value)",cex.lab=1.5)
points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)
dev.off()

install.packages("org.Hlacustris.eg.db",repos = NULL)
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

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         pathway.id = "ko00900",
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

## Carotenoides


barplot.car <- function(carotenoid, car.15mM.mean, car.15mM.sd, car.2mM.mean, car.2mM.sd)
{
  means <- c(car.15mM.mean,car.2mM.mean)
  sds <- c(car.15mM.sd,car.2mM.sd)
  
  png(filename = paste0(carotenoid,".png"),width = 300)
  par(lwd=6,mar=c(7,6,1,1),font.axis=2)
  xpos <- barplot(means,col=c("springgreen4","firebrick2"),names.arg = c("15mM","2mM"),las=2,cex.names = 2.5,ylim=c(0,max(means+sds)*1.2),cex.axis = 2.5,lwd=6)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,code = 3,angle=90,lwd=6)
  dev.off()
}



## Carotenoid enzyme pathway
gene.expression <- read.table(file = "haematococcus_complete_gene_expression.tsv",header = T,sep = "\t")
head(gene.expression)

library("org.Hlacustris.eg.db")
columns(org.Hlacustris.eg.db)
hlacustris.ko <- select(org.Hlacustris.eg.db,columns = c("GID","KO"),
            keys=keys(org.Hlacustris.eg.db,keytype = "GID"))
head(hlacustris.ko)

## Correspondence between genome notation and NCBI notation.
halan.gfh.correspondence <- read.table(file="halaN_GFH_correspondence/gfh_halan_correspondence.tsv",sep="\t",as.is=T)
head(halan.gfh.correspondence)

halan <- halan.gfh.correspondence$V2
gfh <- halan.gfh.correspondence$V1
names(halan) <- gfh
names(gfh) <- halan

barplot.enzyme.ko <- function(ko,enzyme.name,map.gfh.halan,gene.expression)
{
  halan["GFH11097"]
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,KO == ko)$GID],rownames(gene.expression))
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    means <- c(enzyme.mean.15mM,enzyme.mean.2mM)
    
    enzyme.sd.2mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.sd.15mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    sds <- c(enzyme.sd.15mM,enzyme.sd.2mM)
    
    
    png(filename = paste(c(enzyme.name,i,".png"),collapse=""),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("springgreen4","firebrick2"),names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,code = 3,angle=90,lwd=2)
    dev.off()
  }
  
  return(halan.gene.names)
}

gfh.id <- "GFH11097"

barplot.enzyme.gfh <- function(gfh.id,enzyme.name,map.gfh.halan,gene.expression)
{
  halan.gene.names <- halan[gfh.id]
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    means <- c(enzyme.mean.15mM,enzyme.mean.2mM)
    
    enzyme.sd.2mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.sd.15mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    sds <- c(enzyme.sd.15mM,enzyme.sd.2mM)
    
    
    png(filename = paste(c(enzyme.name,i,".png"),collapse=""),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("springgreen4","firebrick2"),names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,code = 3,angle=90,lwd=2)
    dev.off()
  }
  
  return(halan.gene.names)
}

barplot.enzyme.gfh(gfh.id="GFH11097",enzyme.name = "PDH", map.gfh.halan = halan, gene.expression = gene.expression)
barplot.enzyme.gfh(gfh.id="GFH30277",enzyme.name = "SDH", map.gfh.halan = halan, gene.expression = gene.expression)
barplot.enzyme.gfh(gfh.id="GFH33388",enzyme.name = "ACACA", map.gfh.halan = halan, gene.expression = gene.expression)


barplot.enzyme.pfam <- function(pfam,enzyme.name,map.gfh.halan,gene.expression)
{
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,pfam == pfam)$GID],rownames(gene.expression))
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    means <- c(enzyme.mean.15mM,enzyme.mean.2mM)
    
    enzyme.sd.2mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.sd.15mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    sds <- c(enzyme.sd.15mM,enzyme.sd.2mM)
    
    
    png(filename = paste(c(enzyme.name,i,".png"),collapse=""),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("springgreen4","firebrick2"),names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,code = 3,angle=90,lwd=2)
    dev.off()
  }
}

barplot.enzyme.pfam(pfam="PF00067",enzyme.name="carotenoid_epsilon_hydroxylase",map.gfh.halan=halan,gene.expression=gene.expression)



library(gplots)


heatmap.enzyme.ko <- function(ko,enzyme.name,map.gfh.halan,gene.expression,precision)
{
  pos.colfunc <- colorRampPalette(c("white","firebrick2"))
  pos.colors <- pos.colfunc(precision)
  
  neg.colfunc <- colorRampPalette(c("white","springgreen4"))
  neg.colors <- neg.colfunc(precision)
  
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,KO == ko)$GID],rownames(gene.expression))
  
  i <- 1
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    fc <- enzyme.mean.2mM/enzyme.mean.15mM
    if(fc >= 1)
    {
      if(round(fc*10) < precision)
      {
        cols <- pos.colors[round(fc*10)]
      } else
      {
        cols <- pos.colors[precision]
      }
    } else
    {
      if(round(10/fc) < precision)
      {
        cols <- neg.colors[round(10/fc)]
      } else
      {
        cols <- neg.colors[precision]
      }
    }

    png(filename = paste(c("heatmap_",enzyme.name,i,".png"),collapse=""))
    plot(x=0,y=0,col="white",axes=F,xlab="",ylab="",ylim=c(0,10),xlim=c(0,10))
    polygon(x = c(2,8,8,2),y=c(2,2,8,8),lwd=6,col = cols)
    dev.off()
  }
  
  return(list(halan.gene.names,fc,1/fc))
}

heatmap.enzyme.gfh <- function(gfh.id,enzyme.name,map.gfh.halan,gene.expression,precision)
{
  pos.colfunc <- colorRampPalette(c("white","firebrick2"))
  pos.colors <- pos.colfunc(precision)
  
  neg.colfunc <- colorRampPalette(c("white","springgreen4"))
  neg.colors <- neg.colfunc(precision)
  
  halan.gene.names <- halan[gfh.id]
  
  i <- 1
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    fc <- enzyme.mean.2mM/enzyme.mean.15mM
    if(fc >= 1)
    {
      if(round(fc*10) < precision)
      {
        cols <- pos.colors[round(fc*10)]
      } else
      {
        cols <- pos.colors[precision]
      }
    } else
    {
      if(round(10/fc) < precision)
      {
        cols <- neg.colors[round(10/fc)]
      } else
      {
        cols <- neg.colors[precision]
      }
    }
    
    png(filename = paste(c("heatmap_",enzyme.name,i,".png"),collapse=""))
    plot(x=0,y=0,col="white",axes=F,xlab="",ylab="",ylim=c(0,10),xlim=c(0,10))
    polygon(x = c(2,8,8,2),y=c(2,2,8,8),lwd=6,col = cols)
    dev.off()
  }
  
  return(list(halan.gene.names,fc,1/fc))
}
##########################Carotenoids 

barplot.enzyme.gfh(gfh.id="GFH11424",enzyme.name = "ISPH", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.gfh(gfh.id="GFH11424",enzyme.name = "ISPH", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.gfh(gfh.id="GFH09986",enzyme.name = "GGPS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.gfh(gfh.id="GFH09986",enzyme.name = "GGPS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.gfh(gfh.id="GFH18487",enzyme.name = "IDI", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.gfh(gfh.id="GFH18487",enzyme.name = "IDI", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01823",enzyme.name = "IDI", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K01823",enzyme.name = "IDI", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K02291",enzyme.name = "0_PS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K02291",enzyme.name = "0_PS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K02293",enzyme.name = "1_PDS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K02293",enzyme.name = "1_PDS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00514",enzyme.name = "2_ZDS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00514",enzyme.name = "2_ZDS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K06443",enzyme.name = "3_LCYB", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K06443",enzyme.name = "3_LCYB", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

## Rama epsilon
barplot.enzyme.ko(ko="K06444",enzyme.name = "4_LCYE", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K06444",enzyme.name = "4_LCYE", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K09837",enzyme.name = "carotenoid_epsilon_hydroxylase", map.gfh.halan = halan, gene.expression = gene.expression)
## No aparece carotenoid_epsilon_hydroxylase

barplot.enzyme.ko(ko="K15747",enzyme.name = "5_LUT5", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K15747",enzyme.name = "5_LUT5", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

## Rama beta
barplot.enzyme.ko(ko="K09836",enzyme.name = "6_BKT", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K09836",enzyme.name = "6_BKT", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K15746",enzyme.name = "7_crtZ", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K15746",enzyme.name = "7_crtZ", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K09838",enzyme.name = "8_ZEP", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K09838",enzyme.name = "8_ZEP", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

#no está
barplot.enzyme.ko(ko="K09839",enzyme.name = "violaxanthin_de_epoxidase", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K09839",enzyme.name = "violaxanthin_de_epoxidase", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)

#carontenoides cuantificados
barplot.car(carotenoid = "Violaxanthin",car.15mM.mean = 0.045, car.15mM.sd = 0.004, car.2mM.mean = 0.027, car.2mM.sd = 0.002)
barplot.car(carotenoid = "Astaxanthin",car.15mM.mean = 0.014, car.15mM.sd = 0.002, car.2mM.mean = 0.208, car.2mM.sd = 0.036)
barplot.car(carotenoid = "Lutein",car.15mM.mean = 0.129, car.15mM.sd = 0.008, car.2mM.mean = 0.071, car.2mM.sd = 0.003)
barplot.car(carotenoid = "Cantaxanthin",car.15mM.mean = 0.001, car.15mM.sd = 0.000, car.2mM.mean = 0.009, car.2mM.sd = 0.002)
barplot.car(carotenoid = "b-carotene",car.15mM.mean = 0.050, car.15mM.sd = 0.004, car.2mM.mean = 0.022, car.2mM.sd = 0.001)






#####################################################################################
##TCA
barplot.enzyme.ko(ko="K00873",enzyme.name = "1_PK", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00873",enzyme.name = "1_PK", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

#no funciona con el KO, hay que hacerlo con GFH
barplot.enzyme.gfh(gfh.id="GFH11097",enzyme.name = "2_PDH", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.gfh(gfh.id="GFH11097",enzyme.name = "2_PDH", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K01647",enzyme.name = "3_CS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K01647",enzyme.name = "3_CS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K01681",enzyme.name = "4_ACO", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K01681",enzyme.name = "4_ACO", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K00030",enzyme.name = "5_IDH", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00030",enzyme.name = "5_IDH", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K00164",enzyme.name = "6_alpha-KGDH E1", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00164",enzyme.name = "6_alpha-KGDH E1", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K00658",enzyme.name = "7_alpha-KGDH E2", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00658",enzyme.name = "7_alpha-KGDH E2", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K01899",enzyme.name = "8_SCS", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K01899",enzyme.name = "8_SCS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

barplot.enzyme.ko(ko="K00234",enzyme.name = "9_SDH", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K00234",enzyme.name = "9_SDH", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)

#no significativa hacer heatmap blanco
heatmap.enzyme.ko(ko="K01676",enzyme.name = "10_fum", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)
barplot.enzyme.ko(ko="K01676",enzyme.name = "10_fum", map.gfh.halan = halan, gene.expression = gene.expression)

heatmap.enzyme.ko(ko="K00025",enzyme.name = "11_MDH", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)
barplot.enzyme.ko(ko="K00025",enzyme.name = "11_MDH", map.gfh.halan = halan, gene.expression = gene.expression)

#barra de error muy grande
heatmap.enzyme.gfh(gfh.id = "GFH26871",enzyme.name = "12_EM", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)
barplot.enzyme.gfh(gfh.id = "GFH26871",enzyme.name = "12_EM", map.gfh.halan = halan, gene.expression = gene.expression)

#barra de error muy grande
heatmap.enzyme.gfh(gfh.id = "GFH33388",enzyme.name = "13_ACACA", map.gfh.halan = halan, gene.expression = gene.expression,precision = 50)
barplot.enzyme.gfh(gfh.id = "GFH33388",enzyme.name = "13_ACACA", map.gfh.halan = halan, gene.expression = gene.expression)

#ácidos orgánicos
barplot.car(carotenoid = "Glucose",car.15mM.mean = 0.12, car.15mM.sd = 0.01, car.2mM.mean = 0.18, car.2mM.sd = 0.02)
barplot.car(carotenoid = "Piruvate",car.15mM.mean = 0.96, car.15mM.sd = 0.10, car.2mM.mean = 0.55, car.2mM.sd = 0.06)
barplot.car(carotenoid = "Citrate",car.15mM.mean = 1.59, car.15mM.sd = 0.10, car.2mM.mean = 1.98, car.2mM.sd = 0.14)
barplot.car(carotenoid = "cetoglutarato",car.15mM.mean = 0.45, car.15mM.sd = 0.13, car.2mM.mean = 0.17, car.2mM.sd = 0.03)
barplot.car(carotenoid = "succinate",car.15mM.mean = 1.53, car.15mM.sd = 0.07, car.2mM.mean = 1.21, car.2mM.sd = 0.07)
barplot.car(carotenoid = "fumarate",car.15mM.mean = 0.62, car.15mM.sd = 0.21, car.2mM.mean = 0.28, car.2mM.sd = 0.05)
barplot.car(carotenoid = "malate",car.15mM.mean = 5.78, car.15mM.sd = 0.20, car.2mM.mean = 2.88, car.2mM.sd = 0.11)




###################################################################
##Acidos grasos. FAS
#GFH27381 ACP synthase
barplot.enzyme.gfh(gfh.id = "GFH27381",enzyme.name = "ACP synthase", map.gfh.halan = halan, gene.expression = gene.expression)



heatmap.enzyme.ko(ko="K09458",enzyme.name = "14_KS", map.gfh.halan = halan, gene.expression = gene.expression,precision = 40)
barplot.enzyme.ko(ko="K09458",enzyme.name = "14_KS", map.gfh.halan = halan, gene.expression = gene.expression)

heatmap.enzyme.ko(ko="K00059",enzyme.name = "14_MAT", map.gfh.halan = halan, gene.expression = gene.expression,precision = 40)
barplot.enzyme.ko(ko="K00059",enzyme.name = "14_MAT", map.gfh.halan = halan, gene.expression = gene.expression)

barplot.car(carotenoid = "C16:1",car.15mM.mean = 2.25, car.15mM.sd = 0.14, car.2mM.mean = 1.17, car.2mM.sd = 0.04)
barplot.car(carotenoid = "C18:1",car.15mM.mean = 15.81, car.15mM.sd = 0.36, car.2mM.mean = 21.89, car.2mM.sd = 2.33)
barplot.car(carotenoid = "C18:2",car.15mM.mean = 13.03, car.15mM.sd = 0.44, car.2mM.mean = 16.66, car.2mM.sd = 0.68)
barplot.car(carotenoid = "C18:3n4",car.15mM.mean = 20.60, car.15mM.sd = 0.55, car.2mM.mean = 14.47, car.2mM.sd = 0.44)




#Sucrose metabolism, no se expresa
barplot.enzyme.ko(ko="K01193",enzyme.name = "sacA", map.gfh.halan = halan, gene.expression = gene.expression)
heatmap.enzyme.ko(ko="K01193",enzyme.name = "sacA", map.gfh.halan = halan, gene.expression = gene.expression,precision = 78)


##Enzimas sueltas que nos interesan
#GFH21816 aspartato aminotransferasa
barplot.enzyme.gfh(gfh.id="GFH21816",enzyme.name = "Asp aminotrans", map.gfh.halan = halan, gene.expression = gene.expression)
barplot.enzyme.gfh(gfh.id="GFH09772",enzyme.name = "Alan aminotrans", map.gfh.halan = halan, gene.expression = gene.expression)



barplot.enzyme.pfam <- function(pfam,enzyme.name,map.gfh.halan,gene.expression)
{
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,pfam == pfam)$GID],rownames(gene.expression))
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    means <- c(enzyme.mean.15mM,enzyme.mean.2mM)
    
    enzyme.sd.2mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_4mM_1", "NO3_4mM_2")]))
    enzyme.sd.15mM <- sd(unlist(gene.expression[halan.gene.names[i],c("NO3_12mM_1", "NO3_12mM_2")]))
    sds <- c(enzyme.sd.15mM,enzyme.sd.2mM)
    
    
    png(filename = paste(c(enzyme.name,i,".png"),collapse=""),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("darkgreen","darkred"),names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,code = 3,angle=90,lwd=2)
    dev.off()
  }
}

barplot.enzyme.pfam(pfam="PF00067",enzyme.name="carotenoid_epsilon_hydroxylase",map.gfh.halan=halan,gene.expression=gene.expression)
