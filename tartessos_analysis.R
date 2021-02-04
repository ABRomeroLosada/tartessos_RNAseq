## Install and load required libraries
library(clusterProfiler)
install.packages("packages/org.Hlacustris.eg.db",repos = NULL)
install.packages("packages/TxDb.Hlacustris.NCBI",repos = NULL)
library(org.Hlacustris.eg.db)
library(TxDb.Hlacustris.NCBI)
library(pathview)
library(gplots)
library(seqinr)

## Load genome sequence data
hlacustris.info <- read.fasta(file = "haematococcus_lacustris.fa",seqtype = "DNA")
hlacustris.seqs <- getSequence(hlacustris.info)
names(hlacustris.seqs) <- getName(hlacustris.info)

txdb <- TxDb.Hlacustris.NCBI

hlacustris.genes <- as.data.frame(genes(txdb))

## Load gene expression data
gene.expression.data <- read.table(file = "gene_expression.tsv",header = T,sep = "\t")
gene.expression <- as.matrix(gene.expression.data[,2:5])
rownames(gene.expression) <- gene.expression.data$Gene
log.gene.expression <- log2(gene.expression+1)
head(log.gene.expression)

## Differential expression analysis
library(limma)

## Specification of the experimental design
factor.experimental.design <- c(1,1,2,2)
limma.experimental.design <- model.matrix(~ -1+factor(factor.experimental.design))
colnames(limma.experimental.design) <- c("NO3_2mM", "NO3_15mM")

## Linear model fit
linear.fit <- lmFit(log.gene.expression, limma.experimental.design)

## Contrast specification and computation
contrast.matrix <- makeContrasts(NO3_2mM-NO3_15mM,
                                 levels=c("NO3_2mM", "NO3_15mM"))

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Extract differentially expressed genes
degs.results <- topTable(contrast.results, number=nrow(log.gene.expression),coef=1,sort.by="logFC")
head(degs.results)

fold.change <- degs.results$logFC
q.values <- degs.results$adj.P.Val
genes.ids <- rownames(degs.results)

names(fold.change) <- genes.ids
names(q.values) <- genes.ids

fc.threshold <- 2
q.val.threshold <- 0.05

activated.genes <- genes.ids[fold.change > log2(fc.threshold) & q.values < q.val.threshold]
repressed.genes <- genes.ids[fold.change < - log2(fc.threshold) & q.values < q.val.threshold]

length(activated.genes) #414
length(repressed.genes) #5348
length(activated.genes)/nrow(gene.expression)
length(repressed.genes)/nrow(gene.expression)
(length(activated.genes) + length(repressed.genes)) / nrow(gene.expression) # 43.49%

## Generate volcano plot
log10.qval <- -log10(q.values)

png(filename = "volcano_plot.png")
plot(fold.change,log10.qval,pch=19,cex=0.5,col="grey",
     xlim=c(-15,15), ylim=c(0,3),
     xlab="Log2 Fold Change", ylab="-log10(q-value)",cex.lab=1.5)
points(fold.change[activated.genes],log10.qval[activated.genes],cex=0.7,col="red",pch=19)
points(fold.change[repressed.genes],log10.qval[repressed.genes],cex=0.7,col="blue",pch=19)
dev.off()

## Correspondence between genome notation and NCBI notation.
halan.gfh.correspondence <- read.table(file="gfh_halan_correspondence.tsv",sep="\t",as.is=T)
head(halan.gfh.correspondence)

halan <- halan.gfh.correspondence$V2
gfh <- halan.gfh.correspondence$V1
names(halan) <- gfh
names(gfh) <- halan

## KEEG pathway enrichment
hlac.ko <- select(org.Hlacustris.eg.db,columns = c("KO"),keys=keys(org.Hlacustris.eg.db,keytype = "GID"))
ko.universe <- hlac.ko$KO
ko.universe <- ko.universe[!is.na(ko.universe)]

activated.target.ko <- subset(hlac.ko,GID %in% gfh[activated.genes])$KO
activated.target.ko <- activated.target.ko[!is.na(activated.target.ko)]

repressed.target.ko <- subset(hlac.ko,GID %in% gfh[repressed.genes])$KO
repressed.target.ko <- repressed.target.ko[!is.na(repressed.target.ko)]

pathway.enrichment <- enrichKEGG(gene = activated.target.ko, organism = "ko", 
                                 universe = ko.universe,qvalueCutoff = 0.05)
barplot(height = pathway.enrichment,showCategory = 11)

write.table(x = as.data.frame(pathway.enrichment),
            file = "pathway_enrichment/pathway_activated.tsv",
            quote = F,sep = "\t",row.names = F)

head(as.data.frame(pathway.enrichment))

ko.activated <- pathway.enrichment$ID

## Graphical representation of pathways
genes.pathway <- rep(0, length(ko.universe))
names(genes.pathway) <- ko.universe

genes.pathway[activated.target.ko] <- 1
genes.pathway[repressed.target.ko] <- -1

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00906",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/",pathway.id = "ko00020",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00061",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00620",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00640",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00010",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko01212",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00480",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00500",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00710",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko01200",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00720",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00030",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))

pathview(gene.data = sort(genes.pathway,decreasing = TRUE),
         kegg.dir ="pathway_enrichment/", pathway.id = "ko00910",
         species = "ko",
         limit = list(gene=max(abs(genes.pathway)), cpd=1))


## Identification of transcription factors
haematococcus.annotation <- read.table(file = "haematoccocus_functional_annotation.csv",
                                       header = T,as.is = T,fill = T)
head(haematococcus.annotation)

extract.tf <- function(annotation, pfam.id)
{
  return(unique(subset(annotation,PFAM == pfam.id)$GENES))
}

## AP2 
AP2.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00847")
AP2.genes
length(AP2.genes)
ap2.activated <- length(intersect(activated.genes, halan[AP2.genes]))
ap2.repressed <- length(intersect(repressed.genes, halan[AP2.genes]))
ap2.non.diff <- length(AP2.genes) - ap2.activated - ap2.repressed

## B3 
B3.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02362")
b3.activated <- length(intersect(activated.genes, halan[B3.genes]))
b3.repressed <- length(intersect(repressed.genes, halan[B3.genes]))
b3.non.diff <- length(B3.genes) - b3.activated - b3.repressed

##ARR-B
ARR.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00072")
length(ARR.genes)
ARR.activated <- length(intersect(activated.genes, halan[ARR.genes]))
ARR.repressed <- length(intersect(repressed.genes, halan[ARR.genes]))
ARR.non.diff <- length(ARR.genes) - ARR.activated - ARR.repressed

##zf-C2H2 
C2H2.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF12756")
length(C2H2.genes)
C2H2.activated <- length(intersect(activated.genes, halan[C2H2.genes]))
C2H2.repressed <- length(intersect(repressed.genes, halan[C2H2.genes]))
C2H2.non.diff <- length(C2H2.genes) - C2H2.activated - C2H2.repressed

##zf-CCCH
C3H.genes.1 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00642")
C3H.genes.2 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF14608")
C3H.genes.4 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF18044")

C3H.genes <- unique(c(C3H.genes.1, C3H.genes.2, C3H.genes.4))
length(C3H.genes)
C3H.activated <- length(intersect(activated.genes, halan[C3H.genes]))
C3H.repressed <- length(intersect(repressed.genes, halan[C3H.genes]))
C3H.non.diff <- length(C3H.genes) - C3H.activated - C3H.repressed

##CPP
CPP.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03638")
length(CPP.genes)
CPP.activated <- length(intersect(activated.genes, halan[CPP.genes]))
CPP.repressed <- length(intersect(repressed.genes, halan[CPP.genes]))
CPP.non.diff <- length(CPP.genes) - CPP.activated - CPP.repressed

##Dof
Dof.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02701")
length(Dof.genes)
Dof.activated <- length(intersect(activated.genes, halan[Dof.genes]))
Dof.repressed <- length(intersect(repressed.genes, halan[Dof.genes]))
Dof.non.diff <- length(Dof.genes) - Dof.activated - Dof.repressed

##E2F/DP
E2F.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02319")
length(E2F.genes)
E2F.activated <- length(intersect(activated.genes, halan[E2F.genes]))
E2F.repressed <- length(intersect(repressed.genes, halan[E2F.genes]))
E2F.non.diff <- length(E2F.genes) - E2F.activated - E2F.repressed

##GATA
GATA.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00320")
length(GATA.genes)
GATA.activated <- length(intersect(activated.genes, halan[GATA.genes]))
GATA.repressed <- length(intersect(repressed.genes, halan[GATA.genes]))
GATA.non.diff <- length(GATA.genes) - GATA.activated - GATA.repressed

## HB
HB.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00046")
length(HB.genes)
HB.activated <- length(intersect(activated.genes, halan[HB.genes]))
HB.repressed <- length(intersect(repressed.genes, halan[HB.genes]))
HB.non.diff <- length(HB.genes) - HB.activated - HB.repressed

##HSF
HSF.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00447")
length(HSF.genes)
HSF.activated <- length(intersect(activated.genes, halan[HSF.genes]))
HSF.repressed <- length(intersect(repressed.genes, halan[HSF.genes]))
HSF.non.diff <- length(HSF.genes) - HSF.activated - HSF.repressed

##zf-LSD1
LSD.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF06943")
length(LSD.genes)
LSD.activated <- length(intersect(activated.genes, halan[LSD.genes]))
LSD.repressed <- length(intersect(repressed.genes, halan[LSD.genes]))
LSD.non.diff <- length(LSD.genes) - LSD.activated - LSD.repressed

##MYB
myb.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00249")
length(myb.genes)
myb.activated <- length(intersect(activated.genes, halan[myb.genes]))
myb.repressed <- length(intersect(repressed.genes, halan[myb.genes]))
myb.non.diff <- length(myb.genes) - myb.activated - myb.repressed

##NF-Y
NF.Y.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00808")
length(NF.Y.genes)
NF.Y.activated <- length(intersect(activated.genes, halan[NF.Y.genes]))
NF.Y.repressed <- length(intersect(repressed.genes, halan[NF.Y.genes]))
NF.Y.non.diff <- length(NF.Y.genes) - NF.Y.activated - NF.Y.repressed

##Nin-like 
Nin.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02042")
length(Nin.genes)
Nin.activated <- length(intersect(activated.genes, halan[Nin.genes]))
Nin.repressed <- length(intersect(repressed.genes, halan[Nin.genes]))
Nin.non.diff <- length(Nin.genes) - Nin.activated - Nin.repressed

##S1FA
S1FA.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF04689")
length(S1FA.genes)
S1FA.activated <- length(intersect(activated.genes, halan[S1FA.genes]))
S1FA.repressed <- length(intersect(repressed.genes, halan[S1FA.genes]))
S1FA.non.diff <- length(S1FA.genes) - S1FA.activated - S1FA.repressed

##SBP  
SBP.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03110")
length(SBP.genes)
SBP.activated <- length(intersect(activated.genes, halan[SBP.genes]))
SBP.repressed <- length(intersect(repressed.genes, halan[SBP.genes]))
SBP.non.diff <- length(SBP.genes) - SBP.activated - SBP.repressed

##WRKY
WRKY.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03106")
length(WRKY.genes)
WRKY.activated <- length(intersect(activated.genes, halan[WRKY.genes]))
WRKY.repressed <- length(intersect(repressed.genes, halan[WRKY.genes]))
WRKY.non.diff <- length(WRKY.genes) - WRKY.activated - WRKY.repressed

##Whirly
whirly.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF08536")
length(whirly.genes)
whirly.activated <- length(intersect(activated.genes, halan[whirly.genes]))
whirly.repressed <- length(intersect(repressed.genes, halan[whirly.genes]))
whirly.non.diff <- length(whirly.genes) - whirly.activated - whirly.repressed

##BHLH
BHLH.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00010")
length(BHLH.genes)
BHLH.activated <- length(intersect(activated.genes, halan[BHLH.genes]))
BHLH.repressed <- length(intersect(repressed.genes, halan[BHLH.genes]))
BHLH.non.diff <- length(BHLH.genes) - BHLH.activated - BHLH.repressed

##bZIP
Bzip.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00170")
length(Bzip.genes)
Bzip.activated <- length(intersect(activated.genes, halan[Bzip.genes]))
Bzip.repressed <- length(intersect(repressed.genes, halan[Bzip.genes]))
Bzip.non.diff <- length(Bzip.genes) - Bzip.activated - Bzip.repressed

tfs <- c(AP2.genes, B3.genes, ARR.genes, C2H2.genes, C3H.genes, CPP.genes,
  Dof.genes, E2F.genes, GATA.genes, HB.genes, HSF.genes,LSD.genes, myb.genes,
  NF.Y.genes, Nin.genes, S1FA.genes, SBP.genes, WRKY.genes, whirly.genes,
  BHLH.genes, Bzip.genes)

tfs.data <- matrix(c(myb.repressed, myb.non.diff,myb.activated,
                     ARR.repressed, ARR.non.diff, ARR.activated,
                     ap2.repressed,ap2.non.diff, ap2.activated,
                     C3H.repressed, C3H.non.diff, C3H.activated,
                     NF.Y.repressed,NF.Y.non.diff,NF.Y.activated,
                     CPP.repressed, CPP.non.diff, CPP.activated,
                     Nin.repressed,Nin.non.diff,Nin.activated,
                     SBP.repressed,SBP.non.diff,SBP.activated,
                     E2F.repressed, E2F.non.diff, E2F.activated,
                     BHLH.repressed, BHLH.non.diff, BHLH.activated,
                     Bzip.repressed, Bzip.non.diff, Bzip.activated,
                     C2H2.repressed, C2H2.non.diff, C2H2.activated,
                     GATA.repressed, GATA.non.diff, GATA.activated,
                     HSF.repressed, HSF.non.diff, HSF.activated,
                     LSD.repressed, LSD.non.diff,LSD.activated,
                     WRKY.repressed,WRKY.non.diff,WRKY.activated,
                     b3.repressed, b3.non.diff, b3.activated,
                     Dof.repressed, Dof.non.diff, Dof.activated,
                     HB.repressed, HB.non.diff, HB.activated,
                     S1FA.repressed, S1FA.non.diff,S1FA.activated,
                     whirly.repressed,whirly.non.diff,whirly.activated
                     ),ncol=21)
tfs.data <- tfs.data[,21:1]
tfs.names <- rev(c("MYB", "ARR-B", "AP2", "zf-C3H", "NF-Y", "CPP", "Nin-like",
               "SBP", "E2F/DP", "bHLH", "bZIP", "zf-C2H2", "GATA", "HSF",
               "zf-LSD1", "WRKY", "B3", "Dof", "HB", "S1FA", "Whirly"))

png(filename = "tf_barplot.png")
barplot(tfs.data, horiz = T,names.arg = tfs.names, las=2,cex.lab = 1.8, 
        col=c("blue","gray75","red") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        xlab="Number of Genes",
        )
dev.off()


##Barplot gene
barplot.halan.gene <- function(gene.id,gene.expression)
{
    expression.mean.2mM <- mean(unlist(gene.expression[gene.id,
                                                   c("NO3_2mM_1", "NO3_2mM_2")]))
    expression.mean.15mM <- mean(unlist(gene.expression[gene.id,
                                                    c("NO3_15mM_1", "NO3_15mM_2")]))
    means <- c(expression.mean.15mM,expression.mean.2mM)
    
    expression.sd.2mM <- sd(unlist(gene.expression[gene.id,
                                               c("NO3_2mM_1", "NO3_2mM_2")]))
    expression.sd.15mM <- sd(unlist(gene.expression[gene.id,
                                                c("NO3_15mM_1", "NO3_15mM_2")]))
    sds <- c(expression.sd.15mM,expression.sd.2mM)
    
    
    png(filename = paste0(gene.id,".png"),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("springgreen4","firebrick2"),
                    names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,
                    ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                    main=gene.id,
                    cex.main=2)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
           code = 3,angle=90,lwd=2)
    dev.off()
}

barplot.halan.gene(gene.id="HaLaN_03814",gene.expression=gene.expression)
barplot.halan.gene(gene.id="HaLaN_24769",gene.expression=gene.expression)


## Functions for specific pathways analysis
hlacustris.ko <- select(org.Hlacustris.eg.db,columns = c("GID","KO"),
                        keys=keys(org.Hlacustris.eg.db,keytype = "GID"))

barplot.enzyme.ko <- function(ko,enzyme.name,map.gfh.halan,gene.expression,
                              hlacustris.ko,gene.name)
{
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,KO == ko)$GID],
                                rownames(gene.expression))
  
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],
                                                   c("NO3_2mM_1", "NO3_2mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],
                                                    c("NO3_15mM_1", "NO3_15mM_2")]))
    means <- c(enzyme.mean.15mM,enzyme.mean.2mM)
    
    enzyme.sd.2mM <- sd(unlist(gene.expression[halan.gene.names[i],
                                               c("NO3_2mM_1", "NO3_2mM_2")]))
    enzyme.sd.15mM <- sd(unlist(gene.expression[halan.gene.names[i],
                                                c("NO3_15mM_1", "NO3_15mM_2")]))
    sds <- c(enzyme.sd.15mM,enzyme.sd.2mM)
    
    
    png(filename = paste(c(enzyme.name,i,".png"),collapse=""),width = 300)
    par(lwd=3)
    xpos <- barplot(means,col=c("springgreen4","firebrick2"),
                    names.arg = c("15mM","2mM"),las=2,cex.names = 1.5,
                    ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                    main=paste(gene.name,halan.gene.names[i],sep = " - "),
                    cex.main=2)
    arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
           code = 3,angle=90,lwd=2)
    dev.off()
  }
  
  return(halan.gene.names)
}

heatmap.enzyme.ko <- function(ko,enzyme.name,map.gfh.halan,gene.expression,precision)
{
  pos.colfunc <- colorRampPalette(c("white","firebrick2"))
  pos.colors <- pos.colfunc(precision)
  
  neg.colfunc <- colorRampPalette(c("white","springgreen4"))
  neg.colors <- neg.colfunc(precision)
  
  halan.gene.names <- intersect(map.gfh.halan[subset(hlacustris.ko,KO == ko)$GID],rownames(gene.expression))

  fc.max <- 0
  fc.min <- 20
    
  for(i in 1:length(halan.gene.names))
  {
    enzyme.mean.2mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_2mM_1", "NO3_2mM_2")]))
    enzyme.mean.15mM <- mean(unlist(gene.expression[halan.gene.names[i],c("NO3_15mM_1", "NO3_15mM_2")]))
    fc <- enzyme.mean.2mM/enzyme.mean.15mM

    if(fc >= 1)
    {
      if(fc >= fc.max)
      {
        fc.max <- fc
      }
      
      if(round(fc*10) < precision)
      {
        cols <- pos.colors[round(fc*10)]
      } else
      {
        cols <- pos.colors[precision]
      }
      
    } else
    {
      if(fc <= fc.min)
      {
        fc.min <- fc
      }
      
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
  
  return(list(halan.gene.names,fc.max, fc.min, 1/fc.min))
}


## Gradient
colfunc<-colorRampPalette(c("springgreen4","white","firebrick2"))
plot(rep(1,1000),col=(colfunc(1000)), pch=15,cex=20,xlim=c(0,730),axes=F,
     ylab="",xlab="")

## Carotenoids

barplot.enzyme.ko(ko="K03527",enzyme.name = "-3_HDR", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="HDR",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K03527",enzyme.name = "-3_HDR", map.gfh.halan = halan, 
                  gene.expression = gene.expression, 
                  precision = 78)

barplot.enzyme.ko(ko="K01823",enzyme.name = "-2_IPI", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="IPI",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01823",enzyme.name = "-2_IPI", map.gfh.halan = halan, 
                  gene.expression = gene.expression, 
                  precision = 78)

barplot.enzyme.ko(ko="K13789",enzyme.name = "-1_GGPS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GGPS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K13789",enzyme.name = "-1_GGPS", map.gfh.halan = halan, 
                  gene.expression = gene.expression, 
                  precision = 78)

barplot.enzyme.ko(ko="K02291",enzyme.name = "0_PS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PSY",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K02291",enzyme.name = "0_PS", map.gfh.halan = halan, 
                  gene.expression = gene.expression, 
                  precision = 78)

barplot.enzyme.ko(ko="K02293",enzyme.name = "1_PDS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PDS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K02293",enzyme.name = "1_PDS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00514",enzyme.name = "2_ZDS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ZDS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00514",enzyme.name = "2_ZDS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K06443",enzyme.name = "3_LCYB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="LCYB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K06443",enzyme.name = "3_LCYB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

## Rama epsilon
barplot.enzyme.ko(ko="K06444",enzyme.name = "4_LCYE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="LCYE",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K06444",enzyme.name = "4_LCYE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K15747",enzyme.name = "5_BCHB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="BCHB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K15747",enzyme.name = "5_BCHB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

## Rama beta
barplot.enzyme.ko(ko="K09836",enzyme.name = "6_BKT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="BKT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K09836",enzyme.name = "6_BKT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K15746",enzyme.name = "7_CHYB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="CHYB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K15746",enzyme.name = "7_CHYB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K09838",enzyme.name = "8_ZEP", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ZEP",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K09838",enzyme.name = "8_ZEP", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

## Correct the missing annotation for VDE

subset(hlacustris.ko, GID == "GFH24830")
hlacustris.ko[3913,2] <- "K09839"
subset(hlacustris.ko, GID == "GFH24830")

## Correct the missing annotation for glutaminase
subset(hlacustris.ko, GID == "GFH29177")
hlacustris.ko[8350,2] <- "K01425"
subset(hlacustris.ko, GID == "GFH29177")

## Correct the missing annotation for glutaminase
subset(hlacustris.ko, GID == "GFH11184")
hlacustris.ko[8350,2] <- "K01425"
subset(hlacustris.ko, GID == "GFH29177")

gene.expression[map.gfh.halan["GFH11184"],]
gene.expression["HaLaN_04431",]



barplot.enzyme.ko(ko="K09839",enzyme.name = "9_VDE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="VDE",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K09839",enzyme.name = "9_VDE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

## Function to extract promoter sequence

promoter.sequence <- function(genes.names, genes.info, genome.seqs,file.name)
{
  genes.of.interest <- subset(genes.info, gene_id %in% genes.names)
  
  promoter.seqs <- vector(mode="list",length=nrow(genes.of.interest))
  names(promoter.seqs) <- genes.of.interest$gene_id
  
  for(i in 1:nrow(genes.of.interest))
  {
    current.seq.id <- as.character(genes.of.interest$seqnames[i])
    current.seq <- genome.seqs[current.seq.id]
    current.strand <- as.character(genes.of.interest$strand[i])
    
    if(current.strand == "+")
    {
      tss <- genes.of.interest$start[i]
      if(tss > 1000)
      {
        promoter.seqs[[i]] <-current.seq[[1]][(tss-1000):(tss-1)]  
      } else
      {
        promoter.seqs[[i]] <-current.seq[[1]][1:(tss-1)]
      }
      
    } else if (current.strand == "-")
    {
      tss <- genes.of.interest$end[i]
      promoter.seqs[[i]] <- comp(rev(current.seq[[1]][(tss+1):(tss+1000)]))
    }
  }
  
  write.fasta(sequences = promoter.seqs,
              names = names(promoter.seqs),
              file.out = file.name)
  return(promoter.seqs)
}

## Generation of promoters sequence for key enzymes in carotenoids
## biosynthesis to analysis the presence of TFBS 
carotenoids.genes <- c("HaLaN_05222", "HaLaN_06917", "HaLaN_19207", "HaLaN_29333", "HaLaN_01499", "HaLaN_15035",
                       "HaLaN_05151", "HaLaN_15758", "HaLaN_04875")

promoter.sequence(genes.names = carotenoids.genes, 
                  genes.info = hlacustris.genes , 
                  genome.seqs = hlacustris.seqs,
                  file.name ="carotenoids_promoters.fa")
  
## Pyruvate, Fatty acid metabolism
barplot.enzyme.ko(ko="K01960",enzyme.name = "??_PC", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PC",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00873",enzyme.name = "1_PC", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)




barplot.enzyme.ko(ko="K00873",enzyme.name = "1_PK", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PK",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00873",enzyme.name = "1_PK", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01595",enzyme.name = "?_PEPC", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PEPC",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01595",enzyme.name = "?_PEPC", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K01610",enzyme.name = "?_PEPCK", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PEPCK",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01610",enzyme.name = "?_PEPCK", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00161",enzyme.name = "2_PDHE1a", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PDHE1a",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00161",enzyme.name = "2_PDHE1a", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00162",enzyme.name = "2_PDHE1b", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PDHE1b",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00162",enzyme.name = "2_PDHE1b", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01895",enzyme.name = "2_ACS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01895",enzyme.name = "2_ACS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K15231",enzyme.name = "?_ACLY", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACLY",
                  hlacustris.ko=hlacustris.ko)


barplot.enzyme.ko(ko="K00627",enzyme.name = "2_PDHE2", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="PDHE2",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00627",enzyme.name = "2_PDHE2", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)



barplot.enzyme.ko(ko="K01961",enzyme.name = "3_ACCCB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACCCB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01961",enzyme.name = "3_ACCCB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01962",enzyme.name = "3_ACACA", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACACA",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01962",enzyme.name = "3_ACACA", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01963",enzyme.name = "3_ACACB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACACB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01963",enzyme.name = "3_ACACB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00645",enzyme.name = "4_MCAT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="MCAT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00645",enzyme.name = "4_MCAT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00647",enzyme.name = "5_KASI", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="KASI",
                  hlacustris.ko=hlacustris.ko)

barplot.enzyme.ko(ko="K09458",enzyme.name = "5_KASII", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="KASII",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K09458",enzyme.name = "5_KASII", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00648",enzyme.name = "6_KASIII", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="KASIII",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00648",enzyme.name = "6_KASIII", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00059",enzyme.name = "7_KR", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="KR",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00059",enzyme.name = "7_KR", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01897",enzyme.name = "8_LACS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="LACS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01897",enzyme.name = "8_LACS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K03921",enzyme.name = "9_SAD", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="SAD",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K03921",enzyme.name = "9_SAD", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K22849",enzyme.name = "10_DGAT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="DGAT",
                  hlacustris.ko=hlacustris.ko)


 

## Generation of promoters sequence for key enzymes in fatty acid
## biosynthesis to analysis the presence of TFBS 
fatty.acids.genes <- c("HaLaN_04167", "HaLaN_08198", "HaLaN_00479",
                       "HaLaN_32744", "HaLaN_03364", "HaLaN_25691",
                       "HaLaN_16400", "HaLaN_23803")

promoter.sequence(genes.names = fatty.acids.genes, 
                  genes.info = hlacustris.genes , 
                  genome.seqs = hlacustris.seqs,
                  file.name ="fatty_acids_promoters.fa")

## Starch and sucrose metabolism

barplot.enzyme.ko(ko="K00975",enzyme.name = "1_APL", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="APL",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00975",enzyme.name = "1_APL", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00749",enzyme.name = "2_SSS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="SSS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00749",enzyme.name = "2_SSS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K13679",enzyme.name = "3_GBSS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GBSS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K13679",enzyme.name = "3_GBSS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00700",enzyme.name = "4_SBE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="SBE",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00700",enzyme.name = "4_SBE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00688",enzyme.name = "5_GLGP", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GLGP",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00688",enzyme.name = "5_GLGP", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00705",enzyme.name = "6_DPE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="6_DEP",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00705",enzyme.name = "6_DPE", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K16055",enzyme.name = "1_TPS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="TPS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K16055",enzyme.name = "1_TPS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00036",enzyme.name = "G6PDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="G6PDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00036",enzyme.name = "G6PDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00134",enzyme.name = "GAPDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GAPDH",
                  hlacustris.ko=hlacustris.ko)

heatmap.enzyme.ko(ko="K00134",enzyme.name = "GAPDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)



barplot.enzyme.ko(ko="K00033",enzyme.name = "6PGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="6PGDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00033",enzyme.name = "6PGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


## TCA
barplot.enzyme.ko(ko="K01647",enzyme.name = "1_CS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="CS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01647",enzyme.name = "1_CS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01681",enzyme.name = "2_ACO", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ACO",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01681",enzyme.name = "2_ACO", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00030",enzyme.name = "3_IDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="IDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00030",enzyme.name = "3_IDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00031",enzyme.name = "4_IDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="IDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00031",enzyme.name = "4_IDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00164",enzyme.name = "5_OGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="OGDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00164",enzyme.name = "5_OGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00658",enzyme.name = "6_OGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="OGDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00658",enzyme.name = "6_OGDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K01899",enzyme.name = "7_sucA", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="sucA",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01899",enzyme.name = "7_sucA", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01900",enzyme.name = "8_sucB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="sucB",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01900",enzyme.name = "8_sucB", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K00234",enzyme.name = "9_SDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="SDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00234",enzyme.name = "9_SDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01679",enzyme.name = "10_FUM", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="FUM",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01679",enzyme.name = "10_FUM", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00025",enzyme.name = "11_MDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="MDH",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00025",enzyme.name = "11_MDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

## Nitrogen assimiliation
barplot.enzyme.ko(ko="K02575",enzyme.name = "1_NRT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="NRT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K02575",enzyme.name = "1_NRT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K10534",enzyme.name = "2_NR", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="NR",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K10534",enzyme.name = "2_NR", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K01915",enzyme.name = "3_GS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GS",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K01915",enzyme.name = "3_GS", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00284",enzyme.name = "4_GOGAT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GOGAT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K00284",enzyme.name = "4_GOGAT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K00261",enzyme.name = "5_GDH", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GDH",
                  hlacustris.ko=hlacustris.ko)



barplot.enzyme.ko(ko="K14455",enzyme.name = "ATT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ATT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K14455",enzyme.name = "ATT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)

barplot.enzyme.ko(ko="K14272",enzyme.name = "ALT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="ALT",
                  hlacustris.ko=hlacustris.ko)
heatmap.enzyme.ko(ko="K14272",enzyme.name = "ALT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,precision = 78)


barplot.enzyme.ko(ko="K14272",enzyme.name = "GPT", map.gfh.halan = halan, 
                  gene.expression = gene.expression,gene.name="GPT",
                  hlacustris.ko=hlacustris.ko)
