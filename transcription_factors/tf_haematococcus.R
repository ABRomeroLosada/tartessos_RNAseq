##############Read data
haematococcus.annotation <- read.table(file = "haematoccocus_functional_annotation.csv",header = T,as.is = T,fill = T)
head(haematococcus.annotation)
gene.expression <- read.table(file = "./haematococcus_gene_expression.tsv",header = T,sep = "\t")
head(gene.expression)
lowN <- rowMeans(gene.expression[,1:2])
highN <- rowMeans(gene.expression[,3:4])


######### Correspondence between genome notation and NCBI notation.
halan.gfh.correspondence <- read.table(file="./halaN_GFH_correspondence/gfh_halan_correspondence.tsv",sep="\t",as.is=T)
head(halan.gfh.correspondence)

halan <- halan.gfh.correspondence$V2
gfh <- halan.gfh.correspondence$V1
names(halan) <- gfh
names(gfh) <- halan

extract.tf <- function(annotation, pfam.id)
{
  return(subset(annotation,PFAM == pfam.id)$GENES)
}


##############Steps to obtain DEGs list
## Differential expression analysis
library(limma)
library(ballgown)

## Log2 transformation
log.gene.expression <- log2(gene.expression+1)
head(log.gene.expression)
nrow(log.gene.expression)


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

genes.ids <- rownames(degs.results)
names(genes.ids) <- genes.ids


######Manual extraction of genes names corresponding to the different TFs.
## AP2 
AP2.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00847")
genes.ids[halan[AP2.genes]]

## B3 PF02362
B3.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02362")
genes.ids[halan[B3.genes]]


##ARR-B arabidopsis response regulator con MYB motifs PF00249 y PF00072.
##MYB es lo mismo que ARR-B??

ARR.genes.1 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00249")
ARR.genes.2 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00072")
genes.ids[halan[ARR.genes.1]]
genes.ids[halan[ARR.genes.2]]


##C2H2-like zinc fingers PF00096 (no sale ninguno), PF12756 se supone que es lo mismo pero doble
C2H2.genes.1 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00096")
C2H2.genes.2 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF12756")
genes.ids[halan[C2H2.genes.1]] #vacio
genes.ids[halan[C2H2.genes.2]]


##C3H-like zinc finder PF00097 y otro2 dos con 3 y 4 copias (PF13920 y PF15227)
C3H.genes.1 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00097")
C3H.genes.2 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF15227")
C3H.genes.3 <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF13920")
genes.ids[halan[C3H.genes.1]] 
genes.ids[halan[C3H.genes.2]]
genes.ids[halan[C3H.genes.3]]


## Co-like (las que tengan CCT motifs?) PF06203
CO.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF06203")
genes.ids[halan[CO.genes]] 


##CPP-like (family TCR??), lo que encuentro es TSO1-like CXC domain (coincide con lo que dice tair) 
##PF03638
CPP.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03638")
genes.ids[halan[CPP.genes]] 


##Dof zinc fingers PF02701
Dof.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02701")
genes.ids[halan[Dof.genes]] 

##E2F/DP family PF02319
E2F.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02319")
genes.ids[halan[Dof.genes]] 

##AP2/ERF el pfam solo menciona a AP2 PF00847
AP2.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00847")
genes.ids[halan[AP2.genes]] 

##G2-like no lo encuentro

##GATA PF00320
GATA.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00320")
genes.ids[halan[GATA.genes]] 

##HB homeobox PF00046
HB.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00046")
genes.ids[halan[HB.genes]] #vacio

##HSF heat shock transciption factors PF00447
HSF.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00447")
genes.ids[halan[HSF.genes]] #vacio

##LSD1 zinc fingers PF06943
LSD.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF06943")
genes.ids[halan[LSD.genes]] 

##M_type_MADS no lo encuentro.

##NF-X1 zinc fingers PF01422
NF.X1.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF01422")
genes.ids[halan[NF.X1.genes]] #vacio

##NF-Y (incluye a YA, YB e YC aunque chlamy solo tiene B y C) factor PF00808
NF.Y.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00808")
genes.ids[halan[NF.Y.genes]] 

## Nin like proteins: nodule inception proteins. En pfam solo encuentro el dominio que 
##se relaciona con la regulación de genes según la cantidad de nitrogeno. RWP-RK (PF02042)
Nin.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF02042")
genes.ids[halan[Nin.genes]] #vacio?


##S1FA (PF04689)
S1FA.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF04689")
genes.ids[halan[S1FA.genes]] 

## SBP (PF03110)
SBP.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03110")
genes.ids[halan[SBP.genes]] 

##WRKY (PF03106)
WRKY.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF03106")
genes.ids[halan[WRKY.genes]] 

##Whirly (PF08536)
whirly.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF08536")
genes.ids[halan[whirly.genes]] 

##YABBY (PF04690)
YABBY.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF04690")
genes.ids[halan[YABBY.genes]] #VACIO

##BHLH (PF00010)
BHLH.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00010")
genes.ids[halan[BHLH.genes]] 

##bZIP (PF00170)
Bzip.genes <- extract.tf(annotation = haematococcus.annotation,pfam.id = "PF00170")
genes.ids[halan[Bzip.genes]]


##Read the final FT DEGS list
ft.degs <- read.csv(file = "GENES_DEGS_TF.csv", as.is = T, header=F)
head(ft.degs)


#Generate barplots
for(i in 1:nrow(ft.degs))
{
  current.gene <- ft.degs[i,1] #gene.name
  name <- strsplit(ft.degs[i,2],split=",")[[1]] #gene id
   
        sds <- c(sd(gene.expression[current.gene,3:4]),
                 sd(gene.expression[current.gene,1:2]))
        means <- c(highN[current.gene],lowN[current.gene])
        png(filename = paste(name, "_", current.gene,".png",sep = ""))
        xpos <- barplot(means,ylim = c(0,1.1*max(means+sds)),col=c("blue","red"),main=name,names.arg = c("12mM","4mM"),ylab="FPKM")
        arrows(x0 = xpos,y0 =means-sds,x1 = xpos,y1 = means+sds,code=3,angle=90,length=0.1,lwd=2 )
        dev.off()
    
    }
graphics.off() 
