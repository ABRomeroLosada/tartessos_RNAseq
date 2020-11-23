library("ape")

gff <- read.gff(file = "haematococcus_lacustris.gtf")
nrow(gff)
dim(gff)

gff <- subset(gff, (type == "gene"))
nrow(gff)

row.type <- as.vector(gff$type)
row.attributes <- as.vector(gff$attributes)

gene.ids <- vector(mode = "character",length=nrow(gff))
gene.names <- vector(mode = "character",length=nrow(gff))
for(i in 1:nrow(gff))
{
    gene.ids[i] <- strsplit(strsplit(row.attributes[i],split=";")[[1]][1],split=" ")[[1]][2]
    gene.names[i] <- strsplit(strsplit(row.attributes[i],split=";")[[1]][3],split=" ")[[1]][3]
}

write.table(x = data.frame(gene.names,gene.ids),file = "gfh_halan_correspondence.tsv",quote = F,sep = "\t",row.names = F,col.names = F)
