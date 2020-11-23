install.packages("TxDb.Hlacustris.NCBI",repos = NULL)

library(TxDb.Hlacustris.NCBI)
library(seqinr)


hlacustris.info <- read.fasta(file = "haematococcus_lacustris.fa",seqtype = "DNA")
hlacustris.seqs <- getSequence(hlacustris.info)
names(hlacustris.seqs) <- getName(hlacustris.info)

txdb <- TxDb.Hlacustris.NCBI

hlacustris.genes <- as.data.frame(genes(txdb))
head(hlacustris.genes)
carotenoids.genes <- c("HaLaN_05222", "HaLaN_06917", "HaLaN_19207", "HaLaN_29333", "HaLaN_01499", "HaLaN_15035",
                       "HaLaN_05151", "HaLaN_15758", "HaLaN_04875")
gluthatione.genes <- c("HaLaN_11988", "HaLaN_06758", "HaLaN_06760", "HaLaN_22187", "HaLaN_30166")

piruvate.genes <- c("HaLaN_13453", "HaLaN_06731", "HaLaN_03083")

genes.of.interest <- piruvate.genes 

hlacustris.genes.of.interest <- subset(hlacustris.genes, gene_id %in% genes.of.interest)
hlacustris.genes.of.interest

promoter.seqs <- vector(mode="list",length=nrow(hlacustris.genes.of.interest))
names(promoter.seqs) <- hlacustris.genes.of.interest$gene_id

for(i in 1:nrow(hlacustris.genes.of.interest))
{
  current.seq.id <- as.character(hlacustris.genes.of.interest$seqnames[i])
  current.seq <- hlacustris.seqs[current.seq.id]
  current.strand <- as.character(hlacustris.genes.of.interest$strand[i])

  if(current.strand == "+")
  {
    tss <- hlacustris.genes.of.interest$start[i]
    if(tss > 1000)
    {
      promoter.seqs[[i]] <-current.seq[[1]][(tss-1000):(tss-1)]  
    } else
    {
      promoter.seqs[[i]] <-current.seq[[1]][1:(tss-1)]
    }
    
  } else if (current.strand == "-")
  {
    tss <- hlacustris.genes.of.interest$end[i]
    promoter.seqs[[i]] <- comp(rev(current.seq[[1]][(tss+1):(tss+1000)]))
  }
}

write.fasta(sequences = promoter.seqs,names = names(promoter.seqs),file.out = "piruvate_promoters.fa")
