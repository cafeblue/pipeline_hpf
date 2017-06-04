# load libraries
library(edgeR)

# Args[1]    current running folder
# Args[2]    sampleID
# Args[3]    postprocID
# Args[4]    rawcount.txt

# setup directory with data to normalize
setwd(commandArgs(TRUE)[1])
# read raw data
raw.data <- read.table( commandArgs(TRUE)[4] , header = FALSE, sep="\t", row.names =1 )
# remove no features and ambiguous count
remove <- c("__alignment_not_unique", "__not_aligned", "__too_low_aQual", "__ambiguous", "__no_feature", scan('/hpf/largeprojects/pray/wei.wang/RNASeq-pipeline/RNA-Seq_blacklist.txt', character(), quote = ""), scan('/hpf/largeprojects/pray/wei.wang/RNASeq-pipeline/RNA-Seq_shortgene.txt', character(), quote = ""))
rows_to_remove <- which(row.names(raw.data) %in% remove)
raw.data <- raw.data[-rows_to_remove,]
DGE <- DGEList(counts = raw.data)

#### COMPUTE RPKM #####
# extract read lenght from gtf file
library(GenomicRanges)
library(rtracklayer)
GTFfile = "/hpf/largeprojects/pray/wei.wang/RNASeq-pipeline/Homo_sapiens.GRCh37.75.gtf"
GTF <- import.gff(GTFfile, format="gtf", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)
calc_length <- function(x) {
sum(elementMetadata(x)$widths)
}
gene_length <- as.data.frame(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length))
# extract read length for genes selected in DGE
genelist <- rownames(DGE$counts)
list <- match( genelist , rownames(gene_length) )
gene_length <- as.vector(gene_length[list,])
rpkm <- rpkm(DGE, gene.length=gene_length)
colnames(rpkm) <- c("RPKM", "RPKM1")
rpkm <- rpkm[,c("RPKM")]
# add average RPKM per gene

# sub for tpm calculation
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tpm<- counts_to_tpm(DGE$counts, gene_length, c(100,100))
colnames(tpm) <- c("TPM", "TPM1")
tpm <- tpm[,c("TPM")]
#final_table <- cbind(rpkm, tpm)
finaltab <- as.data.frame(cbind(rpkm, tpm))
library(dplyr)
finaltab <- add_rownames(finaltab, var = 'Gene')
finaltab <- mutate(finaltab, rpkm_rank = dense_rank(desc(rpkm)))
finaltab <- mutate(finaltab, postprocID = commandArgs(TRUE)[3])
finaltab <- finaltab[,c(5,1,2,4,3)]

outputfilename <-  paste( commandArgs(TRUE)[2], commandArgs(TRUE)[3], "RPKM_TPM.txt", sep=".") 
write.table( finaltab , file = outputfilename , sep = "\t" , quote = FALSE, row.names = FALSE , col.names = FALSE )
