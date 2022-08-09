# 2022 07 26 created
# processing of PLATE-seq data to get the gene counts

library(tximport)
library(ggplot2)
library(DESeq2)
library(openxlsx)


setwd("D:/pCloud Sync/Columbia/SchwabeLab/HepSC/RSlab/Experimental/DrugScreen")
# *****************************************************cutadapt 
# read in the barcode information
barcodes = read.xlsx("data/20220520_validation/Plateseq/20220616_musHSC_val/Schwabe 1.xlsx")
head(barcodes)
del.seq = "GTTCAGAGTTCTACAGTCCGACGATC"
nchar(del.seq)


# # obtain demultiplexing file for 'sabre'
# tb.new = data.frame(barcode=gsub(del.seq,"",barcodes$Sequence),
#                     out1=paste("plate1_",barcodes$Well.Position,"_R1.fastq",sep=""),
#                     out2=paste("plate1_",barcodes$Well.Position,"_R2.fastq",sep=""))
# write.table(tb.new, file="data/plate1_barcodes_sabre_formatted.txt",sep="\t",quote=F,row.names = F,col.names = F)


# ## Sabre on Linux
# ### using 'sabre' on Linux *** The barcodes will be stripped from the reads and the quality values of the barcode bases will also be removed.
# 
# sabre pe -m 2 -f Card_A_S7_L003_R1_001.fastq.gz -r Card_A_S7_L003_R2_001.fastq.gz -b plateA-barcodes_py-sabre_formatted.txt -u plateA_unknown_barcodes_R1.fastq -w plateA_unknown_barcodes_R2.fastq



# obtain the demultiplexing file in fastq format for 'cutadapt'
wells <- paste0(">",barcodes$Well)
barcodesTrim <- gsub(del.seq,"",barcodes$Sequence)

fastq <- cbind(c(rbind(paste0(">",barcodes$Well),gsub(del.seq,"",barcodes$Sequence))))
head(fastq)

# write.table(fastq, file="data/20220520_validation/Plateseq/20220616_musHSC_val/barcodes_trimmed.fasta",quote=F,row.names = F,col.names = F)

# ## cutadapt command

system(cutadapt --error-rate 0.061 --no-indels --cores 10 --action=trim -g file:barcodes_trimmed.fasta -o trimmed-{name}.R1.fastq.gz -p trimmed-{name}.R2.fastq.gz Schwabe_S7_L003_R1_001.fastq.gz Schwabe_S7_L003_R2_001.fastq.gz)

# explanation:
# --error-rate 0.061 because 0.061*33=2.01; so two mismatches were allowed for a 33nt long trimmed adapter
# -g file:barcodes_trimmed.fasta is the file containing the trimmed adapter in fastq format




# *****************************************************kallisto
# run Kallisto
# get the list of available fastq files
list.dirs(".")
tmp <- list.files("./fastq")
tmp
tmp <- gsub("trimmed-","",tmp)
tmp <- gsub(".R1.fastq.gz","",tmp)
tmp <- gsub(".R2.fastq.gz","",tmp)
tmp <- unique(tmp)
tmp
# tmp <- tmp[-1*c(1:140)] #removing already performed wells
# tmp_wells <- c(unlist(strsplit("A1,B11,B12,B15,B3,C5,D17,G6,H16,I11,I13,I8,K17,L3,N19,N3,O12",split = ",")))
# tmp <- setdiff(tmp,tmp_wells)
# tmp <- c("F7","F8","F9",paste0("G",seq(1:24)))
# tmp <- c(paste0("H",seq(1:24)))
# tmp <- c(paste0("I",seq(1:24)),paste0("j",seq(1:24)))
# tmp <- c(paste0("K",seq(1:24)),paste0("L",seq(1:24)),paste0("M",seq(1:24)))
# tmp <- c(paste0("N",seq(1:24)),paste0("O",seq(1:24)),paste0("P",seq(1:24)))
# tmp <- setdiff(tmp,tmp_wells)

for(i in tmp){
  print(i)
  command <- paste0("kallisto quant -i transcriptome_mus.idx -o counts_R2/output_R2_",i," -b 100 --single -l 180 -s 20 --threads=10 --single-overhang fastq/trimmed-",i,".R2.fastq.gz")
  # if (file.exists(paste0("fastq/trimmed-",i,".R2.fastq.gz"))){print(i)}
  # print(command)
  system(command)
}

# *******************************************************tximport
# use tximport to get the gene level estimates from the transcript level estimates 
# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html


tmp <- list.dirs("./counts_R2")
tmp
dirs <- tmp[-1]
dirs
wells <- gsub("./counts_R2/output_R2_","",dirs)
paste0(dirs,"/abundance.h5")
files <- paste0(dirs,"/abundance.h5")
names(files) <- wells
files
# files <- 'output_A1/abundance.tsv'
# files <- 'output_A1/abundance.h5'

# abundance <- read.delim("D:/kallisto/Plateseq/output_A1/abundance.tsv")
# head(abundance)
transcripts_to_genes <- read.delim("transcripts_to_genes.txt", header=FALSE)
head(transcripts_to_genes)
tx2gene <- data.frame(transcripts=transcripts_to_genes[,1],genes=transcripts_to_genes[,3])
head(tx2gene)
# tmp_pos <- match(abundance$target_id,tx2gene$transcripts)
# length(which(is.na(tmp_pos)))
# abundance[which(is.na(tmp_pos)),]

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi.kallisto.tsv)
head(txi.kallisto.tsv$counts) #the gene counts
sum(txi.kallisto.tsv$counts)
head(txi.kallisto.tsv$abundance) #TPM probably and not correct here for the 3' data
sum(txi.kallisto.tsv$abundance)
head(txi.kallisto.tsv$infReps)
head(txi.kallisto.tsv$length)
head(txi.kallisto.tsv$countsFromAbundance)

# write.table(txi.kallisto.tsv$counts, file = "counts_R2/genecounts_R2.tsv", sep = "\t")
remove(txi.kallisto.tsv)