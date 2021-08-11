# Loading the library
library(dada2)

#Loading library
library(dada2)

# Setting seed
set.seed(20233083)

# Setting filename
raw_files <- file.path("AbsoulteRaw")
filtered_path_f <- file.path("filtered/fwd")
filtered_path_r <- file.path("filtered/rev")

# Reading filename
fns <- sort(list.files(raw_files, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# Sorting file-names
fns <- sort(list.files(raw_files, full.names = TRUE))

# Setting path for filtered reads
filtFs <- file.path(filtered_path_f, basename(fnFs))
filtRs <- file.path(filtered_path_r, basename(fnRs))

# Trimming
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
                c(filtFs[[i]], filtRs[[i]]),
                trimLeft=15,
                truncLen=245,
                maxN = 0, # Removing Ns
                maxEE = 2, # Maximum expected errors, using EE = sum(10^(-Q/10))
                multithread=TRUE, # Parallel compute 
                matchIDs=TRUE,
                truncQ=2,
                compress=TRUE)
}


# Declaring paths
filtpathF <- "filtered/fwd"
filtpathR <- "filtered/rev"

# Sort
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Taking base name
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

# Re-name
names(filtFs) <- sample.names
names(filtRs) <- sample.namesR

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
saveRDS(errF, "RDS_Objects/errF.RDS")

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
saveRDS(errR, "RDS_Objects/errR.RDS")

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

# Sample Inference
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}

saveRDS(derepF,"RDS_Objects/derepF.RDS")
saveRDS(ddF,"RDS_Objects/ddF.RDS")
saveRDS(derepR,"RDS_Objects/derepR.RDS")
saveRDS(ddR,"RDS_Objects/ddR.RDS")
saveRDS(mergers,"RDS_Objects/mergers.RDS")

# Accuracy
#mergers <- readRDS("RDS_Objects/mergers.RDS")
#seqtab.nochim <- removeBimeraDenovo(mergers, multithread=TRUE)
#unqs.mock <- as.data.frame(seqtab.nochim$Mock)$sequence
#unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
#cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
#mock.ref <- getSequences(file.path("mock_16s.fasta"))
#match.ref <- sum(sapply(unqs.mock, function(x) any(grepl(x, mock.ref))))
#cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

# Make Sequence Table
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
saveRDS(seqtab.all ,"RDS_Objects/seqtab.all.RDS")

# Removing Chimera
seqtab <- removeBimeraDenovo(seqtab.all,multithread=TRUE)
saveRDS(seqtab,"RDS_Objects/seqtab.RDS")

# writing downstream files
write.table(t(seqtab), "output/ASV_abundance_raw.tsv",  quote=FALSE, sep='\t', col.names = NA)

asvMap <- as.data.frame(rownames(t(seqtab)))
write.table(asvMap, "output/Map_silva_raw.tsv", quote=FALSE, sep='\t', col.names = NA)

# Silva Taxonomy
asv_silva_tax <- assignTaxonomy(seqtab,"silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
saveRDS(asv_silva_tax, "RDS_Objects/asv_silva_tax.RDS")
write.table(asv_silva_tax, "output/asv_silva_tax_raw.tsv", quote=FALSE, sep='\t', col.names = NA)
