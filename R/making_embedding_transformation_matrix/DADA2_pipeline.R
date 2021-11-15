# DADA2 pipeline

# get libraries

library(dada2); packageVersion("dada2")

path <- "fastq_files" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


### Filter and Trim
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


### Perform Filtering and Trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))



#FILTER FORWARD AND REVERSE READS
print("FILTERING FORWARD AND REVERSE READS")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(140,140),
             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=FALSE, multithread=FALSE)

head(out)


### Learn the Error Rates
print("Learning the Error Rates")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


### Dereplication
print("DEREPLICATE THE FILTERED FASTQ FILES")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


### Sample Inference
print("Sample Inference: infering the sequence variants in each sample")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
# Inspecting the dada-class object returned by dada:
dadaFs[[1]]


### Merge Paired Reads
print("**Merge the denoised forward and reverse reads:**")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


### Construct Sequence Table
print("Construct Sequence Table")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


### Remove Chimeras
print("removing chimeras")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



### Track Reads Through the Pipeline
print("tracking reads through pipeline")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)



### Save Final Sequence Table
print("saving final sequence table")
write.csv(seqtab.nochim, "seqtab_nochim.csv", row.names=TRUE) # change the file name to something of your choosing


