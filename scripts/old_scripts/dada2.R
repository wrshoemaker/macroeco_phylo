
library(dada2); packageVersion("dada2")
#update("dada2")

path <- "/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/fastq"
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# trim primers
#noPrimeFs <- file.path(path, "noprimers", basename(fnFs))
#noPrimeRs <- file.path(path, "noprimers", basename(fnRs))
#Fs_no_primer <- dada2::removePrimers(fnFs, noPrimeFs, primer.fwd='GTGCCAGCMGCCGCGGTAA', primer.rev='GGACTACHVGGGTWTCTAAT', orient=TRUE)
#Rs_no_primer <- dada2::removePrimers(fnRs, noPrimeRs, primer.fwd='GTGCCAGCMGCCGCGGTAA', primer.rev='GGACTACHVGGGTWTCTAAT', orient=TRUE)


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_2_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=200,
#                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,225), trimLeft=20,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


#truncLen=250,

head(out)

# learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

# depreplication
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)



#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)




#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/Users/williamrshoemaker/GitHub/strain_macroecology/data/poyet/poyet_16s/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

# make species level assignments based on exact matching between ASVs and sequenced reference strains
taxa <- addSpecies(taxa, "/Users/williamrshoemaker/GitHub/strain_macroecology/data/poyet/poyet_16s/silva_species_assignment_v138.1.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# export site-by-species
write.table(t(seqtab.nochim), "/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
# export taxonomy
write.table(taxa, "/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim-taxa.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
