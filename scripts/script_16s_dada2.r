library(dada2); packageVersion("dada2")
library(glue)
library(ShortRead)
library(jsonlite)


# simplifyVector = TRUE turns JSON arrays into R vectors automatically
config <- read_json("./config.json", simplifyVector = TRUE)

# Access your parameters using the $ operator
n_threads <- config$threads
max_ee    <- config$max_ee
trunc_l   <- config$trunc_len

# ## providing parameters FROM docker compose
# threads   <- as.numeric(Sys.getenv("THREADS", unset = 4))
# max_ee    <- as.numeric(Sys.getenv("MAX_EE", unset = 2))
# TRUNC_F <- as.numeric(Sys.getenv("TRUNC_F", unset = 220))
# TRUNC_R <- as.numeric(Sys.getenv("TRUNC_R",unset = 220))

db_path   <- Sys.getenv("REF_DB", unset = "/ref_db/silva.fa")


### Function: get_adaptive_seqtab
# seqtab: The sequence table (typically seqtab.nochim)
# buffer: How many base pairs to allow on either side of the peak (default 15)
# get_adaptive_seqtab <- function(seqtab, buffer = 15) {
#   # 1. Get lengths of all unique sequences
#   seq_lens <- nchar(getSequences(seqtab))
#   # 2. Find the most frequent length (the biological peak)
#   len_table <- table(seq_lens)
#   mode_len <- as.numeric(names(len_table)[which.max(len_table)])
#   # 3. Define the adaptive boundaries
#   min_len <- mode_len - buffer
#   max_len <- mode_len + buffer
#   # 4. Filter the table
#   # colnames(seqtab) are the sequences themselves
#   keep_indices <- nchar(colnames(seqtab)) >= min_len & nchar(colnames(seqtab)) <= max_len
#   seqtab.filt <- seqtab[, keep_indices]
#   # 5. Report to console
#   cat(paste0("--- Adaptive Filtering ---\n",
#              "Biological Peak (Mode): ", mode_len, " bp\n",
#              "Range Applied: ", min_len, " to ", max_len, " bp\n",
#              "Sequences Kept: ", ncol(seqtab.filt), " (out of ", ncol(seqtab), ")\n"))
  
#   return(seqtab.filt)
# }


print(paste("Running with", n_threads, "threads"))
start_time <- Sys.time()
out_dir <- "../output"
##Setting up namings and paths
# File parsing
print("Collecting for and rev reads")
pathF <- "../forward_reads" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "../reverse_reads" # CHANGE ME ...

head(file.path(pathF))
head(file.path(pathR))

filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

pattern_fwd <- "_R1_001\\.fastq\\.gz"
pattern_rv <- "_R2_001\\.fastq\\.gz"

# More robust pattern for your specific filenames
fastqFs <- sort(list.files(pathF, pattern=pattern_fwd, full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern=pattern_rv, full.names = TRUE))


if (length(fastqFs) == 0) {
  stop(paste("Errror no file found in:", path_data, "pattern:", pattern_fwd))
} else {
  message(paste("Successo: Trovati", length(fnFs), "file per l'analisi."))
}

print(fastqFs[1])
print('Checkpoint trimming')
plotQualityProfile(fastqFs[1:2]) # Forward sequences
plotQualityProfile(fastqRs[1:2]) # reverse sequences

# Extract sample names, assuming filenames have format: SAMPLENAME_type_XXX.fastq
parts <- strsplit(basename(fastqFs), "_")
# Extract piece 1, extract piece 2, and paste them with "_"
sample.names1 <- sapply(parts, function(x) paste(x[1], x[2], sep="_"))

# Place filtered files in path1/filtered/ subdirectory
filtFs1 <- file.path(filtpathF, paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(filtpathR, paste0(sample.names1, "_R_filt.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
#The V3V4 amplicon with primers sequences is 440-462nts long, 
#so the two truncation lengths must add up to at least 462 + 12 nts to get proper merging.
#######################################################################################
out_check <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(TRUNC_F, TRUNC_R) 
              maxEE=c(2,2),
              trimLeft=c(17,21),   #ASSUMING you are using  Illumina V3-V4
              maxN=0,
              # minQ=20, # Consider increasing to remove low-quality reads
              rm.phix=TRUE,
              compress=TRUE, 
              verbose=TRUE, 
              multithread=n_threads)
########################################################################################
head(out_check)

# File parsing
filtpathF <- "../forward_reads/filtered/" # CHANGE ME to the directory containing your filtered forward fastqs
filtpathR <- "../reverse_reads/filtered/" # CHANGE ME ...

filtFs <- list.files(filtpathF, pattern="_F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="_R_filt.fastq.gz", full.names = TRUE)

parts_f <- strsplit(basename(filtFs), "_")
parts_r <- strsplit(basename(filtRs), "_")

# Extract piece 1, extract piece 2, and paste them with "_"
sample.names_for <- sapply(parts_f, function(x) paste(x[1], x[2], sep="_"))
sample.names_rev <- sapply(parts_r, function(x) paste(x[1], x[2], sep="_"))

print(sample.names_for)
print(sample.names_rev)

# sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
# sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names_for, sample.names_rev)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names_for
names(filtRs) <- sample.names_for
set.seed(100)

# IF ERR models already exists, it load them from output dir, 
# if you wish to re-compute those, be sure to remove the same files 
# before running the pipeline
# Construct file names using the variable

errF_file <- file.path(out_dir, "errF_model.rds")
errR_file <- file.path(out_dir, "errR_model.rds")

if (file.exists(errF_file) && file.exists(errR_file)) {
    print("Loading existing error models...")
    errF <- readRDS(errF_file)
    errR <- readRDS(errR_file)
} else {
    print("Learning Errors...")
    errF <- learnErrors(filtFs, nbases=1e8, multithread=n_threads)
    errR <- learnErrors(filtRs, nbases=1e8, multithread=n_threads)
    saveRDS(errF, errF_file)
    saveRDS(errR, errR_file)
}

# creating Abundance table
seqtab_file <- file.path(out_dir, "seqtab.rds")

if (file.exists(seqtab_file)) {
    cat("Loading existing sequence table:", seqtab_file, "\n")
    seqtab <- readRDS(seqtab_file)
} else {
    # Sample inference and merger of paired-end reads
    mergers <- vector("list", length(sample.names_for))
    names(mergers) <- sample.names_for

    for(sam in sample.names_for) {
        cat("Processing:", sam, "\n")
        # Direct path usage is more memory efficient than derepFastq
        ddF <- dada(filtFs[[sam]], err=errF, multithread=n_threads, verbose=FALSE)
        ddR <- dada(filtRs[[sam]], err=errR, multithread=n_threads, verbose=FALSE)
        
        mergers[[sam]] <- mergePairs(ddF, filtFs[[sam]], ddR, filtRs[[sam]])
    }
    # Construct sequence table
    seqtab <- makeSequenceTable(mergers)
    # Save using the variable
    saveRDS(seqtab, seqtab_file)
    cat("Sequence table saved to:", seqtab_file, "\n")
}

#seqtab.filt <- get_adaptive_seqtab(seqtab, buffer = 15)
#print('How many?',dim(seqtab))
cat("How many samples and ASVs?", dim(seqtab), "\n")
#dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
#setwd(file.path("./", "output"))
#####################################################################
# Merge multiple runs (if necessary)
#st1 <- readRDS("path/to/run1/output/seqtab.rds")
#st2 <- readRDS("path/to/run2/output/seqtab.rds")
#st3 <- readRDS("path/to/run3/output/seqtab.rds")
#st.all <- mergeSequenceTables(st1, st2, st3)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads)
# Assign taxonomy
tax <- assignTaxonomy(seqtab.nochim, "../database/silva_nr99_v138.2_toSpecies_trainset.fa", multithread=threads)
date <- format(Sys.Date(), "%Y-%m-%d")

results <- "../results"
# Write to disk
table_name <-glue("seq_tab_{date}_dada2.rds")
tax_final_ <- glue("tax_final_{date}.rds")

taxa_table <- file.path(results, tax_final_)
seq_table <- file.path(results, table_name)

saveRDS(seqtab.nochim, seq_table) 
saveRDS(tax, taxa_table)

# --- END OF SCRIPT ---
end_time <- Sys.time()
execution_time <- end_time - start_time

cat("\n--- Execution Summary ---\n")
cat("Started at:", as.character(start_time), "\n")
cat("Finished at:", as.character(end_time), "\n")
print(execution_time)