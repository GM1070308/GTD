library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DGEobj.utils, lib='/well/ludwig/users/qhr865/R/4.2/skylake/')

####################### EDIT THE BULK RNA FILE #####################################

#Open the bulk RNA file
gtd <- read.table('GSE135727_rna.processed.GTC_log2.txt', header = T)

#Remove entrez symbols from the Gene column - only symbol remains
gtd$Gene <- sub("[|].*","",gtd$Gene)

#Reverse log transformation
gtd[2:11] <- as.data.frame(lapply(gtd[2:11], function(x) 2^x))
gtd[2:11] <- as.data.frame(lapply(gtd[2:11], function(x) x-1))

#Set gene codes as row names & remove Gene column
gtd <- aggregate(. ~ Gene, gtd, mean)
#rownames(gtd) <- gtd$Gene
#gtd <- gtd[-1]

#De-normalise gtd data
#norm_counts <- as.matrix(gtd)
#upper_quartile_counts <- apply(norm_counts, 2, quantile, probs = 0.75)
#denorm_counts <- norm_counts * upper_quartile_counts

#Convert counts to integers
#integer_counts <- round(denorm_counts)

#Read count table from CHM
chm <- read.table('/well/ludwig/users/qhr865/gtd/public_data/transcriptomics/bulk/chm/GSE138250_read_count_table.txt', header = T)

#Convert symbols to ensembl IDs

symbols <- rownames(gtd)

conversion <- select(org.Hs.eg.db,
                     keys = symbols,
                     column = "ENSEMBL",
                     keytype = "SYMBOL")
colnames(conversion) <- c('Gene', 'ENSEMBL_ID')

#Remove genes that don't have a matching ensembl id and remove duplicated rows
conversion <- conversion[!is.na(conversion$ENSEMBL_ID),]
conversion <- conversion[!duplicated(conversion$ENSEMBL_ID),]

#Convert symbols to ensembl IDs
gtd_ensembl <- merge(x=gtd, y=conversion, by = 'Gene')
gtd_ensembl <- gtd_ensembl[-1]

#Merge gtd & chm dataframes

merged_df <- merge(x=gtd_ensembl, y = chm, by = 'ENSEMBL_ID')

#Load samples from 39 healthy placentas
h39 <- read.csv('/well/ludwig/users/qhr865/gtd/public_data/transcriptomics/bulk/first_trim_healthy/GSE109082_genecounts.txt',
                sep = '\t')
h39$ENSEMBL_ID <- h39$X

#Merge again

merged_df <- merge(x=merged_df, y=h39, by = 'ENSEMBL_ID')
merged_df <- merged_df[,-which(names(merged_df) == 'X')]

#Add clearer sample names

colnames(merged_df) <- c("ENSEMBL_ID",'PSTT2', 'CC3', 'ETT1', 'PSTT8', 'ETT9', 'PSTT5', 'PSTT7', 'ETT10','ETT4', 
                         'CC6', 'CHM6', 'CHM11', 'H1', 'CHM17', 'CHM21', paste0("H",seq(2,40)))

merged_df[2:55] <- round(merged_df[2:55])

#Find gene length of each gene

library(biomaRt)

# Define list of ENSEMBL IDs
ensembl_ids <- rownames(merged_df)

# Connect to the appropriate BioMart database
ensembl <- useMart("ensembl", host = "https://useast.ensembl.org")

# Specify the dataset and attributes to retrieve
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
attributes <- c("ensembl_gene_id", "end_position", "start_position")

# Get the gene length for each ENSEMBL ID
gene_info <- getBM(attributes = attributes, filters = "ensembl_gene_id", values = ensembl_ids, mart = dataset)

gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

#Sort gene_info and df so that ENSEMBL ids are in the same order

merged_df <- sort(merged_df, ENSEMBL_ID)
gene_info <- sort(gene_info, gene_length)
lengths <- as.vector(gene_info$gene_length)

rownames(merged_df) <- merged_df$ENSEMBL_ID
merged_df <- merged_df[-1]

merged_tpm <- convertCounts(as.matrix(merged_df), unit = 'TPM', geneLength = lengths,
                            log = FALSE, normalize = 'none')

#Save as a tab-delimited text file of raw counts
write.table(merged_df, "bk2.txt", quote = F, sep = " ")

########################### EDIT SINGLE CELL FILES ############################

#In the terminal, run 

#awk -F'\t' '{sub(/_.*/, "", $1)}1' raw_data_10x.txt > tmpfile.txt && mv tmpfile.txt raw_data_10x.txt
#awk -F'\t' '{sub(/_.*/, "", $1)}1' raw_data_ss2.txt > tmpfile.txt && mv tmpfile.txt raw_data_ss2.txt

#This will remove ENCODE IDs from the Gene column - only gene codes remain

#Also run

#sed -i '1s/\t/ /g' raw_data_10x.txt
#sed -i '1s/\t/ /g' raw_data_ss2.txt

#This step is required because the header row in raw files is separated by tab, 
                #whereas counts are separated by space
#This step makes sure that the header row is also delimited by space

#---------------------- Get the reference -----------------------------#
library(data.table)
raw_10x <- fread('raw_data_10x.txt', sep = '\t', header = T, fill = TRUE)
raw_10x$Gene <- sub('.*_', '', raw_10x$Gene)

ensembl_ids <- raw_10x$Gene

gene_info <- getBM(attributes = attributes, filters = "ensembl_gene_id", values = ensembl_ids, mart = dataset)
gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

#Only retain ensembl ids that have been mapped
raw_10x <- filter(raw_10x, Gene %in% gene_info$ensembl_gene_id)

#Sort raw_10x and gene_info so that ensembl ids are in the same order
raw_10x <- raw_10x[order(raw_10x[, 'Gene']), ]
gene_info <- gene_info[order(gene_info[, 'ensembl_gene_id']),]
identical(raw_10x$Gene, gene_info$ensembl_gene_id)

lengths <- gene_info$gene_length
raw_10x <- as.matrix(raw_10x, rownames = 1)

#Convert counts to TPM
raw_tpm <- convertCounts(raw_10x, unit = 'TPM', geneLength = lengths,
                         log = FALSE, normalize = 'none')
raw_tpm <- t(raw_tpm)

#Read metadata file
meta <- read.csv('meta_10x.txt', sep = '\t', header = TRUE)
meta <- meta[4]
meta$annotation <- sub(" ", "_", meta$annotation)

#Add annotation column to the big matrix
raw_tpm <- merge(x=raw_tpm, y=meta, by = 'row.names')
raw_tpm <- raw_tpm[-1]

raw_tpm_mean <- aggregate( . ~ annotation, raw_tpm, mean)
raw_tpm_mean <- column_to_rownames(raw_tpm_mean, 'annotation')
raw_tpm_mean <- t(raw_tpm_mean)
saveRDS(raw_tpm_mean, 'mean_tpm2.rds')

raw_tpm_median <- aggregate( . ~ annotation, raw_tpm, median)
raw_tpm_median <- column_to_rownames(raw_tpm_median, 'annotation')
raw_tpm_median <- t(raw_tpm_median)
saveRDS(raw_tpm_median, 'median_tpm.rds')

raw_tpm_sd <- aggregate( . ~ annotation, raw_tpm, sd)
raw_tpm_sd<- column_to_rownames(raw_tpm_sd, 'annotation')
raw_tpm_sd <- t(raw_tpm_sd)
saveRDS(raw_tpm_sd, 'sd_tpm.rds')

raw_tpm_iqr <- aggregate( . ~ annotation, raw_tpm, IQR)
raw_tpm_iqr <- column_to_rownames(raw_tpm_iqr, 'annotation')
raw_tpm_iqr <- t(raw_tpm_iqr)
saveRDS(raw_tpm_iqr, 'iqr_tpm2.rds')

sigGenes <- readRDS('signature_genes.rds')
my_ref01 <- list(refProfiles = raw_tpm_median, sigGenes = sigGenes, refProfiles.var = raw_tpm_iqr)
saveRDS(my_ref01, 'my_ref01.rds')
my_ref02 <- list(refProfiles = raw_tpm_mean, sigGenes = sigGenes, refProfiles.var = raw_tpm_sd)
saveRDS(my_ref02, 'my_ref02.rds')


#Read raw counts files
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)

raw_10x <- fread('raw_data_10x.txt', sep = '\t', header = T, fill = TRUE)
raw_ss2 <- fread('raw_data_ss2.txt', sep = ' ', header = T, fill = TRUE)

#Merge the data frames by genes, for those rows missing in ss2, add 0s
merged_df <- merge(raw_10x, raw_ss2, by = 'Gene')

#Save the merged df

fwrite(merged_df, 'merged_df.txt', sep = ' ', quote = F)

#----------- Get cell.state.labels and cell.type.labels ---------------#

#Read metadata files
meta_10x <- read.csv('meta_10x.txt', sep = '\t')
meta_ss2 <- read.csv('meta_ss2.txt', sep = '\t')

#Remove spaces in annotation column in metadata

meta_10x$annotation <- sub(" ", "_", meta_10x$annotation)
meta_ss2$annotation <- sub(" ", "_", meta_ss2$annotation)

#Join the two data frames
meta <- rbind(meta_10x, meta_ss2)
write.table(meta, 'meta.txt', quote = F, sep = " ",
            row.names = T, col.names = T)

