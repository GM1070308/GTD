library(EPIC)
library(biomaRt)
library(DGEobj.utils)

bulk <- read.csv('~/EPIC/bk2.txt', sep = ' ')
ensembl_ids <- rownames(Z.evt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org", port = 443)

# Specify the dataset and attributes to retrieve
dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
attributes <- c("ensembl_gene_id", "hgnc_symbol")

# Get the gene length for each ENSEMBL ID
gene_info <- getBM(attributes = attributes, filters = "ensembl_gene_id", values = ensembl_ids, mart = dataset)
colnames(gene_info) <- c('gene_id', 'symbol')


gene_info$gene_length <- gene_info$end_position - gene_info$start_position + 1

gene_info <- gene_info[!duplicated(gene_info$hgnc_symbol),]
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id),]
rownames(gene_info) <- NULL
gene_info <- column_to_rownames(gene_info, 'hgnc_symbol')


gene_info <- gene_info[3]

merged <- merge(x=bulk,y=gene_info, by = 0)
merged <- column_to_rownames(merged, 'Row.names')

bulk_tpm <- convertCounts(counts = as.matrix(merged[-ncol(merged)]), unit = 'TPM', geneLength = as.vector(merged$gene_length),
                          log = FALSE, normalize = 'none')

out_qc <- EPIC(bulk = readRDS('~/EPIC/pseudo_ensembl_tpm.rds'), ref = readRDS('~/EPIC/my_ref01.rds'))
                          