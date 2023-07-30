library(BayesPrism)
library(tibble)
print('Loading input files')

#bk.dat <- t(read.csv('~/EPIC/bk2.txt', sep = ' '))
bk.dat <- t(read.csv('~/EPIC/bk2.txt', sep = ' '))
sc.dat <- as.matrix(readRDS('~/BayesPrism/sc.dat.rds'))

# Create pseudo-bulk
#CAFs <- read.csv('~/GSE72056_melanoma_single_cell_revised_v2.txt', sep = '\t')
#CAFs <- CAFs[!duplicated(CAFs$Cell),]
#rownames(CAFs) <- NULL
#CAFs <- column_to_rownames(CAFs, 'Cell')
#CAFs <- CAFs[,CAFs[3,]== 5]
#CAFs <- CAFs[-c(1:3),]
#CAFs <- merge(x=CAFs, y=gene_info, by = 0)


#CAFs <- CAFs[!duplicated(CAFs$ensembl_gene_id),]

#CAFs <- column_to_rownames(CAFs, 'ensembl_gene_id')
#CAFs <- CAFs[-1]
#sums <- rowSums(CAFs)  
#CAFs$cafs <- sums
#CAFs$cafs <- CAFs$cafs*400
#bk.dat <- (merge(x=bk.dat, y = CAFs[ncol(CAFs)], by = 0))
#bk.dat$sum <- bk.dat$sum + bk.dat$cafs
#bk.dat <- column_to_rownames(bk.dat, 'Row.names')
#bk.dat <- bk.dat[1]
#bk.dat <- t(bk.dat)

metadata <- readRDS('~/thesis/BayesPrism/updated_meta.rds')
#metadata <- read.csv('../EPIC/meta_10x.txt', sep = '\t')
#metadata <- metadata[4]

#Chnange T cell annotation
#metadata <- filter(metadata, annotation != 'Tcells')
#metadata <- rbind(metadata, new_annot_tcells)

#Order metadata rows in the same way as sc.dat rows
metadata <- metadata[match(rownames(sc.dat), rownames(metadata)), ]
cell.type.labels <- metadata$cell.type
cell.state.labels <- metadata$cell.state

any(rownames(sc.dat) == rownames(metadata))

print('Filtering sc data')
sc.dat.filtered <- cleanup.genes (input=sc.dat, input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

print('Include only protein-coding genes')
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")


myPrism <- new.prism(
  reference=sc.dat.filtered.pc, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.state.labels,
  key='Extravillous trophoblast',
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

bp.res <- run.prism(prism = myPrism, n.cores=8)


saveRDS(bp.res, 'bp.res')