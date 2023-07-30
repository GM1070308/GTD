library(DESeq2)
library(BayesPrism)
library(dplyr)
library(tidyverse)
library(pheatmap)
library("RColorBrewer")

comp <- as.data.frame(t(get.fraction (bp=readRDS('~/thesis/bp_troph.res.rds'),
                                      which.theta="first",
                                      state.or.type="state")))

############################ Z.evt #############################################

Z.evt <- t(get.exp(bp=readRDS('~/thesis/bp_troph.res.rds'),
                   state.or.type = 'state',
                   cell.name = 'EVT'))

evt <- filter(comp, rownames(comp) == 'EVT')
evt <- as.data.frame(t(evt))
colnames(evt) <- 'value'
evt <- filter(evt,  value > 0.01)

Z.evt <- Z.evt[,colnames(Z.evt) %in% rownames(evt)]

metadata <- data.frame(matrix(ncol = 1, nrow = ncol(Z.evt)))
rownames(metadata) <- colnames(Z.evt)
colnames(metadata) <- 'status'
metadata$status <- rownames(metadata)
metadata$status <- sub('H[0-9].*', 'Healthy placenta',metadata$status)
metadata$status <- sub('CC[0-9].*', 'CC',metadata$status)
metadata$status <- sub('CHM[0-9].*', 'Complete mole',metadata$status)
metadata$status <- sub('PSTT[0-9].*', 'PSTT',metadata$status)
metadata$status <- sub('ETT[0-9].*', 'ETT',metadata$status)

Z.evt <- round(Z.evt)
Z.evt<- Z.evt + 1

input <- Z.evt

################################################################################
res_df <- data.frame("gene_id" = character(), 
                     "baseMean" = numeric(), 
                     "log2FoldChange" = numeric(), 
                     "lfcSE" = numeric(),
                     "stat" = numeric(),
                     "pvalue" = numeric(),
                     "padj" = numeric(),
                     "symbol" = character(),
                     "result_name" = character())

dds <- DESeqDataSetFromMatrix(countData = input,
                              colData = metadata,
                              design = ~ status)

dds$status <- factor(dds$status, levels = c('Healthy placenta', 'Complete mole','CC','PSTT', 'ETT'))

dds <- DESeq(dds)

results_names <- resultsNames(dds)[-1]

for(i in 1:length(results_names)) {
  # grabbing the name of the result file i
  print(i)
  results_name <- results_names[i]
  # populating the res_df with results(dds)
  # x <- results(dds, name = results_name)
  res <- results(dds, name = results_name)
  
  # populating data.frame 1 : temp_res_df
  tmp_res_df <- res %>% as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    merge(gene_info) %>%
    mutate(result_name = results_name)
  
  # Append to full data.frame
  res_df <- bind_rows(res_df, tmp_res_df)
}

metadata <- filter(metadata,status != 'Healthy placenta')
input <- Z.evt[, !startsWith(names(as.data.frame(Z.evt)), "H")]

dds <- DESeqDataSetFromMatrix(countData = input,
                              colData = metadata,
                              design = ~ status)

dds$status <- factor(dds$status, levels = c('Complete mole','CC','PSTT', 'ETT'))

dds <- DESeq(dds)

results_names <- resultsNames(dds)[-1]

for(i in 1:length(results_names)) {
  # grabbing the name of the result file i
  print(i)
  results_name <- results_names[i]
  # populating the res_df with results(dds)
  # x <- results(dds, name = results_name)
  res <- results(dds, name = results_name)
  
  # populating data.frame 1 : temp_res_df
  tmp_res_df <- res %>% as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    merge(gene_info) %>%
    mutate(result_name = results_name)
  
  # Append to full data.frame
  res_df <- bind_rows(res_df, tmp_res_df)
}

metadata <- filter(metadata,status != 'Complete mole')
input <- input[, !startsWith(names(as.data.frame(input)), "CHM")]

dds <- DESeqDataSetFromMatrix(countData = input,
                              colData = metadata,
                              design = ~ status)

dds$status <- factor(dds$status, levels = c('CC','PSTT', 'ETT'))

dds <- DESeq(dds)

results_names <- resultsNames(dds)[-1]

for(i in 1:length(results_names)) {
  # grabbing the name of the result file i
  print(i)
  results_name <- results_names[i]
  # populating the res_df with results(dds)
  # x <- results(dds, name = results_name)
  res <- results(dds, name = results_name)
  
  # populating data.frame 1 : temp_res_df
  tmp_res_df <- res %>% as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    merge(gene_info) %>%
    mutate(result_name = results_name)
  
  # Append to full data.frame
  res_df <- bind_rows(res_df, tmp_res_df)
}

res_df <- res_df[!res_df$symbol == '',]

markers_res <-merge(res_df,markers_df)

markers_res$symbol <- factor(markers_res$symbol, levels = markers)
markers_res <- arrange(markers_res, symbol)

results_order <- c('status_Complete.mole_vs_Healthy.placenta', 'status_CC_vs_Healthy.placenta', 'status_PSTT_vs_Healthy.placenta',
                   'status_ETT_vs_Healthy.placenta', 'status_CC_vs_Complete.mole', 'status_PSTT_vs_Complete.mole', 'status_ETT_vs_Complete.mole',
                   'status_PSTT_vs_CC', 'status_ETT_vs_CC')

markers_res$result_name <- factor(markers_res$result_name, levels = results_order)
markers_res <- arrange(markers_res, result_name)

sig_genes <- markers_res %>%
  filter(padj < 0.05, abs(log2FoldChange) > 2)
 
genes_to_plot <- res_df %>%
  filter(symbol %in% sig_genes$symbol)

genes_to_plot$symbol <- factor(genes_to_plot$symbol, levels = markers)
genes_to_plot <- arrange(genes_to_plot, symbol)

genes_to_plot$result_name <- factor(genes_to_plot$result_name, levels = results_order)
genes_to_plot <- arrange(genes_to_plot, result_name)

genes_to_plot$log2FoldChange[genes_to_plot$log2FoldChange > 0] <- 3

genes_to_plot$log2FoldChange[genes_to_plot$log2FoldChange < 0] <- -3

genes_to_plot$log2FoldChange[genes_to_plot$padj < 0.05] <- 0

lfc_matrix <- genes_to_plot %>% 
  dplyr::select(symbol, log2FoldChange, result_name) %>%
  pivot_wider(names_from = "result_name", values_from = "log2FoldChange") %>%
  column_to_rownames("symbol") %>%
  as.matrix()

annotation <- data.frame(matrix(nrow = ncol(lfc_matrix), ncol = 1))
colnames(annotation) <- 'Contrast'
rownames(annotation) <- colnames(lfc_matrix)
annotation[c(1:4),1] <- rep('Healthy', 4)
annotation[c(5:7),1] <- rep('CHM', 3)
annotation[c(8:9),1] <- rep('CC', 2)

colnames(lfc_matrix) <- sub('status_','',colnames(lfc_matrix))
colnames(lfc_matrix) <- sub('_vs_.*','',colnames(lfc_matrix))
colnames(lfc_matrix) <- sub('Complete.mole','CHM',colnames(lfc_matrix))

my_colours <- c("#91D1C299", "#8491B499", "#F39B7FFF")
names(my_colours) <- unique(annotation$Contrast)
my_colours <- list(Contrast =c(Healthy = "#91D1C299", CHM ="#8491B499", CC = "#F39B7FFF" ))

col <- colorRampPalette(c("#4DBBD5FF", 'white',"#E64B35FF"))
pheatmap::pheatmap((t(lfc_matrix)), show_rownames = T, border_color = 'black',
                   breaks = seq(-3, 3, length.out = 100), 
                   #annotation_row = annotation, annotation_colors = my_colours, annotation_legend = F,
                   cluster_cols=F,
                   legend = F, legend_labels = c('Proliferation', 'Invasion'),
                   cluster_rows=F, color = col(100),
                   fontsize = 10, cellheight = 10, angle_col = 90,
                   legend_breaks = c(-3,3),gaps_row = c(4,4,7,7), gaps_col = c(7,10,13,17,21,33,36))
table(sig_genes$result_name)

ggsave('pheatmap_vs_CHM.png')

########################### Genes of interest from the literature ###############

immune_genes <- c('HLA-G', 'ACKR2', 'PVR', 'CLEC2D', 'CD274', 
                   'IL10RA', 'IL10RB', 'CXCR6', 
                  'PTMA', 'CALML6', 'CCR7', 'CD40LG', 'DNTT', 'IL1RN', 'IRF2BP2', 'ITK')
adh_recr <- c('XCR1', 'EPHB2', 'CSF1R','CCR1', 'PLXNB1', 'FGFR1')
invasion <- c('PAPPA','PAPPA2', 'MMP9')
angio <- c('FLT1','KDR', 'NCAM1', 'B3GAT1','MKI67', 'WNT10B', 'AHSA1', 'ACE2','MYOZ1','OR51E2', 'OR7A5','HEXIM1')
growth <- c('EGFR','OSMR', 'CCND2', 'EGF', 'GATA4', 'IGF2R','BMP7', 'TGFB1','TGFBR1', 'TGFBR2', 'TGFBR3','SMAD3', 'HSP5', 
            'TP53', 'RB1', 'CDKN1A', 'ERBB3')
pi3k <- c('PTEN', 'PIK3CA','MTOR','AKT1','AKT2','AKT3')
meth <- c('DNMT1', 'DNMT3A', 'DNMT3B')
metabolism <- c('ASPSCR1', 'FADS2','PGK2', 'SLC2A1', 'SLC2A2')
markers <- c(immune_genes, adh_recr, invasion, pi3k, angio, growth, meth, metabolism)
markers_df <- data.frame(matrix(ncol = 1, nrow = length(markers)))
colnames(markers_df) <- 'symbol'
markers_df$symbol <- markers
markers_df$symbol <- as.factor(markers_df$symbol)

########################## Run FGSEA ##########################################
library(fgsea)

CHM <- filter(res_df, result_name == 'status_Complete.mole_vs_Healthy.placenta' )
CC <- filter(res_df, result_name == 'status_CC_vs_Healthy.placenta' )
PSTT <- filter(res_df, result_name == 'status_PSTT_vs_Healthy.placenta' )
ETT <- filter(res_df, result_name == 'status_ETT_vs_Healthy.placenta' )

CC <- filter(res_df, result_name == 'status_CC_vs_Complete.mole' )
PSTT <- filter(res_df, result_name == 'status_PSTT_vs_Complete.mole' )
ETT <- filter(res_df, result_name == 'status_ETT_vs_Complete.mole' )

PSTT <- filter(res_df, result_name == 'status_PSTT_vs_CC' )
ETT <- filter(res_df, result_name == 'status_ETT_vs_CC' )

res2 <- ETT %>% 
  dplyr::select(symbol, stat) %>% 
  na.omit() %>% 
  filter(symbol != '') %>%
  distinct() %>% 
  group_by(symbol) %>% 
  summarize(stat=mean(stat))

ranks <- deframe(res2)
head(ranks, 20)

pathways.hallmark <- gmtPathways("~/thesis/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgsea_chm <- as.data.frame(fgseaResTidy)
fgsea_chm$results_name <- rep('CHM', nrow(fgsea_chm))

fgsea_cc <- as.data.frame(fgseaResTidy)
fgsea_cc$results_name <- rep('CC', nrow(fgsea_cc))

fgsea_pstt <- as.data.frame(fgseaResTidy)
fgsea_pstt$results_name <- rep('PSTT', nrow(fgsea_pstt))

fgsea_ett <- as.data.frame(fgseaResTidy)
fgsea_ett$results_name <- rep('ETT', nrow(fgsea_ett))

fgseaResTidy <- rbind(fgsea_chm,fgsea_cc)
fgseaResTidy <- rbind(fgsea_cc, fgsea_pstt)
fgseaResTidy <- rbind(fgsea_pstt, fgsea_ett)

fgseaResTidy$pathway <- sub('HALLMARK_','',fgseaResTidy$pathway)
fgseaResTidy$pathway <- sub('GOBP_','',fgseaResTidy$pathway)
fgseaResTidy$pathway <- sub('KEGG_','',fgseaResTidy$pathway)

fgseaResTidy$pathway <- gsub('_',' ',fgseaResTidy$pathway)
fgseaResTidy <- filter(fgseaResTidy, padj < 0.05)


gradient <- colorRampPalette(c("#4DBBD5FF", "white"))
gradient <- colorRampPalette(c("#3C5488FF", "white"))
gradient <- colorRampPalette(rev(c("#00A087FF", "white")))
gradient <- colorRampPalette(c("#E64B35FF", "white"))

fgsea_vs_control <- readRDS('~/thesis/BayesPrism/hallmarks_vs_control.rds')
fgsea_vs_control$results_name <- factor(fgsea_vs_control$results_name, levels = c('CHM', 'CC', 'PSTT','ETT'))
fgsea_vs_control <- arrange(fgsea_vs_control, results_name)
fgsea_vs_control$pathway <- sub('HALLMARK_','',fgsea_vs_control$pathway)
fgsea_vs_control$pathway <- gsub('_',' ',fgsea_vs_control$pathway)

fgsea_kegg_chm <- readRDS('~/thesis/BayesPrism/keg')
fgsea_kegg_chm$pathway <- sub('KEGG_','',fgsea_kegg_chm$pathway)
fgsea_kegg_chm$pathway <- gsub('_',' ',fgsea_kegg_chm$pathway)


ggplot(fgsea_vs_control[-c(2,5,7,10,14,24,26,34),], aes(reorder(pathway, NES), NES),alpha = -log(padj)) +
  scale_fill_manual(values = c("#4DBBD5FF","#91D1C2FF","#F39B7FFF","#E64B35FF")) +
  scale_alpha_continuous(range = c(0.7,1), limits = c(3.3,3.6), guide = 'none') +
  geom_col(aes(fill = results_name, alpha = -log(padj))) +
  #geom_col(aes(alpha=log_p)) +
  #scale_fill_gradientn(colours = gradient(10), limits = c(3,4)) +
  scale_y_continuous(limits = c(-5,3)) +
  coord_flip() +
  guides(fill = guide_legend(title = 'GTD subtype')) +
  labs(x=NULL, y="Normalized Enrichment Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
        axis.text = element_text(size = 13, color = 'black'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13, face = 'bold'))

ggsave('gsea_evt_ett_vs_cc_kegg.png')
