library(ggplot2)
library(gridExtra)
library(tibble)
library(ggsci)
library(dplyr)
library(scales)

bp54.res <- readRDS('~/thesis/bp54.res.rds')
comp <- as.data.frame(t(get.fraction (bp=bp54.res,
                                      which.theta="final",
                                      state.or.type="type")))
rownames(comp) <- sub('Plasma cells cells' , 'Plasma cells', rownames(comp))

cc <- comp[, startsWith(names(comp), "CC")]
chm <- comp[, startsWith(names(comp), "CHM")]
healthy <- comp[, startsWith(names(comp), "H")]
pstt <- comp[, startsWith(names(comp), "PSTT")]
ett <- comp[, startsWith(names(comp), "ETT")]

grid_layout <- rbind(c(1, 2, 3), c(4, 5, NA))

cc <- rownames_to_column(cc, var = 'cell.type')
chm <- rownames_to_column(chm, var = 'cell.type')
ett <- rownames_to_column(ett, var = 'cell.type')
pstt <- rownames_to_column(pstt, var = 'cell.type')
healthy <- rownames_to_column(healthy, var = 'cell.type') 

boxplot_h <- data.frame(matrix(nrow = 0, ncol = 2))
boxplot_chm <- data.frame(matrix(nrow = 0, ncol = 2))
boxplot_cc <- data.frame(matrix(nrow = 0, ncol = 2))
boxplot_pstt <- data.frame(matrix(nrow = 0, ncol = 2))
boxplot_ett <- data.frame(matrix(nrow = 0, ncol = 2))

colnames(boxplot_h) <- c('cell.type', 'cell.fraction')
colnames(boxplot_chm) <- c('cell.type', 'cell.fraction')
colnames(boxplot_cc) <- c('cell.type', 'cell.fraction')
colnames(boxplot_pstt) <- c('cell.type', 'cell.fraction')
colnames(boxplot_ett) <- c('cell.type', 'cell.fraction')


for (i in (2:ncol(healthy))){
  colnames(healthy)[i] <- 'cell.fraction'
  boxplot_h <- rbind(boxplot_h, healthy[c(1,i)])
}

for (i in (2:ncol(chm))){
  colnames(chm)[i] <- 'cell.fraction'
  boxplot_chm <- rbind(boxplot_chm, chm[c(1,i)])
  print(dim(boxplot_chm))
  print(i)
}

for (i in (2:ncol(cc))){
  colnames(cc)[i] <- 'cell.fraction'
  boxplot_cc <- rbind(boxplot_cc, cc[c(1,i)])
}     

for (i in (2:ncol(pstt))){
  colnames(pstt)[i] <- 'cell.fraction'
  boxplot_pstt <- rbind(boxplot_pstt, pstt[c(1,i)])
}

for (i in (2:ncol(ett))){
  colnames(ett)[i] <- 'cell.fraction'
  boxplot_ett <- rbind(boxplot_ett, ett[c(1,i)])
}

order_list <- c(rev(c('Blood NKs', 'Blood T cells', 'Decidual macrophages', 'Decidual NKs', 'Decidual T cells', 'Dendritic cells', 'Hofbauer cells',
                'Granulocytes', 'Innate lymphoid cells', 'Monocytes', 'Plasma cells')), rev(c('Decidual perivascular cells', 'Decidual stromal cells',
                'Endothelial cells', 'Epithelial glandular cells', 'Extravillous trophoblast', 'Foetal fibroblasts', 'Syncytiotrophoblast',
                'Villous cytotrophoblast')))

boxplot_h$cell.type <- factor(boxplot_h$cell.type, levels = order_list)
boxplot_h <- arrange(boxplot_h, cell.type)

boxplot_chm$cell.type <- factor(boxplot_chm$cell.type, levels = order_list)
boxplot_chm <- arrange(boxplot_chm, cell.type)

boxplot_cc$cell.type <- factor(boxplot_cc$cell.type, levels = order_list)
boxplot_cc <- arrange(boxplot_cc, cell.type)

boxplot_ett$cell.type <- factor(boxplot_ett$cell.type, levels = order_list)
boxplot_ett <- arrange(boxplot_ett, cell.type)

boxplot_pstt$cell.type <- factor(boxplot_pstt$cell.type, levels = order_list)
boxplot_pstt <- arrange(boxplot_pstt, cell.type)

my_colours <- rev(c(pal_npg("nrc")(10), pal_npg("nrc", alpha = 0.6)(10)))

grid.arrange(plot1, plot2, plot3, plot4, plot5, layout_matrix = grid_layout)

plot1 <- ggplot(boxplot_h, aes(x=cell.type, y = cell.fraction ,fill = cell.type)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,1.001)) +
  ylab("Cell fraction") +
  ggtitle('Healthy placenta (n = 40)') +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank()) +
  guides(fill = "none")

plot2 <- ggplot(boxplot_chm, aes(x=cell.type, y = cell.fraction, fill = cell.type)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,1.001)) +
  ylab("Cell fraction") +
  ggtitle('Complete mole (n = 4)') +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank()) +
  guides(fill = "none")

plot3 <- ggplot(boxplot_cc, aes(x=cell.type, y = cell.fraction, fill = cell.type)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,1.001)) +
  ylab("Cell fraction") +
  ggtitle('Choriocarcinoma (n = 2)') +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank()) +
  guides(fill = "none")

plot4 <- ggplot(boxplot_ett, aes(x=cell.type, y = cell.fraction, fill = cell.type)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,1.001)) +
  ylab("Cell fraction") +
  ggtitle('ETT (n = 4)') +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank()) +
  guides(fill = "none")

plot5 <- ggplot(boxplot_pstt, aes(x=cell.type, y = cell.fraction, fill = cell.type)) + 
  geom_boxplot() +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0,1.001)) +
  ylab("Cell fraction") +
  ggtitle('PSTT (n = 4)') +
  coord_flip() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.title.y = element_blank()) +
  guides(fill = "none")

pdf("comp54_boxplots.pdf", width = 15, height = 8)
grid.arrange(plot1, plot2, plot3, plot4, plot5, layout_matrix = grid_layout)
dev.off()

