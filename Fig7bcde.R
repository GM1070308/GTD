library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggsci)

bp54.res <- readRDS('~/thesis/bp54.res.rds')
comp <- as.data.frame(t(get.fraction (bp=readRDS('~/thesis/bp54.res.rds'),
                        which.theta="final",
                        state.or.type="type")))
rownames(comp) <- sub('Plasma cells cells' , 'Plasma cells', rownames(comp))

cc <- comp[, startsWith(names(comp), "CC")]
chm <- comp[, startsWith(names(comp), "CHM")]
healthy <- comp[, startsWith(names(comp), "H")]
pstt <- comp[, startsWith(names(comp), "PSTT")]
ett <- comp[, startsWith(names(comp), "ETT")]

cc$mean <- apply(cc[1:2], 1, mean)
cc$sd <- apply(cc[1:2], 1, sd)
chm$mean <- apply(chm[1:4],1,mean)
chm$sd <- apply(chm[1:4],1,sd)
ett$mean <- apply(ett[1:4],1,mean)
ett$sd <- apply(ett[1:4],1,sd)
pstt$mean <- apply(pstt[1:4],1,mean)
pstt$sd <- apply(pstt[1:4],1,sd)
healthy$mean <- apply(healthy[1:40], 1, mean)
healthy$sd<- apply(healthy[1:40], 1, sd)

grid_layout <- rbind(c(1, 2, 3), c(4, 5, NA))

cc <- rownames_to_column(cc, var = 'cell.type')
chm <- rownames_to_column(chm, var = 'cell.type')
ett <- rownames_to_column(ett, var = 'cell.type')
pstt <- rownames_to_column(pstt, var = 'cell.type')
healthy <- rownames_to_column(healthy, var = 'cell.type') 

order_list <- c(rev(c('Blood NKs', 'Blood T cells', 'Decidual macrophages', 'Decidual NKs', 'Decidual T cells', 'Dendritic cells', 'Hofbauer cells',
                      'Granulocytes', 'Innate lymphoid cells', 'Monocytes', 'Plasma cells')), rev(c('Decidual perivascular cells', 'Decidual stromal cells',
                                                                                                    'Endothelial cells', 'Epithelial glandular cells', 'Extravillous trophoblast', 'Foetal fibroblasts', 'Syncytiotrophoblast',
                                                                                                    'Villous cytotrophoblast')))

healthy$cell.type <- factor(healthy$cell.type, levels = order_list)
healthy <- arrange(healthy, cell.type)

chm$cell.type <- factor(chm$cell.type, levels = order_list)
chm <- arrange(chm, cell.type)

cc$cell.type <- factor(cc$cell.type, levels = order_list)
cc <- arrange(cc, cell.type)

ett$cell.type <- factor(ett$cell.type, levels = order_list)
ett <- arrange(ett, cell.type)

pstt$cell.type <- factor(pstt$cell.type, levels = order_list)
pstt <- arrange(pstt, cell.type)

#Plot total % of immune cells in each GTD subtype

tot_immune_per_sample <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(tot_immune_per_sample) <- c('sample', 'immune_prop')

b <- 1
  
  for (i in (2:(ncol(healthy)))) {
        
    tot_immune_per_sample[b,1] <- colnames(healthy)[i]
    tot_immune_per_sample[b,2] <- sum(healthy[1:11,i])
    b <- b+1
  }
  
  for (i in (2:(ncol(chm)))) {
    
    tot_immune_per_sample[b,1] <- colnames(chm)[i]
    tot_immune_per_sample[b,2] <- sum(chm[1:11,i])
    b <- b+1
  }
  
  for (i in (2:(ncol(cc)))) {
    
    tot_immune_per_sample[b,1] <- colnames(cc)[i]
    tot_immune_per_sample[b,2] <- sum(cc[1:11,i])
    b <- b+1
  }
  
  for (i in (2:(ncol(pstt)))) {
    
    tot_immune_per_sample[b,1] <- colnames(pstt)[i]
    tot_immune_per_sample[b,2] <- sum(pstt[1:11,i])
    b <- b+1
  }
  
  for (i in (2:(ncol(ett)))) {
    
    tot_immune_per_sample[b,1] <- colnames(ett)[i]
    tot_immune_per_sample[b,2] <- sum(ett[1:11,i])
    b <- b+1
  }

df <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(df) <- c('tissue', 'mean', 'sd')
df$tissue <- c('Healthy placenta', 'Complete mole', 'CC', 'PSTT', 'ETT')
df$mean <- c(mean(tot_immune_per_sample$immune_prop[1:40]), mean(tot_immune_per_sample$immune_prop[41:44]), mean(tot_immune_per_sample$immune_prop[45:46]),
             mean(tot_immune_per_sample$immune_prop[47:50]), mean(tot_immune_per_sample$immune_prop[51:54]))
df$sd <- c(sd(tot_immune_per_sample$immune_prop[1:40]), sd(tot_immune_per_sample$immune_prop[41:44]), sd(tot_immune_per_sample$immune_prop[45:46]),
           sd(tot_immune_per_sample$immune_prop[47:50]), sd(tot_immune_per_sample$immune_prop[51:54]))
df$tissue <- factor(df$tissue, levels = df$tissue)
df <- arrange(df, tissue)

tot_immune <- ggplot(df, aes(x = tissue, y = mean, fill = tissue)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = rep("#8491B4FF",5)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "black") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.4)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 15)) +
  ylab('Total immune cell content') +
  guides(fill = "none") +
  theme(plot.margin = unit(c(0.5, 0.5, 2.5, 0.5), "cm"))

#Plot a stacked barplot for immune cells in each group

comp <- rownames_to_column(comp, 'cell.type')

comp$cell.type <- factor(comp$cell.type, levels = order_list)
comp <- arrange(comp, cell.type)
comp <- comp[c(1:11),]
rownames(comp) <- NULL
comp <- column_to_rownames(comp,'cell.type')

cc <- comp[, startsWith(names(comp), "CC")]
chm <- comp[, startsWith(names(comp), "CHM")]
healthy <- comp[, startsWith(names(comp), "H")]
pstt <- comp[, startsWith(names(comp), "PSTT")]
ett <- comp[, startsWith(names(comp), "ETT")]

cc$CC <- apply(cc[1:2], 1, mean)
cc$CC <- cc$CC / sum(cc$CC)
chm$CHM <- apply(chm[1:4],1,mean)
chm$CHM <- chm$CHM / sum(chm$CHM)
ett$ETT<- apply(ett[1:4],1,mean)
ett$ETT <- ett$ETT / sum(ett$ETT)
pstt$PSTT <- apply(pstt[1:4],1,mean)
pstt$PSTT <- pstt$PSTT / sum(pstt$PSTT)
healthy$Healthy_placenta<- apply(healthy[1:40], 1, mean)
healthy$Healthy_placenta <- healthy$Healthy_placenta / sum(healthy$Healthy_placenta)

cc <- cc[ncol(cc)]
chm <- chm[ncol(chm)]
ett <- ett[ncol(ett)]
pstt <- pstt[ncol(pstt)]
healthy <- healthy[ncol(healthy)]

cc <- rownames_to_column(cc, var = 'cell.type')
chm <- rownames_to_column(chm, var = 'cell.type')
ett <- rownames_to_column(ett, var = 'cell.type')
pstt <- rownames_to_column(pstt, var = 'cell.type')
healthy <- rownames_to_column(healthy, var = 'cell.type') 

mean_comp <- merge(cc, chm, 'cell.type')
mean_comp <- merge(mean_comp, ett, 'cell.type')
mean_comp <- merge(mean_comp, pstt, 'cell.type')
mean_comp <- merge(mean_comp, healthy, 'cell.type')

my_colours <- c("#B09C85FF","#E64B35FF", "#00A087FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#3C5488FF",
                "#4DBBD5FF", "#4DBBD599", "#7E6148FF","#8491B499", "#91D1C299", "#7E614899", "#B09C8599")

stacked <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(stacked) <- c('tissue', 'cell.type', 'value')

x <- 1

for (i in (2:ncol(mean_comp))) {
  
  print(i)
  stacked[c(x:(x+10)),1] <- rep(colnames(mean_comp)[i], 11)
  stacked[c(x:(x+10)),2] <- as.character(mean_comp$cell.type)
  stacked[c(x:(x+10)),3] <- mean_comp[,i]
  x <- x + 11
}

stacked$tissue <- sub('Healthy_placenta', 'Healthy placenta', stacked$tissue)
stacked$tissue <- sub('CHM', 'Complete mole', stacked$tissue)
stacked$tissue <- sub('cc', 'CC', stacked$tissue)
stacked$tissue <- sub('pstt', 'PSTT', stacked$tissue)
stacked$tissue <- sub('ett', 'ETT', stacked$tissue)

stacked$tissue <- factor(stacked$tissue, levels = c('Healthy placenta', 'Complete mole' ,'CC', 'PSTT','ETT'))
stacked <- arrange(stacked, tissue)

stacked$cell.type <- as.character(stacked$cell.type)

imm_types <- ggplot(stacked, aes(x = tissue, y = value, fill = cell.type)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = my_colours) +
  scale_y_continuous(name = 'Mean proportion', expand = c(0,0), limits = c(0,1))+
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = 'bold'),
        axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust =1)) +
  guides(fill = guide_legend(title = 'Immune cell type'))

#Plot stacked barplot of non-immune cell types

comp <- rownames_to_column(comp, 'cell.type')

comp$cell.type <- factor(comp$cell.type, levels = order_list)
comp <- arrange(comp, cell.type)
comp <- comp[c(12:19),]
rownames(comp) <- NULL
comp <- column_to_rownames(comp,'cell.type')

cc <- comp[, startsWith(names(comp), "CC")]
chm <- comp[, startsWith(names(comp), "CHM")]
healthy <- comp[, startsWith(names(comp), "H")]
pstt <- comp[, startsWith(names(comp), "PSTT")]
ett <- comp[, startsWith(names(comp), "ETT")]

cc$CC <- apply(cc[1:2], 1, mean)
warngincc$CC <- cc$CC / sum(cc$CC)
chm$CHM <- apply(chm[1:4],1,mean)
chm$CHM <- chm$CHM / sum(chm$CHM)
ett$ETT<- apply(ett[1:4],1,mean)
ett$ETT <- ett$ETT / sum(ett$ETT)
pstt$PSTT <- apply(pstt[1:4],1,mean)
pstt$PSTT <- pstt$PSTT / sum(pstt$PSTT)
healthy$Healthy_placenta<- apply(healthy[1:40], 1, mean)
healthy$Healthy_placenta <- healthy$Healthy_placenta / sum(healthy$Healthy_placenta)

cc <- cc[ncol(cc)]
chm <- chm[ncol(chm)]
ett <- ett[ncol(ett)]
pstt <- pstt[ncol(pstt)]
healthy <- healthy[ncol(healthy)]

cc <- rownames_to_column(cc, var = 'cell.type')
chm <- rownames_to_column(chm, var = 'cell.type')
ett <- rownames_to_column(ett, var = 'cell.type')
pstt <- rownames_to_column(pstt, var = 'cell.type')
healthy <- rownames_to_column(healthy, var = 'cell.type') 

mean_comp <- merge(cc, chm, 'cell.type')
mean_comp <- merge(mean_comp, ett, 'cell.type')
mean_comp <- merge(mean_comp, pstt, 'cell.type')
mean_comp <- merge(mean_comp, healthy, 'cell.type')

stacked2 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(stacked2) <- c('tissue', 'cell.type', 'value')

y <- 1

for (i in (2:ncol(mean_comp))) {
  
  print(i)
  stacked2[c(y:(y+7)),1] <- rep(colnames(mean_comp)[i], 8)
  stacked2[c(y:(y+7)),2] <- as.character(mean_comp$cell.type)
  stacked2[c(y:(y+7)),3] <- mean_comp[,i]
  y <- y + 8
  
}

stacked2$tissue <- sub('Healthy_placenta', 'Healthy placenta', stacked2$tissue)
stacked2$tissue <- sub('CHM', 'Complete mole', stacked2$tissue)
stacked2$tissue <- sub('CC', 'CC', stacked2$tissue)
stacked2$tissue <- sub('PSTT', 'PSTT', stacked2$tissue)
stacked2$tissue <- sub('ETT', 'ETT', stacked2$tissue)

stacked2$tissue <- factor(stacked2$tissue, levels = c('Healthy placenta', 'Complete mole' ,'CC', 'PSTT','ETT'))
stacked2 <- arrange(stacked2, tissue)
stacked2$cell.type <- factor(stacked2$cell.type, levels = c('Decidual perivascular cells', 'Decidual stromal cells',
                                                             'Endothelial cells', 'Epithelial glandular cells', 'Foetal fibroblasts','Extravillous trophoblast', 'Syncytiotrophoblast',
                                                             'Villous cytotrophoblast'))
stacked2 <- arrange(stacked2, cell.type)

pastel_colours <- c("#91D1C299", "#00A08799", "#8491B499", "#3C548899","#4DBBD599", "#F39B7F99", "#E64B3599","#DC000099")

non_imm_stacked <- ggplot(stacked2, aes(x = tissue, y = value, fill = cell.type)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = pastel_colours) +
  scale_y_continuous(name = 'Mean proportion', expand = c(0,0), limits = c(0,1))+
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13, face = 'bold'),
        axis.text.x = element_text(size = 12,angle = 45, vjust = 1, hjust =1)) +
  guides(fill = guide_legend(title = 'Non-immune cell type'))

#Plot boxplots for distribution of epithelial glandular cells in each subgroup

epi <- comp[5,]
boxplot_epi <- t(epi)
boxplot_epi <- rownames_to_column(as.data.frame(boxplot_epi), 'tissue')

boxplot_epi$tissue <- sub('CC[1-9].*', 'CC', boxplot_epi$tissue)
boxplot_epi$tissue <- sub('CHM[1-9].*', 'Complete mole', boxplot_epi$tissue)
boxplot_epi$tissue <- sub('PSTT[1-9].*', 'PSTT', boxplot_epi$tissue)
boxplot_epi$tissue <- sub('ETT[1-9].*', 'ETT', boxplot_epi$tissue)
boxplot_epi$tissue <- sub('H[1-9].*', 'Healthy placenta', boxplot_epi$tissue)
colnames(boxplot_epi) <- c('tissue', 'mean')

boxplot_epi$tissue <- factor(boxplot_epi$tissue, levels = c('Healthy placenta', 'Complete mole' ,'CC', 'PSTT','ETT'))
boxplot_epi <- arrange(boxplot_epi, tissue)

boxplots <- ggplot(boxplot_epi, aes(x = tissue, y = mean)) +
  geom_boxplot(fill = "#8491B499") +
  stat_summary(fun = 'mean', color = "#3C5488FF") +
  #geom_point(size = 3) +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "black") +
  #geom_line(linetype="dashed") +
  theme_classic() +
  scale_color_manual( values = "#8491B4FF") +
  guides(color = "none") +
  scale_y_continuous(name = 'Proportion of epithelial glandular cells', expand = c(0,0), limits = c(0,0.301)) +
  scale_x_discrete(expand = c(0.3,-0.5)) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 11),
        legend.text = element_text(size = 15)) +
  theme(plot.margin = unit(c(0.5, 0.5, 2.5, 0.5), "cm"))

t.test(boxplot_epi[1:40,2], c(boxplot_epi[c(45:46),2]), var.equal = F)

pdf('comp.pdf', width = 15, height = 15)
grid.arrange(non_imm_stacked, boxplots, tot_immune, imm_types, ncol = 2, nrow = 2)
dev.off()

  




