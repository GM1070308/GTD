library(ggrepel)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tibble)
library(dplyr)
library(ggsci)

scaleFUN <- function(x) sprintf("%.2f", x)

epic <- readRDS('~/EPIC/obs_vs_exp_cell_sepc_pseudobulk_ref01_ref02.rds')

prop_og.type$cell.type <- as.factor(prop_og.type$cell.type)

my_colours <- c(pal_npg("nrc")(10),pal_npg("nrc")(10))

#Estimate vs Observed plot for EPIC
epic_qc <- ggplot(epic, aes(x = exp, y = sum, color = cell.type)) +
  geom_point(size = 5, shape = 1, stroke = 2) +
  labs(x = "Expected", y = "Estimated") +
  scale_color_manual(values = my_colours) + 
  theme_minimal() +
  guides(color = "none") +
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              linewidth = 1, 
              linetype = 2) +
  expand_limits(x = 0, y = 0) +
  ggtitle('EPIC') +
  scale_x_continuous(expand = c(0, 0), limits = c(0.00001,0.202)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.202)) +
  annotate('text', x=0.01, y = 0.185, label = paste('p = 0.153'), fontface = 'bold', color = 'black', size = 7, hjust= 'left') +
  annotate('text', x=0.01, y = 0.175, label = paste0('r = 0.351'), fontface = 'bold', color = 'black', size = 7, hjust = 'left') +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.5)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 17)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  

t_epic <- cor.test(epic$exp, epic$sum, method = "pearson")

prop_og.type <- readRDS('~/thesis/prop_og.type.rds')

bp_og <- ggplot(epic, aes(x = exp, y = bp, color = cell.type)) +
  geom_point(size = 5,shape = 1, stroke = 2) +
  labs(x = "Expected", y = "Estimated") +
  scale_color_manual(values = my_colours) +
  theme_minimal() +
  guides(color = "none") +
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              linewidth = 1, 
              linetype = 2) +
  expand_limits(x = 0, y = 0) +
  ggtitle('BayesPrism') +
  scale_x_continuous(expand = c(0.00, 0.00), limits = c(0.00001,0.3015), labels = scaleFUN) + 
  scale_y_continuous(expand = c(0.00, 0.00), limits = c(0.00,0.3015), labels = scaleFUN) +
  annotate('text', x=0.015, y = 0.278, label = paste('p = 0.000312'), fontface = 'bold', color = 'black', size = 6, hjust= 'left') +
  annotate('text', x=0.015, y = 0.262, label = paste0('r = 0.753'), fontface = 'bold', color = 'black', size = 6, hjust = 'left') +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1.5)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 17))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

t_bp_og <- cor.test(epic$exp, epic$bp, method = 'pearson')

#Improved annotation

prop.type <- readRDS('~/thesis/prop.type.rds')

bp_qc <- 
  ggplot(prop.type, aes(x = exp, y = bp, color = cell.type)) +
  geom_text(data = subset(prop.type, exp > 0.05 & exp < 0.19 & bp > 0.03), label=(subset(prop.type, exp > 0.05 & exp < 0.19 & bp > 0.03))$cell.type, nudge_y = 0.01) +
  geom_text(data = subset(prop.type, exp > 0.19 & bp > 0.03), label=(subset(prop.type, exp > 0.19 & bp > 0.03))$cell.type, nudge_y = 0.01, nudge_x = 0.025) +
  geom_text(data = subset(prop.type, exp > 0.1 & bp < 0.05), label=(subset(prop.type, exp > 0.1 & bp < 0.05))$cell.type, nudge_y = 0.01) +
  geom_point(size = 5, shape = 1, stroke = 2) +
  scale_color_manual(values = my_colours) +
  labs(x = "Expected", y = "Estimated") +
  theme_minimal() +
  guides(color = "none") +
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              linewidth = 1, 
              linetype = 2) +
  expand_limits(x = 0, y = 0) +
  ggtitle('BayesPrism') +
  scale_x_continuous(expand = c(0, 0), limits = c(0.00001,0.3015), labels = scaleFUN ) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.3015), labels = scaleFUN) +
  annotate('text', x=0.01, y = 0.283, label = paste('p = 0.000114'), fontface = 'bold', color = 'black', size = 6, hjust= 'left') +
  annotate('text', x=0.01, y = 0.267, label = paste0('r = 0.770'), fontface = 'bold', color = 'black', size = 6, hjust = 'left') +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.5)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 17))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

t_bp <- cor.test(prop.type$exp, prop.type$bp, method = 'pearson')

ggsave('~/thesis/bp_updated_atlas.png')

#Plot contaminated qc
prop.type <- readRDS('~/thesis/prop.type.rds')

cont <- t(get.fraction (bp=readRDS('~/thesis/bp_qc_contaminated.res.rds'),
                        which.theta="final",
                        state.or.type="type"))

cont <- rownames_to_column(as.data.frame(cont), var = 'cell.type')

prop <- as.data.frame(table(cell.type.labels))
prop$exp <- prop$Freq / 89134
colnames(prop) <- c('cell.type', 'Freq', 'exp')
prop <- merge(prop, cont)

bp_cont <- ggplot(prop, aes(x = exp, y = sum, color = cell.type)) +
  geom_text(data = subset(prop, exp > 0.03 & exp < 0.14 & sum > 0.03), label=(subset(prop, exp > 0.03 & exp < 0.14 & sum > 0.03))$cell.type, nudge_y = 0.01) +
  geom_text(data = subset(prop, exp > 0.14 & sum > 0.03), label=(subset(prop, exp > 0.14 & sum > 0.03))$cell.type, nudge_y = 0.01, nudge_x = 0.025) +
  geom_text(data = subset(prop, exp > 0.07 & sum < 0.05), label=(subset(prop, exp > 0.07 & sum < 0.05))$cell.type, nudge_y = 0.01) +
  geom_point(size = 5, shape = 1, stroke = 2) +
  scale_color_manual(values = my_colours) +
  labs(x = "Expected", y = "Estimated") +
  theme_minimal() +
  guides(color = "none") +
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              linewidth = 1, 
              linetype = 2) +
  expand_limits(x = 0, y = 0) +
  ggtitle('BayesPrism') +
  scale_x_continuous(expand = c(0, 0), limits = c(0.00001,0.3015), labels = scaleFUN ) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.3015), labels = scaleFUN) +
  annotate('text', x=0.01, y = 0.283, label = paste('p = 0.000164'), fontface = 'bold', color = 'black', size = 6, hjust= 'left') +
  annotate('text', x=0.01, y = 0.267, label = paste0('R = 0.759'), fontface = 'bold', color = 'black', size = 6, hjust = 'left') +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 1.5)) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 17))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

t <- cor.test(prop$exp, prop$sum, method = 'pearson')

pdf('qc_tools.pdf', width = 15, height = 15)
grid.arrange(bp_og,epic_qc, bp_qc, bp_cont, ncol = 2, nrow = 2)
dev.off()
