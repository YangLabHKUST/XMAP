library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
# library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(ggridges)
library(dplyr)
library(pheatmap)
library(grid)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)

all_traits <- c("Lym", "WBC", "Neutro", "Mono", "Eosino", "Baso", "Plt", "RBC", "MCH", "MCHC", "MCV", "Hb")
# method <- "xmap"
method <- "ukb"
# method <- "bbj"
out_all <- data.frame()
for (i in 1:length(all_traits)) {
  trait <- all_traits[i]
  trait_mat <- fread(paste0("/home/share/mingxuan/fine_mapping/analysis/results/scavenge_output_", trait, "_", method, ".txt"))
  out_all <- rbind(out_all, data.frame(trait_mat, trait = all_traits[i]))
  cat(i, "-th trait finished.\n")
}
out_all$BioClassification <- factor(out_all$BioClassification,
                                    levels = c("Late.Eryth", "Early.Eryth", "Early.Baso", "CD14.Mono.1", "GMP.Neut", "cDC", "GMP", "CD14.Mono.2", "HSC", "CMP.LMPP", "CLP.1", "CLP.2", "Plasma", "B", "Pre.B", "pDC", "CD8.EM", "CD8.CM", "CD4.M", "CD8.N", "NK", "CD4.N1", "CD4.N2"))

################################## UMAP ##################################
class_avg <- trait_mat %>%
  group_by(BioClassification) %>%
  summarise(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2)
  )

colourCount = length(unique(trait_mat$BioClassification))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(colourCount)
p <- ggplot(data = trait_mat, aes(UMAP1, UMAP2, color = BioClassification)) +
  geom_point(size = 1, na.rm = TRUE) +
  geom_text(aes(label = BioClassification), data = class_avg, color = "black") +
  scale_color_manual(values = mycolors) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "none")
p

################################## heatmap ##################################


out_median <- out_all %>%
  group_by(trait, BioClassification,) %>%
  summarise_at("TRS", list(~median(., na.rm = TRUE)))



## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function(coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 75, gp = gpar(...))
  return(res) }

## 'Overwrite' default draw_colnames with your own version
assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
                  ns = asNamespace("pheatmap"))
mat_median <- dcast(out_median, trait ~ BioClassification)
ann_trait <- data.frame(Category = c("Lymphoid white cell",
                                     rep("Myeloid/compound white cell", 5),
                                     "Platelet",
                                     rep("Red cell", 5)))
rownames(mat_median) <- mat_median$trait
rownames(ann_trait) <- mat_median$trait
my_col <- RColorBrewer::brewer.pal(4, "Paired")
names(my_col) <- c("Lymphoid white cell", "Myeloid/compound white cell", "Platelet", "Red cell")
pheatmap(data.matrix(mat_median[, -1]), annotation_row = ann_trait, annotation_names_row = FALSE,
         annotation_colors = list(Category = my_col),
         color = colorRampPalette(c("white", "red"))(50), cluster_rows = FALSE, cluster_cols = FALSE)


viridis = c("#440154FF", "#472D7BFF", "#3B528BFF", "#2C728EFF", "#21908CFF", "#27AD81FF", "#5DC863FF", "#AADC32FF", "#FDE725FF")


# bar plot of prop
out_all$label <- out_all$label_text <- as.character(out_all$BioClassification)
out_prop <- out_all %>%
  group_by(true_cell_top_idx, trait) %>%
  count(BioClassification) %>%
  mutate(cell_prop = n / sum(n)) %>%
  ungroup
for (i in 1:nrow(out_prop)) {
  if (out_prop$cell_prop[i] < 0.03) {
    out_all$label_text[out_all$trait == out_prop$trait[i] &
                         out_all$BioClassification == out_prop$BioClassification[i] &
                         out_all$true_cell_top_idx == out_prop$true_cell_top_idx[i]] <- "Other"
    out_all$label[out_all$trait == out_prop$trait[i] &
                    out_all$BioClassification == out_prop$BioClassification[i] &
                    out_all$true_cell_top_idx == out_prop$true_cell_top_idx[i]] <- "Other"
  } else if (out_prop$cell_prop[i] < 0.05) {
    idx <- which(out_all$trait == out_prop$trait[i] &
                   out_all$BioClassification == out_prop$BioClassification[i] &
                   out_all$true_cell_top_idx == out_prop$true_cell_top_idx[i] &
                   nchar(out_all$label_text) > 5)
    out_all$label_text[idx] <- sub("[.]", ".\n", out_all$label_text[idx])
  }
}
# Define the number of colors you want
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(out_all$label)))
ggplot(out_all[out_all$true_cell_top_idx,], aes(x = trait, fill = label)) +
  geom_bar(position = position_fill()) +
  geom_text(stat = 'count', aes(label = label_text), position = position_fill(vjust = 0.5), size = 2.6, lineheight = 0.8) +
  # geom_text_repel(stat = 'count', aes(label = label), position = position_fill(vjust = 0.5),direction = "x") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.01))) +
  ylab("Trait-enriched cell proportion") +
  xlab(element_blank()) +
  coord_flip() +
  labs(fill = "Cell type")
#1500*400 120


# bar plot of count
out_all$label_text <- as.character(out_all$BioClassification)
out_prop <- out_all %>%
  group_by(true_cell_top_idx, trait) %>%
  count(BioClassification) %>%
  mutate(cell_prop = n / sum(n)) %>%
  ungroup

for (i in 1:nrow(out_prop)) {
  if (out_prop$n[i] / nchar(as.character(out_prop$BioClassification[i])) < 50) {
    out_all$label_text[out_all$trait == out_prop$trait[i] &
                         out_all$BioClassification == out_prop$BioClassification[i] &
                         out_all$true_cell_top_idx == out_prop$true_cell_top_idx[i]] <- NA
  }
}
# Define the number of colors you want
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(out_all$BioClassification)))
ggplot(out_all[out_all$true_cell_top_idx,], aes(x = trait, fill = as.character(BioClassification))) +
  geom_bar(position = position_stack(vjust = 0.5)) +
  geom_text(stat = 'count', aes(label = label_text), position = position_stack(vjust = 0.5), lineheight = 0.8) +
  # geom_text_repel(stat = 'count', aes(label = label), position = position_fill(vjust = 0.5),direction = "x") +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  ylab("Trait-enriched cell count") +
  xlab(element_blank()) +
  coord_flip() +
  labs(fill = "Cell type")
#1500*400 120


for (i in 1:length(all_traits)) {
  trait_name <- all_traits[i]
  p_umap <- ggplot(data = out_all[out_all$trait == trait_name,], aes(UMAP1, UMAP2, color = TRS)) +
    geom_point(size = 1, na.rm = TRUE, alpha = 0.6) +
    scale_color_gradientn(colors = viridis) +
    scale_alpha() +
    theme_classic() +
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.position = c(0.8, 0.8),
          axis.title = element_text(size = 18)) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    annotate("text", x = -9, y = 13, label = trait_name, size = 8)
  # p_umap
  #size: 750*650

  ############# Enriched cell prop boxplot-circle #############
  p_bar <- out_all[out_all$trait == trait_name,] %>%
    group_by(BioClassification, trait) %>%
    summarise(enriched_cell = mean(true_cell_top_idx)) %>%
    ggplot(aes(x = BioClassification, y = enriched_cell, fill = BioClassification)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(y = enriched_cell, label = scales::percent(enriched_cell, accuracy = 0.1)), position = position_dodge(width = 0.8)) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(-1, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ylab("Trait-enriched cell proportion") +
    xlab(element_blank()) +
    theme(
      panel.border = element_blank(),
      plot.margin = unit(rep(0, 4), "cm"),
      legend.position = "none",
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      axis.title.y = element_text(size = 18)
    ) +
    labs(fill = "Cell type") +
    coord_polar(start = 0) +
    annotate("text", x = 0, y = -1, label = trait_name, size = 8)
  #size: 1000*950
  # p_bar




  ggsave(paste0("/home/share/mingxuan/fine_mapping/figure/SCAVENGE_", method, "_blood_", trait_name, ".pdf"),
         ggarrange(p_umap, p_bar,
                   ncol = 2, nrow = 1),
         width = 15, height = 7.5)

  # pdf(paste0("/home/share/mingxuan/fine_mapping/figure/SCAVENGE_",method,"_blood_",trait_name,".pdf"),
  #     width = 19,height=9.5)
  # print(ggarrange(p_umap, p_bar,
  #         ncol = 2, nrow = 1))
  # dev.off()
  #size: 1900*950
}
