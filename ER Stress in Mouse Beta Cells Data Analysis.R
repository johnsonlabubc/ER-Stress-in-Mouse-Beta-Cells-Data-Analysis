library(ComplexHeatmap)
library(dplyr)
library(furrr)
library(ggbeeswarm)
library(ggforce)
library(ggplot2)
library(grid)
library(gridSVG)
library(pcaMethods)
library(scales)
library(stringr)
library(tidyr)
library(VennDiagram)
library(zoo)


### Figure 1 ###

## Fig 1A ##

#clear environment
rm(list = ls())

#import files
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set undefined column name as gene
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#get list of common genes
df_NDVsT2DBeta_common <- df_NDVsT2DBeta_male[df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene,]

#get lists of unique genes for each sex
df_NDVsT2DBeta_male_unique <- subset(df_NDVsT2DBeta_male, !df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene)
df_NDVsT2DBeta_female_unique <- subset(df_NDVsT2DBeta_female, !df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene)

#create Venn diagram
venn <- VennDiagram::draw.pairwise.venn(nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_male_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_common),
                                        category = c("Female", "Male"), 
                                        inverted = TRUE,
                                        ext.text = FALSE,
                                        lwd = rep(8, 2),
                                        fill = c("white", "white"),
                                        col = c("#E7B81E", "#84B429"),
                                        alpha = rep(0.5, 2), 
                                        cat.pos = c(47, -30),
                                        fontfamily = rep("sans", 3),
                                        cat.fontfamily = rep("sans", 2),
                                        cex = c(2, 2, 2),
                                        cat.cex = rep(2, 2),
                                        cat.dist = c(0.075, 0.075),
                                        margin = 0)

#clear existing graphics
grid::grid.newpage()

#draw Venn diagram
grid::grid.draw(venn)

#add percent annotations
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_female_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[5]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_common)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[7]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_male_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[6]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)


## Fig 1B ##

#clear environment
rm(list = ls())

#import files
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set undefined column name as gene
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance and fold change
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05 & avg_logFC < 0)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05 & avg_logFC < 0)

#get list of common genes
df_NDVsT2DBeta_common <- df_NDVsT2DBeta_male[df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene,]

#get lists of unique genes for each sex
df_NDVsT2DBeta_male_unique <- subset(df_NDVsT2DBeta_male, !df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene)
df_NDVsT2DBeta_female_unique <- subset(df_NDVsT2DBeta_female, !df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene)

#create Venn diagram
venn <- VennDiagram::draw.pairwise.venn(nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_male_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_common),
                                        category = c("Female", "Male"), 
                                        inverted = TRUE,
                                        ext.text = FALSE,
                                        # lty = rep("1", 2), 
                                        lwd = rep(8, 2),
                                        fill = c("white", "white"),
                                        col = c("#E7B81E", "#84B429"),
                                        alpha = rep(0.5, 2), 
                                        cat.pos = c(42, -31),
                                        fontfamily = rep("sans", 3),
                                        cat.fontfamily = rep("sans", 2),
                                        cex = c(2, 2, 2),
                                        cat.cex = rep(2, 2),
                                        cat.dist = c(0.065, 0.065),
                                        margin = 0)

#clear existing graphics
grid::grid.newpage()

#draw Venn diagram
grid::grid.draw(venn)

#add percent annotations
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_female_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[5]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_common)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[7]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_male_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[6]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)


## Fig 1C ##

#clear environment
rm(list = ls())

#import files
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set undefined column name as gene
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance and fold change
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05 & avg_logFC > 0)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05 & avg_logFC > 0)

#get list of common genes
df_NDVsT2DBeta_common <- df_NDVsT2DBeta_male[df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene,]

#get lists of unique genes for each sex
df_NDVsT2DBeta_male_unique <- subset(df_NDVsT2DBeta_male, !df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene)
df_NDVsT2DBeta_female_unique <- subset(df_NDVsT2DBeta_female, !df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene)

#create Venn diagram
venn <- VennDiagram::draw.pairwise.venn(nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_male_unique) + nrow(df_NDVsT2DBeta_common), 
                                        nrow(df_NDVsT2DBeta_common),
                                        category = c("Female", "Male"), 
                                        inverted = TRUE,
                                        ext.text = FALSE,
                                        lwd = rep(8, 2),
                                        fill = c("white", "white"),
                                        col = c("#E7B81E", "#84B429"),
                                        alpha = rep(0.5, 2), 
                                        cat.pos = c(48, -28),
                                        fontfamily = rep("sans", 3),
                                        cat.fontfamily = rep("sans", 2),
                                        cex = c(2, 2, 2),
                                        cat.cex = rep(2, 2),
                                        cat.dist = c(0.08, 0.08),
                                        margin = 0)

#clear existing graphics
grid::grid.newpage()

#draw Venn diagram
grid::grid.draw(venn)

#add percent annotations
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_female_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[5]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_common)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[7]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)
grid::grid.text(paste0("(", round(nrow(df_NDVsT2DBeta_male_unique)/(nrow(df_NDVsT2DBeta_common) + nrow(df_NDVsT2DBeta_female_unique) + nrow(df_NDVsT2DBeta_male_unique)) * 100, 0), "%)"), x = as.numeric(venn[[6]][2]), y = as.numeric(venn[[7]][3]) - 0.06, gp = grid::gpar(col = "black", fontsize = 16), hjust = 0.5)


## Fig 1D ##

#clear environment
rm(list = ls())

#import files
df_NDvsT2D_common_T2D <- read.csv(".\\Data\\ND_vs_T2D - Common_T2D_0.05_Reactome_Results.csv")
df_NDvsT2D_common_ND <- read.csv(".\\Data\\ND_vs_T2D - Common_ND_0.05_Reactome_Results.csv")

#fix column names
colnames(df_NDvsT2D_common_T2D) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")
colnames(df_NDvsT2D_common_ND) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")

#select top 10 pathways
df_NDvsT2D_common_T2D <- df_NDvsT2D_common_T2D[c(1:10),]
df_NDvsT2D_common_ND <- df_NDvsT2D_common_ND[c(1:10),]

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_NDvsT2D_common_T2D$Gene.ratio <- df_NDvsT2D_common_T2D$Entities.found/45
df_NDvsT2D_common_ND$Gene.ratio <- df_NDvsT2D_common_ND$Entities.found/9

#get -log10 of FDR
df_NDvsT2D_common_T2D$transf.log.1.FDR <- -log10(df_NDvsT2D_common_T2D$Entities.FDR)
df_NDvsT2D_common_ND$transf.log.1.FDR <- -log10(df_NDvsT2D_common_ND$Entities.FDR)

#order by gene ratio
df_NDvsT2D_common_T2D_ordered <- dplyr::arrange(df_NDvsT2D_common_T2D, Gene.ratio)
df_NDvsT2D_common_ND_ordered <- dplyr::arrange(df_NDvsT2D_common_ND, desc(Gene.ratio))

#order by -log10(FDR)
df_NDvsT2D_common_T2D_ordered <- dplyr::arrange(df_NDvsT2D_common_T2D_ordered, transf.log.1.FDR)
df_NDvsT2D_common_ND_ordered <- dplyr::arrange(df_NDvsT2D_common_ND_ordered, desc(transf.log.1.FDR))

#flip ND to negative side
df_NDvsT2D_common_ND_ordered$transf.log.1.FDR <- -df_NDvsT2D_common_ND_ordered$transf.log.1.FDR

#combine ND and T2D
df_NDvsT2D_common_ordered_combined <- rbind(df_NDvsT2D_common_ND_ordered, df_NDvsT2D_common_T2D_ordered)

#categorize significance
df_NDvsT2D_common_ordered_combined$significance <- ifelse(df_NDvsT2D_common_ordered_combined$Entities.FDR < 0.01, "< 0.01", ifelse(df_NDvsT2D_common_ordered_combined$Entities.FDR > 0.01 & df_NDvsT2D_common_ordered_combined$Entities.FDR < 0.05, "0.01-0.05", "> 0.05"))

#set factor levels
df_NDvsT2D_common_ordered_combined$Pathway.name <- factor(df_NDvsT2D_common_ordered_combined$Pathway.name, levels = df_NDvsT2D_common_ordered_combined$Pathway.name)

#set x-axis limits
left_limit <- -5
right_limit <- 7.5

#create plot
ggplot2::ggplot(df_NDvsT2D_common_ordered_combined, ggplot2::aes(x = transf.log.1.FDR, y = Pathway.name)) +
  
  ggplot2::geom_point(alpha = 0) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1", 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1,
                      fill = "#848484",
                      data = df_NDvsT2D_common_ordered_combined[df_NDvsT2D_common_ordered_combined$significance == "< 0.01",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1",
                                   shape = significance, 
                                   size = Gene.ratio),
                      stroke = 1,
                      data = df_NDvsT2D_common_ordered_combined[df_NDvsT2D_common_ordered_combined$significance == "0.01-0.05",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(shape = significance, 
                                   size = Gene.ratio), 
                      stroke = 1,
                      colour = "#848484", 
                      data = df_NDvsT2D_common_ordered_combined[df_NDvsT2D_common_ordered_combined$significance == "> 0.05",]
  ) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, 20.5), expand = FALSE, clip = "off") +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = 20.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = 10.5, yend = 10.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = 21, yend = 21, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = 21, yend = 21, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -0.06*2, xend = 0.06*2, y = 21, yend = 21, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = 21, label = c("ND")) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = 21, label = c("T2D")) +
  ggplot2::annotate(geom = "text", x = left_limit-0.65, y = 21, label = c("Condition"), hjust = 1) +
  
  ggplot2::scale_color_manual(values = c("#848484")) +
  
  ggplot2::scale_size_continuous(limits = c(0, 0.58),
                                 breaks = seq(0.08, 0.56, by = 0.16)) +
  
  ggplot2::scale_shape_manual(values = c(21, 21)) +
  
  ggplot2::scale_y_discrete(labels = function(y) stringr::stringr::str_wrap(y, width = 40)) +
  
  ggplot2::labs(y = "Pathway", x = bquote(-log[10]*"(FDR)"), colour = "Gene Ratio", size = "Gene Ratio", shape = "Significance") +
  
  ggplot2::guides(
    colour = "none",
    size = ggplot2::guide_legend(override.aes = list(colour = "#848484"), reverse = T),
    shape = ggplot2::guide_legend(override.aes = list(colour = "black", fill = c("black", "white"), size = 4))
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
                 axis.text.y = ggplot2::element_text(size = 11, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 21/5)


## Fig 1E ##

#clear environment
rm(list = ls())

#import files
df_NDvsT2D_female_T2D <- read.csv(".\\Data\\ND_vs_T2D - Female_T2D_0.05_Reactome_Results.csv")
df_NDvsT2D_female_ND <- read.csv(".\\Data\\ND_vs_T2D - Female_ND_0.05_Reactome_Results.csv")

#fix column names
colnames(df_NDvsT2D_female_T2D) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")
colnames(df_NDvsT2D_female_ND) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")

#select top 10 pathways
df_NDvsT2D_female_T2D <- df_NDvsT2D_female_T2D[c(1:10),]
df_NDvsT2D_female_ND <- df_NDvsT2D_female_ND[c(1:10),]

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_NDvsT2D_female_T2D$Gene.ratio <- df_NDvsT2D_female_T2D$Entities.found/81
df_NDvsT2D_female_ND$Gene.ratio <- df_NDvsT2D_female_ND$Entities.found/66

#get -log10 of FDR
df_NDvsT2D_female_T2D$transf.log.1.FDR <- -log10(df_NDvsT2D_female_T2D$Entities.FDR)
df_NDvsT2D_female_ND$transf.log.1.FDR <- -log10(df_NDvsT2D_female_ND$Entities.FDR)

#order by gene ratio
df_NDvsT2D_female_T2D_ordered <- dplyr::arrange(df_NDvsT2D_female_T2D, Gene.ratio)
df_NDvsT2D_female_ND_ordered <- dplyr::arrange(df_NDvsT2D_female_ND, desc(Gene.ratio))

#order by -log10(FDR)
df_NDvsT2D_female_T2D_ordered <- dplyr::arrange(df_NDvsT2D_female_T2D_ordered, transf.log.1.FDR)
df_NDvsT2D_female_ND_ordered <- dplyr::arrange(df_NDvsT2D_female_ND_ordered, desc(transf.log.1.FDR))

#flip ND to negative side
df_NDvsT2D_female_ND_ordered$transf.log.1.FDR <- -df_NDvsT2D_female_ND_ordered$transf.log.1.FDR

#combine ND and T2D
df_NDvsT2D_female_ordered_combined <- rbind(df_NDvsT2D_female_ND_ordered, df_NDvsT2D_female_T2D_ordered)

#categorize significance
df_NDvsT2D_female_ordered_combined$significance <- ifelse(df_NDvsT2D_female_ordered_combined$Entities.FDR < 0.01, "< 0.01", ifelse(df_NDvsT2D_female_ordered_combined$Entities.FDR > 0.01 & df_NDvsT2D_female_ordered_combined$Entities.FDR < 0.05, "0.01-0.05", "> 0.05"))

#set factor levels
df_NDvsT2D_female_ordered_combined$Pathway.name <- factor(df_NDvsT2D_female_ordered_combined$Pathway.name, levels = df_NDvsT2D_female_ordered_combined$Pathway.name)

#set x-axis limits
left_limit <- -15
right_limit <- 5

#create plot
ggplot2::ggplot(df_NDvsT2D_female_ordered_combined, ggplot2::aes(x = transf.log.1.FDR, y = Pathway.name)) +
  
  ggplot2::geom_point(alpha = 0) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1", 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1,
                      fill = "#E7B91E",
                      data = df_NDvsT2D_female_ordered_combined[df_NDvsT2D_female_ordered_combined$significance == "< 0.01",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1",
                                   shape = significance, 
                                   size = Gene.ratio), 
                      stroke = 1,
                      data = df_NDvsT2D_female_ordered_combined[df_NDvsT2D_female_ordered_combined$significance == "0.01-0.05",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(shape = significance, 
                                   size = Gene.ratio), 
                      stroke = 1,
                      colour = "grey", 
                      data = df_NDvsT2D_female_ordered_combined[df_NDvsT2D_female_ordered_combined$significance == "> 0.05",]
  ) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, 20.5), expand = FALSE, clip = "off") +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = 20.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = 10.5, yend = 10.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = 21, yend = 21, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = 21, yend = 21, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -0.06*4, xend = 0.06*4, y = 21, yend = 21, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = 21, label = c("ND")) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = 21, label = c("T2D")) +
  ggplot2::annotate(geom = "text", x = left_limit-0.65, y = 21, label = c("Condition"), hjust = 1) +
  
  ggplot2::scale_color_manual(values = c("#E7B91E")) +
  
  ggplot2::scale_size_continuous(limits = c(0, 0.58),
                                 breaks = seq(0.08, 0.56, by = 0.16)) +
  
  ggplot2::scale_shape_manual(values = c(21, 21)) +
  
  ggplot2::scale_y_discrete(labels = function(y) stringr::stringr::str_wrap(y, width = 40)) +
  
  ggplot2::labs(y = "Pathway", x = bquote(-log[10]*"(FDR)"), colour = "Gene Ratio", size = "Gene Ratio", shape = "Significance") +
  
  ggplot2::guides(
    colour = "none",
    size = ggplot2::guide_legend(override.aes = list(colour = "#E7B91E"), reverse = T),
    shape = ggplot2::guide_legend(override.aes = list(colour = "black", fill = c("black", "white"), size = 4))
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size = 12, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 21/5)


## Fig 1F ##

#clear environment
rm(list = ls())

#import files
df_NDvsT2D_male_T2D <- read.csv(".\\Data\\ND_vs_T2D - Male_T2D_0.05_Reactome_Results.csv")
df_NDvsT2D_male_ND <- read.csv(".\\Data\\ND_vs_T2D - Male_ND_0.05_Reactome_Results.csv")

#fix column names
colnames(df_NDvsT2D_male_T2D) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")
colnames(df_NDvsT2D_male_ND) <- c("Pathway.identifier",	"Pathway.name",	"Entities.found",	"Entities.total",	"Entities.ratio",	"Entities.pValue",	"Entities.FDR",	"Reactions.found",	"Reactions.total",	"Reactions.ratio",	"Species.identifier",	"Species.name",	"Submitted.entities.found",	"Mapped.entities",	"Found.reaction.identifiers")


#select top 10 pathways
df_NDvsT2D_male_T2D <- df_NDvsT2D_male_T2D[c(1:10),]
df_NDvsT2D_male_ND <- df_NDvsT2D_male_ND[c(1:10),]


#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_NDvsT2D_male_T2D$Gene.ratio <- df_NDvsT2D_male_T2D$Entities.found/264
df_NDvsT2D_male_ND$Gene.ratio <- df_NDvsT2D_male_ND$Entities.found/104


#get -log10 of FDR
df_NDvsT2D_male_T2D$transf.log.1.FDR <- -log10(df_NDvsT2D_male_T2D$Entities.FDR)
df_NDvsT2D_male_ND$transf.log.1.FDR <- -log10(df_NDvsT2D_male_ND$Entities.FDR)

#order by gene ratio
df_NDvsT2D_male_T2D_ordered <- dplyr::arrange(df_NDvsT2D_male_T2D, Gene.ratio)
df_NDvsT2D_male_ND_ordered <- dplyr::arrange(df_NDvsT2D_male_ND, desc(Gene.ratio))

#order by -log10(FDR)
df_NDvsT2D_male_T2D_ordered <- dplyr::arrange(df_NDvsT2D_male_T2D_ordered, transf.log.1.FDR)
df_NDvsT2D_male_ND_ordered <- dplyr::arrange(df_NDvsT2D_male_ND_ordered, desc(transf.log.1.FDR))

#flip ND to negative side
df_NDvsT2D_male_ND_ordered$transf.log.1.FDR <- -df_NDvsT2D_male_ND_ordered$transf.log.1.FDR

#combine ND and T2D
df_NDvsT2D_male_ordered_combined <- rbind(df_NDvsT2D_male_ND_ordered, df_NDvsT2D_male_T2D_ordered)

#categorize significance
df_NDvsT2D_male_ordered_combined$significance <- ifelse(df_NDvsT2D_male_ordered_combined$Entities.FDR < 0.01, "< 0.01", ifelse(df_NDvsT2D_male_ordered_combined$Entities.FDR > 0.01 & df_NDvsT2D_male_ordered_combined$Entities.FDR < 0.05, "0.01-0.05", "> 0.05"))

#set factor levels
df_NDvsT2D_male_ordered_combined$Pathway.name <- factor(df_NDvsT2D_male_ordered_combined$Pathway.name, levels = df_NDvsT2D_male_ordered_combined$Pathway.name)

#set x-axis limits
left_limit <- -5
right_limit <- 15

#create plot
ggplot2::ggplot(df_NDvsT2D_male_ordered_combined, ggplot2::aes(x = transf.log.1.FDR, y = Pathway.name)) +
  
  ggplot2::geom_point(alpha = 0) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1", 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1,
                      fill = "#84B429",
                      data = df_NDvsT2D_male_ordered_combined[df_NDvsT2D_male_ordered_combined$significance == "< 0.01",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(colour = "1",
                                   shape = significance, 
                                   size = Gene.ratio), 
                      stroke = 1,
                      data = df_NDvsT2D_male_ordered_combined[df_NDvsT2D_male_ordered_combined$significance == "0.01-0.05",]
  ) +
  
  ggplot2::geom_point(ggplot2::aes(shape = significance, 
                                   size = Gene.ratio), 
                      colour = "grey", 
                      stroke = 1,
                      data = df_NDvsT2D_male_ordered_combined[df_NDvsT2D_male_ordered_combined$significance == "> 0.05",]
  ) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, 20.5), expand = FALSE, clip = "off") +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = 20.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = 10.5, yend = 10.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = 21, yend = 21, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = 21, yend = 21, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -0.06*4, xend = 0.06*4, y = 21, yend = 21, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = 21, label = c("ND")) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = 21, label = c("T2D")) +
  ggplot2::annotate(geom = "text", x = left_limit-0.65, y = 21, label = c("Condition"), hjust = 1) +
  
  ggplot2::scale_color_manual(values = c("#84B429")) +
  
  ggplot2::scale_size_continuous(limits = c(0, 0.58),
                                 breaks = seq(0.08, 0.56, by = 0.16)) +
  
  ggplot2::scale_shape_manual(values = c(21, 21)) +
  
  ggplot2::scale_y_discrete(labels = function(y) stringr::stringr::str_wrap(y, width = 40)) +
  
  ggplot2::labs(y = "Pathway", x = bquote(-log[10]*"(FDR)"), colour = "Gene Ratio", size = "Gene Ratio", shape = "Significance") +
  
  ggplot2::guides(
    colour = "none",
    size = ggplot2::guide_legend(override.aes = list(colour = "#84B429"), reverse = T),
    shape = ggplot2::guide_legend(override.aes = list(colour = "black", fill = c("black", "white"),size = 4))
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size = 12, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 21/5)


### Figure 1 - Supplement 1 ###

## Fig 1 S1A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set gene column name
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#get common gene list
df_NDVsT2DBeta_common <- subset(df_NDVsT2DBeta_male, df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene)

#get logFC for common genes in each sex
df_NDVsT2DBeta_male_common <- subset(df_NDVsT2DBeta_male, df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_common$gene)
df_NDVsT2DBeta_female_common <- subset(df_NDVsT2DBeta_female, df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_common$gene)

#unique column name for each sex
colnames(df_NDVsT2DBeta_male_common)[which(colnames(df_NDVsT2DBeta_male_common) == "avg_logFC")] <- "avg_logFC_M"
colnames(df_NDVsT2DBeta_female_common)[which(colnames(df_NDVsT2DBeta_female_common) == "avg_logFC")] <- "avg_logFC_F"

colnames(df_NDVsT2DBeta_male_common)[which(colnames(df_NDVsT2DBeta_male_common) == "p_val_adj")] <- "p_val_adj_M"
colnames(df_NDVsT2DBeta_female_common)[which(colnames(df_NDVsT2DBeta_female_common) == "p_val_adj")] <- "p_val_adj_F"

#merge common female and male
df_merge_common <- merge(df_NDVsT2DBeta_female_common[, c("gene", "avg_logFC_F", "p_val_adj_F")], df_NDVsT2DBeta_male_common[, c("gene", "avg_logFC_M", "p_val_adj_M")], by = "gene")

#make matrix of logFC
matrix_logFC_common <- data.matrix(df_merge_common[, c("avg_logFC_F", "avg_logFC_M")])

#set sex column names
colnames(matrix_logFC_common) <- c("Male", "Female")

#set gene as rownames
rownames(matrix_logFC_common) <- df_merge_common$gene

#create heatmap
heatmap <- ComplexHeatmap::Heatmap(matrix_logFC_common,
                                   name = "matrix_logFC_heatmap",
                                   row_names_side = "right",
                                   row_names_gp = grid::gpar(fontface = "italic", fontsize = 8),
                                   show_row_dend = FALSE,
                                   top_annotation = ComplexHeatmap::columnAnnotation(
                                     empty = ComplexHeatmap::anno_empty(border = FALSE)),
                                   column_order = 1:ncol(matrix_logFC_common),
                                   show_column_names =  FALSE,
                                   width = ncol(matrix_logFC_common)*grid::unit(40, "mm"),
                                   height = nrow(matrix_logFC_common)*grid::unit(3.5, "mm"),
                                   border_gp = grid::gpar(col = "black"),
                                   col = circlize::colorRamp2(c(-1.25, 0, 1.25), c("blue", "white", "red")),
                                   na_col = "black",
                                   heatmap_legend_param = list(title = NULL,
                                                               title = "logFC",
                                                               title_position = "lefttop",
                                                               legend_height = grid::unit(8, "cm"),
                                                               direction = "horizontal",
                                                               legend_width = grid::unit(4, "cm")
                                   )
)

#draw heatmap
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.015,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x)*0.5,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E", 
                                col = "#E7B81E", 
                                lwd = 0))

grid::grid.text("Female",
                x = 0.25,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

grid::grid.rect(x = 0.5,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x)*0.5,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0))

grid::grid.text("Male",
                x = 0.75,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

grid::grid.rect(x = 0.5,
                y = 0,
                width = grid::unit(0.03, "inches"),
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("center", "bottom"),
                gp = grid::gpar(fill = "white", 
                                col = "white", 
                                lwd = 1))

#select legend viewport
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("logFC",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


## Fig 1 S1B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set gene column name
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#filter for unique genes
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, !df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene)

#star unique genes
df_NDVsT2DBeta_female$gene <- ifelse(!df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene, paste(df_NDVsT2DBeta_female$gene, "*"), df_NDVsT2DBeta_female$gene)

#select gene, logFC, and colour index columns
df_NDVsT2DBeta_female <- df_NDVsT2DBeta_female[, c("gene", "avg_logFC")]

#sort be abs(logFC)
df_NDVsT2DBeta_female <- dplyr::arrange(df_NDVsT2DBeta_female, desc(abs(df_NDVsT2DBeta_female$avg_logFC)))

#select top 60
df_NDVsT2DBeta_female <- df_NDVsT2DBeta_female[1:60,]

#make matrix of logFC
matrix_logFC_female <- data.matrix(df_NDVsT2DBeta_female[, "avg_logFC"])

#set sex column names
colnames(matrix_logFC_female) <- c("Female")

#set gene as rownames
rownames(matrix_logFC_female) <- df_NDVsT2DBeta_female$gene  

#create heatmap
heatmap <- ComplexHeatmap::Heatmap(matrix_logFC_female,
                                   name = "matrix_logFC_heatmap",
                                   row_names_side = "right",
                                   row_names_gp = grid::gpar(fontface = "italic", 
                                                             fontsize = 8),
                                   show_row_dend = FALSE,
                                   top_annotation = ComplexHeatmap::columnAnnotation(
                                     empty = ComplexHeatmap::anno_empty(border = FALSE)
                                   ),
                                   column_order = 1:ncol(matrix_logFC_female),
                                   show_column_names =  FALSE,
                                   width = ncol(matrix_logFC_female)*grid::unit(40, "mm"),
                                   height = 71*grid::unit(3.5, "mm"),
                                   border_gp = grid::gpar(col = "black"),
                                   col = circlize::colorRamp2(c(-1.25, 0, 1.25), c("blue", "white", "red")),
                                   na_col = "black",
                                   heatmap_legend_param = list(title = NULL,
                                                               title_position = "lefttop",
                                                               legend_height = grid::unit(8, "cm"),
                                                               direction = "horizontal",
                                                               legend_width = grid::unit(4, "cm")
                                   )
)

#draw heatmap
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.03,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E", 
                                col = "#E7B81E", 
                                lwd = 0))

grid::grid.text("Female",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

#select legend viewport
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("logFC",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


## Fig 1 S1C ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set gene column name
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#filter for unique genes
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, !df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene)

#star unique genes
df_NDVsT2DBeta_male$gene <- ifelse(!df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene, paste(df_NDVsT2DBeta_male$gene, "*"), df_NDVsT2DBeta_male$gene)

#select gene, logFC, and colour index columns
df_NDVsT2DBeta_male <- df_NDVsT2DBeta_male[, c("gene", "avg_logFC")]

#sort be abs(logFC)
df_NDVsT2DBeta_male <- dplyr::arrange(df_NDVsT2DBeta_male, desc(abs(df_NDVsT2DBeta_male$avg_logFC)))

#select top 60
df_NDVsT2DBeta_male <- df_NDVsT2DBeta_male[1:60,]

#create label colour list
male_unique_colour_index <- df_NDVsT2DBeta_male$unique

#make matrix of logFC
matrix_logFC_male <- data.matrix(df_NDVsT2DBeta_male[, "avg_logFC"])

#set sex column names
colnames(matrix_logFC_male) <- c("Male")

#set gene as rownames
rownames(matrix_logFC_male) <- df_NDVsT2DBeta_male$gene  

#create heatmap
heatmap <- ComplexHeatmap::Heatmap(matrix_logFC_male,
                                   name = "matrix_logFC_heatmap",
                                   row_names_side = "right",
                                   row_names_gp = grid::gpar(fontface = "italic", 
                                                             fontsize = 8,
                                                             col = male_unique_colour_index),
                                   show_row_dend = FALSE,
                                   top_annotation = ComplexHeatmap::columnAnnotation(
                                     empty = ComplexHeatmap::anno_empty(border = FALSE)),
                                   column_order = 1:ncol(matrix_logFC_male),
                                   show_column_names =  FALSE,
                                   width = ncol(matrix_logFC_male)*grid::unit(40, "mm"),
                                   height = 71*grid::unit(3.5, "mm"),
                                   border_gp = grid::gpar(col = "black"),
                                   col = circlize::colorRamp2(c(-1.25, 0, 1.25), c("blue", "white", "red")),
                                   na_col = "black",
                                   heatmap_legend_param = list(title = NULL,
                                                               title_position = "lefttop",
                                                               legend_height = grid::unit(8, "cm"),
                                                               direction = "horizontal",
                                                               legend_width = grid::unit(4, "cm")
                                   )
)

#draw heatmap
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.03,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0
                ))

grid::grid.text("Male",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

#select legend viewport
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("logFC",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


### Figure 1 - Supplement 2 ###

## Fig 1 S2A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set gene column name
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#index colour of unique genes
df_NDVsT2DBeta_female$unique <- ifelse(!df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene, "red", "black")

#star unique genes
df_NDVsT2DBeta_female$gene <- ifelse(!df_NDVsT2DBeta_female$gene %in% df_NDVsT2DBeta_male$gene, paste(df_NDVsT2DBeta_female$gene, "*"), df_NDVsT2DBeta_female$gene)

#select gene, logFC, and colour index columns
df_NDVsT2DBeta_female <- df_NDVsT2DBeta_female[, c("gene", "avg_logFC", "unique")]

#sort be abs(logFC)
df_NDVsT2DBeta_female <- dplyr::arrange(df_NDVsT2DBeta_female, desc(abs(df_NDVsT2DBeta_female$avg_logFC)))

#select top 60
df_NDVsT2DBeta_female <- df_NDVsT2DBeta_female[1:60,]

#create label colour list
female_unique_colour_index <- df_NDVsT2DBeta_female$unique

#make matrix of logFC
matrix_logFC_female <- data.matrix(df_NDVsT2DBeta_female[, "avg_logFC"])

#set sex column names
colnames(matrix_logFC_female) <- c("Female")

#set gene as rownames
rownames(matrix_logFC_female) <- df_NDVsT2DBeta_female$gene  

#create heatmap
heatmap <- ComplexHeatmap::Heatmap(matrix_logFC_female,
                                   name = "matrix_logFC_heatmap",
                                   row_names_side = "right",
                                   row_names_gp = grid::gpar(fontface = "italic", 
                                                             fontsize = 8,
                                                             col = female_unique_colour_index),
                                   show_row_dend = FALSE,
                                   top_annotation = ComplexHeatmap::columnAnnotation(
                                     empty = ComplexHeatmap::anno_empty(border = FALSE)),
                                   column_order = 1:ncol(matrix_logFC_female),
                                   show_column_names =  FALSE,
                                   width = ncol(matrix_logFC_female)*grid::unit(40, "mm"),
                                   height = 71*grid::unit(3.5, "mm"),
                                   border_gp = grid::gpar(col = "black"),
                                   col = circlize::colorRamp2(c(-1.25, 0, 1.25), c("blue", "white", "red")),
                                   na_col = "black",
                                   heatmap_legend_param = list(title = NULL,
                                                               title_position = "lefttop",
                                                               legend_height = grid::unit(8, "cm"),
                                                               direction = "horizontal",
                                                               legend_width = grid::unit(4, "cm")
                                   )
)

#draw heatmap
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.03,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E", 
                                col = "#E7B81E", 
                                lwd = 0))

grid::grid.text("Female",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

#select legend viewport
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("logFC",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


## Fig 1 S2B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_NDVsT2DBeta_male <- read.csv(".\\Data\\NDVsT2DBeta_male.csv")
df_NDVsT2DBeta_female <- read.csv(".\\Data\\NDVsT2DBeta_female.csv")

#set gene column name
colnames(df_NDVsT2DBeta_male)[1] <- "gene"
colnames(df_NDVsT2DBeta_female)[1] <- "gene"

#filter for significance
df_NDVsT2DBeta_male <- subset(df_NDVsT2DBeta_male, p_val_adj < 0.05)
df_NDVsT2DBeta_female <- subset(df_NDVsT2DBeta_female, p_val_adj < 0.05)

#index colour of unique genes
df_NDVsT2DBeta_male$unique <- ifelse(!df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene, "red", "black")

#star unique genes
df_NDVsT2DBeta_male$gene <- ifelse(!df_NDVsT2DBeta_male$gene %in% df_NDVsT2DBeta_female$gene, paste(df_NDVsT2DBeta_male$gene, "*"), df_NDVsT2DBeta_male$gene)

#select gene, logFC, and colour index columns
df_NDVsT2DBeta_male <- df_NDVsT2DBeta_male[, c("gene", "avg_logFC", "unique")]

#sort be abs(logFC)
df_NDVsT2DBeta_male <- dplyr::arrange(df_NDVsT2DBeta_male, desc(abs(df_NDVsT2DBeta_male$avg_logFC)))

#select top 60
df_NDVsT2DBeta_male <- df_NDVsT2DBeta_male[1:60,]

#create label colour list
male_unique_colour_index <- df_NDVsT2DBeta_male$unique

#make matrix of logFC
matrix_logFC_male <- data.matrix(df_NDVsT2DBeta_male[, "avg_logFC"])

#set sex column names
colnames(matrix_logFC_male) <- c("Male")

#set gene as rownames
rownames(matrix_logFC_male) <- df_NDVsT2DBeta_male$gene  

#create heatmap
heatmap <- ComplexHeatmap::Heatmap(matrix_logFC_male,
                                   name = "matrix_logFC_heatmap",
                                   row_names_side = "right",
                                   row_names_gp = grid::gpar(fontface = "italic", 
                                                             fontsize = 8,
                                                             col = male_unique_colour_index),
                                   show_row_dend = FALSE,
                                   top_annotation = ComplexHeatmap::columnAnnotation(
                                     empty = ComplexHeatmap::anno_empty(border = FALSE)),
                                   column_order = 1:ncol(matrix_logFC_male),
                                   show_column_names =  FALSE,
                                   width = ncol(matrix_logFC_male)*grid::unit(40, "mm"),
                                   height = 71*grid::unit(3.5, "mm"),
                                   border_gp = grid::gpar(col = "black"),
                                   col = circlize::colorRamp2(c(-1.25, 0, 1.25), c("blue", "white", "red")),
                                   na_col = "black",
                                   heatmap_legend_param = list(title = NULL,
                                                               title_position = "lefttop",
                                                               legend_height = grid::unit(8, "cm"),
                                                               direction = "horizontal",
                                                               legend_width = grid::unit(4, "cm")
                                   )
)

#draw heatmap
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.03,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0))

grid::grid.text("Male",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 10))

#select legend viewport
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("logFC",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


### Figure 2 ###

## Fig 2A ##

#clear environment
rm(list = ls())

#import file
df_rlog <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.rlog.csv")

#number replicated genes
df_rlog[df_rlog$Gene == "44256", "Gene"] <- c("44256 (1)", "44256 (2)")
df_rlog[df_rlog$Gene == "44257", "Gene"] <- c("44257 (1)", "44257 (2)")

#set rownames as genes
rownames(df_rlog) <- df_rlog[, 1]
df_rlog <- df_rlog[, -1]

#transform data
df_rlog_transformed <- as.data.frame(t(df_rlog))

#annotate sex
df_rlog_transformed$sex <- c(rep("Female", 8), rep("Male", 9))

#create PCA data
pc <- pcaMethods::pca(df_rlog_transformed)

#add annotations to PCA data
pc_plot <- merge(pcaMethods::scores(pc), df_rlog_transformed, by = 0)

#plot PCA data
ggplot2::ggplot(pc_plot, ggplot2::aes(PC1, PC2, colour = sex)) +
  
  ggplot2::geom_hline(yintercept = 0, color = "black") +
  ggplot2::geom_vline(xintercept = 0, color = "black") +
  
  ggplot2::geom_point(size = 5) +
  ggplot2::stat_ellipse(geom = "polygon", ggplot2::aes(fill = sex), alpha = 0.3, color = NA) +
  
  ggplot2::scale_color_manual("Sex", values = c("#E7B81E", "#84B429")) +  
  ggplot2::scale_fill_manual("Sex", values = c("#E7B81E", "#84B429")) + 
  ggplot2::xlab(paste0("PC1: ", format(round(pc@R2[1] * 100, 0), nsmall = 0), "% of the variance")) +
  ggplot2::ylab(paste0("PC2: ", format(round(pc@R2[2] * 100, 0), nsmall = 0), "% of the variance")) +
  ggplot2::theme_bw() +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = 0, fill = NA))) + 
  ggplot2::theme(
    legend.title = ggplot2::element_text(size = 30),
    legend.spacing.y = grid::unit(5, "mm"),
    legend.key.height = grid::unit(1.25, "cm"),
    panel.grid = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 28),
    axis.title = ggplot2::element_text(size = 30),
    axis.text = ggplot2::element_text(size = 28, color = "black"),
    axis.ticks = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 28, family = "sans"),
    aspect.ratio = 1/1)


## Fig 2C ##

#clear environment
rm(list = ls())

#import files
df_pathway_genes_female <- read.csv(".\\Data\\F_vs_M - WT_0.01_Top_1000_Female_Reactome_Results.csv")[, c(2, 3:5, 8)]
df_pathway_genes_male <- read.csv(".\\Data\\F_vs_M - WT_0.01_Top_1000_Male_Reactome_Results.csv")[, c(2:5, 8)]

#add column for sex
df_pathway_genes_female$sex <- "Female"
df_pathway_genes_male$sex <- "Male"

#fix column names
colnames(df_pathway_genes_female) <- c("Pathway", "Entities.found", "Entities.total", "Entities.ratio","FDR", "sex")
colnames(df_pathway_genes_male) <- c("Pathway", "Entities.found", "Entities.total", "Entities.ratio","FDR", "sex")

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_pathway_genes_female$Gene.ratio <- df_pathway_genes_female$Entities.found/611
df_pathway_genes_male$Gene.ratio <- df_pathway_genes_male$Entities.found/550

#get -log10 of FDR
df_pathway_genes_female$transf.log.10.FDR <- -log10(df_pathway_genes_female$FDR)
df_pathway_genes_male$transf.log.10.FDR <- -log10(df_pathway_genes_male$FDR)

#order by gene ratio
df_pathway_genes_female <- dplyr::arrange(df_pathway_genes_female, Gene.ratio)
df_pathway_genes_male <- dplyr::arrange(df_pathway_genes_male, desc(Gene.ratio))

#order by -log10(FDR)
df_pathway_genes_female <- dplyr::arrange(df_pathway_genes_female, desc(transf.log.10.FDR))
df_pathway_genes_male <- dplyr::arrange(df_pathway_genes_male, desc(transf.log.10.FDR))

#flip males to negative side
df_pathway_genes_male$transf.log.10.FDR <- -df_pathway_genes_male$transf.log.10.FDR

#select top 20 pathways
df_pathway_genes_female <- head(df_pathway_genes_female, 20)
df_pathway_genes_male <- head(df_pathway_genes_male, 20)

#filter for FDR < 0.05
df_pathway_genes_female <- df_pathway_genes_female[df_pathway_genes_female$FDR < 0.05,] %>% na.omit()
df_pathway_genes_male <- df_pathway_genes_male[df_pathway_genes_male$FDR < 0.05,] %>% na.omit()

#order by -log10(FDR)
df_pathway_genes_female <- dplyr::arrange(df_pathway_genes_female, transf.log.10.FDR)

#combine males and females
df_separated <- rbind(df_pathway_genes_male, df_pathway_genes_female)

#categorize significance
df_separated$significance <- ifelse(df_separated$FDR < 0.01, "< 0.01", ifelse(df_separated$FDR > 0.01 & df_separated$FDR < 0.05, "0.01-0.05", "> 0.05"))

#set pathway factor levels
df_separated$Pathway <- factor(df_separated$Pathway, levels = df_separated$Pathway)

#set significance factor levels
df_separated$significance <- factor(df_separated$significance, levels = sort(unique(df_separated$significance))[c(1, 3, 2)])

#set axis limits and label position
left_limit <- -5
right_limit <- 10

top_label_height_multiplier <- 1.015

num_row <- nrow(df_separated)
num_row_bottom <- length(df_separated$sex[df_separated$sex == "Male"])

#create plot
ggplot2::ggplot(df_separated, ggplot2::aes(x = transf.log.10.FDR, y = Pathway)) +
  
  ggplot2::geom_point(alpha = 0) +
  
  ggplot2::geom_point(ggplot2::aes(colour = sex, 
                                   fill = sex,
                                   size = Gene.ratio, 
                                   shape = significance),
                      stroke = 1,
                      data = df_separated[df_separated$significance == "< 0.01",]) +
  
  ggplot2::geom_point(aes(size = Gene.ratio, 
                          shape = significance),
                      stroke = 1,
                      colour = "grey", 
                      data = df_separated[df_separated$significance == "> 0.05",]) +
  
  ggplot2::geom_point(ggplot2::aes(colour = sex, 
                                   size = Gene.ratio, 
                                   shape = significance),
                      stroke = 1,
                      data = df_separated[df_separated$significance == "0.01-0.05",]) +
  
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, nrow(df_separated) + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = nrow(df_separated) + 0.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = num_row_bottom + 0.5, yend = num_row_bottom + 0.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = num_row * top_label_height_multiplier + 0.75, yend = num_row * top_label_height_multiplier + 0.75, colour = "#84B429", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = num_row * top_label_height_multiplier + 0.75, yend = num_row * top_label_height_multiplier + 0.75, colour = "#E7B81E", size = 8) +
  ggplot2::annotate(geom = "segment", x = -(right_limit - left_limit) * 0.0075, xend = (right_limit - left_limit) * 0.0075, y = num_row * top_label_height_multiplier + 0.75, yend = num_row * top_label_height_multiplier + 0.75, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = num_row * top_label_height_multiplier + 0.75, label = c("Male"), size = 4) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = num_row * top_label_height_multiplier + 0.75, label = c("Female")) +
  ggplot2::annotate(geom = "text", x = left_limit - (right_limit - left_limit) * 0.035, y = num_row * top_label_height_multiplier + 0.75, label = c("Sex"), hjust = 1) +
  
  ggplot2::scale_color_manual(values = c("#E7B81E", "#84B429")) +
  
  ggplot2::scale_fill_manual(values = c("#E7B81E", "#84B429")) +
  
  ggplot2::scale_size_continuous(limits = c(min(df_separated$Gene.ratio) - 0.05, max(df_separated$Gene.ratio) + 0.05),
                                 breaks = seq(0, 0.3, by = 0.05)) +
  
  ggplot2::scale_shape_manual(values = 21) +
  
  ggplot2::scale_y_discrete(labels = gsub("\\\\n", "\n", as.character(df_separated$Pathway))) +
  
  ggplot2::labs(x = bquote(-log[10]*"(FDR)"), colour = "Sex", size = "Gene Ratio", shape = "Significance") +
  
  ggplot2::guides(
    fill = "none",
    shape = guide_legend(override.aes = list(size = 4, 
                                             fill = c("black", "white")), 
                         order = 3),
    colour = guide_legend(reverse = T, override.aes = list(size = 4), 
                          order = 1),
    size = guide_legend(reverse = T, 
                        order = 2)
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.text = ggplot2::element_text(size = 12, colour = "black"),
                 axis.title.x = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = 16/5)


## Fig 2D ##

#clear environment
rm(list = ls())

#import GO database files

#response to endoplasmic reticulum stress
df_GO_ER_stress_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211029_200313_response to endoplasmic reticulum stress_GO_0034976.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_ER_stress_pathway_genes)[1:11] <- colnames(df_GO_ER_stress_pathway_genes)[2:12]
df_GO_ER_stress_pathway_genes <- df_GO_ER_stress_pathway_genes[, -12]

#transcription, DNA-templated
df_GO_transcription_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_184307_transcription, DNA-templated_GO_0006351.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_transcription_pathway_genes)[1:11] <- colnames(df_GO_transcription_pathway_genes)[2:12]
df_GO_transcription_pathway_genes <- df_GO_transcription_pathway_genes[, -12]

#translation
df_GO_translation_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_184449_translation_GO_0006412.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_translation_pathway_genes)[1:11] <- colnames(df_GO_translation_pathway_genes)[2:12]
df_GO_translation_pathway_genes <- df_GO_translation_pathway_genes[, -12]

#protein processing
df_GO_protein_processing_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_185002_protein processing_GO_0016485.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_protein_processing_pathway_genes)[1:11] <- colnames(df_GO_protein_processing_pathway_genes)[2:12]
df_GO_protein_processing_pathway_genes <- df_GO_protein_processing_pathway_genes[, -12]

#protein folding
df_GO_protein_folding_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_030539_protein folding_GO_0006457.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_protein_folding_pathway_genes)[1:11] <- colnames(df_GO_protein_folding_pathway_genes)[2:12]
df_GO_protein_folding_pathway_genes <- df_GO_protein_folding_pathway_genes[, -12]

#secretion
df_GO_secretion_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_185237_secretion_GO_0046903.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_secretion_pathway_genes)[1:11] <- colnames(df_GO_secretion_pathway_genes)[2:12]
df_GO_secretion_pathway_genes <- df_GO_secretion_pathway_genes[, -12]

#protein quality control
df_GO_protein_quality_control_pathway_genes <- read.csv(".\\Data\\GO_term_summary_20211102_225240_protein quality control_GO_0006515.tsv", check.names = F, header = TRUE, sep = "\t", row.names = NULL)
colnames(df_GO_protein_quality_control_pathway_genes)[1:11] <- colnames(df_GO_protein_quality_control_pathway_genes)[2:12]
df_GO_protein_quality_control_pathway_genes <- df_GO_protein_quality_control_pathway_genes[, -12]

#import files
data_rlog <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.rlog.csv")
data_res <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.res.csv")

#get logFC
data_rlog$log2FoldChange <- data_res$log2FoldChange

#get padj
data_rlog$padj <- data_res$padj

#filter for ER stress genes
data_rlog <- subset(data_rlog, tolower(data_rlog$Gene) %in% tolower(df_GO_ER_stress_pathway_genes$Symbol))

#arrange by descending logFC
data_rlog <- dplyr::arrange(data_rlog, desc(log2FoldChange))

#filter for significance
data_rlog <- subset(data_rlog, padj < 0.01)

#add pathway annotation columns
data_rlog$pathway_1 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_transcription_pathway_genes$Symbol), "red", NA)
data_rlog$pathway_2 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_translation_pathway_genes$Symbol), "red", NA)
data_rlog$pathway_3 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_protein_processing_pathway_genes$Symbol), "red", NA)
data_rlog$pathway_4 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_protein_folding_pathway_genes$Symbol), "red", NA)
data_rlog$pathway_5 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_secretion_pathway_genes$Symbol), "red", NA)
data_rlog$pathway_6 <- ifelse(tolower(data_rlog$Gene) %in% tolower(df_GO_protein_quality_control_pathway_genes$Symbol), "red", NA)

#combined pathway annotations
data_rlog$pathway_x <- paste0(data_rlog$pathway_1, 
                              data_rlog$pathway_2,
                              data_rlog$pathway_3,
                              data_rlog$pathway_4,
                              data_rlog$pathway_5,
                              data_rlog$pathway_6)

#separate by sex
data_rlog_M <- data_rlog[data_rlog$log2FoldChange > 0,]
data_rlog_F <- data_rlog[data_rlog$log2FoldChange < 0,]

#manually ordering heatmap
data_rlog_F <- data_rlog_F %>%
  dplyr::mutate(pathway_x = dplyr::case_when(data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[1] ~ 11,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[2] ~ 9,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[3] ~ 6,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[4] ~ 7,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[5] ~ 8,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[6] ~ 10,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[7] ~ 2,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[8] ~ 1,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[9] ~ 5,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[10] ~ 3,
                                             data_rlog_F$pathway_x == sort(unique(data_rlog_F$pathway_x))[11] ~ 4))

data_rlog_M <- data_rlog_M %>% 
  dplyr::mutate(pathway_x = dplyr::case_when(data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[1] ~ 8,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[2] ~ 7,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[3] ~ 6,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[4] ~ 5,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[5] ~ 4,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[6] ~ 1,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[7] ~ 2,
                                             data_rlog_M$pathway_x == sort(unique(data_rlog_M$pathway_x))[8] ~ 3))

#arrange per sex
data_rlog_M <- dplyr::arrange(data_rlog_M, pathway_x)
data_rlog_F <- dplyr::arrange(data_rlog_F, pathway_x)

#combine sexes
data_rlog <- rbind(data_rlog_F, data_rlog_M)

#convert to Z-score matrix
data_rlog$mean_all <- rowMeans(data_rlog[, grepl("RB", colnames(data_rlog))])

data_rlog$sd_all <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 1, sd)

Z_score <- function(x) with(data_rlog, (x - mean_all)/sd_all)

matrix_Z_score_total <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 2, Z_score)

rownames(matrix_Z_score_total) <- data_rlog$Gene

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total",
                                              column_title_gp = grid::gpar(fontsize = 12),
                                              row_names_gp = grid::gpar(fontface = "italic", 
                                                                        fontsize = 9,
                                                                        col = data_rlog$in_UPR_colour), 
                                              row_names_side = "left",
                                              show_column_names = FALSE,
                                              col = circlize::colorRamp2(c(-ceiling(max(na.omit(reshape2::melt(matrix_Z_score_total)$value))), 0, ceiling(max(na.omit(reshape2::melt(matrix_Z_score_total)$value)))),
                                                                         c("blue", "white", "red")),
                                              show_row_dend = FALSE,
                                              layer_fun = function(j, i, x, y, width, height, fill) {
                                                mat = restore_matrix(j, i, x, y)
                                                ind = unique(c(mat[, c(8)]))
                                                grid::grid.rect(x = x[ind] + grid::unit(0.5/ncol(matrix_Z_score_total), "npc"), 
                                                                y = y[ind], 
                                                                width = grid::unit(0.03, "inches"), 
                                                                height = grid::unit(1/nrow(matrix_Z_score_total), "npc"),
                                                                gp = grid::gpar(col = "white")
                                                )
                                              },
                                              
                                              top_annotation = ComplexHeatmap::columnAnnotation(empty2 = ComplexHeatmap::anno_empty(border = FALSE,
                                                                                                                                    width = grid::unit(20, "mm"))),
                                              
                                              
                                              right_annotation = ComplexHeatmap::rowAnnotation(`Transcription` = data_rlog$pathway_1,
                                                                                               
                                                                                               `Translation` = data_rlog$pathway_2,
                                                                                               
                                                                                               `Protein Processing` = data_rlog$pathway_3,
                                                                                               
                                                                                               `Protein Folding` = data_rlog$pathway_4,
                                                                                               
                                                                                               `Secretion` = data_rlog$pathway_5,
                                                                                               
                                                                                               `Protein Quality Control` = data_rlog$pathway_6,
                                                                                               
                                                                                               col = list(`Transcription` = c("red" = scales::hue_pal()(6)[1]),
                                                                                                          `Translation` = c("red" = scales::hue_pal()(6)[2]),
                                                                                                          `Protein Processing` = c("red" = scales::hue_pal()(6)[3]),
                                                                                                          `Protein Folding` = c("red" = scales::hue_pal()(6)[4]),
                                                                                                          `Secretion` = c("red" = scales::hue_pal()(6)[5]),
                                                                                                          `Protein Quality Control` = c("red" = scales::hue_pal()(6)[6])
                                                                                               ),
                                                                                               
                                                                                               annotation_legend_param = list(
                                                                                                 
                                                                                                 `Transcription` = list(
                                                                                                   title = c(""), 
                                                                                                   labels = c("Transcription")),
                                                                                                 `Translation` = list(
                                                                                                   title = c(""), 
                                                                                                   labels = c("Translation")),
                                                                                                 `Protein Processing` = list(
                                                                                                   title = c(""),
                                                                                                   labels = c("Protein Processing")),
                                                                                                 `Protein Folding` = list(
                                                                                                   title = c(""),
                                                                                                   labels = c("Protein Folding")),
                                                                                                 `Secretion` = list(
                                                                                                   title = c(""),
                                                                                                   labels = c("Secretion")),
                                                                                                 `Protein Quality Control` = list(
                                                                                                   title = c(""),
                                                                                                   labels = c("Protein Quality Control"))),
                                                                                               
                                                                                               show_annotation_name = FALSE,
                                                                                               na_col = "white"
                                              ),
                                              
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              row_order = 1:nrow(matrix_Z_score_total),
                                              height = nrow(matrix_Z_score_total)*grid::unit(3.5*0.85, "mm"),
                                              width = ncol(matrix_Z_score_total)*grid::unit(3*1.07, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(title = NULL,
                                                                          title_position = "topcenter",
                                                                          direction = "horizontal",
                                                                          legend_height = grid::unit(6.4, "cm")
                                              )
)

#draw heatmap
ComplexHeatmap::draw(heatmap_list_total, heatmap_legend_side = "bottom", annotation_legend_side = "right", legend_grouping = "original")

#select viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend title
grid::grid.text("Row Z-score",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))

#select viewport and define dimensions
grid::seekViewport("annotation_empty2_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.text("Sex",
                x = -0.015,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 12),
                hjust = 1)

grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x)*8/17,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E", 
                                col = "#E7B81E", 
                                lwd = 0))

grid::grid.text("Female",
                x = (8/2)/17,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 12))

grid::grid.rect(x = 8/17,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x)*9/17,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0))

grid::grid.text("Male",
                x = (8 + 9/2)/17,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 12))

grid::grid.rect(x = 8/17,
                y = 0,
                width = grid::unit(0.03, "inches"),
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("center", "bottom"),
                gp = grid::gpar(fill = "white", 
                                col = "white", 
                                lwd = 1))


## Fig 2E ##

#clear environment
rm(list = ls())

#import files
data_rlog <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.rlog.csv")
data_res <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.res.csv")

#filter for significance
data_res <- subset(data_res, data_res$padj < 0.01)

#filter for ribosomal genes
data_res <- subset(data_res, substr(tolower(data_res$GeneName), 1, 3) == "rpl" | substr(tolower(data_res$GeneName), 1, 3) == "rps")

#filter rlog for same genes
data_rlog <- subset(data_rlog, data_rlog$Gene %in% data_res$GeneName)

#set gene column name
colnames(data_rlog)[1] <- "GeneName"

#add logFC column
data_rlog <- merge(data_rlog, data_res[, c("GeneName", "log2FoldChange")], by = c("GeneName"))

#annotate large and small ribosomal subunits
data_rlog$Ribosome_Size <- ifelse(grepl("Rps", data_rlog$GeneName), "Rps", "Rpl")

#sort by logFC
data_rlog <- data_rlog[with(data_rlog, rev(order(Ribosome_Size, log2FoldChange))),]

#store sample names
sample_names <- colnames(data_rlog)[grepl("RB", colnames(data_rlog))]

#convert to Z-score matrix
data_rlog$mean_all <- rowMeans(data_rlog[, grepl("RB", colnames(data_rlog))])

data_rlog$sd_all <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 1, sd)

Z_score <- function(x) with(data_rlog, (x - mean_all)/sd_all)

df_Z_score_combined <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 2, Z_score)

#transform to horizontal layout
matrix_Z_score_total <- t(df_Z_score_combined)

rownames(matrix_Z_score_total) <- sample_names

colnames(matrix_Z_score_total) <- data_rlog$GeneName

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total",
                                              column_title_gp = grid::gpar(fontsize = 12),
                                              column_names_gp = grid::gpar(fontface = "italic",
                                                                           fontsize = 9),
                                              show_row_names = FALSE,
                                              col = circlize::colorRamp2(c(-ceiling(max(na.omit(reshape2::melt(matrix_Z_score_total)$value))), 0, ceiling(max(na.omit(reshape2::melt(matrix_Z_score_total)$value)))),
                                                                         c("blue", "white", "red")),
                                              show_row_dend = FALSE,
                                              layer_fun = function(j, i, x, y, width, height, fill) {
                                                mat = restore_matrix(j, i, x, y)
                                                ind = unique(c(mat[c(9), ]))
                                                grid::grid.rect(x = x[ind], 
                                                                y = y[ind] + grid::unit(0.5/nrow(matrix_Z_score_total), "npc"),
                                                                width = grid::unit(1/ncol(matrix_Z_score_total), "npc"),
                                                                height = grid::unit(0.03, "inches"), 
                                                                gp = grid::gpar(col = "white")
                                                )
                                              },
                                              left_annotation = ComplexHeatmap::rowAnnotation(empty1 = ComplexHeatmap::anno_empty(border = FALSE,
                                                                                                                                  width = grid::unit(10, "mm"))),
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              row_order = 1:nrow(matrix_Z_score_total),
                                              height = nrow(matrix_Z_score_total)*grid::unit(3, "mm"),
                                              width = ncol(matrix_Z_score_total)*grid::unit(4, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(title = NULL,
                                                                          title_position = "topcenter",
                                                                          direction = "vertical",
                                                                          legend_height = grid::unit(3*17, "mm")
                                              )
)

#draw heatmap
ComplexHeatmap::draw(heatmap_list_total, heatmap_legend_side = "right", legend_grouping = "original")

#select viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend label
grid::grid.text("Column Z-score",
                x = 1.2,
                y = 0.5,
                rot = 90,
                gp = grid::gpar(fontsize = 11))

#select viewport and define dimensions
grid::seekViewport("annotation_empty1_1")
anno_empty_1_loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
anno_empty_1_loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create left labels
grid::grid.rect(x = 0.35,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 0.65,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y)*9/17,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0))

grid::grid.text("Male",
                x = 0.35 + 0.65/2,
                y = (9/2)/17,
                rot = 90,
                gp = grid::gpar(fontsize = 12))

grid::grid.rect(x = 0.35,
                y = 9/17,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 0.65,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y)*8/17,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E"
                                , col = "#E7B81E"
                                , lwd = 0))
grid::grid.text("Female",
                x = 0.35 + 0.65/2,
                y = (9 + 8/2)/17,
                rot = 90,
                gp = grid::gpar(fontsize = 12))

grid::grid.rect(x = 0.35,
                y = 9/17,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 0.65,
                height = grid::unit(0.03, "inches"),
                just = c("left", "center"),
                gp = grid::gpar(fill = "white", 
                                col = "white", 
                                lwd = 1))


### Figure 2 - Supplement 2 ###

## Fig 2 S2A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
data_res <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.res.csv")
data_rlog <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.rlog.csv")
data_counts <- read.csv(".\\Data\\txt.Female_vs_Male.deseq.counts.csv")

#add logFC column
data_rlog$log2FoldChange <- data_res$log2FoldChange

#add padj column
data_rlog$padj <- data_res$padj

#sort by logFC
data_rlog <- dplyr::arrange(data_rlog, desc(log2FoldChange))

#remove blanks
data_rlog <- na.omit(data_rlog)

#remove rows with only zeros
data_rlog <- subset(data_rlog, RB01 != 0)

#filter counts for rows with more than 10 total counts
data_counts <- data_counts[rowSums(data_counts[, c(2:18)]) > 10, ]

#filter rlog for only genes left in counts
data_rlog <- subset(data_rlog, data_rlog$Gene %in% data_counts$Gene)

#set column names
sample_names <- c(paste("Female", sprintf("%02d", 1:8)),  paste("Male", sprintf("%02d", 1:9)))

#convert to Z-score matrix
data_rlog$mean_all <- rowMeans(data_rlog[, grepl("RB", colnames(data_rlog))])

data_rlog$sd_all <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 1, sd)

Z_score <- function(x) with(data_rlog, (x - mean_all)/sd_all)

df_Z_score_total <- apply(data_rlog[, grepl("RB", colnames(data_rlog))], 2, Z_score)

matrix_Z_score_total <- t(df_Z_score_total)

rownames(matrix_Z_score_total) <- sample_names

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total", 
                                              row_names_side = "left", 
                                              show_column_names = FALSE,
                                              row_names_gp = grid::gpar(col = c(rep("#E7B81E", 8), rep("#84B429", 9))),
                                              right_annotation = rowAnnotation(empty2 = anno_empty(border = FALSE,
                                                                                                   width = grid::unit(4, "mm"))),
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              width = grid::unit(180, "mm"), 
                                              height = nrow(matrix_Z_score_total)*grid::unit(6, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(title = "Column\nZ-score",
                                                                          title_position = "topcenter",
                                                                          legend_height = grid::unit(4, "cm")))

#draw heatmap
ComplexHeatmap::draw(heatmap_list_total)


### Figure 3 - Supplement 1 ###

## Fig 3 S1B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import metadata
OPP_Run3_Metadata <- read.csv(".\\Data\\OPP Run3 Metadata.csv", stringsAsFactors = F, check.names = T, sep = ",")[, -9]
OPP_Run4_Metadata <- read.csv(".\\Data\\OPP Run4 Metadata.csv", stringsAsFactors = F, check.names = T, sep = ",")[, -9]

#filter out empty treatment rows in metadata 
OPP_Run3_Metadata <- subset(OPP_Run3_Metadata, !is.na(Treatment))
OPP_Run4_Metadata <- subset(OPP_Run4_Metadata, !is.na(Treatment))

#create group column
OPP_Run3_Metadata$Group <- paste0(OPP_Run3_Metadata$Sex, "_", OPP_Run3_Metadata$Treatment)
OPP_Run4_Metadata$Group <- paste0(OPP_Run4_Metadata$Sex, "_", OPP_Run4_Metadata$Treatment)

#import files
OPP_Run3 <- read.csv(".\\Data\\01292021_OPP_B6_MI_Run3.1_preprocessed.csv", stringsAsFactors = F, check.names = T, sep = ",")
OPP_Run4 <- read.csv(".\\Data\\02122021_OPP_B6_MI_Run4_v3.1_red cell size_preprocessed.csv", stringsAsFactors = F, check.names = T, sep = ",")

#name well column
colnames(OPP_Run3_Metadata)[1] <- "Well_Name"
colnames(OPP_Run4_Metadata)[1] <- "Well_Name"

colnames(OPP_Run3)[2] <- "Well_Name"
colnames(OPP_Run4)[2] <- "Well_Name"

#merge data and metadata
OPP_Run3 <- merge(OPP_Run3, OPP_Run3_Metadata[, c("Well_Name", "Sex", "Group")], by = "Well_Name", all.x = TRUE)
OPP_Run4 <- merge(OPP_Run4, OPP_Run4_Metadata[, c("Well_Name", "Sex", "Group")], by = "Well_Name", all.x = TRUE)

#remove rows with sex missing
OPP_Run3 <- subset(OPP_Run3, Sex != "#N/A")

#get mean per run and sex of OPP-
OPP_Run3_mean_M <- mean(OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Group == "Male_OPP-"])
OPP_Run3_mean_F <- mean(OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Group == "Female_OPP-"])

OPP_Run4_mean_M <- mean(OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Group == "Male_OPP-"])
OPP_Run4_mean_F <- mean(OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Group == "Female_OPP-"])

#correct intensity for OPP-
OPP_Run3$Intensity_Corrected[OPP_Run3$Sex == "Male"] <- OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Sex == "Male"] - OPP_Run3_mean_M
OPP_Run3$Intensity_Corrected[OPP_Run3$Sex == "Female"] <- OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Sex == "Female"] - OPP_Run3_mean_F

OPP_Run4$Intensity_Corrected[OPP_Run4$Sex == "Male"] <- OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Sex == "Male"] - OPP_Run4_mean_M
OPP_Run4$Intensity_Corrected[OPP_Run4$Sex == "Female"] <- OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Sex == "Female"] - OPP_Run4_mean_F

#combine runs 3 and 4
OPP_Run3_4 <- rbind(OPP_Run3, OPP_Run4)

#select only for rows with positive corrected intensity
OPP_Run3_4 <- subset(OPP_Run3_4, OPP_Run3_4$Intensity_Corrected > 0)

#select for FBS+ rows
OPP_Run3_4_FBS_positive <- OPP_Run3_4[grep("[:+:]FBS", OPP_Run3_4$Group), ]

#make rows unique
Group_List_Run3_4_FBS_positive <- unique(OPP_Run3_4_FBS_positive$Group)

#set FBS+ group factor levels
OPP_Run3_4_FBS_positive$Group <- factor(OPP_Run3_4_FBS_positive$Group, levels = c("Female_DMSO_+FBS",
                                                                                  "Male_DMSO_+FBS"))
#calculate stats for FBS+
OPP_Run3_4_FBS_positive_norm <- OPP_Run3_4_FBS_positive %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    sd = sd(Intensity_Corrected),
    mean = mean(Intensity_Corrected),
    n = length(Intensity_Corrected),
    interquartile_bar = IQR(Intensity_Corrected)/2,
    median = median(Intensity_Corrected)
  )

OPP_Run3_4_FBS_positive_norm$SEM <- OPP_Run3_4_FBS_positive_norm$sd/sqrt(OPP_Run3_4_FBS_positive_norm$n)

#set axis limits and label position
left_limit <- 0.5
right_limit <- length(Group_List_Run3_4_FBS_positive) + 0.5

#create plot
ggplot2::ggplot(OPP_Run3_4_FBS_positive_norm, aes(y = Intensity_Corrected, x = Group)) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0, 700000), expand = FALSE, clip = "off") +
  
  ggbeeswarm::geom_quasirandom(
    aes(
      y = Intensity_Corrected,
      colour = Sex
    ),
    size = 0.7,
    data = subset(OPP_Run3_4_FBS_positive, OPP_Run3_4_FBS_positive$Intensity_Corrected < 700000)) +
  
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(x) mean(x),
                        fun.max = function(x) mean(x),
                        geom = "errorbar",
                        ggplot2::aes(group = Group),
                        size = 0.5,
                        colour = "black",
                        width = 0.4,
                        data = OPP_Run3_4_FBS_positive) +
  
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(x) mean(x) - std.error(x)/2,
                        fun.max = function(x) mean(x) + std.error(x)/2,
                        geom = "errorbar",
                        ggplot2::aes(group = Group),
                        colour = "black",
                        width = 0.1,
                        data = OPP_Run3_4_FBS_positive) +
  
  ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE), 
                              limits = c(NA, 700000), 
                              breaks = seq(0, 700000, by = 100000)) +
  
  ggplot2::labs(y = "Integrated Fluorescence Intensity",
                x = "DMSO (+FBS)") +
  
  ggplot2::guides(
    colour = ggplot2::guide_legend(reverse = FALSE, override.aes = list(size = 5))
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(label = round(mean, 1), y = -40000, x = Group),
    colour = "black",
    hjust = 0.5,
    vjust = 0) +
  
  ggplot2::scale_color_manual(values = c("#E7B81E", "#84B429")) +
  
  ggplot2::theme_bw() + 
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(t = 0.25, r = 0, b = 1, l = 0, "cm"),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(vjust = -8),
    axis.ticks.x.bottom = ggplot2::element_blank(),
    axis.ticks.length.y.left = grid::unit(2, "mm"),
    axis.ticks.y.left = ggplot2::element_line(size = 1, colour = "black"),
    panel.border = ggplot2::element_blank(),
    axis.line.y.left = ggplot2::element_line(colour = "black", size = 1, lineend = "square"),
    axis.title.y.left = ggplot2::element_text(size = 12),
    axis.line.x.bottom = ggplot2::element_line(colour = "black", size = 1),
    axis.text.y = ggplot2::element_text(colour = "black", size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11)
  )

## Fig 3 S1C ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import metadata
OPP_Run3_Metadata <- read.csv(".\\Data\\OPP Run3 Metadata.csv", stringsAsFactors = F, check.names = T, sep = ",")[, -9]
OPP_Run4_Metadata <- read.csv(".\\Data\\OPP Run4 Metadata.csv", stringsAsFactors = F, check.names = T, sep = ",")[, -9]

#filter out empty treatment rows in metadata 
OPP_Run3_Metadata <- subset(OPP_Run3_Metadata, !is.na(Treatment))
OPP_Run4_Metadata <- subset(OPP_Run4_Metadata, !is.na(Treatment))

#create group column
OPP_Run3_Metadata$Group <- paste0(OPP_Run3_Metadata$Sex, "_", OPP_Run3_Metadata$Treatment)
OPP_Run4_Metadata$Group <- paste0(OPP_Run4_Metadata$Sex, "_", OPP_Run4_Metadata$Treatment)

#import files
OPP_Run3 <- read.csv(".\\Data\\01292021_OPP_B6_MI_Run3.1_preprocessed.csv", stringsAsFactors = F, check.names = T, sep = ",")
OPP_Run4 <- read.csv(".\\Data\\02122021_OPP_B6_MI_Run4_v3.1_red cell size_preprocessed.csv", stringsAsFactors = F, check.names = T, sep = ",")

#name well column
colnames(OPP_Run3_Metadata)[1] <- "Well_Name"
colnames(OPP_Run4_Metadata)[1] <- "Well_Name"

colnames(OPP_Run3)[2] <- "Well_Name"
colnames(OPP_Run4)[2] <- "Well_Name"

#merge data and metadata
OPP_Run3 <- merge(OPP_Run3, OPP_Run3_Metadata[, c("Well_Name", "Sex", "Group")], by = "Well_Name", all.x = TRUE)
OPP_Run4 <- merge(OPP_Run4, OPP_Run4_Metadata[, c("Well_Name", "Sex", "Group")], by = "Well_Name", all.x = TRUE)

#remove rows with sex missing
OPP_Run3 <- subset(OPP_Run3, Sex != "#N/A")

#get mean per run and sex of OPP-
OPP_Run3_mean_M <- mean(OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Group == "Male_OPP-"])
OPP_Run3_mean_F <- mean(OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Group == "Female_OPP-"])

OPP_Run4_mean_M <- mean(OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Group == "Male_OPP-"])
OPP_Run4_mean_F <- mean(OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Group == "Female_OPP-"])

#correct intensity for OPP-
OPP_Run3$Intensity_Corrected[OPP_Run3$Sex == "Male"] <- OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Sex == "Male"] - OPP_Run3_mean_M
OPP_Run3$Intensity_Corrected[OPP_Run3$Sex == "Female"] <- OPP_Run3$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run3$Sex == "Female"] - OPP_Run3_mean_F

OPP_Run4$Intensity_Corrected[OPP_Run4$Sex == "Male"] <- OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Sex == "Male"] - OPP_Run4_mean_M
OPP_Run4$Intensity_Corrected[OPP_Run4$Sex == "Female"] <- OPP_Run4$Cell..W2.Cell.Integr.Intensity..MultiWaveScoring.[OPP_Run4$Sex == "Female"] - OPP_Run4_mean_F

#combine runs 3 and 4
OPP_Run3_4 <- rbind(OPP_Run3, OPP_Run4)

#select only for rows with positive corrected intensity
OPP_Run3_4 <- subset(OPP_Run3_4, OPP_Run3_4$Intensity_Corrected > 0)

#select for FBS- rows
OPP_Run3_4_FBS_negative <- OPP_Run3_4[grep("-FBS", OPP_Run3_4$Group), ]

#select for FBS+ rows
OPP_Run3_4_FBS_positive <- OPP_Run3_4[grep("[:+:]FBS", OPP_Run3_4$Group), ]

#make rows unique
Group_List_Run3_4_FBS_negative <- unique(OPP_Run3_4_FBS_negative$Group)
Group_List_Run3_4_FBS_positive <- unique(OPP_Run3_4_FBS_positive$Group)

#set FBS- group factor levels
OPP_Run3_4_FBS_negative$Group <- factor(OPP_Run3_4_FBS_negative$Group, levels = c("Female_DMSO_-FBS",
                                                                                  "Female_2hrTg_-FBS",
                                                                                  "Female_24hrTg_-FBS",
                                                                                  "Male_DMSO_-FBS",
                                                                                  "Male_2hrTg_-FBS",
                                                                                  "Male_24hrTg_-FBS"))

#calculate stats for FBS-
OPP_Run3_4_FBS_negative_norm <- OPP_Run3_4_FBS_negative %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    sd = sd(Intensity_Corrected),
    mean = mean(Intensity_Corrected),
    n = length(Intensity_Corrected),
    interquartile_bar = IQR(Intensity_Corrected)/2,
    median = median(Intensity_Corrected)
  )

OPP_Run3_4_FBS_negative_norm$SEM <- OPP_Run3_4_FBS_negative_norm$sd/sqrt(OPP_Run3_4_FBS_negative_norm$n)

#set axis limits and label position
left_limit <- 0.5
right_limit <- length(Group_List_Run3_4_FBS_negative) + 0.5

#create plot
ggplot2::ggplot(OPP_Run3_4_FBS_negative_norm, ggplot2::aes(y = Intensity_Corrected, x = Group)) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0, 700000), expand = FALSE, clip = "off") +
  
  ggbeeswarm::geom_quasirandom(
    ggplot2::aes(
      y = Intensity_Corrected,
      colour = Sex
    ),
    size = 0.7,
    data = subset(OPP_Run3_4_FBS_negative, OPP_Run3_4_FBS_negative$Intensity_Corrected < 700000)
  ) +
  
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(x) mean(x),
                        fun.max = function(x) mean(x),
                        geom = "errorbar",
                        ggplot2::aes(group = Group),
                        size = 0.5,
                        colour = "black",
                        width = 0.4,
                        data = OPP_Run3_4_FBS_negative) +
  
  ggplot2::stat_summary(fun = mean,
                        fun.min = function(x) mean(x) - std.error(x)/2,
                        fun.max = function(x) mean(x) + std.error(x)/2,
                        geom = "errorbar",
                        ggplot2::aes(group = Group),
                        colour = "black",
                        width = 0.1,
                        data = OPP_Run3_4_FBS_negative) +
  
  ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE), 
                              limits = c(NA, 700000), 
                              breaks = seq(0, 700000, by = 100000)) +
  
  ggplot2::labs(y = "Integrated Fluorescence Intensity", x = "Group") +
  
  ggplot2::guides(
    colour = ggplot2::guide_legend(reverse = FALSE, override.aes = list(size = 5))
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(label = round(mean, 1), y = -40000, x = Group),
    colour = "black",
    hjust = 0.5,
    vjust = 0) +
  
  ggplot2::annotate(geom = "segment", x = 3.5, xend = 3.5, y = 0, yend = 700000, colour = "Black", size = 0.5, linetype = 2) +
  
  ggplot2::scale_x_discrete(labels = c("DMSO", "2 hr Tg", "24 hr Tg", "DMSO", "2 hr Tg", "24 hr Tg")) +
  
  ggplot2::scale_color_manual(values = c("#E7B81E", "#84B429")) +
  
  ggplot2::theme_bw() + 
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(t = 0.25, r = 0, b = 1, l = 0, "cm"),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x.bottom = ggplot2::element_blank(),
    axis.ticks.length.y.left = grid::unit(2, "mm"),
    axis.ticks.y.left = ggplot2::element_line(size = 1, colour = "black"),
    panel.border = ggplot2::element_blank(),
    axis.line.y.left = ggplot2::element_line(colour = "black", size = 1, lineend = "square"),
    axis.title.y.left = ggplot2::element_text(size = 12),
    axis.line.x.bottom = ggplot2::element_line(colour = "black", size = 1),
    axis.text.x = ggplot2::element_text(angle = 0, vjust = -8, hjust = 0.5, size = 11, colour = "black"),
    axis.text.y = ggplot2::element_text(colour = "black", size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11))


### Figure 3 - Supplement 2 ###

## Fig 3 S2A ##

#clear environment
rm(list = ls())

#import run 1 file (too large to upload to GitHub)
df_Ins2_GFP_Run_1 <-  read.csv(
  ".\\Data\\Fullsheet_raw.tsv",
  check.names = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

#add time column
df_Ins2_GFP_Run_1$`Time (Hours)` <- (df_Ins2_GFP_Run_1$Frame - 1)/2

#order data by well, run, cell, and frame
df_Ins2_GFP_Run_1 <-  df_Ins2_GFP_Run_1[with(df_Ins2_GFP_Run_1, order(Well, ObjectNumber, Frame)),]

#manually filtering out bad wells
bad_well_list <- c("G02", "H02", "H03", "H04", "H05", "J10", "J16", "M10", "M11", "M15")

df_Ins2_GFP_Run_1 <- df_Ins2_GFP_Run_1[!df_Ins2_GFP_Run_1$Well %in% bad_well_list,]

#manually smoothing out blips
list_of_wells <- c("B04", "G07", "G13")

list_of_time_point_lists <- list(list(18),
                                 list(3),
                                 list(2, 10, 11, 15, 40, 97, 121, 132, 137)
)

for (i in 1:3){
  
  well <- list_of_wells[i]
  
  time_point_list <- list_of_time_point_lists[[i]]
  
  for (time_point in time_point_list) {
    
    df_Ins2_GFP_Run_1$FITC[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    df_Ins2_GFP_Run_1$TEXASRED[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    df_Ins2_GFP_Run_1$DAPI[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    
  }
  
}

#smoothing via interpolation
df_Ins2_GFP_Run_1$FITC <- zoo::na.approx(df_Ins2_GFP_Run_1$FITC)
df_Ins2_GFP_Run_1$TEXASRED <- zoo::na.approx(df_Ins2_GFP_Run_1$TEXASRED)
df_Ins2_GFP_Run_1$DAPI <- zoo::na.approx(df_Ins2_GFP_Run_1$DAPI)


#import run 2 file (too large to upload to GitHub)
df_Ins2_GFP_Run_2 <-  read.csv(
  ".\\Data\\Fullsheet_raw2.tsv",
  check.names = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

#add time column
df_Ins2_GFP_Run_2$`Time (Hours)` <- (df_Ins2_GFP_Run_2$Frame - 1)/2

#order data by well, run, cell, and frame
df_Ins2_GFP_Run_2 <-  df_Ins2_GFP_Run_2[with(df_Ins2_GFP_Run_2, order(Well, ObjectNumber, Frame)),]


#manually filtering out bad wells
bad_well_list <- c("D10", "D12", "D15", "E02", "E04", "E11", "F09", "G10", "H02", "H03", "H04", "H05", "H10", "H12", "H13", "H14", "H16", "I09", "I10", "I16", "I18", "J02", "K13", "L04", "L10", "L13", "L15", "L21")

df_Ins2_GFP_Run_2 <- df_Ins2_GFP_Run_2[!df_Ins2_GFP_Run_2$Well %in% bad_well_list,]

#blip removal filter not needed as no imaging artifacts of intensity dropping were detected

#3 filter function
filter_data <- function(dataset){
  
  unique_CellID_list <- unique(dataset$CellID)
  
  for (Cell_ID in unique_CellID_list) {
    
    #frame must start at 1 or 2 minimum filter
    if (min(dataset$Frame[dataset$CellID == Cell_ID] > 2)) {
      
      dataset$FITC[dataset$CellID == Cell_ID] <- NA
      dataset$TEXASRED[dataset$CellID == Cell_ID] <- NA
      dataset$DAPI[dataset$CellID == Cell_ID] <- NA
      
    }
    
    #number of frames per cell must be greater than 3 filter
    if (length(dataset$Frame[dataset$CellID == Cell_ID]) < 3) {
      
      dataset$FITC[dataset$CellID == Cell_ID] <- NA
      dataset$TEXASRED[dataset$CellID == Cell_ID] <- NA
      dataset$DAPI[dataset$CellID == Cell_ID] <- NA
      
    }
    
    #all GFP deleted per cell after Texas Red hitting 2000 filter
    if_cell_dead <- FALSE
    
    for (Frame in sort(dataset$Frame[dataset$CellID == Cell_ID])) {
      
      if (if_cell_dead == FALSE) {
        
        if (is.na(dataset$TEXASRED[dataset$CellID == Cell_ID & 
                                   dataset$Frame == Frame])) {
          
        } else if (dataset$TEXASRED[dataset$CellID == Cell_ID & 
                                    dataset$Frame == Frame] > 2000) {
          
          if_cell_dead <- TRUE
          
        }
        
      }
      
      if (if_cell_dead == TRUE) {
        
        dataset$FITC[dataset$CellID == Cell_ID & 
                       dataset$Frame == Frame] <- NA
        
      }
      
    }
    
  }
  
  dataset
  
}

#filter datasets in parallel by multicore processing
future::plan(cluster, workers = future::availableCores() - 1)

data_split_Run_1 <- split(df_Ins2_GFP_Run_1, df_Ins2_GFP_Run_1$Well)

df_Ins2_GFP_Run_1 <- furrr::future_map_dfr(data_split_Run_1, function(.data){
  
  filter_data(.data)
  
})

data_split_Run_2 <- split(df_Ins2_GFP_Run_2, df_Ins2_GFP_Run_2$Well)

df_Ins2_GFP_Run_2 <- furrr::future_map_dfr(data_split_Run_2, function(.data){
  
  filter_data(.data)
  
})

plan(cluster, workers = 1)

#label each run
df_Ins2_GFP_Run_1$Run <- "CH1_R1"
df_Ins2_GFP_Run_2$Run <- "CH1_R2"

#combine runs
df_Ins2_GFP_Run_1_and_2 <- rbind(df_Ins2_GFP_Run_1, df_Ins2_GFP_Run_2)[, c(1:7, 10:12)]

#import metadata
df_metadata <-
  read.csv(
    ".\\Data\\Ins2GFP+ CH1_R1+R2 metadata.csv",
    stringsAsFactors = F,
    check.names = FALSE,
    sep = ","
  ) %>% na.omit()

#combine sex and treatment columns
df_metadata$Sex_w_Treatment <- paste(df_metadata$Sex, df_metadata$Treatment, sep = "_")

#format Sex_w_Treatment column
df_metadata$Sex_w_Treatment <- gsub("_", " ", gsub("_percent", "%", gsub("_positive", "+", gsub("_negative", "-", df_metadata$Sex_w_Treatment))))

#define treatments needed
treatments_of_interest <- c(
  "0.1_percent_DMSO_FBS_positive", 
  "0.1_uM_Tg_FBS_positive", 
  "1_uM_Tg_FBS_positive")

#filter metadata for treatments needed
df_metadata <- df_metadata[df_metadata$Treatment %in% treatments_of_interest,]

#format treatment column
df_metadata$Treatment <- gsub("_", " ", gsub("_percent", "%", gsub("_positive", "+", gsub("_negative", "-", df_metadata$Treatment))))

#merge data and metadata
df_merged <- merge(df_Ins2_GFP_Run_1_and_2, df_metadata, by = c("Well", "Run"), all.y = TRUE) %>% na.omit()

#order data by well, run, cell, and frame
df_merged <- df_merged[with(df_merged, order(Well, Run, ObjectNumber, Frame)),]

#calculate average intensity by sex and treatment for first 2 hours
df_2_hr_avg <- df_merged[df_merged$`Time (Hours)` <= 2,] %>%
  dplyr::group_by(Sex_w_Treatment) %>%
  dplyr::summarise(
    mean_FITC_2_hr = mean(FITC, na.rm = TRUE),
  )

#normalize intensity to average of first two hours by sex and treatment
for (i in df_2_hr_avg$Sex_w_Treatment) {
  
  df_merged$FITC_norm[df_merged$Sex_w_Treatment == i] <- df_merged$FITC[df_merged$Sex_w_Treatment == i]/df_2_hr_avg$mean_FITC_2_hr[df_2_hr_avg$Sex_w_Treatment == i]
  
}

#calculate stats for intensity for each time point, sex, treatment, and run
df_merged_mean <- df_merged %>%
  dplyr::group_by(`Time (Hours)`, Frame, Treatment, Sex, Sex_w_Treatment, Run) %>%
  dplyr::summarise(
    sd_norm = sd(FITC_norm, na.rm = TRUE),
    n = length(FITC),
    mean_FITC_norm = mean(FITC_norm, na.rm = TRUE),
    mean_TEXASRED = mean(TEXASRED, na.rm = TRUE),
    mean_DAPI = mean(DAPI, na.rm = TRUE)
  )

#calculate sum x and sum x^2
df_merged_mean$ex <- df_merged_mean$n * df_merged_mean$mean_FITC_norm
df_merged_mean$exx <- df_merged_mean$sd_norm^2 * (df_merged_mean$n - 1) + df_merged_mean$ex^2/df_merged_mean$n

#calculate total n, sum x, and sum x^2
df.summary2 <- df_merged_mean %>%
  dplyr::group_by(`Time (Hours)`, Frame, Treatment, Sex, Sex_w_Treatment) %>%
  dplyr::summarise(
    tn = sum(n),
    tx = sum(ex),
    txx = sum(exx)
  )

#calculate combined mean and sd from each run together
df.summary2$mean_FITC_norm <- df.summary2$tx/df.summary2$tn
df.summary2$sd_norm_combined <- sqrt((df.summary2$txx - df.summary2$tx^2/df.summary2$tn)/(df.summary2$tn - 1))

#calculate sem
df.summary2$sem_norm <- df.summary2$sd_norm_combined/sqrt(df.summary2$tn)

#calculate mean at 0 hr for each sex and treatment
df_0_hr <- df.summary2[df.summary2$`Time (Hours)` == 0,] %>%
  dplyr::group_by(Sex_w_Treatment) %>%
  dplyr::summarise(
    mean_FITC_norm_0_hr = mean_FITC_norm,
  )

#initialize new minus 0 hr columns
df.summary2$mean_FITC_norm_minus_0_hr <- NA
df.summary2$mean_FITC_norm_minus_0_hr_plus_sem <- NA
df.summary2$mean_FITC_norm_minus_0_hr_minus_sem <- NA

#subtract by mean at 0 hr for each sex and treatment
for (i in df_0_hr$Sex_w_Treatment) {
  
  df.summary2$mean_FITC_norm_minus_0_hr[df.summary2$Sex_w_Treatment == i] <- df.summary2$mean_FITC_norm[df.summary2$Sex_w_Treatment == i] - df_0_hr$mean_FITC_norm_0_hr[df_0_hr$Sex_w_Treatment == i]
  
}

#calculate new plus and minus sem values
df.summary2$mean_FITC_norm_minus_0_hr_plus_sem <- df.summary2$mean_FITC_norm_minus_0_hr + df.summary2$sem_norm
df.summary2$mean_FITC_norm_minus_0_hr_minus_sem <- df.summary2$mean_FITC_norm_minus_0_hr - df.summary2$sem_norm

#format treatment names
df.summary2$Treatment <- gsub("0.1% DMSO FBS+", "DMSO", df.summary2$Treatment, fixed = TRUE)
df.summary2$Treatment <- gsub("FBS+", "", df.summary2$Treatment, fixed = TRUE)
df.summary2$Treatment <- gsub("u", "μ", df.summary2$Treatment, fixed = TRUE)

#filter for sex and range from 0 to 60 hr
df.summary3 <- df.summary2[df.summary2$Sex == "Female" & df.summary2$`Time (Hours)` <= 60,]

#order treatment factor levels
df.summary3$Treatment <- factor(df.summary3$Treatment, levels = sort(unique(df.summary3$Treatment))[c(3, 1, 2)])

#create plot
ggplot2::ggplot(df.summary3, aes(x = `Time (Hours)`)) +
  
  ggplot2::geom_line(ggplot2::aes(y = mean_FITC_norm_minus_0_hr, 
                                  group = Treatment, 
                                  colour = Treatment),
                     data = df.summary3,
                     size = 1) +
  
  ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_FITC_norm_minus_0_hr_minus_sem, 
                                    ymax = mean_FITC_norm_minus_0_hr_plus_sem, 
                                    group = Treatment, 
                                    fill = Treatment, 
                                    colour = Treatment),
                       alpha = 0.2,
                       linetype = 0,
                       outline.type = "both",
                       show.legend = FALSE) +
  
  ggplot2::scale_color_manual(values = c("black", "#E7821E", "#E7B91E")) +
  
  ggplot2::scale_fill_manual(values = c("black", "#E7821E", "#E7B91E")) +
  
  ggplot2::coord_cartesian(xlim = c(0, 60),
                           ylim = c(-0.2, 0.2),
                           expand = FALSE,
                           clip = "on") +
  
  ggplot2::scale_x_continuous(breaks = c(seq(0, 96, 12))) +
  
  ggplot2::scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  
  ggplot2::labs(y = "Normalized GFP Intensity", x = "Time (hours)") +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0.25,
      r = 0.25,
      b = 0,
      l = 0.25,
      "cm"
    ),
    axis.ticks.length.y.left = grid::unit(2, "mm"),
    axis.ticks.y.left = ggplot2::element_line(size = 1, colour = "black"),
    axis.ticks.length.x.bottom = grid::unit(2, "mm"),
    axis.ticks.x.bottom = ggplot2::element_line(size = 1, colour = "black"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", size = 1, lineend = "square"),
    axis.title.y.left = ggplot2::element_text(size = 12),
    axis.line.x.bottom = ggplot2::element_line(colour = "black", size = 1),
    axis.text.x = ggplot2::element_text(size = 11, colour = "black"),
    axis.text.y = ggplot2::element_text(colour = "black", size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = grid::unit(1, "cm"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11),
    aspect.ratio = 1 / 2)


## Fig 3 S2B ##

#clear environment
rm(list = ls())

#import run 1 file
df_Ins2_GFP_Run_1 <-  read.csv(
  ".\\Data\\Fullsheet_raw.tsv",
  check.names = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

#add time column
df_Ins2_GFP_Run_1$`Time (Hours)` <- (df_Ins2_GFP_Run_1$Frame - 1)/2

#order data by well, run, cell, and frame
df_Ins2_GFP_Run_1 <-  df_Ins2_GFP_Run_1[with(df_Ins2_GFP_Run_1, order(Well, ObjectNumber, Frame)),]

#manually filtering out bad wells
bad_well_list <- c("G02", "H02", "H03", "H04", "H05", "J10", "J16", "M10", "M11", "M15")

df_Ins2_GFP_Run_1 <- df_Ins2_GFP_Run_1[!df_Ins2_GFP_Run_1$Well %in% bad_well_list,]

#manually smoothing out blips
list_of_wells <- c("B04", "G07", "G13")

list_of_time_point_lists <- list(list(18),
                                 list(3),
                                 list(2, 10, 11, 15, 40, 97, 121, 132, 137)
)

for (i in 1:3){
  
  well <- list_of_wells[i]
  
  time_point_list <- list_of_time_point_lists[[i]]
  
  for (time_point in time_point_list) {
    
    df_Ins2_GFP_Run_1$FITC[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    df_Ins2_GFP_Run_1$TEXASRED[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    df_Ins2_GFP_Run_1$DAPI[df_Ins2_GFP_Run_1$Frame == time_point & df_Ins2_GFP_Run_1$Well == well] <- NA
    
  }
  
}

#smoothing via interpolation
df_Ins2_GFP_Run_1$FITC <- zoo::na.approx(df_Ins2_GFP_Run_1$FITC)
df_Ins2_GFP_Run_1$TEXASRED <- zoo::na.approx(df_Ins2_GFP_Run_1$TEXASRED)
df_Ins2_GFP_Run_1$DAPI <- zoo::na.approx(df_Ins2_GFP_Run_1$DAPI)

#import run 2 file
df_Ins2_GFP_Run_2 <-  read.csv(
  ".\\Data\\Fullsheet_raw2.tsv",
  check.names = TRUE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

#add time column
df_Ins2_GFP_Run_2$`Time (Hours)` <- (df_Ins2_GFP_Run_2$Frame - 1)/2

#order data by well, run, cell, and frame
df_Ins2_GFP_Run_2 <-  df_Ins2_GFP_Run_2[with(df_Ins2_GFP_Run_2, order(Well, ObjectNumber, Frame)),]

#manually filtering out bad wells
bad_well_list <- c("D10", "D12", "D15", "E02", "E04", "E11", "F09", "G10", "H02", "H03", "H04", "H05", "H10", "H12", "H13", "H14", "H16", "I09", "I10", "I16", "I18", "J02", "K13", "L04", "L10", "L13", "L15", "L21")

df_Ins2_GFP_Run_2 <- df_Ins2_GFP_Run_2[!df_Ins2_GFP_Run_2$Well %in% bad_well_list,]

#blip removal filter not needed as no imaging artifacts of intensity dropping were detected

#3 filter function
filter_data <- function(dataset){
  
  unique_CellID_list <- unique(dataset$CellID)
  
  
  for (Cell_ID in unique_CellID_list) {
    
    #frame must start at 1 or 2 minimum filter
    if (min(dataset$Frame[dataset$CellID == Cell_ID] > 2)) {
      
      dataset$FITC[dataset$CellID == Cell_ID] <- NA
      dataset$TEXASRED[dataset$CellID == Cell_ID] <- NA
      dataset$DAPI[dataset$CellID == Cell_ID] <- NA
      
    }
    
    #number of frames per cell must be greater than 3 filter
    if (length(dataset$Frame[dataset$CellID == Cell_ID]) < 3) {
      
      dataset$FITC[dataset$CellID == Cell_ID] <- NA
      dataset$TEXASRED[dataset$CellID == Cell_ID] <- NA
      dataset$DAPI[dataset$CellID == Cell_ID] <- NA
      
    }
    
    #all GFP deleted per cell after Texas Red hitting 2000 filter
    if_cell_dead <- FALSE
    
    for (Frame in sort(dataset$Frame[dataset$CellID == Cell_ID])) {
      
      if (if_cell_dead == FALSE) {
        
        if (is.na(dataset$TEXASRED[dataset$CellID == Cell_ID & 
                                   dataset$Frame == Frame])) {
          
        } else if (dataset$TEXASRED[dataset$CellID == Cell_ID & 
                                    dataset$Frame == Frame] > 2000) {
          
          if_cell_dead <- TRUE
          
        }
        
      }
      
      if (if_cell_dead == TRUE) {
        
        dataset$FITC[dataset$CellID == Cell_ID & 
                       dataset$Frame == Frame] <- NA
        
      }
      
    }
    
  }
  
  dataset
  
}

#filter datasets in parallel by multicore processing
future::plan(cluster, workers = future::availableCores() - 1)

data_split_Run_1 <- split(df_Ins2_GFP_Run_1, df_Ins2_GFP_Run_1$Well)

df_Ins2_GFP_Run_1 <- furrr::future_map_dfr(data_split_Run_1, function(.data){
  
  filter_data(.data)
  
})

data_split_Run_2 <- split(df_Ins2_GFP_Run_2, df_Ins2_GFP_Run_2$Well)

df_Ins2_GFP_Run_2 <- furrr::future_map_dfr(data_split_Run_2, function(.data){
  
  filter_data(.data)
  
})

plan(cluster, workers = 1)

#label each run
df_Ins2_GFP_Run_1$Run <- "CH1_R1"
df_Ins2_GFP_Run_2$Run <- "CH1_R2"

#combine runs
df_Ins2_GFP_Run_1_and_2 <- rbind(df_Ins2_GFP_Run_1, df_Ins2_GFP_Run_2)[, c(1:7, 10:12)]

#import metadata
df_metadata <-
  read.csv(
    ".\\Data\\Ins2GFP+ CH1_R1+R2 metadata.csv",
    stringsAsFactors = F,
    check.names = FALSE,
    sep = ","
  ) %>% na.omit()

#combine sex and treatment columns
df_metadata$Sex_w_Treatment <- paste(df_metadata$Sex, df_metadata$Treatment, sep = "_")

#format Sex_w_Treatment column
df_metadata$Sex_w_Treatment <- gsub("_", " ", gsub("_percent", "%", gsub("_positive", "+", gsub("_negative", "-", df_metadata$Sex_w_Treatment))))

#define treatments needed
treatments_of_interest <- c(
  "0.1_percent_DMSO_FBS_positive", 
  "0.1_uM_Tg_FBS_positive", 
  "1_uM_Tg_FBS_positive")

#filter metadata for treatments needed
df_metadata <- df_metadata[df_metadata$Treatment %in% treatments_of_interest,]

#format treatment column
df_metadata$Treatment <- gsub("_", " ", gsub("_percent", "%", gsub("_positive", "+", gsub("_negative", "-", df_metadata$Treatment))))

#merge data and metadata
df_merged <- merge(df_Ins2_GFP_Run_1_and_2, df_metadata, by = c("Well", "Run"), all.y = TRUE) %>% na.omit()

#order data by well, run, cell, and frame
df_merged <- df_merged[with(df_merged, order(Well, Run, ObjectNumber, Frame)),]

#calculate average intensity by sex and treatment for first 2 hours
df_2_hr_avg <- df_merged[df_merged$`Time (Hours)` <= 2,] %>%
  dplyr::group_by(Sex_w_Treatment) %>%
  dplyr::summarise(
    mean_FITC_2_hr = mean(FITC, na.rm = TRUE),
  )

#normalize intensity to average of first two hours by sex and treatment
for (i in df_2_hr_avg$Sex_w_Treatment) {
  
  df_merged$FITC_norm[df_merged$Sex_w_Treatment == i] <- df_merged$FITC[df_merged$Sex_w_Treatment == i]/df_2_hr_avg$mean_FITC_2_hr[df_2_hr_avg$Sex_w_Treatment == i]
  
}

#calculate stats for intensity for each time point, sex, treatment, and run
df_merged_mean <- df_merged %>%
  dplyr::group_by(`Time (Hours)`, Frame, Treatment, Sex, Sex_w_Treatment, Run) %>%
  dplyr::summarise(
    sd_norm = sd(FITC_norm, na.rm = TRUE),
    n = length(FITC),
    mean_FITC_norm = mean(FITC_norm, na.rm = TRUE),
    mean_TEXASRED = mean(TEXASRED, na.rm = TRUE),
    mean_DAPI = mean(DAPI, na.rm = TRUE)
  )

#calculate sum x and sum x^2
df_merged_mean$ex <- df_merged_mean$n * df_merged_mean$mean_FITC_norm
df_merged_mean$exx <- df_merged_mean$sd_norm^2 * (df_merged_mean$n - 1) + df_merged_mean$ex^2/df_merged_mean$n

#calculate total n, sum x, and sum x^2
df.summary2 <- df_merged_mean %>%
  dplyr::group_by(`Time (Hours)`, Frame, Treatment, Sex, Sex_w_Treatment) %>%
  dplyr::summarise(
    tn = sum(n),
    tx = sum(ex),
    txx = sum(exx)
  )

#calculate combined mean and sd from each run together
df.summary2$mean_FITC_norm <- df.summary2$tx/df.summary2$tn
df.summary2$sd_norm_combined <- sqrt((df.summary2$txx - df.summary2$tx^2/df.summary2$tn)/(df.summary2$tn - 1))

#calculate sem
df.summary2$sem_norm <- df.summary2$sd_norm_combined/sqrt(df.summary2$tn)

#calculate mean at 0 hr for each sex and treatment
df_0_hr <- df.summary2[df.summary2$`Time (Hours)` == 0,] %>%
  dplyr::group_by(Sex_w_Treatment) %>%
  dplyr::summarise(
    mean_FITC_norm_0_hr = mean_FITC_norm,
  )

#initialize new minus 0 hr columns
df.summary2$mean_FITC_norm_minus_0_hr <- NA
df.summary2$mean_FITC_norm_minus_0_hr_plus_sem <- NA
df.summary2$mean_FITC_norm_minus_0_hr_minus_sem <- NA

#subtract by mean at 0 hr for each sex and treatment
for (i in df_0_hr$Sex_w_Treatment) {
  
  df.summary2$mean_FITC_norm_minus_0_hr[df.summary2$Sex_w_Treatment == i] <- df.summary2$mean_FITC_norm[df.summary2$Sex_w_Treatment == i] - df_0_hr$mean_FITC_norm_0_hr[df_0_hr$Sex_w_Treatment == i]
  
}

#calculate new plus and minus sem values
df.summary2$mean_FITC_norm_minus_0_hr_plus_sem <- df.summary2$mean_FITC_norm_minus_0_hr + df.summary2$sem_norm
df.summary2$mean_FITC_norm_minus_0_hr_minus_sem <- df.summary2$mean_FITC_norm_minus_0_hr - df.summary2$sem_norm

#format treatment names
df.summary2$Treatment <- gsub("0.1% DMSO FBS+", "DMSO", df.summary2$Treatment, fixed = TRUE)
df.summary2$Treatment <- gsub("FBS+", "", df.summary2$Treatment, fixed = TRUE)
df.summary2$Treatment <- gsub("u", "μ", df.summary2$Treatment, fixed = TRUE)

#filter for sex and range from 0 to 60 hr
df.summary3 <- df.summary2[df.summary2$Sex == "Male" & df.summary2$`Time (Hours)` <= 60,]

#order treatment factor levels
df.summary3$Treatment <- factor(df.summary3$Treatment, levels = sort(unique(df.summary3$Treatment))[c(3, 1, 2)])

#create plot
ggplot2::ggplot(df.summary3, aes(x = `Time (Hours)`)) +
  
  ggplot2::geom_line(ggplot2::aes(y = mean_FITC_norm_minus_0_hr, 
                                  group = Treatment, 
                                  colour = Treatment),
                     data = df.summary3,
                     size = 1) +
  
  ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_FITC_norm_minus_0_hr_minus_sem, 
                                    ymax = mean_FITC_norm_minus_0_hr_plus_sem, 
                                    group = Treatment, 
                                    fill = Treatment, 
                                    colour = Treatment),
                       alpha = 0.2,
                       linetype = 0,
                       outline.type = "both",
                       show.legend = FALSE) +
  
  ggplot2::scale_color_manual(values = c("black", "#3E5514", "#84B42A")) +
  
  ggplot2::scale_fill_manual(values = c("black", "#3E5514", "#84B42A")) +
  
  ggplot2::coord_cartesian(xlim = c(0, 60),
                           ylim = c(-0.2, 0.2),
                           expand = FALSE,
                           clip = "on") +
  
  ggplot2::scale_x_continuous(breaks = c(seq(0, 96, 12))) +
  
  ggplot2::scale_y_continuous(breaks = seq(-0.2, 0.2, 0.1)) +
  
  ggplot2::labs(y = "Normalized GFP Intensity", x = "Time (hours)") +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0.25,
      r = 0.25,
      b = 0,
      l = 0.25,
      "cm"
    ),
    axis.ticks.length.y.left = grid::unit(2, "mm"),
    axis.ticks.y.left = ggplot2::element_line(size = 1, colour = "black"),
    axis.ticks.length.x.bottom = grid::unit(2, "mm"),
    axis.ticks.x.bottom = ggplot2::element_line(size = 1, colour = "black"),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", size = 1, lineend = "square"),
    axis.title.y.left = ggplot2::element_text(size = 12),
    axis.line.x.bottom = ggplot2::element_line(colour = "black", size = 1),
    axis.text.x = ggplot2::element_text(size = 11, colour = "black"),
    axis.text.y = ggplot2::element_text(colour = "black", size = 11),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = grid::unit(1, "cm"),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 11),
    aspect.ratio = 1 / 2)


### Figure 5 ###

## Fig 5A ##

#clear environment
rm(list = ls())

#import file
df_rlog <- read.csv(".\\Data\\MouseIslet.TgTreatment_all_genes.rlog.csv")

#number replicated genes
df_rlog[df_rlog$X == "44256", 1] <- c("44256 (1)", "44256 (2)")
df_rlog[df_rlog$X == "44257", 1] <- c("44257 (1)", "44257 (2)")

#filter out rows with only 0
df_rlog_no_empty <- filter(df_rlog, df_rlog[, 2] != 0)

#remove gene column
df_rlog_no_empty_no_names <- df_rlog_no_empty[,-1]

#create data frame of info for each sample
colData <- data.frame(Sample = colnames(df_rlog_no_empty_no_names), 
                      Sex = as.factor(ifelse(grepl("F...", colnames(df_rlog_no_empty_no_names)), "Female", "Male")),
                      Time = as.factor(ifelse(grepl("6hr", colnames(df_rlog_no_empty_no_names)), "6 hr", "12 hr")),
                      Treatment = as.factor(ifelse(grepl("Tg", colnames(df_rlog_no_empty_no_names)), "Tg", "DMSO")))

#add condition column
colData$Condition <- paste(colData$Sex, colData$Time, colData$Treatment)

#set row names as gene names
row.names(df_rlog_no_empty_no_names) <- df_rlog_no_empty$X

#transform data
df_rlog_no_empty_no_names_transformed <- as.data.frame(t(df_rlog_no_empty_no_names))

#add condition column
df_rlog_no_empty_no_names_transformed$Condition <- colData$Condition

#create PCA data
pc <- pcaMethods::pca(df_rlog_no_empty_no_names_transformed, nPcs = 2)

#add annotations to PCA data
pc_plot <- merge(pcaMethods::scores(pc), df_rlog_no_empty_no_names_transformed, by = 0)

#remove unneeded columns
pc_plot <- pc_plot[, -c(4:18925)]

#flip PC1
pc_plot$PC1 <- -pc_plot$PC1

#order PCA data by condition
pc_plot <- dplyr::arrange(pc_plot, Condition)

#set condition factor level order
pc_plot$Condition <- factor(pc_plot$Condition, levels = sort(unique(pc_plot$Condition))[c(3, 4, 1, 2, 7, 8, 5, 6)])

#plot PCA data
ggplot2::ggplot(pc_plot, ggplot2::aes(PC1, PC2, colour = Condition)) +
  
  ggplot2::geom_point(
    ggplot2::aes(group = Condition, 
                 shape = Condition,
                 fill = Condition),
    data = pc_plot, size = c(rep(4.5, 3),
                             rep(4.5, 3),
                             rep(4.5, 4),
                             rep(4.5, 4),
                             rep(4.5, 3),
                             rep(4.5, 3),
                             rep(4.5, 4),
                             rep(4.5, 3))) +
  
  ggplot2::annotate(geom = "text", 
                    x = c(pc_plot$PC1[1] + 3.5, 
                          pc_plot$PC1[4] - 1.8,
                          pc_plot$PC1[7] + 2,
                          pc_plot$PC1[11] + 2.2,
                          pc_plot$PC1[15] + 5,
                          pc_plot$PC1[18] + 4.5,
                          pc_plot$PC1[21] + 3,
                          pc_plot$PC1[25] - 3), 
                    y = c(pc_plot$PC2[1] - 3,
                          pc_plot$PC2[4],
                          pc_plot$PC2[7] + 0.2,
                          pc_plot$PC2[11] - 0.8,
                          pc_plot$PC2[15] + 1.4,
                          pc_plot$PC2[18] + 1.4,
                          pc_plot$PC2[21],
                          pc_plot$PC2[25] - 2.5),
                    label = unique(pc_plot$Condition)[1:8], 
                    colour = c("#606060", 
                               "#E7B91E", 
                               "#000000", 
                               "#E7821E", 
                               "#606060", 
                               "#84B42A", 
                               "#000000", 
                               "#3E5514"),
                    size = 6, 
                    hjust = c(1, 1, 0, 0, 1, 1, 0, 1)) +
  
  ggforce::geom_mark_ellipse(x = pc_plot$PC1, 
                             y = pc_plot$PC2, 
                             fill =  c(rep("#606060", 3),
                                       rep("#E7B91E", 3),
                                       rep("#000000", 4),
                                       rep("#E7821E", 4),
                                       rep("#606060", 3),
                                       rep("#84B42A", 3),
                                       rep("#000000", 4),
                                       rep("#3E5514", 3)),
                             group = pc_plot$Condition, 
                             color = NA,
                             data = pc_plot) +
  
  ggplot2::scale_shape_manual(values = c(21, 22, 25, 24, 21, 22, 25, 24)) +
  
  ggplot2::scale_color_manual("Condition", values = c(rep("Black", 8))) +
  
  ggplot2::scale_fill_manual("Condition", 
                             values = c("#000000", 
                                        "#E7821E", 
                                        "#606060", 
                                        "#E7B91E", 
                                        "#000000", 
                                        "#3E5514", 
                                        "#606060", 
                                        "#84B42A")) +
  
  ggplot2::xlab(paste0("PC1: ", format(round(pc@R2[1] * 100, 0), nsmall = 0), "% of the variance")) +
  
  ggplot2::ylab(paste0("PC2: ", format(round(pc@R2[2] * 100, 0), nsmall = 0), "% of the variance")) +
  
  ggplot2::guides(color = "none") +
  
  ggplot2::theme_bw() +
  
  ggplot2::coord_cartesian(
    xlim = c(-31, 31),
    ylim = c(-16, 16),
    expand = FALSE,
    clip = "off") +
  
  ggplot2::scale_y_continuous(breaks = seq(-15, 15, by = 5)) +
  
  ggplot2::geom_hline(yintercept = 0, color = "black") +
  
  ggplot2::geom_vline(xintercept = 0, color = "black") +
  
  ggplot2::guides(colour = guide_legend(override.aes = list(size = 5, 
                                                            linetype = 0,
                                                            fill = c("#000000", 
                                                                     "#E7821E", 
                                                                     "#606060", 
                                                                     "#E7B91E", 
                                                                     "#000000", 
                                                                     "#3E5514", 
                                                                     "#606060", 
                                                                     "#84B42A"), 
                                                            colour = "black"))) + 
  
  ggplot2::theme(
    legend.position = "none",
    legend.title = ggplot2::element_text(size = 30),
    legend.spacing.y = grid::unit(5, "mm"),
    legend.key.height = grid::unit(1.25, "cm"),
    legend.text = ggplot2::element_text(size = 28),
    panel.grid = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 30),
    axis.text = ggplot2::element_text(size = 28, colour = "black"),
    axis.ticks = ggplot2::element_blank(),
    text = ggplot2::element_text(size = 28),
    aspect.ratio = 1/1)


## Fig 5B ##

#clear environment
rm(list = ls())

#import file
df_rlog <- read.csv(".\\Data\\MouseIslet.TgTreatment_all_genes.rlog.csv")

#number replicated genes
df_rlog[df_rlog$X == "44256", 1] <- c("44256 (1)", "44256 (2)")
df_rlog[df_rlog$X == "44257", 1] <- c("44257 (1)", "44257 (2)")

#filter out rows with only 0
df_rlog_no_empty <- filter(df_rlog, df_rlog[, 2] != 0)

#remove gene column
df_rlog_no_empty_no_names <- df_rlog_no_empty[,-1]

#transform data
df_rlog_no_empty_no_names_transformed <- as.data.frame(t(df_rlog_no_empty_no_names))

#create PCA data
pc <- pcaMethods::pca(df_rlog_no_empty_no_names_transformed, nPcs = 5)

#get PCA scores
pc_scores <- as.data.frame(pc@scores)

#add sex, time, and treatment numbers
pc_scores$Sex <- ifelse(substr(rownames(pc_scores), 1, 1) == "F", 1, 0)
pc_scores$Time <- ifelse(grepl("6hr", rownames(pc_scores)), 1, 0)
pc_scores$Treatment <- ifelse(grepl("DMSO", rownames(pc_scores)), 1, 0)

#perform Spearman's correlation on PCA scores
cormatrix <- cor(pc_scores, method = "spearman")

#get correlations needed
cormatrix <- cormatrix[c(1:5), c("Sex", "Time", "Treatment")]

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(cormatrix,
                                              name = "cormatrix",
                                              column_title_gp = grid::gpar(fontsize = 12),
                                              row_names_side = "left", 
                                              column_names_gp = grid::gpar(fontsize = 12),
                                              column_names_side = "top",
                                              col = circlize::colorRamp2(c(-ceiling(max(na.omit(reshape2::melt(cormatrix)$value))), 0, ceiling(max(na.omit(reshape2::melt(cormatrix)$value)))),
                                                                         c("blue", "white", "red")),
                                              row_dend_side = "right",
                                              right_annotation = ComplexHeatmap::rowAnnotation(empty2 = ComplexHeatmap::anno_empty(border = FALSE, width = grid::unit(10, "mm"))),
                                              cell_fun = function(j, i, x, y, width, height, fill) {
                                                grid::grid.text(sprintf("%.1f", cormatrix[i, j]), x, y, gp = grid::gpar(fontsize = 10))
                                              },
                                              column_order = 1:ncol(cormatrix),
                                              row_order = 1:nrow(cormatrix),
                                              height = nrow(cormatrix)*grid::unit(10, "mm"),
                                              width = ncol(cormatrix)*grid::unit(10, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(title = NULL,
                                                                          title_position = "topcenter",
                                                                          direction = "horizontal",
                                                                          legend_height = grid::unit(3, "cm")
                                              )
)

#draw heatmap
ComplexHeatmap::draw(heatmap_list_total, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend label
grid::grid.text("Spearman's Correlation",
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 11))


## Fig 5C ##

#clear environment
rm(list = ls())

#import files
df_F_DMSO_vs_Tg_Tg_6hr <- read.csv(".\\Data\\Female DMSO vs Tg Time by treatment interaction_higher at 6hr_cat 367 0.01 Reactome Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_F_DMSO_vs_Tg_Tg_12hr <- read.csv(".\\Data\\Female DMSO vs Tg Time by treatment interaction_higher at 12hrs_cat 458 0.01 Reactome Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_DMSO_vs_Tg_Tg_6hr <- read.csv(".\\Data\\Male DMSO vs Tg Time by treatment interaction_higher at 6hrs_cat 367 0.01 Reactome Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_DMSO_vs_Tg_Tg_12hr <- read.csv(".\\Data\\Male DMSO vs Tg Time by treatment interaction_higher at 12hrs_cat 458 0.01 Reactome Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]

#annotate time
df_F_DMSO_vs_Tg_Tg_6hr$Time <- "6 hr"
df_F_DMSO_vs_Tg_Tg_12hr$Time <- "12 hr"
df_M_DMSO_vs_Tg_Tg_6hr$Time <- "6 hr"
df_M_DMSO_vs_Tg_Tg_12hr$Time <- "12 hr"

#fix column names
colnames(df_F_DMSO_vs_Tg_Tg_6hr) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Time")
colnames(df_F_DMSO_vs_Tg_Tg_12hr) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Time")
colnames(df_M_DMSO_vs_Tg_Tg_6hr) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Time")
colnames(df_M_DMSO_vs_Tg_Tg_12hr) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Time")

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_F_DMSO_vs_Tg_Tg_6hr$Gene.ratio <- df_F_DMSO_vs_Tg_Tg_6hr$Entities.found/369
df_F_DMSO_vs_Tg_Tg_12hr$Gene.ratio <- df_F_DMSO_vs_Tg_Tg_12hr$Entities.found/255
df_M_DMSO_vs_Tg_Tg_6hr$Gene.ratio <- df_M_DMSO_vs_Tg_Tg_6hr$Entities.found/509
df_M_DMSO_vs_Tg_Tg_12hr$Gene.ratio <- df_M_DMSO_vs_Tg_Tg_12hr$Entities.found/530

#get -log10 of FDR
df_F_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR <- -log10(df_F_DMSO_vs_Tg_Tg_6hr$FDR)
df_F_DMSO_vs_Tg_Tg_12hr$transf.log.10.FDR <- -log10(df_F_DMSO_vs_Tg_Tg_12hr$FDR)
df_M_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR <- -log10(df_M_DMSO_vs_Tg_Tg_6hr$FDR)
df_M_DMSO_vs_Tg_Tg_12hr$transf.log.10.FDR <- -log10(df_M_DMSO_vs_Tg_Tg_12hr$FDR)

#order by gene ratio
df_F_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_F_DMSO_vs_Tg_Tg_6hr, desc(Gene.ratio))
df_F_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_F_DMSO_vs_Tg_Tg_12hr, desc(Gene.ratio))
df_M_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_M_DMSO_vs_Tg_Tg_6hr, desc(Gene.ratio))
df_M_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_M_DMSO_vs_Tg_Tg_12hr, desc(Gene.ratio))

#order by -log10(FDR)
df_F_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_F_DMSO_vs_Tg_Tg_6hr, desc(transf.log.10.FDR))
df_F_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_F_DMSO_vs_Tg_Tg_12hr, desc(transf.log.10.FDR))
df_M_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_M_DMSO_vs_Tg_Tg_6hr, desc(transf.log.10.FDR))
df_M_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_M_DMSO_vs_Tg_Tg_12hr, desc(transf.log.10.FDR))

#categorize significance
df_F_DMSO_vs_Tg_Tg_6hr$significance <- ifelse(df_F_DMSO_vs_Tg_Tg_6hr$FDR < 0.01, "< 0.01", ifelse(df_F_DMSO_vs_Tg_Tg_6hr$FDR > 0.01 & df_F_DMSO_vs_Tg_Tg_6hr$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_F_DMSO_vs_Tg_Tg_12hr$significance <- ifelse(df_F_DMSO_vs_Tg_Tg_12hr$FDR < 0.01, "< 0.01", ifelse(df_F_DMSO_vs_Tg_Tg_12hr$FDR > 0.01 & df_F_DMSO_vs_Tg_Tg_12hr$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_DMSO_vs_Tg_Tg_6hr$significance <- ifelse(df_M_DMSO_vs_Tg_Tg_6hr$FDR < 0.01, "< 0.01", ifelse(df_M_DMSO_vs_Tg_Tg_6hr$FDR > 0.01 & df_M_DMSO_vs_Tg_Tg_6hr$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_DMSO_vs_Tg_Tg_12hr$significance <- ifelse(df_M_DMSO_vs_Tg_Tg_12hr$FDR < 0.01, "< 0.01", ifelse(df_M_DMSO_vs_Tg_Tg_12hr$FDR > 0.01 & df_M_DMSO_vs_Tg_Tg_12hr$FDR < 0.05, "0.01-0.05", "> 0.05"))

#annotate sex
colnames(df_F_DMSO_vs_Tg_Tg_6hr)[c(2:4, 6:8)] <- paste0(colnames(df_F_DMSO_vs_Tg_Tg_6hr)[c(2:4, 6:8)], "_F")
colnames(df_F_DMSO_vs_Tg_Tg_12hr)[c(2:4, 6:8)] <- paste0(colnames(df_F_DMSO_vs_Tg_Tg_12hr)[c(2:4, 6:8)], "_F")
colnames(df_M_DMSO_vs_Tg_Tg_6hr)[c(2:4, 6:8)] <- paste0(colnames(df_M_DMSO_vs_Tg_Tg_6hr)[c(2:4, 6:8)], "_M")
colnames(df_M_DMSO_vs_Tg_Tg_12hr)[c(2:4, 6:8)] <- paste0(colnames(df_M_DMSO_vs_Tg_Tg_12hr)[c(2:4, 6:8)], "_M")

#merge sexes for each time point
df_DMSO_vs_Tg_Tg_6hr <- merge(df_F_DMSO_vs_Tg_Tg_6hr, df_M_DMSO_vs_Tg_Tg_6hr, by = c("Pathway", "Time"), all.x = TRUE, all.y = TRUE)
df_DMSO_vs_Tg_Tg_12hr <- merge(df_F_DMSO_vs_Tg_Tg_12hr, df_M_DMSO_vs_Tg_Tg_12hr, by = c("Pathway", "Time"), all.x = TRUE, all.y = TRUE)

#remove rows that neither sex is significant
df_DMSO_vs_Tg_Tg_6hr <- df_DMSO_vs_Tg_Tg_6hr[!(df_DMSO_vs_Tg_Tg_6hr$significance_F == "No" & df_DMSO_vs_Tg_Tg_6hr$significance_M == "No"),] %>% na.omit()
df_DMSO_vs_Tg_Tg_12hr <- df_DMSO_vs_Tg_Tg_12hr[!(df_DMSO_vs_Tg_Tg_12hr$significance_F == "No" & df_DMSO_vs_Tg_Tg_12hr$significance_M == "No"),] %>% na.omit()

#get larger -log10 FDR per row
df_DMSO_vs_Tg_Tg_6hr$order_transf.log.10.FDR_max <- do.call(pmax, c(df_DMSO_vs_Tg_Tg_6hr[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))
df_DMSO_vs_Tg_Tg_12hr$order_transf.log.10.FDR_max <- do.call(pmax, c(df_DMSO_vs_Tg_Tg_12hr[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))

#get larger gene ratio per row
df_DMSO_vs_Tg_Tg_6hr$order_Gene.ratio <- do.call(pmax, c(df_DMSO_vs_Tg_Tg_6hr[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))
df_DMSO_vs_Tg_Tg_12hr$order_Gene.ratio <- do.call(pmax, c(df_DMSO_vs_Tg_Tg_12hr[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))

#get top 5 pathways for each sex based on respective FDR
df_DMSO_vs_Tg_Tg_6hr_F <- dplyr::arrange(df_DMSO_vs_Tg_Tg_6hr, FDR_F)[1:5,]
df_DMSO_vs_Tg_Tg_6hr_M <- dplyr::arrange(df_DMSO_vs_Tg_Tg_6hr, FDR_M)[1:5,]
df_DMSO_vs_Tg_Tg_12hr_F <- dplyr::arrange(df_DMSO_vs_Tg_Tg_12hr, FDR_F)[1:5,]
df_DMSO_vs_Tg_Tg_12hr_M <- dplyr::arrange(df_DMSO_vs_Tg_Tg_12hr, FDR_M)[1:5,]

#combine sexes for each time point
df_DMSO_vs_Tg_Tg_6hr <- rbind(df_DMSO_vs_Tg_Tg_6hr_F, df_DMSO_vs_Tg_Tg_6hr_M)
df_DMSO_vs_Tg_Tg_12hr <- rbind(df_DMSO_vs_Tg_Tg_12hr_F, df_DMSO_vs_Tg_Tg_12hr_M)

#remove duplicated rows
df_DMSO_vs_Tg_Tg_6hr <- df_DMSO_vs_Tg_Tg_6hr[!duplicated(df_DMSO_vs_Tg_Tg_6hr$Pathway),]
df_DMSO_vs_Tg_Tg_12hr <- df_DMSO_vs_Tg_Tg_12hr[!duplicated(df_DMSO_vs_Tg_Tg_12hr$Pathway),]

#order by gene ratio
df_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_DMSO_vs_Tg_Tg_6hr, desc(order_Gene.ratio))
df_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_DMSO_vs_Tg_Tg_12hr, order_Gene.ratio)

#order by -log10 FDR
df_DMSO_vs_Tg_Tg_6hr <- dplyr::arrange(df_DMSO_vs_Tg_Tg_6hr, desc(order_transf.log.10.FDR_max))
df_DMSO_vs_Tg_Tg_12hr <- dplyr::arrange(df_DMSO_vs_Tg_Tg_12hr, order_transf.log.10.FDR_max)

#flip -log10 FDR for 6 hr
df_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR_F <- -df_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR_F
df_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR_M <- -df_DMSO_vs_Tg_Tg_6hr$transf.log.10.FDR_M

#combine time points
df_DMSO_vs_Tg_Tg <- rbind(df_DMSO_vs_Tg_Tg_6hr, df_DMSO_vs_Tg_Tg_12hr)

#remove rows for non-significance for both sexes
df_DMSO_vs_Tg_Tg <- df_DMSO_vs_Tg_Tg[!(df_DMSO_vs_Tg_Tg$significance_F == "> 0.05" & df_DMSO_vs_Tg_Tg$significance_M == "> 0.05"),] %>% na.omit()

#format pathway names
df_DMSO_vs_Tg_Tg$Pathway[df_DMSO_vs_Tg_Tg$Pathway == "Activation of the mRNA upon binding of the cap-binding complex and eIFs, and subsequent binding to 43S"] <- "mRNA activation upon binding of the cap-binding\ncomplex and eIFs and binding to 43S"
df_DMSO_vs_Tg_Tg$Pathway[df_DMSO_vs_Tg_Tg$Pathway == "Response of EIF2AK4 (GCN2) to amino acid deficiency"] <- "EIF2AK4 (GCN2) response to amino acid deficiency"
df_DMSO_vs_Tg_Tg$Pathway[df_DMSO_vs_Tg_Tg$Pathway == "Assembly of collagen fibrils and other multimeric structures"] <- "Assembly of collagen fibrils\nand other multimeric structures"

df_DMSO_vs_Tg_Tg$Pathway <- gsub("the ", "", df_DMSO_vs_Tg_Tg$Pathway)

#set pathway factor level order
df_DMSO_vs_Tg_Tg$Pathway <- factor(df_DMSO_vs_Tg_Tg$Pathway, levels = df_DMSO_vs_Tg_Tg$Pathway)

#set significance factor level order
df_DMSO_vs_Tg_Tg$significance_F <- factor(df_DMSO_vs_Tg_Tg$significance_F, levels = sort(unique(df_DMSO_vs_Tg_Tg$significance_F))[c(1, 3, 2)])

#set x-axis limits, y-axis annotations, and label position
top_label_height_multiplier <- 1.0325

left_limit <- -5
right_limit <- 20

num_row <- nrow(df_DMSO_vs_Tg_Tg)
num_row_bottom <- nrow(df_DMSO_vs_Tg_Tg_6hr[!(df_DMSO_vs_Tg_Tg_6hr$significance_F == "> 0.05" & df_DMSO_vs_Tg_Tg_6hr$significance_M == "> 0.05"),] %>% na.omit())

#create plot
ggplot2::ggplot(df_DMSO_vs_Tg_Tg, ggplot2::aes(y = Pathway)) +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = num_row + 0.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = num_row_bottom + 0.5, yend = num_row_bottom + 0.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -(right_limit - left_limit) * 0.0075, xend = (right_limit - left_limit) * 0.0075, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("Down"), size = 2.5) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("Up")) +
  ggplot2::annotate(geom = "text", x = left_limit - (right_limit - left_limit) * 0.035, y = num_row * top_label_height_multiplier + 0.5, label = c("Regulation at 12 hr"), hjust = 1) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_F), 
                                   shape = significance_F), 
                      alpha = 0) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_F), 
                                   colour = "Female", 
                                   shape = significance_F, 
                                   size = Gene.ratio_F), 
                      alpha = 0.8, 
                      stroke = 1, 
                      fill = adjustcolor("#E7B81E", alpha.f = 0), 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_F == "< 0.01",]) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_F), 
                                   colour = "Female", 
                                   shape = significance_F, 
                                   size = Gene.ratio_F), 
                      alpha = 0.8, 
                      stroke = 1, 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_F == "0.01-0.05",]) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_F), 
                                   shape = significance_F, 
                                   size = Gene.ratio_F), 
                      colour = "grey", 
                      stroke = 1, 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_F == "> 0.05",]) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_M), 
                                   colour = "Male", 
                                   shape = significance_M, 
                                   size = Gene.ratio_M), 
                      alpha = 0.8, 
                      stroke = 1, 
                      fill = adjustcolor("#84B429", alpha.f = 0), 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_M == "< 0.01",]) +  
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_M), 
                                   colour = "Male", 
                                   shape = significance_M, 
                                   size = Gene.ratio_M), 
                      alpha = 0.8, 
                      stroke = 1, 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_M == "0.01-0.05",]) +  
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR_M), 
                                   shape = significance_M, 
                                   size = Gene.ratio_M), 
                      alpha = 0.8, 
                      stroke = 1, 
                      colour = "grey", 
                      data = df_DMSO_vs_Tg_Tg[df_DMSO_vs_Tg_Tg$significance_M == "> 0.05",]) +  
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, num_row + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::scale_colour_manual("Sex",
                               values = c("#E7B81E", "#84B429"),
                               labels = c("Female", "Male")) +
  
  ggplot2::scale_y_discrete(labels = gsub("\\\\n", "\n", as.character(df_DMSO_vs_Tg_Tg$Pathway))) +
  
  ggplot2::scale_alpha(guide = "none") +
  ggplot2::labs(x = bquote(-log[10]*"(FDR)"), 
                colour = "Sex", 
                fill = "Sex", 
                size = "Gene Ratio", 
                shape = "Significance") +
  
  ggplot2::scale_shape_manual(values = 21) +
  
  ggplot2::guides(
    fill = "none", 
    colour = ggplot2::guide_legend(reverse = F, 
                                   override.aes = list(size = 5), 
                                   order = 1), 
    size = ggplot2::guide_legend(reverse = T, 
                                 override.aes = list(colour = "black"),
                                 order = 2), 
    shape = ggplot2::guide_legend(reverse = F, 
                                  override.aes = list(colour = c("black", "black", "grey"),
                                                      fill = c("black", "white", "white"),
                                                      shape = 21,
                                                      size = 4),
                                  order = 3)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 12),
                 axis.text = ggplot2::element_text(colour = "black"),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0.5, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = (nrow(df_DMSO_vs_Tg_Tg) + 1)/5)


## Fig 5D ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
df_6hr_Tg_Proteomics <- read.csv(".\\Data\\6hr Tg serum free_20wk B6 Mouse Islet Proteomics log2_removed_sample12_ imputed_norm to median 0.05.csv", check.names = F, header = TRUE, sep = ",", row.names = NULL)

#format table
colnames(df_6hr_Tg_Proteomics) <- df_6hr_Tg_Proteomics[1,]
df_6hr_Tg_Proteomics <- df_6hr_Tg_Proteomics[-1,]

#get data for sex
df_6hr_Tg_Proteomics_F <- df_6hr_Tg_Proteomics[, c(1:3)]
df_6hr_Tg_Proteomics_M <- df_6hr_Tg_Proteomics[, c(5:7)]

#remove empty rows
df_6hr_Tg_Proteomics_F <- df_6hr_Tg_Proteomics_F[!df_6hr_Tg_Proteomics_F$`Gene Name` == "",]
df_6hr_Tg_Proteomics_M <- df_6hr_Tg_Proteomics_M[!df_6hr_Tg_Proteomics_M$`Gene Name` == "",]

#set gene column name
colnames(df_6hr_Tg_Proteomics_F)[1] <- "gene"
colnames(df_6hr_Tg_Proteomics_M)[1] <- "gene"

#order by Pvalue
df_6hr_Tg_Proteomics_F <- df_6hr_Tg_Proteomics_F[with(df_6hr_Tg_Proteomics_F, order(Pvalue)),]
df_6hr_Tg_Proteomics_M <- df_6hr_Tg_Proteomics_M[with(df_6hr_Tg_Proteomics_M, order(Pvalue)),]

#find common genes for each sex
df_6hr_Tg_Proteomics_common_F <- subset(df_6hr_Tg_Proteomics_F, df_6hr_Tg_Proteomics_F$gene %in% df_6hr_Tg_Proteomics_M$gene)
df_6hr_Tg_Proteomics_common_M <- subset(df_6hr_Tg_Proteomics_M, df_6hr_Tg_Proteomics_M$gene %in% df_6hr_Tg_Proteomics_F$gene)

#annotate column names with sex
colnames(df_6hr_Tg_Proteomics_common_F)[2:3] <- paste0(colnames(df_6hr_Tg_Proteomics_common_F)[2:3], "_F")
colnames(df_6hr_Tg_Proteomics_common_M)[2:3] <- paste0(colnames(df_6hr_Tg_Proteomics_common_M)[2:3], "_M")

#make male and female graphs unique
df_6hr_Tg_Proteomics_F_unique_list <- subset(df_6hr_Tg_Proteomics_F, !df_6hr_Tg_Proteomics_F$gene %in% df_6hr_Tg_Proteomics_M$gene)
df_6hr_Tg_Proteomics_M_unique_list <- subset(df_6hr_Tg_Proteomics_M, !df_6hr_Tg_Proteomics_M$gene %in% df_6hr_Tg_Proteomics_F$gene)

#index colour of unique genes
df_6hr_Tg_Proteomics_F$unique <- ifelse(df_6hr_Tg_Proteomics_F$gene %in% df_6hr_Tg_Proteomics_F_unique_list$gene, "red", "black")
df_6hr_Tg_Proteomics_M$unique <- ifelse(df_6hr_Tg_Proteomics_M$gene %in% df_6hr_Tg_Proteomics_M_unique_list$gene, "red", "black")

#star unique genes
df_6hr_Tg_Proteomics_F$gene <- ifelse(df_6hr_Tg_Proteomics_F$gene %in% df_6hr_Tg_Proteomics_F_unique_list$gene, paste(df_6hr_Tg_Proteomics_F$gene, "*"), df_6hr_Tg_Proteomics_F$gene)
df_6hr_Tg_Proteomics_M$gene <- ifelse(df_6hr_Tg_Proteomics_M$gene %in% df_6hr_Tg_Proteomics_M_unique_list$gene, paste(df_6hr_Tg_Proteomics_M$gene, "*"), df_6hr_Tg_Proteomics_M$gene)

#select logFC and colour annoation columns
df_6hr_Tg_Proteomics_F_specific_columns <- df_6hr_Tg_Proteomics_F[, c("gene", "log2FC", "unique")]
df_6hr_Tg_Proteomics_M_specific_columns <- df_6hr_Tg_Proteomics_M[, c("gene", "log2FC", "unique")]

#rename columns
colnames(df_6hr_Tg_Proteomics_F_specific_columns)[1:2] <- c("gene", "Female")
colnames(df_6hr_Tg_Proteomics_M_specific_columns)[1:2] <- c("gene", "Male")

#order by abs(logFC)
df_6hr_Tg_Proteomics_F_specific_columns <- dplyr::arrange(df_6hr_Tg_Proteomics_F_specific_columns, desc(abs(as.numeric(df_6hr_Tg_Proteomics_F_specific_columns$Female))))
df_6hr_Tg_Proteomics_M_specific_columns <- dplyr::arrange(df_6hr_Tg_Proteomics_M_specific_columns, desc(abs(as.numeric(df_6hr_Tg_Proteomics_M_specific_columns$Male))))

#select top 45 proteins
df_6hr_Tg_Proteomics_F_specific_columns <- df_6hr_Tg_Proteomics_F_specific_columns[c(1:45),]
df_6hr_Tg_Proteomics_M_specific_columns <- df_6hr_Tg_Proteomics_M_specific_columns[c(1:45),]

#create label colour index
female_unique_colour_index <- df_6hr_Tg_Proteomics_F_specific_columns$unique
male_unique_colour_index <- df_6hr_Tg_Proteomics_M_specific_columns$unique

#select logFC column
df_logFC_female <- df_6hr_Tg_Proteomics_F_specific_columns[, "Female"]
df_logFC_male <- df_6hr_Tg_Proteomics_M_specific_columns[, "Male"]

#create matrix
matrix_logFC_male <- as.matrix(as.numeric(df_logFC_male))
matrix_logFC_female <- as.matrix(as.numeric(df_logFC_female))

#set sex as column name
colnames(matrix_logFC_male) <- c("Male")
colnames(matrix_logFC_female) <- c("Female")

#set genes as row names
rownames(matrix_logFC_male) <- df_6hr_Tg_Proteomics_M_specific_columns$gene  
rownames(matrix_logFC_female) <- df_6hr_Tg_Proteomics_F_specific_columns$gene 

#create female heatmap
heatmap_female <- ComplexHeatmap::Heatmap(matrix_logFC_female,
                                          name = "matrix_logFC_heatmap",
                                          row_names_side = "right",
                                          row_names_gp = grid::gpar(fontface = "italic", 
                                                                    fontsize = 12,
                                                                    col = female_unique_colour_index),
                                          show_row_dend = FALSE,
                                          top_annotation = ComplexHeatmap::columnAnnotation(
                                            empty = ComplexHeatmap::anno_empty(border = FALSE)
                                          ),
                                          column_order = 1:ncol(matrix_logFC_female),
                                          show_column_names =  FALSE,
                                          width = ncol(matrix_logFC_female) * grid::unit(40, "mm"),
                                          height = 45 * grid::unit(7, "mm"),
                                          border_gp = grid::gpar(col = "black"),
                                          col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                                          na_col = "black",
                                          heatmap_legend_param = list(title = NULL,
                                                                      title = "logFC",
                                                                      title_position = "lefttop",
                                                                      legend_height = grid::unit(8, "cm"),
                                                                      direction = "horizontal",
                                                                      legend_width = grid::unit(4, "cm"),
                                                                      labels_gp = grid::gpar(fontsize = 12)))

#create male heatmap
heatmap_male <- ComplexHeatmap::Heatmap(matrix_logFC_male,
                                        name = "matrix_logFC_heatmap",
                                        row_names_side = "right",
                                        row_names_gp = grid::gpar(fontface = "italic", 
                                                                  fontsize = 12,
                                                                  col = male_unique_colour_index),
                                        show_row_dend = FALSE,
                                        top_annotation = ComplexHeatmap::columnAnnotation(
                                          empty = ComplexHeatmap::anno_empty(border = FALSE)
                                        ),
                                        column_order = 1:ncol(matrix_logFC_male),
                                        show_column_names =  FALSE,
                                        width = ncol(matrix_logFC_male) * grid::unit(40, "mm"),
                                        height = 45 * grid::unit(7, "mm"),
                                        border_gp = grid::gpar(col = "black"),
                                        col = circlize::colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red")),
                                        na_col = "black",
                                        heatmap_legend_param = list(title = NULL,
                                                                    title = "logFC",
                                                                    title_position = "lefttop",
                                                                    legend_height = grid::unit(8, "cm"),
                                                                    direction = "horizontal",
                                                                    legend_width = grid::unit(4, "cm"),
                                                                    labels_gp = grid::gpar(fontsize = 12)))

#draw female heatmap
ComplexHeatmap::draw(heatmap_female, heatmap_legend_side = "bottom")

#select viewport and define dimensions
grid::seekViewport("annotation_empty_1")
grid::anno_empty_1_loc1 = deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
grid::anno_empty_1_loc2 = deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#E7B81E", 
                                col = "#E7B81E", 
                                lwd = 0))

grid::grid.text("Female",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 12))

#select legend viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend label
grid::grid.text(bquote(log[2]*"(Fold Change)"),
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 12))

#draw male heatmap
ComplexHeatmap::draw(heatmap_male, heatmap_legend_side = "bottom")

grid::seekViewport("annotation_empty_1")
grid::anno_empty_1_loc1 = deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
grid::anno_empty_1_loc2 = deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create top labels
grid::grid.rect(x = 0,
                y = 0,
                width = (anno_empty_1_loc2$x - anno_empty_1_loc1$x) * 1,
                height = (anno_empty_1_loc2$y - anno_empty_1_loc1$y) * 0.65,
                just = c("left", "bottom"),
                gp = grid::gpar(fill = "#84B429", 
                                col = "#84B429", 
                                lwd = 0))

grid::grid.text("Male",
                x = 0.5,
                y = 0.65/2,
                gp = grid::gpar(fontsize = 12))

#select legend viewport and define dimensions
grid::seekViewport("heatmap_legend")
loc1 = grid::deviceLoc(x = grid::unit(0, "npc"), y = grid::unit(0, "npc"))
loc2 = grid::deviceLoc(x = grid::unit(1, "npc"), y = grid::unit(1, "npc"))

#create legend label
grid::grid.text(bquote(log[2]*"(Fold Change)"),
                x = 0.5,
                y = -0.5,
                gp = grid::gpar(fontsize = 12))


### Figure 5 - Supplement 1 ###

## Fig 5 S1A ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import file
data_rlog <- read.csv(".\\Data\\global_genes.rlog.txt", stringsAsFactors = F, check.names = T, sep = ",")
data_counts <- read.csv(".\\Data\\Tg_global_genes.counts.csv")

#filter out empty and zero-only rows
data_rlog <- na.omit(data_rlog)
data_rlog <- subset(data_rlog, F1.12hr.Tg.and.F8.12hr.Tg != 0)

#filter for only rows with more than 10 counts total
data_counts <- data_counts[rowSums(data_counts[, c(2:28)]) > 10,]

#filter rlog data for same as filtered count data genes
data_rlog <- subset(data_rlog, data_rlog$X %in% data_counts$X)

#get Fold Change of Male DMSO vs Female DMSO
data_rlog$FoldChange <- rowMeans(as.matrix(data_rlog[, grepl("DMSO", colnames(data_rlog)) & grepl("M\\d", colnames(data_rlog))]))/rowMeans(as.matrix(data_rlog[, grepl("DMSO", colnames(data_rlog)) & grepl("F\\d", colnames(data_rlog))]))

#order by Fold Change
data_rlog <- data_rlog[with(data_rlog, rev(order(FoldChange))),][, -ncol(data_rlog)]

#convert to Z-scores
data_rlog$mean_all <- rowMeans(data_rlog[, c(2:28)])
data_rlog$sd_all <- apply(data_rlog[, c(2:28)], 1, sd)
Z_score <- function(x) with(data_rlog, (x - mean_all)/sd_all)
df_Z_score_total <- cbind(apply(data_rlog[, c(2:28)], 2, Z_score))

#flip rows and columns
df_Z_score_total_flipped <- as.data.frame(t(df_Z_score_total))                                                                                         #flipped columns and rows

#add old sample name column
df_Z_score_total_flipped$Old <- rownames(df_Z_score_total_flipped)

#extract sex, time, and treatment
df_Z_score_total_flipped$Sex <- ifelse(grepl("^F", df_Z_score_total_flipped$Old), "Female", "Male")
df_Z_score_total_flipped$Time <- ifelse(grepl("6hr", df_Z_score_total_flipped$Old), "6 hr", "12 hr")
df_Z_score_total_flipped$Treatment <- ifelse(grepl("Tg", df_Z_score_total_flipped$Old), "Tg", "DMSO")

#create condition column
df_Z_score_total_flipped$Condition <- paste(df_Z_score_total_flipped$Sex, df_Z_score_total_flipped$Time, df_Z_score_total_flipped$Treatment)

#import sample renaming table
Tg_rename <- read.csv(".\\Data\\Tg RNA Seq Sample Renaming.csv", stringsAsFactors = F, check.names = T, sep = ",")

#combine data with renaming table
df_Z_score_total_flipped <- merge(df_Z_score_total_flipped, Tg_rename, by = "Old", all.x = TRUE)

#rename row names
rownames(df_Z_score_total_flipped) <- df_Z_score_total_flipped$New

#order based on condition
df_Z_score_total_flipped <- df_Z_score_total_flipped[order(df_Z_score_total_flipped$Condition),]

#select data columns
df_Z_score_total_flipped <- df_Z_score_total_flipped[, grepl("\\d", colnames(df_Z_score_total_flipped))]

#convert to matrix
matrix_Z_score_total <- data.matrix(df_Z_score_total_flipped)

#remove column names
colnames(matrix_Z_score_total) <- NULL

#format row names
rownames(matrix_Z_score_total) <- gsub("\\.", " ", rownames(matrix_Z_score_total))
rownames(matrix_Z_score_total) <- gsub("6hr", "6 hr", rownames(matrix_Z_score_total))   
rownames(matrix_Z_score_total) <- gsub("12hr", "12 hr", rownames(matrix_Z_score_total))

#create heatmap
heatmap_list_total <- ComplexHeatmap::Heatmap(matrix_Z_score_total,
                                              name = "matrix_Z_score_total",
                                              row_names_side = "left", 
                                              row_names_gp = grid::gpar(
                                                col = c(
                                                  rep("#606060", 3),
                                                  rep("#E7B91E", 3),
                                                  rep("#000000", 4),
                                                  rep("#E7821E", 4),
                                                  rep("#606060", 3),
                                                  rep("#316219", 3),
                                                  rep("#000000", 4),
                                                  rep("#84B42A", 3)
                                                ),
                                                fontface = c(rep("italic", 14),
                                                             rep("plain", 13))
                                              ), 
                                              row_dend_width = grid::unit(2, "cm"),
                                              right_annotation = ComplexHeatmap::rowAnnotation(empty2 = ComplexHeatmap::anno_empty(border = FALSE,
                                                                                                                                   width = grid::unit(4, "mm"))),
                                              column_order = 1:ncol(matrix_Z_score_total),
                                              width = grid::unit(240, "mm"), 
                                              height = nrow(matrix_Z_score_total) * grid::unit(6, "mm"),
                                              border_gp = grid::gpar(col = "black"), 
                                              heatmap_legend_param = list(title = "Column\nZ-score",
                                                                          title_position = "topcenter",
                                                                          legend_height = grid::unit(4, "cm")))
#draw heatmap
ComplexHeatmap::draw(heatmap_list_total)


### Figure 5 - Supplement 2 ###

## Fig 5 S2A ##

#clear environment
rm(list = ls())

#import files
df_F_6hr_DMSOvsTg_DMSO <- read.csv(".\\Data\\DMSO_vs_Tg - FemaleControl6hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_F_6hr_DMSOvsTg_Tg <- read.csv(".\\Data\\DMSO_vs_Tg - FemaleTg6hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_6hr_DMSOvsTg_DMSO <- read.csv(".\\Data\\DMSO_vs_Tg - MaleControl6hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_6hr_DMSOvsTg_Tg <- read.csv(".\\Data\\DMSO_vs_Tg - MaleTg6hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]

#annotate treatment
df_F_6hr_DMSOvsTg_DMSO$Treatment <- "DMSO"
df_F_6hr_DMSOvsTg_Tg$Treatment <- "Tg"
df_M_6hr_DMSOvsTg_DMSO$Treatment <- "DMSO"
df_M_6hr_DMSOvsTg_Tg$Treatment <- "Tg"

#annotate sex
df_F_6hr_DMSOvsTg_DMSO$Sex <- "Female"
df_F_6hr_DMSOvsTg_Tg$Sex <- "Female"
df_M_6hr_DMSOvsTg_DMSO$Sex <- "Male"
df_M_6hr_DMSOvsTg_Tg$Sex <- "Male"

#fix column names
colnames(df_F_6hr_DMSOvsTg_DMSO) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_F_6hr_DMSOvsTg_Tg) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_M_6hr_DMSOvsTg_DMSO) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_M_6hr_DMSOvsTg_Tg) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_F_6hr_DMSOvsTg_DMSO$Gene.ratio <- df_F_6hr_DMSOvsTg_DMSO$Entities.found/592
df_F_6hr_DMSOvsTg_Tg$Gene.ratio <- df_F_6hr_DMSOvsTg_Tg$Entities.found/565
df_M_6hr_DMSOvsTg_DMSO$Gene.ratio <- df_M_6hr_DMSOvsTg_DMSO$Entities.found/606
df_M_6hr_DMSOvsTg_Tg$Gene.ratio <- df_M_6hr_DMSOvsTg_Tg$Entities.found/560

#get -log10 of FDR
df_F_6hr_DMSOvsTg_DMSO$transf.log.10.FDR <- -log10(df_F_6hr_DMSOvsTg_DMSO$FDR)
df_F_6hr_DMSOvsTg_Tg$transf.log.10.FDR <- -log10(df_F_6hr_DMSOvsTg_Tg$FDR)
df_M_6hr_DMSOvsTg_DMSO$transf.log.10.FDR <- -log10(df_M_6hr_DMSOvsTg_DMSO$FDR)
df_M_6hr_DMSOvsTg_Tg$transf.log.10.FDR <- -log10(df_M_6hr_DMSOvsTg_Tg$FDR)

#order by gene ratio
df_F_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_F_6hr_DMSOvsTg_DMSO, desc(Gene.ratio))
df_F_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_F_6hr_DMSOvsTg_Tg, desc(Gene.ratio))
df_M_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_M_6hr_DMSOvsTg_DMSO, desc(Gene.ratio))
df_M_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_M_6hr_DMSOvsTg_Tg, desc(Gene.ratio))

#order by -log10(FDR)
df_F_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_F_6hr_DMSOvsTg_DMSO, desc(transf.log.10.FDR))
df_F_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_F_6hr_DMSOvsTg_Tg, desc(transf.log.10.FDR))
df_M_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_M_6hr_DMSOvsTg_DMSO, desc(transf.log.10.FDR))
df_M_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_M_6hr_DMSOvsTg_Tg, desc(transf.log.10.FDR))

#categorize significance
df_F_6hr_DMSOvsTg_DMSO$significance <- ifelse(df_F_6hr_DMSOvsTg_DMSO$FDR < 0.01, "< 0.01", ifelse(df_F_6hr_DMSOvsTg_DMSO$FDR > 0.01 & df_F_6hr_DMSOvsTg_DMSO$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_F_6hr_DMSOvsTg_Tg$significance <- ifelse(df_F_6hr_DMSOvsTg_Tg$FDR < 0.01, "< 0.01", ifelse(df_F_6hr_DMSOvsTg_Tg$FDR > 0.01 & df_F_6hr_DMSOvsTg_Tg$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_6hr_DMSOvsTg_DMSO$significance <- ifelse(df_M_6hr_DMSOvsTg_DMSO$FDR < 0.01, "< 0.01", ifelse(df_M_6hr_DMSOvsTg_DMSO$FDR > 0.01 & df_M_6hr_DMSOvsTg_DMSO$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_6hr_DMSOvsTg_Tg$significance <- ifelse(df_M_6hr_DMSOvsTg_Tg$FDR < 0.01, "< 0.01", ifelse(df_M_6hr_DMSOvsTg_Tg$FDR > 0.01 & df_M_6hr_DMSOvsTg_Tg$FDR < 0.05, "0.01-0.05", "> 0.05"))

#annotate sex
colnames(df_F_6hr_DMSOvsTg_DMSO)[c(2:4, 6:9)] <- paste0(colnames(df_F_6hr_DMSOvsTg_DMSO)[c(2:4, 6:9)], "_F")
colnames(df_F_6hr_DMSOvsTg_Tg)[c(2:4, 6:9)] <- paste0(colnames(df_F_6hr_DMSOvsTg_Tg)[c(2:4, 6:9)], "_F")
colnames(df_M_6hr_DMSOvsTg_DMSO)[c(2:4, 6:9)] <- paste0(colnames(df_M_6hr_DMSOvsTg_DMSO)[c(2:4, 6:9)], "_M")
colnames(df_M_6hr_DMSOvsTg_Tg)[c(2:4, 6:9)] <- paste0(colnames(df_M_6hr_DMSOvsTg_Tg)[c(2:4, 6:9)], "_M")

#merge sexes for each treatment
df_6hr_DMSOvsTg_DMSO <- merge(df_F_6hr_DMSOvsTg_DMSO, df_M_6hr_DMSOvsTg_DMSO, by = c("Pathway", "Treatment"), all.x = TRUE, all.y = TRUE)
df_6hr_DMSOvsTg_Tg <- merge(df_F_6hr_DMSOvsTg_Tg, df_M_6hr_DMSOvsTg_Tg, by = c("Pathway", "Treatment"), all.x = TRUE, all.y = TRUE)

#get larger -log10 FDR per row
df_6hr_DMSOvsTg_DMSO$order_transf.log.10.FDR_max <- do.call(pmax, c(df_6hr_DMSOvsTg_DMSO[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))
df_6hr_DMSOvsTg_Tg$order_transf.log.10.FDR_max <- do.call(pmax, c(df_6hr_DMSOvsTg_Tg[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))

#get larger gene ratio per row
df_6hr_DMSOvsTg_DMSO$order_Gene.ratio <- do.call(pmax, c(df_6hr_DMSOvsTg_DMSO[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))
df_6hr_DMSOvsTg_Tg$order_Gene.ratio <- do.call(pmax, c(df_6hr_DMSOvsTg_Tg[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))

#get top 10 pathways for each sex based on respective FDR
df_6hr_DMSOvsTg_DMSO_F <- dplyr::arrange(df_6hr_DMSOvsTg_DMSO, FDR_F)[1:10,]
df_6hr_DMSOvsTg_DMSO_M <- dplyr::arrange(df_6hr_DMSOvsTg_DMSO, FDR_M)[1:10,]
df_6hr_DMSOvsTg_Tg_F <- dplyr::arrange(df_6hr_DMSOvsTg_Tg, FDR_F)[1:10,]
df_6hr_DMSOvsTg_Tg_M <- dplyr::arrange(df_6hr_DMSOvsTg_Tg, FDR_M)[1:10,]

#combine sexes for each treatment
df_6hr_DMSOvsTg_DMSO <- rbind(df_6hr_DMSOvsTg_DMSO_F, df_6hr_DMSOvsTg_DMSO_M)
df_6hr_DMSOvsTg_Tg <- rbind(df_6hr_DMSOvsTg_Tg_F, df_6hr_DMSOvsTg_Tg_M)

#remove duplicated rows
df_6hr_DMSOvsTg_DMSO <- df_6hr_DMSOvsTg_DMSO[!duplicated(df_6hr_DMSOvsTg_DMSO$Pathway),]
df_6hr_DMSOvsTg_Tg <- df_6hr_DMSOvsTg_Tg[!duplicated(df_6hr_DMSOvsTg_Tg$Pathway),]

#order by gene ratio
df_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_6hr_DMSOvsTg_DMSO, desc(order_Gene.ratio))
df_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_6hr_DMSOvsTg_Tg, order_Gene.ratio)

#order by -log10 FDR
df_6hr_DMSOvsTg_DMSO <- dplyr::arrange(df_6hr_DMSOvsTg_DMSO, desc(order_transf.log.10.FDR_max))
df_6hr_DMSOvsTg_Tg <- dplyr::arrange(df_6hr_DMSOvsTg_Tg, order_transf.log.10.FDR_max)

#flip -log10 FDR for DMSO
df_6hr_DMSOvsTg_DMSO$transf.log.10.FDR_F <- -df_6hr_DMSOvsTg_DMSO$transf.log.10.FDR_F
df_6hr_DMSOvsTg_DMSO$transf.log.10.FDR_M <- -df_6hr_DMSOvsTg_DMSO$transf.log.10.FDR_M

#combine treatments
df_6hr_DMSOvsTg <- rbind(df_6hr_DMSOvsTg_DMSO, df_6hr_DMSOvsTg_Tg)

#remove rows for non-significance for both sexes
df_6hr_DMSOvsTg <- df_6hr_DMSOvsTg[!(df_6hr_DMSOvsTg$significance_F == "> 0.05" & df_6hr_DMSOvsTg$significance_M == "> 0.05"),] %>% na.omit()

#set pathway factor level order
df_6hr_DMSOvsTg$Pathway <- factor(df_6hr_DMSOvsTg$Pathway, levels = df_6hr_DMSOvsTg$Pathway)

#separate based on treatment
df_6hr_DMSOvsTg_plotted_Tg <- df_6hr_DMSOvsTg[df_6hr_DMSOvsTg$Treatment == "Tg",]
df_6hr_DMSOvsTg_plotted_DMSO <- df_6hr_DMSOvsTg[df_6hr_DMSOvsTg$Treatment == "DMSO",]

#select columns for each sex and treatment
df_6hr_DMSOvsTg_plotted_Tg_F <- df_6hr_DMSOvsTg_plotted_Tg[, c("Pathway", "Treatment", "Entities.found_F", "Entities.total_F", "FDR_F", "Sex_F", "Gene.ratio_F", "transf.log.10.FDR_F", "significance_F", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_6hr_DMSOvsTg_plotted_Tg_M <- df_6hr_DMSOvsTg_plotted_Tg[, c("Pathway", "Treatment", "Entities.found_M", "Entities.total_M", "FDR_M", "Sex_M", "Gene.ratio_M", "transf.log.10.FDR_M", "significance_M", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_6hr_DMSOvsTg_plotted_DMSO_F <- df_6hr_DMSOvsTg_plotted_DMSO[, c("Pathway", "Treatment", "Entities.found_F", "Entities.total_F", "FDR_F", "Sex_F", "Gene.ratio_F", "transf.log.10.FDR_F", "significance_F", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_6hr_DMSOvsTg_plotted_DMSO_M <- df_6hr_DMSOvsTg_plotted_DMSO[, c("Pathway", "Treatment", "Entities.found_M", "Entities.total_M", "FDR_M", "Sex_M", "Gene.ratio_M", "transf.log.10.FDR_M", "significance_M", "order_transf.log.10.FDR_max", "order_Gene.ratio")]

#rename columns
colnames(df_6hr_DMSOvsTg_plotted_Tg_F) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_6hr_DMSOvsTg_plotted_Tg_M) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_6hr_DMSOvsTg_plotted_DMSO_F) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_6hr_DMSOvsTg_plotted_DMSO_M) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")

#combine sexes for each treatment
df_6hr_DMSOvsTg_plotted_Tg <- rbind(df_6hr_DMSOvsTg_plotted_Tg_F, df_6hr_DMSOvsTg_plotted_Tg_M)
df_6hr_DMSOvsTg_plotted_DMSO <- rbind(df_6hr_DMSOvsTg_plotted_DMSO_F, df_6hr_DMSOvsTg_plotted_DMSO_M)

#create dot colour column
df_6hr_DMSOvsTg_plotted_Tg$colour_index <- paste(df_6hr_DMSOvsTg_plotted_Tg$Sex, df_6hr_DMSOvsTg_plotted_Tg$significance)
df_6hr_DMSOvsTg_plotted_DMSO$colour_index <- paste(df_6hr_DMSOvsTg_plotted_DMSO$Sex, df_6hr_DMSOvsTg_plotted_DMSO$significance)

#set significance factor level order
df_6hr_DMSOvsTg_plotted_Tg$significance <- factor(df_6hr_DMSOvsTg_plotted_Tg$significance, levels = sort(unique(df_6hr_DMSOvsTg_plotted_Tg$significance))[c(1, 3, 2)])

#set x-axis limits, y-axis annotations, and label position
top_label_height_multiplier <- 1.0305

left_limit <- -10
right_limit <- 15

num_row <- nrow(df_6hr_DMSOvsTg)
num_row_bottom <- nrow(df_6hr_DMSOvsTg_DMSO[!(df_6hr_DMSOvsTg_DMSO$significance_F == "> 0.05" & df_6hr_DMSOvsTg_DMSO$significance_M == "> 0.05"),] %>% na.omit())

#create plot
ggplot2::ggplot(df_6hr_DMSOvsTg, ggplot2::aes(y = Pathway)) +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = num_row + 0.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = num_row_bottom + 0.5, yend = num_row_bottom + 0.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -(right_limit - left_limit) * 0.0075, xend = (right_limit - left_limit) * 0.0075, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("DMSO"), size = 3) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("Tg")) +
  ggplot2::annotate(geom = "text", x = left_limit - (right_limit - left_limit) * 0.035, y = num_row * top_label_height_multiplier + 0.5, label = c("Treatment"), hjust = 1) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(order_transf.log.10.FDR_max)), 
                      alpha = 0, 
                      data = df_6hr_DMSOvsTg) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR), 
                                   colour = colour_index, 
                                   fill = Sex, 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1, 
                      position = ggstance::position_dodgev(height = 0.2), 
                      data = df_6hr_DMSOvsTg_plotted_Tg) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR), 
                                   colour = colour_index, 
                                   fill = Sex, 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1, 
                      position = ggstance::position_dodgev(height = 0.2), 
                      data = df_6hr_DMSOvsTg_plotted_DMSO) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, num_row + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::scale_size_continuous(range = c(0.5, 5)) +
  
  ggplot2::scale_fill_manual("Sex",
                             values = c("#E7B91E", "#84B429"),
                             labels = c("Female", "Male")) +
  
  ggplot2::scale_colour_manual("Significance", values = c("#E7B91E", "#E7B91E", "#84B429", "grey", "#84B429")) +
  
  ggplot2::scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 45), drop = FALSE) +
  
  ggplot2::scale_alpha(guide = "none") +
  
  ggplot2::labs(x = bquote(-log[10]*"(FDR)"), 
                colour = "Sex", 
                fill = "Sex", 
                size = "Gene Ratio", 
                shape = "Significance") +
  
  ggplot2::scale_shape_manual(values = c(21, 1, 1)) +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(reverse = F, 
                                 override.aes = list(size = 5, 
                                                     colour = c("#E7B91E", "#84B429")), 
                                 order = 1), 
    colour = "none", 
    size = ggplot2::guide_legend(reverse = T, 
                                 override.aes = list(colour = "black"),
                                 order = 2), 
    shape = ggplot2::guide_legend(reverse = F, 
                                  override.aes = list(colour = c("black", "black", "grey"), 
                                                      size = 4, 
                                                      fill = c("black", "white", "white")),
                                  order = 3)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 axis.text = ggplot2::element_text(colour = "black"),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0.5, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = (nrow(df_6hr_DMSOvsTg) + 1)/5)


## Fig 5 S2B ##

#clear environment
rm(list = ls())

#set seed
set.seed(123)

#import files
df_F_12hr_DMSOvsTg_DMSO <- read.csv(".\\Data\\DMSO_vs_Tg - FemaleControl12hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_F_12hr_DMSOvsTg_Tg <- read.csv(".\\Data\\DMSO_vs_Tg - FemaleTg12hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_12hr_DMSOvsTg_DMSO <- read.csv(".\\Data\\DMSO_vs_Tg - MaleControl12hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]
df_M_12hr_DMSOvsTg_Tg <- read.csv(".\\Data\\DMSO_vs_Tg - MaleTg12hr_0.01_top_1000_Reactome_Results.csv", stringsAsFactors = F, check.names = T)[, c(2:4, 7)]

#annotate treatment
df_F_12hr_DMSOvsTg_DMSO$Treatment <- "DMSO"
df_F_12hr_DMSOvsTg_Tg$Treatment <- "Tg"
df_M_12hr_DMSOvsTg_DMSO$Treatment <- "DMSO"
df_M_12hr_DMSOvsTg_Tg$Treatment <- "Tg"

#annotate sex
df_F_12hr_DMSOvsTg_DMSO$Sex <- "Female"
df_F_12hr_DMSOvsTg_Tg$Sex <- "Female"
df_M_12hr_DMSOvsTg_DMSO$Sex <- "Male"
df_M_12hr_DMSOvsTg_Tg$Sex <- "Male"

#fix column names
colnames(df_F_12hr_DMSOvsTg_DMSO) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_F_12hr_DMSOvsTg_Tg) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_M_12hr_DMSOvsTg_DMSO) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")
colnames(df_M_12hr_DMSOvsTg_Tg) <- c("Pathway", "Entities.found", "Entities.total", "FDR", "Treatment", "Sex")

#calculate ratio of number of genes from data in pathway to total number of recognized genes from data
df_F_12hr_DMSOvsTg_DMSO$Gene.ratio <- df_F_12hr_DMSOvsTg_DMSO$Entities.found/600
df_F_12hr_DMSOvsTg_Tg$Gene.ratio <- df_F_12hr_DMSOvsTg_Tg$Entities.found/588
df_M_12hr_DMSOvsTg_DMSO$Gene.ratio <- df_M_12hr_DMSOvsTg_DMSO$Entities.found/596
df_M_12hr_DMSOvsTg_Tg$Gene.ratio <- df_M_12hr_DMSOvsTg_Tg$Entities.found/572

#get -log10 of FDR
df_F_12hr_DMSOvsTg_DMSO$transf.log.10.FDR <- -log10(df_F_12hr_DMSOvsTg_DMSO$FDR)
df_F_12hr_DMSOvsTg_Tg$transf.log.10.FDR <- -log10(df_F_12hr_DMSOvsTg_Tg$FDR)
df_M_12hr_DMSOvsTg_DMSO$transf.log.10.FDR <- -log10(df_M_12hr_DMSOvsTg_DMSO$FDR)
df_M_12hr_DMSOvsTg_Tg$transf.log.10.FDR <- -log10(df_M_12hr_DMSOvsTg_Tg$FDR)

#order by gene ratio
df_F_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_F_12hr_DMSOvsTg_DMSO, desc(Gene.ratio))
df_F_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_F_12hr_DMSOvsTg_Tg, desc(Gene.ratio))
df_M_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_M_12hr_DMSOvsTg_DMSO, desc(Gene.ratio))
df_M_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_M_12hr_DMSOvsTg_Tg, desc(Gene.ratio))

#order by -log10(FDR)
df_F_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_F_12hr_DMSOvsTg_DMSO, desc(transf.log.10.FDR))
df_F_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_F_12hr_DMSOvsTg_Tg, desc(transf.log.10.FDR))
df_M_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_M_12hr_DMSOvsTg_DMSO, desc(transf.log.10.FDR))
df_M_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_M_12hr_DMSOvsTg_Tg, desc(transf.log.10.FDR))

#categorize significance
df_F_12hr_DMSOvsTg_DMSO$significance <- ifelse(df_F_12hr_DMSOvsTg_DMSO$FDR < 0.01, "< 0.01", ifelse(df_F_12hr_DMSOvsTg_DMSO$FDR > 0.01 & df_F_12hr_DMSOvsTg_DMSO$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_F_12hr_DMSOvsTg_Tg$significance <- ifelse(df_F_12hr_DMSOvsTg_Tg$FDR < 0.01, "< 0.01", ifelse(df_F_12hr_DMSOvsTg_Tg$FDR > 0.01 & df_F_12hr_DMSOvsTg_Tg$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_12hr_DMSOvsTg_DMSO$significance <- ifelse(df_M_12hr_DMSOvsTg_DMSO$FDR < 0.01, "< 0.01", ifelse(df_M_12hr_DMSOvsTg_DMSO$FDR > 0.01 & df_M_12hr_DMSOvsTg_DMSO$FDR < 0.05, "0.01-0.05", "> 0.05"))
df_M_12hr_DMSOvsTg_Tg$significance <- ifelse(df_M_12hr_DMSOvsTg_Tg$FDR < 0.01, "< 0.01", ifelse(df_M_12hr_DMSOvsTg_Tg$FDR > 0.01 & df_M_12hr_DMSOvsTg_Tg$FDR < 0.05, "0.01-0.05", "> 0.05"))

#annotate sex
colnames(df_F_12hr_DMSOvsTg_DMSO)[c(2:4, 6:9)] <- paste0(colnames(df_F_12hr_DMSOvsTg_DMSO)[c(2:4, 6:9)], "_F")
colnames(df_F_12hr_DMSOvsTg_Tg)[c(2:4, 6:9)] <- paste0(colnames(df_F_12hr_DMSOvsTg_Tg)[c(2:4, 6:9)], "_F")
colnames(df_M_12hr_DMSOvsTg_DMSO)[c(2:4, 6:9)] <- paste0(colnames(df_M_12hr_DMSOvsTg_DMSO)[c(2:4, 6:9)], "_M")
colnames(df_M_12hr_DMSOvsTg_Tg)[c(2:4, 6:9)] <- paste0(colnames(df_M_12hr_DMSOvsTg_Tg)[c(2:4, 6:9)], "_M")

#merge sexes for each treatment
df_12hr_DMSOvsTg_DMSO <- merge(df_F_12hr_DMSOvsTg_DMSO, df_M_12hr_DMSOvsTg_DMSO, by = c("Pathway", "Treatment"), all.x = TRUE, all.y = TRUE)
df_12hr_DMSOvsTg_Tg <- merge(df_F_12hr_DMSOvsTg_Tg, df_M_12hr_DMSOvsTg_Tg, by = c("Pathway", "Treatment"), all.x = TRUE, all.y = TRUE)

#get larger -log10 FDR per row
df_12hr_DMSOvsTg_DMSO$order_transf.log.10.FDR_max <- do.call(pmax, c(df_12hr_DMSOvsTg_DMSO[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))
df_12hr_DMSOvsTg_Tg$order_transf.log.10.FDR_max <- do.call(pmax, c(df_12hr_DMSOvsTg_Tg[c("transf.log.10.FDR_F", "transf.log.10.FDR_M")], list(na.rm = TRUE)))

#get larger gene ratio per row
df_12hr_DMSOvsTg_DMSO$order_Gene.ratio <- do.call(pmax, c(df_12hr_DMSOvsTg_DMSO[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))
df_12hr_DMSOvsTg_Tg$order_Gene.ratio <- do.call(pmax, c(df_12hr_DMSOvsTg_Tg[c("Gene.ratio_F", "Gene.ratio_M")], list(na.rm = TRUE)))

#get top 10 pathways for each sex based on respective FDR
df_12hr_DMSOvsTg_DMSO_F <- dplyr::arrange(df_12hr_DMSOvsTg_DMSO, FDR_F)[1:10,]
df_12hr_DMSOvsTg_DMSO_M <- dplyr::arrange(df_12hr_DMSOvsTg_DMSO, FDR_M)[1:10,]
df_12hr_DMSOvsTg_Tg_F <- dplyr::arrange(df_12hr_DMSOvsTg_Tg, FDR_F)[1:10,]
df_12hr_DMSOvsTg_Tg_M <- dplyr::arrange(df_12hr_DMSOvsTg_Tg, FDR_M)[1:10,]

#combine sexes for each treatment
df_12hr_DMSOvsTg_DMSO <- rbind(df_12hr_DMSOvsTg_DMSO_F, df_12hr_DMSOvsTg_DMSO_M)
df_12hr_DMSOvsTg_Tg <- rbind(df_12hr_DMSOvsTg_Tg_F, df_12hr_DMSOvsTg_Tg_M)

#remove duplicated rows
df_12hr_DMSOvsTg_DMSO <- df_12hr_DMSOvsTg_DMSO[!duplicated(df_12hr_DMSOvsTg_DMSO$Pathway),]
df_12hr_DMSOvsTg_Tg <- df_12hr_DMSOvsTg_Tg[!duplicated(df_12hr_DMSOvsTg_Tg$Pathway),]

#order by gene ratio
df_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_12hr_DMSOvsTg_DMSO, desc(order_Gene.ratio))
df_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_12hr_DMSOvsTg_Tg, order_Gene.ratio)

#order by -log10 FDR
df_12hr_DMSOvsTg_DMSO <- dplyr::arrange(df_12hr_DMSOvsTg_DMSO, desc(order_transf.log.10.FDR_max))
df_12hr_DMSOvsTg_Tg <- dplyr::arrange(df_12hr_DMSOvsTg_Tg, order_transf.log.10.FDR_max)

#flip -log10 FDR for DMSO
df_12hr_DMSOvsTg_DMSO$transf.log.10.FDR_F <- -df_12hr_DMSOvsTg_DMSO$transf.log.10.FDR_F
df_12hr_DMSOvsTg_DMSO$transf.log.10.FDR_M <- -df_12hr_DMSOvsTg_DMSO$transf.log.10.FDR_M

#combine treatments
df_12hr_DMSOvsTg <- rbind(df_12hr_DMSOvsTg_DMSO, df_12hr_DMSOvsTg_Tg)

#remove rows for non-significance for both sexes
df_12hr_DMSOvsTg <- df_12hr_DMSOvsTg[!(df_12hr_DMSOvsTg$significance_F == "> 0.05" & df_12hr_DMSOvsTg$significance_M == "> 0.05"),] %>% na.omit()

#set pathway factor level order
df_12hr_DMSOvsTg$Pathway <- factor(df_12hr_DMSOvsTg$Pathway, levels = df_12hr_DMSOvsTg$Pathway)

#separate based on treatment
df_12hr_DMSOvsTg_plotted_Tg <- df_12hr_DMSOvsTg[df_12hr_DMSOvsTg$Treatment == "Tg",]
df_12hr_DMSOvsTg_plotted_DMSO <- df_12hr_DMSOvsTg[df_12hr_DMSOvsTg$Treatment == "DMSO",]

#select columns for each sex and treatment
df_12hr_DMSOvsTg_plotted_Tg_F <- df_12hr_DMSOvsTg_plotted_Tg[, c("Pathway", "Treatment", "Entities.found_F", "Entities.total_F", "FDR_F", "Sex_F", "Gene.ratio_F", "transf.log.10.FDR_F", "significance_F", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_12hr_DMSOvsTg_plotted_Tg_M <- df_12hr_DMSOvsTg_plotted_Tg[, c("Pathway", "Treatment", "Entities.found_M", "Entities.total_M", "FDR_M", "Sex_M", "Gene.ratio_M", "transf.log.10.FDR_M", "significance_M", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_12hr_DMSOvsTg_plotted_DMSO_F <- df_12hr_DMSOvsTg_plotted_DMSO[, c("Pathway", "Treatment", "Entities.found_F", "Entities.total_F", "FDR_F", "Sex_F", "Gene.ratio_F", "transf.log.10.FDR_F", "significance_F", "order_transf.log.10.FDR_max", "order_Gene.ratio")]
df_12hr_DMSOvsTg_plotted_DMSO_M <- df_12hr_DMSOvsTg_plotted_DMSO[, c("Pathway", "Treatment", "Entities.found_M", "Entities.total_M", "FDR_M", "Sex_M", "Gene.ratio_M", "transf.log.10.FDR_M", "significance_M", "order_transf.log.10.FDR_max", "order_Gene.ratio")]

#rename columns
colnames(df_12hr_DMSOvsTg_plotted_Tg_F) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_12hr_DMSOvsTg_plotted_Tg_M) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_12hr_DMSOvsTg_plotted_DMSO_F) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")
colnames(df_12hr_DMSOvsTg_plotted_DMSO_M) <- c("Pathway", "Treatment", "Entities.found", "Entities.total", "FDR", "Sex", "Gene.ratio", "transf.log.10.FDR", "significance", "order_transf.log.10.FDR_max", "order_Gene.ratio")

#combine sexes for each treatment
df_12hr_DMSOvsTg_plotted_Tg <- rbind(df_12hr_DMSOvsTg_plotted_Tg_F, df_12hr_DMSOvsTg_plotted_Tg_M)
df_12hr_DMSOvsTg_plotted_DMSO <- rbind(df_12hr_DMSOvsTg_plotted_DMSO_F, df_12hr_DMSOvsTg_plotted_DMSO_M)

#create dot colour column
df_12hr_DMSOvsTg_plotted_Tg$colour_index <- paste(df_12hr_DMSOvsTg_plotted_Tg$Sex, df_12hr_DMSOvsTg_plotted_Tg$significance)
df_12hr_DMSOvsTg_plotted_DMSO$colour_index <- paste(df_12hr_DMSOvsTg_plotted_DMSO$Sex, df_12hr_DMSOvsTg_plotted_DMSO$significance)

#set significance factor level order
df_12hr_DMSOvsTg_plotted_Tg$significance <- factor(df_12hr_DMSOvsTg_plotted_Tg$significance, levels = sort(unique(df_12hr_DMSOvsTg_plotted_Tg$significance))[c(1, 3, 2)])

#set x-axis limits, y-axis annotations, and label position
top_label_height_multiplier <- 1.0305

left_limit <- -10
right_limit <- 15

num_row <- nrow(df_12hr_DMSOvsTg)
num_row_bottom <- nrow(df_12hr_DMSOvsTg_DMSO[!(df_12hr_DMSOvsTg_DMSO$significance_F == "> 0.05" & df_12hr_DMSOvsTg_DMSO$significance_M == "> 0.05"),] %>% na.omit())

#create plot
ggplot2::ggplot(df_12hr_DMSOvsTg, ggplot2::aes(y = Pathway)) +
  
  ggplot2::annotate(geom = "segment", x = 0, xend = 0, y = 0.5, yend = num_row + 0.5, size = 0.5) +
  ggplot2::annotate(geom = "segment", x = left_limit, xend = right_limit, y = num_row_bottom + 0.5, yend = num_row_bottom + 0.5, size = 0.5) +
  
  ggplot2::annotate(geom = "segment", x = left_limit, xend = 0, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#DDDDDD", size = 8) +
  ggplot2::annotate(geom = "segment", x = 0, xend = right_limit, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#AAAAAA", size = 8) +
  ggplot2::annotate(geom = "segment", x = -(right_limit - left_limit) * 0.0075, xend = (right_limit - left_limit) * 0.0075, y = num_row * top_label_height_multiplier + 0.5, yend = num_row * top_label_height_multiplier + 0.5, colour = "#FFFFFF", size = 8) +
  ggplot2::annotate(geom = "text", x = left_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("DMSO"), size = 3) +
  ggplot2::annotate(geom = "text", x = right_limit/2, y = num_row * top_label_height_multiplier + 0.5, label = c("Tg")) +
  ggplot2::annotate(geom = "text", x = left_limit - (right_limit - left_limit) * 0.035, y = num_row * top_label_height_multiplier + 0.5, label = c("Treatment"), hjust = 1) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(order_transf.log.10.FDR_max)), 
                      alpha = 0, 
                      data = df_12hr_DMSOvsTg) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR), 
                                   colour = colour_index, 
                                   fill = Sex, 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1, 
                      position = ggstance::position_dodgev(height = 0.2), 
                      data = df_12hr_DMSOvsTg_plotted_Tg) +
  
  ggplot2::geom_point(ggplot2::aes(x = as.numeric(transf.log.10.FDR), 
                                   colour = colour_index, 
                                   fill = Sex, 
                                   shape = significance, 
                                   size = Gene.ratio), 
                      alpha = 0.8, 
                      stroke = 1, 
                      position = ggstance::position_dodgev(height = 0.2), 
                      data = df_12hr_DMSOvsTg_plotted_DMSO) +
  
  ggplot2::coord_cartesian(xlim = c(left_limit, right_limit), ylim = c(0.5, num_row + 0.5), expand = FALSE, clip = "off") +
  
  ggplot2::scale_size_continuous(range = c(0.5, 5)) +
  
  ggplot2::scale_fill_manual("Sex",
                             values = c("#E7B91E", "#84B429"),
                             labels = c("Female", "Male")) +
  
  ggplot2::scale_colour_manual("Significance", values = c("#E7B91E", "grey", "#E7B91E", "#84B429", "grey", "#84B429")) +
  
  ggplot2::scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 45), drop = FALSE) +
  
  ggplot2::scale_alpha(guide = "none") +
  
  ggplot2::labs(x = bquote(-log[10]*"(FDR)"), 
                colour = "Sex", 
                fill = "Sex", 
                size = "Gene Ratio", 
                shape = "Significance"
  ) +
  
  ggplot2::scale_shape_manual(values = c(21, 1, 1)) +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(reverse = F, 
                                 override.aes = list(size = 5, 
                                                     colour = c("#E7B91E", "#84B429")), 
                                 order = 1), 
    colour = "none", 
    size = ggplot2::guide_legend(reverse = T, 
                                 override.aes = list(colour = "black"),
                                 order = 2), 
    shape = ggplot2::guide_legend(reverse = F, 
                                  override.aes = list(colour = c("black", "black", "grey"), 
                                                      size = 4, 
                                                      fill = c("black", "white", "white")),
                                  order = 3)) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 0.5, unit = "cm")),
                 plot.margin = ggplot2::margin(t = 2, r = 1, b = 1, l = 1, "cm"),
                 axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(size = 8),
                 axis.text = ggplot2::element_text(colour = "black"),
                 panel.grid.minor = ggplot2::element_blank(),
                 legend.title.align = 0.5,
                 legend.direction = "vertical",
                 legend.box.just = "center",
                 legend.key.width = grid::unit(0.75, "cm"),
                 legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0, l = 0.5, "cm"),
                 legend.box.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, "cm"),
                 aspect.ratio = (nrow(df_12hr_DMSOvsTg) + 1)/5
  )


### Figure 5 - Supplement 3 ###

## Fig 5 S3A ##

#clear environment
rm(list = ls())

#import file
data_interaction_F <- read.csv(".\\Data\\deltaFC_overtime_low genes filtered_DeSeq2 interaction by_time_sex_F.tsv", stringsAsFactors = F, check.names = F, sep = "\t")

#remove rows without padj
data_interaction_F <- subset(data_interaction_F, !is.na(interaction_padj_Female))

#categorize interaction change
data_interaction_F <- data_interaction_F %>%
  dplyr::mutate(
    category = dplyr::case_when(
      interaction_padj_Female > 0.05 & FC_Female_6hr_TGvsDMSO > 0 & FC_Female_12hr_TGvsDMSO > 0 ~ 1,
      interaction_padj_Female > 0.05 & FC_Female_6hr_TGvsDMSO < 0 & FC_Female_12hr_TGvsDMSO < 0 ~ 2,
      interaction_padj_Female < 0.05 & deltaFC_Female < 0 & FC_Female_6hr_TGvsDMSO > 0 & FC_Female_12hr_TGvsDMSO > 0 ~ 3,
      interaction_padj_Female < 0.05 & deltaFC_Female > 0 & FC_Female_6hr_TGvsDMSO > 0 & FC_Female_12hr_TGvsDMSO > 0 ~ 4,
      interaction_padj_Female < 0.05 & deltaFC_Female > 0 & FC_Female_6hr_TGvsDMSO < 0 & FC_Female_12hr_TGvsDMSO < 0 ~ 5,
      interaction_padj_Female < 0.05 & deltaFC_Female < 0 & FC_Female_6hr_TGvsDMSO < 0 & FC_Female_12hr_TGvsDMSO < 0 ~ 6,
      FC_Female_6hr_TGvsDMSO > 0 & FC_Female_12hr_TGvsDMSO < 0 ~ 7,
      FC_Female_6hr_TGvsDMSO < 0 & FC_Female_12hr_TGvsDMSO > 0 ~ 8,
    )
  )

#import beta cell gene list
df_beta_cell_genes <- read.csv(".\\Data\\beta cell genes.csv", stringsAsFactors = F, check.names = F, sep = ",")

#categorize gene type
data_interaction_F$GeneType <-  ifelse(tolower(data_interaction_F$`Gene Name`) %in% tolower(df_beta_cell_genes$genes), "Beta Cell Genes", "Other")

#categorize interaction significance
data_interaction_F$Sig_Interaction <-  ifelse(data_interaction_F$interaction_padj_Female < 0.05 & !is.na(data_interaction_F$interaction_padj_Female), "Significant Interaction", "Non-Significant Interaction")

#categorize fold change significance
data_interaction_F$Sig_FC_6hr <-  ifelse(data_interaction_F$padj_Female_6hr_TGvsDMSO < 0.05 & !is.na(data_interaction_F$padj_Female_6hr_TGvsDMSO), "Significant Fold Change", "Non-Significant Fold Change")
data_interaction_F$Sig_FC_12hr <-  ifelse(data_interaction_F$padj_Female_12hr_TGvsDMSO < 0.05 & !is.na(data_interaction_F$padj_Female_12hr_TGvsDMSO), "Significant Fold Change", "Non-Significant Fold Change")

#get -log10 of padj
data_interaction_F$interaction_padj_Female <- -log10(data_interaction_F$interaction_padj_Female)

#change padj column name
colnames(data_interaction_F)[which(colnames(data_interaction_F) == "interaction_padj_Female")] <- "interaction_transf.log.10.padj_Female"

#flip -log10 padj when deltaFC less than 0
data_interaction_F$interaction_transf.log.10.padj_Female[data_interaction_F$deltaFC_Female < 0] <- -data_interaction_F$interaction_transf.log.10.padj_Female[data_interaction_F$deltaFC_Female < 0]

#change fold change column names
colnames(data_interaction_F)[which(colnames(data_interaction_F) %in% c("FC_Female_6hr_TGvsDMSO", "FC_Female_12hr_TGvsDMSO"))] <- c("at 6hr (DMSO vs Tg)", "at 12hr (DMSO vs Tg)")

#lengthens data with categorized data points by type 
data_interaction_F_plotted <- data_interaction_F %>%
  tidyr::pivot_longer(
    cols = c(deltaFC_Female, 
             `at 6hr (DMSO vs Tg)`, 
             `at 12hr (DMSO vs Tg)`
    ),
    names_to = "Variable", 
    values_to = "Data")

#add significance annotation to delta FC rows
data_interaction_F_plotted$Variable <- ifelse(data_interaction_F_plotted$Variable == "deltaFC_Female", paste(data_interaction_F_plotted$Variable, data_interaction_F_plotted$Sig_Interaction), data_interaction_F_plotted$Variable)

#order by -log10 padj
data_interaction_F <- data_interaction_F[with(data_interaction_F, order(interaction_transf.log.10.padj_Female)),]

#set gene label positions
data_interaction_F$Label_Position[data_interaction_F$GeneType == "Beta Cell Genes"] <- c(1, 1, 2, 2, 1, 1.1, 1, 1, 2, 1, 2, 1.1, 1, 2.1, 2, 1, 2, 2, 1, 4, 3, 4)

#create delta FC table
data_interaction_F_plotted_deltaFC_Female <- data.frame(`Gene Name` = data_interaction_F_plotted$`Gene Name`[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Non-Significant Interaction", "deltaFC_Female Significant Interaction")],
                                                        deltaFC_Female = data_interaction_F_plotted$Data[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Non-Significant Interaction", "deltaFC_Female Significant Interaction")],
                                                        check.names = FALSE)

#annotate time point with significance per gene
data_interaction_F_plotted <- data_interaction_F_plotted %>%
  dplyr::mutate(
    Sig_FC = dplyr::case_when(
      Variable == "at 6hr (DMSO vs Tg)"  ~ paste(Sig_FC_6hr, Variable),
      Variable == "at 12hr (DMSO vs Tg)" ~ paste(Sig_FC_12hr, Variable)
    )
  )

#set delta FC point colours
data_interaction_F_plotted$delta_FC_colour <- ifelse(data_interaction_F_plotted$Sig_Interaction == "Significant Interaction", "blue", "grey70")

#merge delta FC colours with data
data_interaction_F_plotted <- merge(data_interaction_F_plotted, data_interaction_F_plotted_deltaFC_Female, by = "Gene Name", all.x = TRUE)

#order by -log10 padj
data_interaction_F_plotted <- data_interaction_F_plotted[with(data_interaction_F_plotted, order(interaction_transf.log.10.padj_Female)),]

#create data frame with only significant fold changes
data_interaction_F_plotted_Sig_FC <- data_interaction_F_plotted[!is.na(data_interaction_F_plotted$Sig_FC),]

#set fold change significance factor level order
data_interaction_F_plotted_Sig_FC$Sig_FC <- factor(data_interaction_F_plotted_Sig_FC$Sig_FC, levels = sort(unique(data_interaction_F_plotted_Sig_FC$Sig_FC))[c(2, 4, 1, 3)])

#number ordered genes with index 
data_interaction_F <- within(data_interaction_F, Index <- match(`Gene Name`, unique(`Gene Name`)))
data_interaction_F_plotted <- within(data_interaction_F_plotted, Index <- match(`Gene Name`, unique(`Gene Name`)))
data_interaction_F_plotted_Sig_FC <- within(data_interaction_F_plotted_Sig_FC, Index <- match(`Gene Name`, unique(`Gene Name`)))

#set label text hjust
bottom_left_hjust_F <- c(0, 0, 0.6, 0, 0, 0.4, 0, 0, 0)
bottom_left_hjust_F_1.1 <- c(0, 0)
top_left_hjust_F <- c(0.6, 0.3, 0, 1, 0, 0, 0)
top_left_hjust_F_2.1 <- c(0)
bottom_right_hjust_F <- 0
top_right_hjust_F <- 0

#create plot
ggplot2::ggplot(data_interaction_F_plotted_Sig_FC, ggplot2::aes(x = Index)) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, shape = Sig_FC, colour = Sig_FC),
             alpha = 0,
             data = data_interaction_F_plotted_Sig_FC) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, fill = Sig_Interaction),
             alpha = 0.05,
             shape = 16,
             size = 2, 
             colour = data_interaction_F_plotted$delta_FC_colour[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Significant Interaction")],
             data = data_interaction_F_plotted[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Significant Interaction"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, fill = Sig_Interaction),
             alpha = 0.1,
             shape = 16,
             size = 2, 
             colour = data_interaction_F_plotted$delta_FC_colour[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Non-Significant Interaction")],
             data = data_interaction_F_plotted[data_interaction_F_plotted$Variable %in% c("deltaFC_Female Non-Significant Interaction"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC),
             alpha = 0.1,
             size = 1, 
             data = data_interaction_F_plotted[data_interaction_F_plotted$Sig_FC %in% c("Non-Significant Fold Change at 6hr (DMSO vs Tg)", "Non-Significant Fold Change at 12hr (DMSO vs Tg)"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC, alpha = GeneType, shape = Sig_FC),
             stroke = 1.5,
             size = 4, 
             data = data_interaction_F_plotted[data_interaction_F_plotted$GeneType == "Beta Cell Genes" & data_interaction_F_plotted$Variable == "at 6hr (DMSO vs Tg)" & !is.na(data_interaction_F_plotted$Sig_FC),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC, alpha = GeneType, shape = Sig_FC),
             stroke = 1.5,
             size = 4, 
             data = data_interaction_F_plotted[data_interaction_F_plotted$GeneType == "Beta Cell Genes" & data_interaction_F_plotted$Variable == "at 12hr (DMSO vs Tg)" & !is.na(data_interaction_F_plotted$Sig_FC),]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)`, 
                                     yend = `at 12hr (DMSO vs Tg)`,
                                     x = Index,
                                     xend = Index,
                                     linetype = Sig_Interaction),
                        colour = "black",
                        size = 1,
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Sig_Interaction == "Significant Interaction",]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` - 0.15,
                                     yend = `at 12hr (DMSO vs Tg)` - 1,
                                     x = Index + 32,
                                     xend = Index + 200),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` - 0.15,
                                     yend = `at 12hr (DMSO vs Tg)` - 2,
                                     x = Index + 32,
                                     xend = Index + 400),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 1.1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` + 0.15,
                                     yend = `at 6hr (DMSO vs Tg)` + 1,
                                     x = Index + 32,
                                     xend = Index + 200),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 2,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` + 0.15,
                                     yend = `at 6hr (DMSO vs Tg)` + 2,
                                     x = Index + 32,
                                     xend = Index + 400),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 2.1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` - 0.15,
                                     yend = `at 6hr (DMSO vs Tg)` - 1,
                                     x = Index + 32,
                                     xend = Index + 200),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 3,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` + 0.15,
                                     yend = `at 12hr (DMSO vs Tg)` + 1,
                                     x = Index + 32,
                                     xend = Index + 200),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 4,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` - 1,
                                  x = Index + 200),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_left_hjust_F,
                     vjust = 1.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` - 2,
                                  x = Index + 400),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_left_hjust_F_1.1,
                     vjust = 1.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 1.1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` + 1,
                                  x = Index + 200),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_left_hjust_F,
                     vjust = -0.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 2,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` + 2,
                                  x = Index + 400),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_left_hjust_F_2.1,
                     vjust = -0.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 2.1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` - 1,
                                  x = Index + 200),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_right_hjust_F,
                     vjust = 1.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 3,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` + 1,
                                  x = Index + 200),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_right_hjust_F,
                     vjust = -0.2,
                     data = data_interaction_F[data_interaction_F$GeneType == "Beta Cell Genes" & data_interaction_F$Label_Position == 4,]) +
  
  ggplot2::scale_color_manual(values = c("#FFA505", "#FFA505", "#D676FF", "#D676FF")) +
  
  ggplot2::scale_shape_manual(values = c(21, 16, 21, 16)) +
  
  ggplot2::scale_alpha_manual(values = c(1, 0.01)) +
  
  ggplot2::scale_fill_manual(values = c("grey70", "blue")) +
  
  ggplot2::coord_cartesian(xlim = c(-100.5, nrow(data_interaction_F) + 101.5),
                           ylim = c(-6.5, 6.5),
                           expand = FALSE,
                           clip = "on") +
  
  ggplot2::scale_linetype_manual(values = c("solid")) +
  
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
  
  ggplot2::geom_vline(xintercept = data_interaction_F$Index[length(data_interaction_F$Index[data_interaction_F$Sig_Interaction == "Significant Interaction" & data_interaction_F$interaction_transf.log.10.padj_Female < 0])], color = "grey70", size = 1, linetype = "dashed") +
  
  ggplot2::geom_vline(xintercept = data_interaction_F$Index[length(data_interaction_F$Index[(data_interaction_F$Sig_Interaction == "Significant Interaction" & data_interaction_F$interaction_transf.log.10.padj_Female < 0) | data_interaction_F$Sig_Interaction == "Non-Significant Interaction"])], color = "grey70", size = 1, linetype = "dashed") +
  
  ggplot2::labs(y = bquote(log[2]*"(Fold Change) (DMSO vs Thapsigargin)"), 
                alpha = "Cell Type", 
                shape = "Fold Change Significance", 
                colour = "Fold Change Significance", 
                linetype = "Interaction Significance",
                fill = "Interaction Points") +
  
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 4, fill = "white")),
    fill = ggplot2::guide_legend(override.aes = list(size = 4, shape = 21, alpha = 1, colour = c("grey70", "blue"))),
    alpha = "none",
    size = "none"
  ) +
  
  ggplot2::theme_bw() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0.25,
      r = 0.25,
      b = 0,
      l = 0.25,
      "cm"
    ),
    axis.ticks.length.y = grid::unit(2, "mm"),
    axis.ticks.y = ggplot2::element_line(size = 1, colour = "black"),
    axis.ticks.length.x.bottom = grid::unit(2, "mm"),
    axis.ticks.x.bottom = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(size = 2, colour = "black"),
    axis.title.y = ggplot2::element_text(size = 20),
    axis.line.x.bottom = ggplot2::element_blank(),
    axis.line.y.left = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(colour = "black", size = 20),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 20),
    legend.text = ggplot2::element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "bottom",
    legend.key.height = grid::unit(1, "cm"),
    legend.spacing.x = grid::unit(0.2, "cm"),
    legend.margin = ggplot2::margin(0, 0.63, 0, 0.63, unit = "cm"),
    legend.background = ggplot2::element_rect(fill = alpha("white", 0)),
    aspect.ratio = 1.5 / 2)


## Fig 5 S3B ##

#clear environment
rm(list = ls())

#import file
data_interaction_M <- read.csv(".\\Data\\deltaFC_overtime_low genes filtered_DeSeq2 interaction by_time_sex_M.tsv", stringsAsFactors = F, check.names = F, sep = "\t")

#remove rows without padj
data_interaction_M <- subset(data_interaction_M, !is.na(interaction_padj_Male))

#categorize interaction change
data_interaction_M <- data_interaction_M %>%
  dplyr::mutate(
    category = dplyr::case_when(
      interaction_padj_Male > 0.05 & FC_Male_6hr_TGvsDMSO > 0 & FC_Male_12hr_TGvsDMSO > 0 ~ 1,
      interaction_padj_Male > 0.05 & FC_Male_6hr_TGvsDMSO < 0 & FC_Male_12hr_TGvsDMSO < 0 ~ 2,
      interaction_padj_Male < 0.05 & deltaFC_Male < 0 & FC_Male_6hr_TGvsDMSO > 0 & FC_Male_12hr_TGvsDMSO > 0 ~ 3,
      interaction_padj_Male < 0.05 & deltaFC_Male > 0 & FC_Male_6hr_TGvsDMSO > 0 & FC_Male_12hr_TGvsDMSO > 0 ~ 4,
      interaction_padj_Male < 0.05 & deltaFC_Male > 0 & FC_Male_6hr_TGvsDMSO < 0 & FC_Male_12hr_TGvsDMSO < 0 ~ 5,
      interaction_padj_Male < 0.05 & deltaFC_Male < 0 & FC_Male_6hr_TGvsDMSO < 0 & FC_Male_12hr_TGvsDMSO < 0 ~ 6,
      FC_Male_6hr_TGvsDMSO > 0 & FC_Male_12hr_TGvsDMSO < 0 ~ 7,
      FC_Male_6hr_TGvsDMSO < 0 & FC_Male_12hr_TGvsDMSO > 0 ~ 8,
    )
  )

#import beta cell gene list
df_beta_cell_genes <- read.csv(".\\Data\\beta cell genes.csv", stringsAsFactors = F, check.names = F, sep = ",")

#categorize gene type
data_interaction_M$GeneType <-  ifelse(tolower(data_interaction_M$`Gene Name`) %in% tolower(df_beta_cell_genes$genes), "Beta Cell Genes", "Other")

#categorize interaction significance
data_interaction_M$Sig_Interaction <-  ifelse(data_interaction_M$interaction_padj_Male < 0.05 & !is.na(data_interaction_M$interaction_padj_Male), "Significant Interaction", "Non-Significant Interaction")

#categorize fold change significance
data_interaction_M$Sig_FC_6hr <-  ifelse(data_interaction_M$padj_Male_6hr_TGvsDMSO < 0.05 & !is.na(data_interaction_M$padj_Male_6hr_TGvsDMSO), "Significant Fold Change", "Non-Significant Fold Change")
data_interaction_M$Sig_FC_12hr <-  ifelse(data_interaction_M$padj_Male_12hr_TGvsDMSO < 0.05 & !is.na(data_interaction_M$padj_Male_12hr_TGvsDMSO), "Significant Fold Change", "Non-Significant Fold Change")

#get -log10 of padj
data_interaction_M$interaction_padj_Male <- -log10(data_interaction_M$interaction_padj_Male)

#change padj column name
colnames(data_interaction_M)[which(colnames(data_interaction_M) == "interaction_padj_Male")] <- "interaction_transf.log.10.padj_Male"

#flip -log10 padj when deltaFC less than 0
data_interaction_M$interaction_transf.log.10.padj_Male[data_interaction_M$deltaFC_Male < 0] <- -data_interaction_M$interaction_transf.log.10.padj_Male[data_interaction_M$deltaFC_Male < 0]

#change fold change column names
colnames(data_interaction_M)[which(colnames(data_interaction_M) %in% c("FC_Male_6hr_TGvsDMSO", "FC_Male_12hr_TGvsDMSO"))] <- c("at 6hr (DMSO vs Tg)", "at 12hr (DMSO vs Tg)")

#lengthens data with categorized data points by type 
data_interaction_M_plotted <- data_interaction_M %>%
  tidyr::pivot_longer(
    cols = c(deltaFC_Male, 
             `at 6hr (DMSO vs Tg)`, 
             `at 12hr (DMSO vs Tg)`
    ),
    names_to = "Variable", 
    values_to = "Data")

#add significance annotation to delta FC rows
data_interaction_M_plotted$Variable <- ifelse(data_interaction_M_plotted$Variable == "deltaFC_Male", paste(data_interaction_M_plotted$Variable, data_interaction_M_plotted$Sig_Interaction), data_interaction_M_plotted$Variable)

#order by -log10 padj
data_interaction_M <- data_interaction_M[with(data_interaction_M, order(interaction_transf.log.10.padj_Male)),]

#set gene label positions
data_interaction_M$Label_Position[data_interaction_M$GeneType == "Beta Cell Genes"] <- c(1, 2, 2, 1.1, 1.1, 2, 1, 1, 1, 2, 1, 2.1, 2, 2, 1, 1, 2, 1, 4, 3, 4, 3)

#create delta FC table
data_interaction_M_plotted_deltaFC_Male <- data.frame(`Gene Name` = data_interaction_M_plotted$`Gene Name`[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Non-Significant Interaction", "deltaFC_Male Significant Interaction")],
                                                      deltaFC_Male = data_interaction_M_plotted$Data[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Non-Significant Interaction", "deltaFC_Male Significant Interaction")],
                                                      check.names = FALSE)

#annotate time point with significance per gene
data_interaction_M_plotted <- data_interaction_M_plotted %>%
  dplyr::mutate(
    Sig_FC = dplyr::case_when(
      Variable == "at 6hr (DMSO vs Tg)"  ~ paste(Sig_FC_6hr, Variable),
      Variable == "at 12hr (DMSO vs Tg)" ~ paste(Sig_FC_12hr, Variable)
    )
  )

#set delta FC point colours
data_interaction_M_plotted$delta_FC_colour <- ifelse(data_interaction_M_plotted$Sig_Interaction == "Significant Interaction", "blue", "grey70")

#merge delta FC colours with data
data_interaction_M_plotted <- merge(data_interaction_M_plotted, data_interaction_M_plotted_deltaFC_Male, by = "Gene Name", all.x = TRUE)

#order by -log10 padj
data_interaction_M_plotted <- data_interaction_M_plotted[with(data_interaction_M_plotted, order(interaction_transf.log.10.padj_Male)),]

#create data frame with only significant fold changes
data_interaction_M_plotted_Sig_FC <- data_interaction_M_plotted[!is.na(data_interaction_M_plotted$Sig_FC),]

#set fold change significance factor level order
data_interaction_M_plotted_Sig_FC$Sig_FC <- factor(data_interaction_M_plotted_Sig_FC$Sig_FC, levels = sort(unique(data_interaction_M_plotted_Sig_FC$Sig_FC))[c(2, 4, 1, 3)])

#number ordered genes with index 
data_interaction_M <- within(data_interaction_M, Index <- match(`Gene Name`, unique(`Gene Name`)))
data_interaction_M_plotted <- within(data_interaction_M_plotted, Index <- match(`Gene Name`, unique(`Gene Name`)))
data_interaction_M_plotted_Sig_FC <- within(data_interaction_M_plotted_Sig_FC, Index <- match(`Gene Name`, unique(`Gene Name`)))

#set label text hjust
bottom_left_hjust_M <- c(0, 0, 1, 0, 0, 0, 0.2, 0)
bottom_left_hjust_M_1.1 <- c(1, 0.2)
top_left_hjust_M <- c(0, 1.1, 0.8, 0.7, 0, 0, 0)
top_left_hjust_M_2.1 <- 0
bottom_right_hjust_M <- 0
top_right_hjust_M <- 0

#create plot
ggplot(data_interaction_M_plotted_Sig_FC, ggplot2::aes(x = Index)) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, shape = Sig_FC, colour = Sig_FC),
                      alpha = 0,
                      data = data_interaction_M_plotted_Sig_FC) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, fill = Sig_Interaction),
                      alpha = 0.05,
                      shape = 16,
                      size = 2, 
                      colour = data_interaction_M_plotted$delta_FC_colour[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Significant Interaction")],
                      data = data_interaction_M_plotted[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Significant Interaction"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, fill = Sig_Interaction),
                      alpha = 0.1,
                      shape = 16,
                      size = 2, 
                      colour = data_interaction_M_plotted$delta_FC_colour[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Non-Significant Interaction")],
                      data = data_interaction_M_plotted[data_interaction_M_plotted$Variable %in% c("deltaFC_Male Non-Significant Interaction"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC),
                      alpha = 0.1,
                      size = 1, 
                      data = data_interaction_M_plotted[data_interaction_M_plotted$Sig_FC %in% c("Non-Significant Fold Change at 6hr (DMSO vs Tg)", "Non-Significant Fold Change at 12hr (DMSO vs Tg)"),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC, alpha = GeneType, shape = Sig_FC),
                      stroke = 1.5,
                      size = 4, 
                      data = data_interaction_M_plotted[data_interaction_M_plotted$GeneType == "Beta Cell Genes" & data_interaction_M_plotted$Variable == "at 6hr (DMSO vs Tg)" & !is.na(data_interaction_M_plotted$Sig_FC),]) +
  
  ggplot2::geom_point(ggplot2::aes(y = Data, colour = Sig_FC, alpha = GeneType, shape = Sig_FC),
                      stroke = 1.5,
                      size = 4, 
                      data = data_interaction_M_plotted[data_interaction_M_plotted$GeneType == "Beta Cell Genes" & data_interaction_M_plotted$Variable == "at 12hr (DMSO vs Tg)" & !is.na(data_interaction_M_plotted$Sig_FC),]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)`, 
                                     yend = `at 12hr (DMSO vs Tg)`,
                                     x = Index,
                                     xend = Index,
                                     linetype = Sig_Interaction),
                        colour = "black",
                        size = 1,
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Sig_Interaction == "Significant Interaction",]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` - 0.12,
                                     yend = `at 12hr (DMSO vs Tg)` - 1,
                                     x = Index + 50,
                                     xend = Index + 490),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` - 0.12,
                                     yend = `at 12hr (DMSO vs Tg)` - 1,
                                     x = Index - 50,
                                     xend = Index - 490),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 1.1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` + 0.12,
                                     yend = `at 6hr (DMSO vs Tg)` + 1,
                                     x = Index + 50,
                                     xend = Index + 490),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 2,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` + 0.12,
                                     yend = `at 6hr (DMSO vs Tg)` + 1.25,
                                     x = Index + 50,
                                     xend = Index + 613),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 2.1,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 6hr (DMSO vs Tg)` - 0.12,
                                     yend = `at 6hr (DMSO vs Tg)` - 1,
                                     x = Index + 50,
                                     xend = Index + 490),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 3,]) +
  
  ggplot2::geom_segment(ggplot2::aes(y = `at 12hr (DMSO vs Tg)` + 0.12,
                                     yend = `at 12hr (DMSO vs Tg)` + 1,
                                     x = Index + 50,
                                     xend = Index + 490),
                        colour = "black",
                        alpha = 1,
                        data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 4,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` - 1,
                                  x = Index + 490),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_left_hjust_M,
                     vjust = 1.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` - 1,
                                  x = Index - 490),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_left_hjust_M_1.1,
                     vjust = 1.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 1.1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` + 1,
                                  x = Index + 490),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_left_hjust_M,
                     vjust = -0.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 2,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` + 1.25,
                                  x = Index + 613),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_left_hjust_M_2.1,
                     vjust = -0.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 2.1,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 6hr (DMSO vs Tg)` - 1,
                                  x = Index + 490),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = bottom_right_hjust_M,
                     vjust = 1.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 3,]) +
  
  ggplot2::geom_text(ggplot2::aes(label = `Gene Name`,
                                  y = `at 12hr (DMSO vs Tg)` + 1,
                                  x = Index + 490),
                     colour = "black",
                     alpha = 1,
                     size = 6,
                     hjust = top_right_hjust_M,
                     vjust = -0.2,
                     data = data_interaction_M[data_interaction_M$GeneType == "Beta Cell Genes" & data_interaction_M$Label_Position == 4,]) +
  
  ggplot2::scale_color_manual(values = c("#75CF4A", "#75CF4A", "#469FD2", "#469FD2")) +
  
  ggplot2::scale_shape_manual(values = c(21, 16, 21, 16)) +
  
  ggplot2::scale_alpha_manual(values = c(1, 0.01)) +
  
  ggplot2::scale_fill_manual(values = c("grey70", "blue")) +
  
  ggplot2::coord_cartesian(xlim = c(-100.5, nrow(data_interaction_M) + 101.5),
                           ylim = c(-6.5, 6.5),
                           expand = FALSE,
                           clip = "on") +
  
  ggplot2::scale_linetype_manual(values = c("solid")) +
  
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 1, linetype = "dashed") +
  
  ggplot2::geom_vline(xintercept = data_interaction_M$Index[length(data_interaction_M$Index[data_interaction_M$Sig_Interaction == "Significant Interaction" & data_interaction_M$interaction_transf.log.10.padj_Male < 0])], color = "grey70", size = 1, linetype = "dashed") +
  
  ggplot2::geom_vline(xintercept = data_interaction_M$Index[length(data_interaction_M$Index[(data_interaction_M$Sig_Interaction == "Significant Interaction" & data_interaction_M$interaction_transf.log.10.padj_Male < 0) | data_interaction_M$Sig_Interaction == "Non-Significant Interaction"])], color = "grey70", size = 1, linetype = "dashed") +
  
  ggplot2::labs(y = bquote(log[2]*"(Fold Change) (DMSO vs Thapsigargin)"), 
                alpha = "Cell Type", 
                shape = "Fold Change Significance", 
                colour = "Fold Change Significance", 
                linetype = "Interaction Significance",
                fill = "Interaction Points") +
  
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 4, fill = "white")),
    fill = ggplot2::guide_legend(override.aes = list(size = 4, shape = 21, alpha = 1, colour = c("grey70", "blue"))),
    alpha = "none",
    size = "none"
  ) +
  
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0.25,
      r = 0.25,
      b = 0,
      l = 0.25,
      "cm"
    ),
    axis.ticks.length.y = grid::unit(2, "mm"),
    axis.ticks.y = ggplot2::element_line(size = 1, colour = "black"),
    axis.ticks.length.x.bottom = grid::unit(2, "mm"),
    axis.ticks.x.bottom = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(size = 2, colour = "black"),
    axis.title.y = ggplot2::element_text(size = 20),
    axis.line.x.bottom = ggplot2::element_blank(),
    axis.line.y.left = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(colour = "black", size = 20),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(size = 20),
    legend.text = ggplot2::element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "bottom",
    legend.key.height = grid::unit(1, "cm"),
    legend.spacing.x = grid::unit(0.2, "cm"),
    legend.margin = ggplot2::margin(0, 0.63, 0, 0.63, unit = "cm"),
    legend.background = ggplot2::element_rect(fill = alpha("white", 0)),
    aspect.ratio = 1.5 / 2)


