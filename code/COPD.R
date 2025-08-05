rm(list=ls())
graphics.off()

# load packages -----------------------------------------------------------

library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(microbiome)
library(readr)
library(vegan)
library(ape)
library(readr)
library(tidyr)
library(circlize)
library(igraph)
library(corrplot)
library(correlation)
library(Hmisc)

# load data ---------------------------------------------------------------

abundance <- read.table("../data/abundance-genus.txt", sep = "\t", header = T)
rownames(abundance) <- abundance$taxonomy
abundance$taxonomy <- NULL

# metadata
metadata_ori <- read_tsv("../data/metadata.tsv")

# ae
ae <- read.table("../data/AE.txt", sep = "\t", header = T)




# process -----------------------------------------------------------------

metadata <- subset(metadata_ori, metadata_ori$Midgroup != "categorical")
colnames(metadata)[1] <- "SampleID"
metadata_plot <- dplyr::select(metadata, SampleID, Midgroup, Group)

metadata_copd <- subset(metadata, Group == "COPD")
metadata_copd$group <- metadata_copd$Genotype
metadata_copd_com <- dplyr::select(metadata_copd, SampleID, group)
metadata_copd_com$SampleID <- paste0(substr(metadata_copd_com$SampleID, 16, 17), "_", str_split_fixed(metadata_copd_com$SampleID, "\\.", 3)[,2])

# ae
ae$SampleID <- paste0(substr(ae$SampleID, 15, 16), "_", str_split_fixed(ae$SampleID, "\\.", 3)[,2])
ae_com <- dplyr::select(ae, SampleID, HospitalizedAEV2, TotalAEV2, FE)

metadata_plot_ae <- left_join(metadata_copd_com, ae_com, by = "SampleID")
metadata_plot_ae <- drop_na(metadata_plot_ae)
metadata_plot_ae$secre <- substr(metadata_plot_ae$group, 1,5)
metadata_plot_ae$class <- paste0(metadata_plot_ae$secre,"_", metadata_plot_ae$TotalAEV2)
metadata_plot_ae$AE <- ifelse(metadata_plot_ae$TotalAEV2 == 0, "No-AE", "AE")
metadata_plot_ae$class_2 <- paste0(metadata_plot_ae$secre,"_", metadata_plot_ae$AE)

# a diversity -------------------------------------------------------------

shan_index <- vegan::diversity(t(abundance), "shannon")
simpson_index <- vegan::diversity(t(abundance), "simpson")
invsimpson_index <- vegan::diversity(t(abundance), "invsimpson")


a_diversity <- data.frame(colnames(abundance), shan_index, simpson_index, invsimpson_index)
colnames(a_diversity)[1] <- "SampleID"

adiversity <- left_join(a_diversity, metadata_plot, by = "SampleID")
colnames(adiversity)[2] <-c("Shannon")
colnames(adiversity)[3] <-c("Simpson")
colnames(adiversity)[4] <-c("Invsimpson")

plot_shannon <- ggplot(adiversity, aes(x = factor(Midgroup, levels = c("HCon", "COPDSecre", "COPDNonse")),  y = Shannon))+ 
  geom_violin(aes(fill = Midgroup)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586" ,  "#BFC096","#89c3eb")) + 
  theme(legend.position = "none")

plot_Simpson <- ggplot(adiversity, aes(x = factor(Midgroup, levels = c("HCon", "COPDSecre", "COPDNonse")),  y = Simpson))+ 
  geom_violin(aes(fill = Midgroup)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586", "#BFC096", "#89c3eb")) + 
  theme(legend.position = "none")


plot_Invsimpson <- ggplot(adiversity, aes(x = factor(Midgroup, levels = c("HCon", "COPDSecre", "COPDNonse")),  y = Invsimpson))+ 
  geom_violin(aes(fill = Midgroup)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586", "#BFC096", "#89c3eb")) + 
  theme(legend.position = "none")


pdf("../version1/figures/shannon_diversity.pdf", width = 6/2.54, height = 4.6/2.54)
plot_shannon
graphics.off()

pdf("../version1/figures/simpson_diversity.pdf", width = 6/2.54, height = 4.6/2.54)
plot_Simpson
graphics.off()

pdf("../version1/figures/insimpson_diversity.pdf", width = 6/2.54, height = 4.6/2.54)
plot_Invsimpson
graphics.off()


# test
adiversity_HCon <- subset(adiversity, Midgroup == "HCon")
adiversity_COPDSecre <- subset(adiversity, Midgroup == "COPDSecre")
adiversity_COPDNonse <- subset(adiversity, Midgroup == "COPDNonse")

wilcox.test(adiversity_HCon$Shannon, adiversity_COPDNonse$Shannon, alternative = "two.sided")
wilcox.test(adiversity_HCon$Shannon, adiversity_COPDSecre$Shannon, alternative = "two.sided")
wilcox.test(adiversity_COPDNonse$Shannon, adiversity_COPDSecre$Shannon, alternative = "two.sided")



wilcox.test(adiversity_HCon$Simpson, adiversity_COPDNonse$Simpson, alternative = "two.sided")
wilcox.test(adiversity_HCon$Simpson, adiversity_COPDSecre$Simpson, alternative = "two.sided")
wilcox.test(adiversity_COPDNonse$Simpson, adiversity_COPDSecre$Simpson, alternative = "two.sided")


wilcox.test(adiversity_HCon$Invsimpson, adiversity_COPDNonse$Invsimpson, alternative = "two.sided")
wilcox.test(adiversity_HCon$Invsimpson, adiversity_COPDSecre$Invsimpson, alternative = "two.sided")
wilcox.test(adiversity_COPDNonse$Invsimpson, adiversity_COPDSecre$Invsimpson, alternative = "two.sided")

write.table(adiversity, file = "../version1/results/a_diversity.txt", sep = "\t", col.names = T, 
            quote = FALSE, row.names = F)




# PCoA --------------------------------------------------------------------

dist <- vegdist(t(abundance), method = "bray")

# Perform PCoA using pcoa() function
pcoa <- pcoa(dist)

# Compute eigenvalues 
eig <- pcoa$values
pct_var <- round(100 * eig$Eigenvalues / sum(eig$Eigenvalues), 2)

# Extract the PCoA scores for the first two axes and create a data frame
pcoa_plot <- pcoa$vectors[,1:2] %>% as.data.frame()
colnames(pcoa_plot) <- c("PcoA1", "PcoA2")
pcoa_plot$SampleID <- rownames(pcoa_plot)

# Merge PCoA data with sample metadata
pcoa_plot <- left_join(pcoa_plot, metadata_plot, by = "SampleID")

# Create PCoA plot using ggplot2
plot_pcoa <- ggplot(pcoa_plot, aes(x = PcoA1, y = PcoA2, 
                                   color = factor(Midgroup, levels = c("HCon", "COPDSecre", "COPDNonse")))) + 
  geom_point(size = 0.5, alpha = 0.7) +
  theme_test(base_size = 8) + stat_ellipse(level = 0.5) +
  scale_color_manual("Group", values = c("#89c3eb", "#BFC096","#c97586"  )) + 
  labs(x = paste0("PcoA1 (", pct_var[1], "%)"),
       y = paste0("PcoA2 (", pct_var[2], "%)")) + theme(legend.position = "none")

pdf("../version1/figures/pcoa.pdf", width = 4.5/2.54, height = 4.6/2.54)
plot_pcoa
graphics.off()


# permanova
metadata_pcoa <- data.frame(metadata_plot)
rownames(metadata_pcoa) <- metadata_pcoa$SampleID
metadata_pcoa <- metadata_pcoa[colnames(abundance),]

adonis <- adonis2(data.frame(t(abundance)) ~ Midgroup, metadata_pcoa, permutations = 999, na.rm = T)
adonis 



# test
pcoa_HCon <- subset(pcoa_plot, Midgroup == "HCon")
pcoa_COPDSecre <- subset(pcoa_plot, Midgroup == "COPDSecre")
pcoa_COPDNonse <- subset(pcoa_plot, Midgroup == "COPDNonse")


wilcox.test(pcoa_COPDSecre$PcoA1, pcoa_COPDNonse$PcoA1, alternative = "two.sided")




# pairwise permanova
# no hc
metadata_pcoa_pairwise_no_hc <- subset(metadata_pcoa, Midgroup != "HCon")
abundance_pcoa_pairwise_no_hc <- abundance[, colnames(abundance) %in% metadata_pcoa_pairwise_no_hc$SampleID]

adonis_no_hc <- adonis2(data.frame(t(abundance_pcoa_pairwise_no_hc)) ~ Midgroup, metadata_pcoa_pairwise_no_hc, 
                  permutations = 999, na.rm = T)
adonis_no_hc 


# no COPDSecre
metadata_pcoa_pairwise_no_COPDSecre <- subset(metadata_pcoa, Midgroup != "COPDSecre")
abundance_pcoa_pairwise_no_COPDSecre <- abundance[, colnames(abundance) %in% metadata_pcoa_pairwise_no_COPDSecre$SampleID]

adonis_no_COPDSecre <- adonis2(data.frame(t(abundance_pcoa_pairwise_no_COPDSecre)) ~ Midgroup, metadata_pcoa_pairwise_no_COPDSecre, 
                               permutations = 999, na.rm = T)
adonis_no_COPDSecre


# no COPDNonse
metadata_pcoa_pairwise_no_COPDNonse <- subset(metadata_pcoa, Midgroup != "COPDNonse")
abundance_pcoa_pairwise_no_COPDNonse <- abundance[, colnames(abundance) %in% metadata_pcoa_pairwise_no_COPDNonse$SampleID]

adonis_no_COPDNonse <- adonis2(data.frame(t(abundance_pcoa_pairwise_no_COPDNonse)) ~ Midgroup, metadata_pcoa_pairwise_no_COPDNonse, 
                               permutations = 999, na.rm = T)
adonis_no_COPDNonse




# barplot -----------------------------------------------------------------

get_barplot_genus <- function(MICRO_DATA, METADATA, GROUP){
  MICRO_DATA = MICRO_DATA[(order(-rowSums(MICRO_DATA))), ]
  metadata <- select(METADATA, SampleID, GROUP)
  colnames(metadata)[2] <- "Group"
  
  # combine low abundance taxonomy into Low abundance
  other = colSums(MICRO_DATA[22 : dim(MICRO_DATA)[1], ])
  MICRO_DATA = MICRO_DATA[1:(22 - 1), ]
  MICRO_DATA = rbind(MICRO_DATA, other)
  rownames(MICRO_DATA)[22] = c("Others")
  
  sample_order <- data.frame(t(MICRO_DATA))
  MICRO_DATA$Taxonomy <- rownames(MICRO_DATA)
  
  data2plot <- data.frame(melt(MICRO_DATA, id.vars = c("Taxonomy")))
  data2plot <- merge(data2plot, metadata, by.x = "variable", by.y = "SampleID")
  
  
  sample_order$variable <- as.factor(rownames(sample_order))
  data4plot <- left_join(sample_order, data2plot, by = "variable")
  data4plot <- data4plot[order(data4plot[,1], decreasing = T), ]
  data4plot$variable <- factor(data4plot$variable, levels = unique(data4plot$variable))
  
  p_indi <- ggplot(data4plot, aes(x = factor(variable), 
                                  y = value, fill = factor(Taxonomy, levels = rownames(MICRO_DATA)[22 : 1]))) + 
    geom_bar(stat = "identity",position="fill", width=1) + theme_test(base_size = 8) + 
    facet_grid(~ factor(Group, levels = c("HCon", "COPDSecre", "COPDNonse")), scales = "free_x", switch = "x" ) +
    scale_y_continuous(labels = scales::percent) + 
    theme(strip.background = element_blank()) +
    theme(axis.text.x=element_blank(), legend.key.size = unit(0.3, "cm"),
          axis.ticks.x=element_blank()) + 
    ylab("Percentage") + xlab("") +
    labs(fill="Taxonomy") + scale_fill_manual(values = c("Others" = "#C0C0C0", 
                                                         "k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Brucellaceae;g__Pseudochrobactrum" = "#F5DEB3",
                                                         "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Selenomonas" = "#1B9E77", 
                                                         "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas" = "#CE93BF",
                                                         "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Carnobacteriaceae;g__Granulicatella" = "#89c3eb", 
                                                         "k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus" = "#c97586", 
                                                         "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Micrococcaceae;g__Rothia" = "#BFC096", 
                                                         "k__Bacteria;p__Spirochaetes;c__Spirochaetes;o__Spirochaetales;f__Spirochaetaceae;g__Treponema" = "#FFE4E1",
                                                         "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Actinobacillus" = "#98FB98", 
                                                         "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Moraxellaceae;g__Moraxella" = "#FF00FF", 
                                                         "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__.Paraprevotellaceae.;g__.Prevotella." = "#F0CFE3", 
                                                         "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Corynebacteriaceae;g__Corynebacterium" = "#DAA520", 
                                                         "k__Bacteria;p__TM7;c__TM7.3;o__;f__;g__" = "#FFB6C1", 
                                                         "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Leptotrichiaceae;g__Leptotrichia" = "#4dbbd5",
                                                         "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium" = "#f39b7f",
                                                         "k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces" = "#3c54bb", 
                                                         "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Porphyromonadaceae;g__Porphyromonas" = "#00a087",
                                                         "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Veillonellaceae;g__Veillonella" = "#91d1c2",
                                                         "k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__Neisseria" = "#8491b4",
                                                         "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus" = "#DC143C", 
                                                         "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella" = "#b09c85",
                                                         "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus" = "#e64b35")) +
    guides(fill = guide_legend(ncol=3, bycol=TRUE))
  # average in group
  data2group <- aggregate(data4plot[1:22], by = data4plot[22+4], FUN = mean)
  data4group = as.data.frame(melt(data2group, id.vars="Group"))
  data4group$variable <- gsub("Low.abundance", "Low abundance", data4group$variable)
  
  p_group = ggplot(data4group, aes(x= factor(Group, levels = c("HCon", "COPDSecre", "COPDNonse")),
                                   y = value, fill = factor(variable, levels = colnames(data2group)[ncol(data2group):2]))) + 
    geom_bar(stat = "identity",position="fill", width=0.7)+ 
    scale_y_continuous(labels = scales::percent) + xlab("") +
    ylab("Percentage") + theme_test(base_size = 8) + 
    theme(legend.key.size = unit(0.3, "cm"), legend.position = "none", axis.text.x = element_text(angle=60, vjust = 0.5)) +
    labs(fill="Taxonomy") + scale_fill_manual(values = c("Others" = "#C0C0C0", 
                                                         "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Brucellaceae.g__Pseudochrobactrum" = "#F5DEB3",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Selenomonas" = "#1B9E77", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Pseudomonadaceae.g__Pseudomonas" = "#CE93BF",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Carnobacteriaceae.g__Granulicatella" = "#89c3eb", 
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Staphylococcaceae.g__Staphylococcus" = "#c97586", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Micrococcaceae.g__Rothia" = "#BFC096", 
                                                         "k__Bacteria.p__Spirochaetes.c__Spirochaetes.o__Spirochaetales.f__Spirochaetaceae.g__Treponema" = "#FFE4E1",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Actinobacillus" = "#98FB98", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__Moraxella" = "#FF00FF", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__.Paraprevotellaceae..g__.Prevotella." = "#F0CFE3", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Corynebacteriaceae.g__Corynebacterium" = "#DAA520", 
                                                         "k__Bacteria.p__TM7.c__TM7.3.o__.f__.g__" = "#FFB6C1", 
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia" = "#4dbbd5",
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium" = "#f39b7f",
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces" = "#3c54bb", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Porphyromonadaceae.g__Porphyromonas" = "#00a087",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Veillonella" = "#91d1c2",
                                                         "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Neisseriales.f__Neisseriaceae.g__Neisseria" = "#8491b4",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus" = "#DC143C", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella" = "#b09c85",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus" = "#e64b35")) +
    guides(fill = guide_legend(ncol=3, bycol=TRUE))
  return(list(p_indi, p_group))
} 

MICRO_DATA = abundance
METADATA = metadata_plot
GROUP = "Midgroup"

plot_indi <- get_barplot_genus(MICRO_DATA = abundance, METADATA = metadata_plot, GROUP = "Midgroup")[1]
pdf("../version1/figures/genus_indi.pdf", width = 50/2.54, height = 4.6/2.54)
plot_indi
graphics.off()

plot_group <- get_barplot_genus(MICRO_DATA = abundance, METADATA = metadata_plot, GROUP = "Midgroup")[2]
pdf("../version1/figures/genus_group.pdf", width = 3/2.54, height = 4.8/2.54)
plot_group
graphics.off()


# differential bacteria ---------------------------------------------------

abundance_diff <- abundance[apply(abundance == 0, 1, sum) <= (ncol(abundance) * 0.9), ]
abundance_t <- data.frame(t(abundance_diff))
abundance_t$SampleID <- rownames(abundance_t)
adiversity_com <- dplyr::select(adiversity, SampleID, Midgroup)
data_abundance_t <- left_join(abundance_t, adiversity_com, by = "SampleID")


data_results <- data.frame()
microbiome_names <- colnames(data_abundance_t)[1:127]

# i <- "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus"
for (i in microbiome_names){
  df <- dplyr::select(data_abundance_t, i, Midgroup, SampleID)
  colnames(df)[1] <- "microbiome"
  df_COPDSecre <- subset(df, Midgroup == "COPDSecre")
  df_TT <- subset(df, Midgroup == "COPDNonse")
  df_HC <- subset(df, Midgroup == "HCon")
  
  p_COPDSecre_TT <- wilcox.test(df_COPDSecre$microbiome, df_TT$microbiome, na.action = na.omit)[["p.value"]]
  p_TT_HC <- wilcox.test(df_TT$microbiome, df_HC$microbiome, na.action = na.omit)[["p.value"]]
  p_COPDSecre_HC <- wilcox.test(df_COPDSecre$microbiome, df_HC$microbiome, na.action = na.omit)[["p.value"]]
  
  if (mean(df_COPDSecre$microbiome) > mean(df_TT$microbiome) & mean(df_COPDSecre$microbiome) > mean(df_HC$microbiome)){
    enrich <- "COPDSecre" } else if (mean(df_HC$microbiome) > mean(df_COPDSecre$microbiome) & mean(df_HC$microbiome) > mean(df_TT$microbiome)){
        enrich <- "HC" } else {
          enrich <- "TT"}
  
  onereport = data.frame(microbiome = i,
                         p_COPDSecre_TT = p_COPDSecre_TT,
                         p_TT_HC = p_TT_HC, 
                         p_COPDSecre_HC = p_COPDSecre_HC,
                         enrich = enrich)
  data_results <- rbind(data_results, onereport)
  

    p_microbiome <- ggplot(data = df) + 
      geom_violin(aes(x = factor(Midgroup, levels =  c("HCon", "COPDSecre", "COPDNonse")), y = microbiome, 
                      fill = factor(Midgroup, levels =  c("HCon", "COPDSecre", "COPDNonse")))) + 
      geom_boxplot(aes(x = Midgroup, y = microbiome), width = 0.1, outlier.shape = NA) +
      xlab(i) + ylab("Abundance") +
      theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + 
      scale_fill_manual(values = c("#89c3eb", "#BFC096", "#c97586")) +
      theme(legend.position = "none")
    
    pdf_file <- paste0("../version1/figures/bacteria_diff/", i, ".pdf")
    pdf(pdf_file, width = 6/2.54, height = 4/2.54)
    print(p_microbiome) 
    dev.off()

}
data_results$FDR_COPDSecre_TT <- p.adjust(data_results$p_COPDSecre_TT, method = "fdr")
data_results$FDR_TT_HC <- p.adjust(data_results$p_TT_HC, method = "fdr")
data_results$FDR_COPDSecre_HC <- p.adjust(data_results$p_COPDSecre_HC, method = "fdr")

write.table(data_results, "../version1/results/differential_bacteria_Midgroup_ra.txt", 
            sep = '\t', col.names = T, quote = FALSE, row.names = F)


# AE ----------------------------------------------------------------------

# a diversity

adiversity_ae <- adiversity
adiversity_ae$SampleID <- paste0(substr(adiversity_ae$SampleID, 16, 17), "_", str_split_fixed(adiversity_ae$SampleID, "\\.", 3)[,2])

adiversity_ae <- left_join(metadata_plot_ae, adiversity_ae, by = "SampleID") %>% drop_na()




write.table(adiversity_ae, "../version1/results/a_diversity_ae.txt", 
            sep = '\t', col.names = T, quote = FALSE, row.names = F)

plot_Shannon_ae <- ggplot(adiversity_ae , aes(x = AE,  y = Shannon))+ 
  geom_violin(aes(fill = AE)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586" , "#BFC096", "#89c3eb")) + 
  theme(legend.position = "none")

pdf("../version1/figures/AE/shannon_diversity.pdf", width = 6/2.54, height = 4.6/2.54)
plot_Shannon_ae
graphics.off()

plot_Simpson_ae <- ggplot(adiversity_ae , aes(x = AE,  y = Simpson))+ 
  geom_violin(aes(fill = AE)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586" , "#BFC096", "#89c3eb")) + 
  theme(legend.position = "none")

pdf("../version1/figures/AE/Simpson_ae.pdf", width = 6/2.54, height = 4.6/2.54)
plot_Simpson_ae
graphics.off()

plot_Invsimpson_ae <- ggplot(adiversity_ae , aes(x = AE,  y = Invsimpson))+ 
  geom_violin(aes(fill = AE)) + geom_boxplot(width = 0.1, outlier.shape = NA) + 
  ylab("Bacterial diversity") + xlab("")  +
  theme(plot.title=element_text(hjust = 0.5)) + theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#c97586" , "#BFC096", "#89c3eb")) + 
  theme(legend.position = "none")

pdf("../version1/figures/AE/Invsimpson_ae.pdf", width = 6/2.54, height = 4.6/2.54)
plot_Invsimpson_ae
graphics.off()

# test
adiversity_ae_value <- subset(adiversity_ae, AE == "AE")
adiversity_not_ae_value <- subset(adiversity_ae, AE == "No-AE")

wilcox.test(adiversity_ae_value$Shannon, adiversity_not_ae_value$Shannon, alternative = "two.sided")
wilcox.test(adiversity_ae_value$Simpson, adiversity_not_ae_value$Simpson, alternative = "two.sided")
wilcox.test(adiversity_ae_value$Invsimpson, adiversity_not_ae_value$Invsimpson, alternative = "two.sided")


abundance_ae <- abundance
colnames(abundance_ae) <- paste0(substr(colnames(abundance_ae), 16, 17), "_", str_split_fixed(colnames(abundance_ae), "\\.", 3)[,2])
abundance_ae <- abundance_ae[, colnames(abundance_ae) %in% metadata_plot_ae$SampleID]

rownames(abundance_ae) <- gsub(";", ".", rownames(abundance_ae))


MICRO_DATA = abundance_ae
METADATA = metadata_plot_ae
GROUP = "AE"

get_barplot_genus_ae <- function(MICRO_DATA, METADATA, GROUP){
  MICRO_DATA = MICRO_DATA[(order(-rowSums(MICRO_DATA))), ]
  metadata <- select(METADATA, SampleID, GROUP)
  colnames(metadata)[2] <- "Group"
  
  # combine low abundance taxonomy into Low abundance
  other = colSums(MICRO_DATA[22 : dim(MICRO_DATA)[1], ])
  MICRO_DATA = MICRO_DATA[1:(22 - 1), ]
  MICRO_DATA = rbind(MICRO_DATA, other)
  rownames(MICRO_DATA)[22] = c("Others")
  
  sample_order <- data.frame(t(MICRO_DATA))
  MICRO_DATA$Taxonomy <- rownames(MICRO_DATA)
  
  data2plot <- data.frame(melt(MICRO_DATA, id.vars = c("Taxonomy")))
  data2plot <- merge(data2plot, metadata, by.x = "variable", by.y = "SampleID")
  
  
  sample_order$variable <- as.factor(rownames(sample_order))
  data4plot <- left_join(sample_order, data2plot, by = "variable")
  data4plot <- data4plot[order(data4plot[,1], decreasing = T), ]
  data4plot$variable <- factor(data4plot$variable, levels = unique(data4plot$variable))
  
  p_indi <- ggplot(data4plot, aes(x = factor(variable), 
                                  y = value, fill = factor(Taxonomy, levels = rownames(MICRO_DATA)[22 : 1]))) + 
    geom_bar(stat = "identity",position="fill", width=1) + theme_test(base_size = 8) + 
    facet_grid(~ factor(Group), scales = "free_x", switch = "x" ) +
    scale_y_continuous(labels = scales::percent) + 
    theme(strip.background = element_blank()) +
    theme(axis.text.x=element_blank(), legend.key.size = unit(0.3, "cm"),
          axis.ticks.x=element_blank()) + 
    ylab("Percentage") + xlab("") +
    labs(fill="Taxonomy") + scale_fill_manual(values = c("Others" = "#C0C0C0", 
                                                         "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Brucellaceae.g__Pseudochrobactrum" = "#F5DEB3",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Selenomonas" = "#1B9E77", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Pseudomonadaceae.g__Pseudomonas" = "#CE93BF",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Carnobacteriaceae.g__Granulicatella" = "#89c3eb", 
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Staphylococcaceae.g__Staphylococcus" = "#c97586", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Micrococcaceae.g__Rothia" = "#BFC096", 
                                                         "k__Bacteria.p__Spirochaetes.c__Spirochaetes.o__Spirochaetales.f__Spirochaetaceae.g__Treponema" = "#FFE4E1",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Actinobacillus" = "#98FB98", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__Moraxella" = "#FF00FF", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__.Paraprevotellaceae..g__.Prevotella." = "#F0CFE3", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Corynebacteriaceae.g__Corynebacterium" = "#DAA520", 
                                                         "k__Bacteria.p__TM7.c__TM7.3.o__.f__.g__" = "#FFB6C1", 
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia" = "#4dbbd5",
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium" = "#f39b7f",
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces" = "#3c54bb", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Porphyromonadaceae.g__Porphyromonas" = "#00a087",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Veillonella" = "#91d1c2",
                                                         "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Neisseriales.f__Neisseriaceae.g__Neisseria" = "#8491b4",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus" = "#DC143C", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella" = "#b09c85",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus" = "#e64b35")) +
    guides(fill = guide_legend(ncol=3, bycol=TRUE))
  # average in group
  data2group <- aggregate(data4plot[1:22], by = data4plot[22+4], FUN = mean)
  data4group = as.data.frame(melt(data2group, id.vars="Group"))
  data4group$variable <- gsub("Low.abundance", "Low abundance", data4group$variable)
  
  p_group = ggplot(data4group, aes(x= factor(Group),
                                   y = value, fill = factor(variable, levels = colnames(data2group)[ncol(data2group):2]))) + 
    geom_bar(stat = "identity",position="fill", width=0.7)+ 
    scale_y_continuous(labels = scales::percent) + xlab("") +
    ylab("Percentage") + theme_test(base_size = 8) + 
    theme(legend.key.size = unit(0.3, "cm"), legend.position = "none", axis.text.x = element_text(angle=60, vjust = 0.5)) +
    labs(fill="Taxonomy") + scale_fill_manual(values = c("Others" = "#C0C0C0", 
                                                         "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales.f__Brucellaceae.g__Pseudochrobactrum" = "#F5DEB3",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Selenomonas" = "#1B9E77", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Pseudomonadaceae.g__Pseudomonas" = "#CE93BF",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Carnobacteriaceae.g__Granulicatella" = "#89c3eb", 
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Staphylococcaceae.g__Staphylococcus" = "#c97586", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Micrococcaceae.g__Rothia" = "#BFC096", 
                                                         "k__Bacteria.p__Spirochaetes.c__Spirochaetes.o__Spirochaetales.f__Spirochaetaceae.g__Treponema" = "#FFE4E1",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Actinobacillus" = "#98FB98", 
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__Moraxella" = "#FF00FF", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__.Paraprevotellaceae..g__.Prevotella." = "#F0CFE3", 
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Corynebacteriaceae.g__Corynebacterium" = "#DAA520", 
                                                         "k__Bacteria.p__TM7.c__TM7.3.o__.f__.g__" = "#FFB6C1", 
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Leptotrichiaceae.g__Leptotrichia" = "#4dbbd5",
                                                         "k__Bacteria.p__Fusobacteria.c__Fusobacteriia.o__Fusobacteriales.f__Fusobacteriaceae.g__Fusobacterium" = "#f39b7f",
                                                         "k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Actinomycetaceae.g__Actinomyces" = "#3c54bb", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Porphyromonadaceae.g__Porphyromonas" = "#00a087",
                                                         "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Veillonellaceae.g__Veillonella" = "#91d1c2",
                                                         "k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Neisseriales.f__Neisseriaceae.g__Neisseria" = "#8491b4",
                                                         "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus" = "#DC143C", 
                                                         "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella" = "#b09c85",
                                                         "k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__Streptococcus" = "#e64b35")) +
    guides(fill = guide_legend(ncol=3, bycol=TRUE))
  return(list(p_indi, p_group))
} 

p <- get_barplot_genus_ae(MICRO_DATA = abundance_ae, METADATA = metadata_plot_ae, GROUP = "AE")
pdf("../version1/figures/AE/genus_indi.pdf", width = 50/2.54, height = 4.6/2.54)
p[1]
graphics.off()

pdf("../version1/figures/AE/genus_group.pdf", width = 2.6/2.54, height = 4.8/2.54)
p[2]
graphics.off()



# differential bacteria
abundance_diff_ae <- abundance_ae[apply(abundance_ae == 0, 1, sum) <= (ncol(abundance_ae) * 0.9), ]

abundance_t <- data.frame(t(abundance_diff_ae))
abundance_t$SampleID <- rownames(abundance_t)
adiversity_com <- dplyr::select(metadata_plot_ae, SampleID, AE)
data_abundance_t <- left_join(abundance_t, adiversity_com, by = "SampleID")


data_results <- data.frame()
microbiome_names <- colnames(data_abundance_t)[1:128]

i <- "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus"

# i <- "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus"
for (i in microbiome_names){
  df <- dplyr::select(data_abundance_t, i, AE, SampleID)
  colnames(df)[1] <- "microbiome"
  df_AE <- subset(df, AE == "AE")
  df_No_AE <- subset(df, AE == "No-AE")

  
  p <- wilcox.test(df_AE$microbiome, df_No_AE$microbiome, na.action = na.omit, alternative = "greater")[["p.value"]]

  if (mean(df_AE$microbiome) > mean(df_No_AE$microbiome)){
    enrich <- "AE" }  else {
        enrich <- "No_AE"}
  
  onereport = data.frame(microbiome = i,
                         p = p,
                         enrich = enrich)
  data_results <- rbind(data_results, onereport)
  
    p_microbiome <- ggplot(data = df) + 
      geom_violin(aes(x = factor(AE), y = microbiome, 
                     fill = factor(AE))) + 
     geom_boxplot(aes(x = AE, y = microbiome), width = 0.1, outlier.shape = NA) +
     xlab(i) + ylab("Abundance") +
     theme_classic(base_size = 8) + theme(panel.border = element_rect( fill = NA)) + 
     scale_fill_manual(values = c("#89c3eb", "#BFC096", "#c97586")) +
     theme(legend.position = "none")
    
   pdf_file <- paste0("../version1/figures/AE/AE_diff/", i, ".pdf")
   pdf(pdf_file, width = 6/2.54, height = 4/2.54)
   print(p_microbiome) 
   dev.off()
}

data_results$FDR <- p.adjust(data_results$p, method = "fdr")

write.table(data_results, "../version1/results/differential_bacteria_AE_ra.txt", 
            sep = '\t', col.names = T, quote = FALSE, row.names = F)

# haem and shannon --------------------------------------------------------------

abundance_ae_t <- data.frame(t(abundance_ae))
abundance_correlation <- dplyr::select(abundance_ae_t, k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus)
abundance_correlation$SampleID <- rownames(abundance_correlation)

a_diversity_ae <- a_diversity
a_diversity_ae$SampleID <- paste0(substr(a_diversity_ae$SampleID, 16, 17), "_", str_split_fixed(a_diversity_ae$SampleID, "\\.", 3)[,2])
abundance_correlation <- left_join(abundance_correlation, a_diversity_ae, by = "SampleID")


# heam correlation

plot_haem <- ggplot(data = abundance_correlation, aes(x = k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus, 
                                                      y = shan_index)) + xlab("Haemophilus")  + xlab("Shannon index") + 
  geom_point(alpha = 1, size = 0.5, color = "#5F9EA0") + geom_smooth(method = "lm", color = "#008080") + 
  theme_test(base_size = 8 ) 


pdf("../version1/figures/Haemophilus_shannon.pdf", width = 4.6/2.54, height = 4.6/2.54)
plot_haem
graphics.off()

fit_1 <- lm(shan_index ~ k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus, 
            data = abundance_correlation)
summary(fit_1)


spearman_result_1 <- cor.test(abundance_correlation$shan_index, abundance_correlation$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus)
spearman_result_1 


# hea
abundance_hae <- left_join(abundance_correlation, metadata_plot_ae, by = "SampleID")
p_hae <- ggplot(data = abundance_hae) + 
  geom_violin(aes(x = factor(AE), y = log2(k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus * 100 ), 
                  fill = factor(AE))) + 
  geom_boxplot(aes(x = AE, y = log2(k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus * 100)), width = 0.1, outlier.shape = NA) +
  xlab("Haemophilus") + ylab("Log2(Relative abundance)") +
  theme_classic(base_size = 8) + 
  scale_fill_manual(values = c("#7bccc4", "#A6559D")) +
  theme(legend.position = "none")


pdf("../version1/figures/AE/Haemophilus.pdf", width = 4.5/2.54, height = 4.6/2.54)
print(p_hae) 
dev.off()

abundance_ae_hae <- subset(abundance_hae, AE == "AE")
abundance_no_ae_hae <- subset(abundance_hae, AE == "No-AE")

wilcox.test(abundance_ae_hae$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus, 
            abundance_no_ae_hae$k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus, alternative = "two.sided")






