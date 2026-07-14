# Prepare random subsets of variants from VQSR files
# We used the 4 million randomly selected variant positions to create the subsets for VQSR evaluation

### R Script 
### Variant Calling in a mini-core collection of Brassica napus
### Author: Thomas Bergmann

### Load libraries
library(VariantAnnotation)
library(dplyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(tibble)
library(ClusterR)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(purrr)
library(scales)
library(ggrepel)
library(mclust)
library(circlize)
library(paletteer)
library(nnet)
library(ggalluvial)
library(viridisLite)  # load package

### Set path
setwd("C:/Users/00115569/OneDrive - UWA/GATK_Pipeline/Minicore_Collection/VQSR")

# Note: genotype information was removed with bcftools prior to reduce VCF size

# Prepare data

prep_VCF <- function(inputVCF) {

  vcf <- readVcf(inputVCF)

  df <- info(vcf) %>%
    as.data.frame() %>%
    rownames_to_column("VariantID") %>%
    mutate(
      FILTER = fixed(vcf)$FILTER,
      PASS   = FILTER == "PASS"
    ) %>%
    select(
      VariantID,
      FILTER,
      PASS,
      QD,
      FS,
      SOR,
      MQ,
      MQRankSum,
      ReadPosRankSum
    ) %>%
	mutate(
		log_FS = log(FS + 1)
	)

  df
}

Tranche999 <- prep_VCF("Minicore.tranche.99.9.random.4M.subset.annotations.vcf.gz")
Tranche995 <- prep_VCF("Minicore.tranche.99.5.random.4M.subset.annotations.vcf.gz")
Tranche990 <- prep_VCF("Minicore.tranche.99.0.random.4M.subset.annotations.vcf.gz")
Tranche985 <- prep_VCF("Minicore.tranche.98.5.random.4M.subset.annotations.vcf.gz")
Tranche980 <- prep_VCF("Minicore.tranche.98.0.random.4M.subset.annotations.vcf.gz")
Hardfiltered <- prep_VCF("Minicore.tranche.hardfiltered.random.4M.subset.annotations.vcf.gz")

# 
table(Tranche999$PASS)
#  FALSE    TRUE 
#1708497 2286350

#Passing rate: 57

table(Tranche995$PASS)
#  FALSE    TRUE 
#1991651 2003196

#Passing rate: 50

table(Tranche990$PASS)
#  FALSE    TRUE 
#2152660 1842187

#Passing rate: 46

table(Tranche985$PASS)
#  FALSE    TRUE 
#2277046 1717801

#Passing rate: 43

table(Tranche980$PASS)
#  FALSE    TRUE 
#2399994 1594853

#Passing rate: 40

table(Hardfiltered$PASS)
#  FALSE    TRUE 
#1454988 2539859

#Passing rate:63

# Plot distribution

# Set different scales 
QD_breaks <- c(0,10,20,30,40)
QD_limits <- c(0,50)
logDP_breaks <- c(0,2,4,6,8,10,12)
logDP_limits <- c(0,13)
logFS_breaks <- c(-1,0,2,4,6,8)
logFS_limits <- c(-2,9)
MQ_breaks <- c(0,10,20,30,40,50,60)
MQ_limits <- c(0,70)
SOR_breaks <- c(0,2,4,6,8,10,12,14,16)
SOR_limits <- c(0,18)
MQRS_breaks <- c(-60,-40,-20,0,20,40,60)
MQRS_limits <- c(-80,80)
RPRS_breaks <- c(-30,-20,-10,0,10,20,30)
RPRS_limits <- c(-35,35)

## Function to make density plots
plot_density <- function(data, var, threshold, breaks, limits, label, show_legend = FALSE) {
  ggplot(data, aes(x = .data[[var]], fill = PASS, color = PASS)) + #fill = chrom_var)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = threshold, color = "red", linewidth = 1, linetype = "dashed") +
	labs(title = paste(label), x = "", y = "Density") +
	scale_y_continuous(expand=c(0,0), breaks=c(0,1,2,3), limits=c(0,3.5))+
	scale_x_continuous(expand=c(0,0), breaks=breaks, limits=limits)+
	#scale_color_manual(
	#	values = 
	#	)+
	theme(
		plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
		axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
		legend.position = if (show_legend) c(0.9,0.9) else "none",
		legend.justification = "center",
		legend.background = element_blank(),
		legend.key = element_blank(),
		legend.key.width  = unit(1.2, "lines"),
		legend.key.height = unit(1.0, "lines"),
		legend.title = element_text(face = "bold", size = 12),
		legend.text  = element_text(size = 12),
		plot.title.position = "plot") +
	guides(
		color = guide_legend(
		direction = "horizontal",
		nrow = 1,
		byrow = TRUE,
		title.position = "left",
		title.hjust = 0.5,
		override.aes = list(linewidth = 1)
  )
)
}

# Plot QD
QD_Tranche999 <- plot_density(Tranche999, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Tranche 99.9]", show_legend = TRUE)
QD_Tranche995 <- plot_density(Tranche995, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Tranche 99.5]", show_legend = TRUE)
QD_Tranche990 <- plot_density(Tranche990, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Tranche 99.0]", show_legend = TRUE)
QD_Tranche985 <- plot_density(Tranche985, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Tranche 98.5]", show_legend = TRUE)
QD_Tranche980 <- plot_density(Tranche980, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Tranche 98.0]", show_legend = TRUE)
QD_HF <- plot_density(Hardfiltered, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
QD_Plots <- cowplot::plot_grid(QD_Tranche999, QD_Tranche995, QD_Tranche990, 
							   QD_Tranche985, QD_Tranche980, QD_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./QD_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
QD_Plots
#dev.off()

## Tranche 99.9 vs. HF
QD_Plots <- cowplot::plot_grid(QD_Tranche999, QD_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./QD_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
QD_Plots
#dev.off()

# Plot FS
FS_Tranche999 <- plot_density(Tranche999, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 99.9]", show_legend = TRUE)
FS_Tranche995 <- plot_density(Tranche995, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 99.5]", show_legend = TRUE)
FS_Tranche990 <- plot_density(Tranche990, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 99.0]", show_legend = TRUE)
FS_Tranche985 <- plot_density(Tranche985, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 98.5]", show_legend = TRUE)
FS_Tranche980 <- plot_density(Tranche980, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 98.0]", show_legend = TRUE)
FS_HF <- plot_density(Hardfiltered, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
FS_Plots <- cowplot::plot_grid(FS_Tranche999, FS_Tranche995, FS_Tranche990, 
							   FS_Tranche985, FS_Tranche980, FS_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
png("./FS_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
FS_Plots
dev.off()

## Tranche 99.9 vs. HF
FS_Plots <- cowplot::plot_grid(FS_Tranche999, FS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./FS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
FS_Plots
#dev.off()

# Plot SOR
SOR_Tranche999 <- plot_density(Tranche999, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 99.9]", show_legend = TRUE)
SOR_Tranche995 <- plot_density(Tranche995, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 99.5]", show_legend = TRUE)
SOR_Tranche990 <- plot_density(Tranche990, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 99.0]", show_legend = TRUE)
SOR_Tranche985 <- plot_density(Tranche985, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 98.5]", show_legend = TRUE)
SOR_Tranche980 <- plot_density(Tranche980, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 98.0]", show_legend = TRUE)
SOR_HF <- plot_density(Hardfiltered, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
SOR_Plots <- cowplot::plot_grid(SOR_Tranche999, SOR_Tranche995, SOR_Tranche990, 
							   SOR_Tranche985, SOR_Tranche980, SOR_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./SOR_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
SOR_Plots
#dev.off()

## Tranche 99.9 vs. HF
SOR_Plots <- cowplot::plot_grid(SOR_Tranche999, SOR_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./SOR_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
SOR_Plots
#dev.off()


# Plot MQ
MQ_Tranche999 <- plot_density(Tranche999, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 99.9]", show_legend = TRUE)
MQ_Tranche995 <- plot_density(Tranche995, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 99.5]", show_legend = TRUE)
MQ_Tranche990 <- plot_density(Tranche990, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 99.0]", show_legend = TRUE)
MQ_Tranche985 <- plot_density(Tranche985, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 98.5]", show_legend = TRUE)
MQ_Tranche980 <- plot_density(Tranche980, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 98.0]", show_legend = TRUE)
MQ_HF <- plot_density(Hardfiltered, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
MQ_Plots <- cowplot::plot_grid(MQ_Tranche999, MQ_Tranche995, MQ_Tranche990, 
							   MQ_Tranche985, MQ_Tranche980, MQ_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./MQ_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
MQ_Plots
#dev.off()

## Tranche 99.9 vs. HF
MQ_Plots <- cowplot::plot_grid(MQ_Tranche999, MQ_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./MQ_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
MQ_Plots
#dev.off()


# Plot MQRS
MQRS_Tranche999 <- plot_density(Tranche999, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 99.9]", show_legend = TRUE)
MQRS_Tranche995 <- plot_density(Tranche995, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 99.5]", show_legend = TRUE)
MQRS_Tranche990 <- plot_density(Tranche990, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 99.0]", show_legend = TRUE)
MQRS_Tranche985 <- plot_density(Tranche985, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 98.5]", show_legend = TRUE)
MQRS_Tranche980 <- plot_density(Tranche980, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 98.0]", show_legend = TRUE)
MQRS_HF <- plot_density(Hardfiltered, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
MQRS_Plots <- cowplot::plot_grid(MQRS_Tranche999, MQRS_Tranche995, MQRS_Tranche990, 
							   MQRS_Tranche985, MQRS_Tranche980, MQRS_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./MQRS_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
MQRS_Plots
#dev.off()

## Tranche 99.9 vs. HF
MQRS_Plots <- cowplot::plot_grid(MQRS_Tranche999, MQRS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./MQRS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
MQRS_Plots
#dev.off()


# Plot RPRS
RPRS_Tranche999 <- plot_density(Tranche999, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.9]", show_legend = TRUE)
RPRS_Tranche995 <- plot_density(Tranche995, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.5]", show_legend = TRUE)
RPRS_Tranche990 <- plot_density(Tranche990, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.0]", show_legend = TRUE)
RPRS_Tranche985 <- plot_density(Tranche985, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 98.5]", show_legend = TRUE)
RPRS_Tranche980 <- plot_density(Tranche980, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 98.0]", show_legend = TRUE)
RPRS_HF <- plot_density(Hardfiltered, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
RPRS_Plots <- cowplot::plot_grid(RPRS_Tranche999, RPRS_Tranche995, RPRS_Tranche990, 
							   RPRS_Tranche985, RPRS_Tranche980, RPRS_HF,
							   nrow=2, ncol=3, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./RPRS_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=6000, height=3000, res=300)
RPRS_Plots
#dev.off()

## Tranche 99.9 vs. HF
RPRS_Plots <- cowplot::plot_grid(RPRS_Tranche999, RPRS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./RPRS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=3000, height=3000, res=300)
RPRS_Plots
#dev.off()

## Plot metrics
add_tranche_info <- function(df, tranche_name) {
  df %>%
    mutate(
      Tranche = tranche_name,
      Status  = if_else(PASS, "PASS", "FAIL")
    )
} 

t999 <- add_tranche_info(Tranche999, "Tranche_99.9")
hardfiltered <- add_tranche_info(Hardfiltered, "Hardfiltered")

master_tranche <- bind_rows(t999,hardfiltered)
master_tranche$Tranche <- factor(master_tranche$Tranche, levels = c("Hardfiltered","Tranche_99.9"))
master_tranche$Status <- factor(master_tranche$Status, levels = c("PASS","FAIL"))

make_boxplots <- function(data, var, threshold, breaks, limits, label, show_legend = FALSE) {
  ggplot(
    data,
    aes(
      x = Tranche,
      y = .data[[var]],
      fill = Status
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      outlier.size = 0.3
    ) +
    scale_fill_manual(values = c("PASS" = "#00BFC4", "FAIL" = "#F8766D")) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = breaks,
      limits = limits
    ) +
    geom_hline(
      yintercept = threshold,
      color = "red",
      linewidth = 1,
      linetype = "dashed"
    ) +
	labs(title = paste(label), x = "") +
		
    theme(
      plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
      axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
      axis.text.x  = element_text(colour = "black", size = 14, face = "bold"),
      axis.text.y  = element_text(colour = "black", size = 14, face = "bold"),
      legend.position = if (show_legend) c(0.9, 0.9) else "none",
      legend.justification = "center",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.width  = unit(1.2, "lines"),
      legend.key.height = unit(1.0, "lines"),
      legend.title = element_blank(),
      legend.text  = element_text(size = 12),
      plot.title.position = "plot"
    )
}

# QD plot
QD_Boxplots <- make_boxplots(master_tranche, "QD", 2.0, QD_breaks, QD_limits, "QualByDepth", show_legend=TRUE)
#png("./QD_Boxplot_Tranche999_vs_Hardfiltered_Darmor.png", width=2000, height=3000, res=300)
QD_Boxplots
#dev.off()

# FS plot
logFS_limits <- c(-0.5,8)
FS_Boxplots <- make_boxplots(master_tranche, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand)", show_legend=TRUE)
png("./FS_Boxplot_Tranche999_vs_Hardfiltered_Darmor.png", width=2000, height=3000, res=300)
FS_Boxplots
dev.off()
