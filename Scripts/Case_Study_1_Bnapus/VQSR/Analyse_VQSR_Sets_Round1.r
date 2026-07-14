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
setwd("PATH/TO/VQSR/Round_1")

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

Tranche999 <- prep_VCF("Samples_Bna_Minicore_1-135.subset.SNPs.sites.VQSR.99.9.tranche.vcf.gz")

# Label according to hard filtering
Tranche999 <- Tranche999 %>%
						mutate(PASS_HF = !(
								(QD < 2.0 & !is.na(QD)) | 
								(FS > 60.0 & !is.na(FS)) |
								(SOR > 3.0 & !is.na(SOR)) | 
								(MQ < 40.0 & !is.na(MQ)) |
								(MQRankSum < -12.5 & !is.na(MQRankSum)) |
								(ReadPosRankSum < -8 & !is.na(ReadPosRankSum))
									)
								)

# Compare passing rates
table(Tranche999$PASS)
# FALSE   TRUE 
#374123 625877

#Passing rate: 62.6

table(Tranche999$PASS_HF)
# FALSE   TRUE 
#358727 641273

#Passing rate: 64.1

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
MQRS_breaks <- c(-16,-12,-8,-4,0,4,8,12,16)
MQRS_limits <- c(-20,20)
RPRS_breaks <- c(-10,-8,-6,-4,-2,0,2,4,6,8,10)
RPRS_limits <- c(-12,12)

## Function to make density plots (adjust scales according to annotations!)
plot_density <- function(data, var, threshold, breaks, limits, label, show_legend = FALSE) {
  ggplot(data, aes(x = .data[[var]], fill = PASS, color = PASS)) + #fill = chrom_var)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = threshold, color = "red", linewidth = 1, linetype = "dashed") +
	labs(title = paste(label), x = "", y = "Density") +
	scale_y_continuous(expand=c(0,0), breaks=c(0,.2,.4,.6,.8), limits=c(0,.825))+
	scale_x_continuous(expand=c(0,0), breaks=breaks, limits=limits)+
	#scale_color_manual(
	#	values = 
	#	)+
	theme(
		plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
		axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
		legend.position = if (show_legend) c(0.85,0.85) else "none",
		legend.justification = "center",
		legend.background = element_blank(),
		legend.key = element_blank(),
		legend.key.width  = unit(1.2, "lines"),
		legend.key.height = unit(1.0, "lines"),
		legend.title = element_text(face = "bold", size = 12),
		legend.text  = element_text(size = 12)) +
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
QD_HF <- plot_density(Tranche999, "QD", 2, QD_breaks, QD_limits, "QualByDepth [Hardfiltered]", show_legend = TRUE)

## Combine in one plot
QD_Plots <- cowplot::plot_grid(QD_Tranche999, 
							   QD_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/QD_Plot_Tranches_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
QD_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/QD_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
QD_Tranche999
#dev.off()

# Plot FS
FS_Tranche999 <- plot_density(Tranche999, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 99.9]", show_legend = FALSE)
FS_HF <- plot_density(Tranche999, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [Hardfiltered]", show_legend = TRUE)

## Tranche 99.9 vs. HF
FS_Plots <- cowplot::plot_grid(FS_Tranche999, FS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/FS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
FS_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/FS_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
FS_Tranche999
#dev.off()

# Plot SOR
SOR_Tranche999 <- plot_density(Tranche999, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 99.9]", show_legend = FALSE)
SOR_HF <- plot_density(Tranche999, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [Hardfiltered]", show_legend = TRUE)

## Tranche 99.9 vs. HF
SOR_Plots <- cowplot::plot_grid(SOR_Tranche999, SOR_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/SOR_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
SOR_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/SOR_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
SOR_Tranche999
#dev.off()


# Plot MQ
MQ_Tranche999 <- plot_density(Tranche999, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Tranche 99.9]", show_legend = FALSE)
MQ_HF <- plot_density(Tranche999, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [Hardfiltered]", show_legend = TRUE)

## Tranche 99.9 vs. HF
MQ_Plots <- cowplot::plot_grid(MQ_Tranche999, MQ_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/MQ_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
MQ_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/MQ_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
MQ_Tranche999
#dev.off()

# Plot MQRS
MQRS_Tranche999 <- plot_density(Tranche999, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 99.9]", show_legend = FALSE)
MQRS_HF <- plot_density(Tranche999, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [Hardfiltered]", show_legend = TRUE)

## Tranche 99.9 vs. HF
MQRS_Plots <- cowplot::plot_grid(MQRS_Tranche999, MQRS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/MQRS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
MQRS_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/MQRS_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
MQRS_Tranche999
#dev.off()

# Plot RPRS
RPRS_Tranche999 <- plot_density(Tranche999, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.9]", show_legend = FALSE)
RPRS_HF <- plot_density(Tranche999, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [Hardfiltered]", show_legend = TRUE)

## Tranche 99.9 vs. HF
RPRS_Plots <- cowplot::plot_grid(RPRS_Tranche999, RPRS_HF,
							   nrow=2, ncol=1, align="hv", 
							   labels = "AUTO", label_size = 26)
#png("./Figures/RPRS_Plot_Tranche999_vs_Hardfiltered_Darmor.png", width=2500, height=3000, res=300)
RPRS_Plots
#dev.off()

## Plot only for tranche 99.9
#png("./Figures/RPRS_Plot_Tranche999_Darmor.png", width=3000, height=2000, res=300)
RPRS_Tranche999
#dev.off()

## Plot together
All_Plots <- cowplot::plot_grid(QD_Tranche999, FS_Tranche999,
								SOR_Tranche999, MQ_Tranche999,
								MQRS_Tranche999, RPRS_Tranche999,
							    nrow=3, ncol=2, align="hv", scale = 0.95,
							    labels = "AUTO", label_size = 28)
#png("./Figures/All_Plots_Tranche999_Darmor.png", width=4000, height=4000, res=300)
All_Plots
#dev.off()
