# Prepare random subsets of variants from VQSR files
# We used the one million randomly selected variant positions to create the subsets for VQSR evaluation

### R Script 
### Variant Calling in a wild collection of Brassica oleracea
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
setwd("PATH/TO/VCF/")

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

Tranche999 <- prep_VCF("Samples_Bol_Wild_1-7.INDELs.VQSR.99.90.tranche.round1.subset.annotations.vcf.gz")

# Label according to hard filtering
Tranche999 <- Tranche999 %>%
						mutate(PASS_HF = !(
								(QD < 2.0 & !is.na(QD)) | 
								(FS > 200.0 & !is.na(FS)) |
								(ReadPosRankSum < -20 & !is.na(ReadPosRankSum))
									)
								)

# Compare passing rates
table(Tranche999$PASS)
# FALSE   TRUE 
#491965 508035 

#Passing rate: 50.8

table(Tranche999$PASS_HF)
# FALSE   TRUE 
# 44083 955917

#Passing rate: 95.6

# Plot distribution

# Set different scales 
QD_breaks <- c(0,10,20,30,40)
QD_limits <- c(0,45)
logFS_breaks <- c(0,2,4)
logFS_limits <- c(-1,5)
RPRS_breaks <- c(-4,0,4)
RPRS_limits <- c(-8,8)

## Function to make density plots (adjust scales according to annotations!)
plot_density <- function(data, var, breaks, limits, label, show_legend = FALSE) {
  ggplot(data, aes(x = .data[[var]], fill = PASS, color = PASS)) + #fill = chrom_var)) +
    geom_density(alpha = 0.6) +
	labs(title = paste(label), x = "", y = "Density") +
	#scale_y_continuous(expand=c(0,0), breaks=c(0,0.1,0.2,0.3), limits=c(0,0.4))+ # adjust to annotation!
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
QD_Tranche999 <- plot_density(Tranche999, "QD", QD_breaks, QD_limits, "QualByDepth [Tranche 99.9]", show_legend = TRUE)
# Plot FS
FS_Tranche999 <- plot_density(Tranche999, "log_FS", logFS_breaks, logFS_limits, "log(FisherStrand) [Tranche 99.9]", show_legend = FALSE)
# Plot RPRS
RPRS_Tranche999 <- plot_density(Tranche999, "ReadPosRankSum", RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.9]", show_legend = FALSE)

## Plot together
All_Plots <- cowplot::plot_grid(QD_Tranche999, 
								FS_Tranche999, 
								RPRS_Tranche999,
							    nrow=3, ncol=1, align="hv", scale = 0.95,
							    labels = "AUTO", label_size = 28)
#png("./Figures/All_Plots_INDEL_Tranche999_W03.png", width=3500, height=4500, res=300)
All_Plots
#dev.off()

table(Tranche999$PASS, Tranche999$PASS_HF)
#         FALSE   TRUE
#  FALSE  36115 455850
#  TRUE    7968 500067

# VQSR passing rate: 7968 + 500067 ~ 50.8 - 49.2 fail
# HF 455850 + 500067 ~ 95.6 - 4.4 fail
# Agreement 36115 + 500067 --> 536182 / 1000000 - 53.6% --> 46.4% disagreement!

