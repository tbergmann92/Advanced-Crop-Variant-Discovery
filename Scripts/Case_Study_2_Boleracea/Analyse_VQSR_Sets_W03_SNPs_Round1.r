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

Tranche999 <- prep_VCF("Samples_Bol_Wild_1-7.SNPs.VQSR.99.90.tranche.round1.subset.annotations.vcf.gz")

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
#266998 733002

#Passing rate: 73.3

table(Tranche999$PASS_HF)
# FALSE   TRUE 
#272981 727019

#Passing rate: 72.7

# Plot distribution

# Set different scales 
QD_breaks <- c(0,10,20,30,40)
QD_limits <- c(0,45)
logFS_breaks <- c(-1,0,2,4,6)
logFS_limits <- c(-2,6.5)
MQ_breaks <- c(20,30,40,50,60)
MQ_limits <- c(18,65)
SOR_breaks <- c(0,2,4,6)
SOR_limits <- c(0,8)
MQRS_breaks <- c(-8,-4,0,4)
MQRS_limits <- c(-10,6)
RPRS_breaks <- c(-4,-2,0,2,4)
RPRS_limits <- c(-5,5)

## Function to make density plots (adjust scales according to annotations!)
plot_density <- function(data, var, breaks, limits, label, show_legend = FALSE) {
  ggplot(data, aes(x = .data[[var]], fill = PASS, color = PASS)) + #fill = chrom_var)) +
    geom_density(alpha = 0.6) +
	labs(title = paste(label), x = "", y = "Density") +
	#scale_y_continuous(expand=c(0,0), breaks=c(0,0.2,0.4), limits=c(0,0.4))+ # adjust to annotation!
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
# Plot SOR
SOR_Tranche999 <- plot_density(Tranche999, "SOR", SOR_breaks, SOR_limits, "StrandOddsRatio [Tranche 99.9]", show_legend = FALSE)
# Plot MQ
MQ_Tranche999 <- plot_density(Tranche999, "MQ", MQ_breaks, MQ_limits, "MappingQuality [Tranche 99.9]", show_legend = FALSE)
# Plot MQRS
MQRS_Tranche999 <- plot_density(Tranche999, "MQRankSum", MQRS_breaks, MQRS_limits, "MQRankSum [Tranche 99.9]", show_legend = FALSE)
# Plot RPRS
RPRS_Tranche999 <- plot_density(Tranche999, "ReadPosRankSum", RPRS_breaks, RPRS_limits, "ReadPosRankSum [Tranche 99.9]", show_legend = FALSE)

## Plot together
All_Plots <- cowplot::plot_grid(QD_Tranche999, FS_Tranche999,
								SOR_Tranche999, MQ_Tranche999,
								MQRS_Tranche999, RPRS_Tranche999,
							    nrow=3, ncol=2, align="hv", scale = 0.95,
							    labels = "AUTO", label_size = 28)
#png("./Figures/All_Plots_Tranche999_W03.png", width=4000, height=4000, res=300)
All_Plots
#dev.off()

table(Tranche999$PASS, Tranche999$PASS_HF)
#       
#         FALSE   TRUE
#  FALSE 157471 109527
#  TRUE  115510 617492

# VQSR passing rate: 617492 + 115510 ~ 73.3 - 26.7 fail
# HF 617492 + 109527 ~ 72.7 - 27.3 fail
# Agreement 157471 + 617492 --> 774963 / 1000000 - 77.5% --> 22.5% disagreement!

