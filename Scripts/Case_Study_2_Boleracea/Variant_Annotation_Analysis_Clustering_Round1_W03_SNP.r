### R Script 
### Variant Calling in a wild population  of Brassica oleracea
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
setwd("C:/Users/00115569/OneDrive - UWA/GATK_Pipeline/Wild_Brassica/VCF/Round_1")

#####################################################################
############################ PART 1 #################################
#####################################################################

### Load and prepare data

# Function to load and prepare VCFs into R
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

# Read the two VCF sets
Clustering_Subset <- prep_VCF("./SUBSET/Samples_Bol_Wild_1-7.merged.raw.SNP.SUBSET.annotations.vcf.gz")
Truth_Set <- prep_VCF("./TRUTH_SET/Samples_Bol_Wild_1-7.TRUTH_SET.annotations.vcf.gz")

### Check distribution of both sets

## Convert Data
# Prepare Annos for Clustering Subset
Clustering_Subset_Annos <- Clustering_Subset %>%
	separate(VariantID, into = c("CHROM", "POS_REF"), sep = ":", remove = TRUE) %>%
	separate(POS_REF, into = c("POS", "REF_ALT"), sep = "_") %>%
	separate(REF_ALT, into = c("REF", "ALT"), sep = "/") %>%
	mutate(CHROM = as.factor(CHROM), POS = as.integer(POS), log_FS = log(FS + 1)) %>%
	arrange(CHROM, POS)

# Prepare Annos for Truth_Set
Truth_Set_Annos <- Truth_Set %>%
	separate(VariantID, into = c("CHROM", "POS_REF"), sep = ":", remove = TRUE) %>%
	separate(POS_REF, into = c("POS", "REF_ALT"), sep = "_") %>%
	separate(REF_ALT, into = c("REF", "ALT"), sep = "/") %>%
	mutate(CHROM = as.factor(CHROM), POS = as.integer(POS), log_FS = log(FS + 1))
		
## Create data frame for plotting SNP distribution

# Clustering Subset
Subset_SNP_Distribution <- data.frame(
  Chromosome = names(table(Clustering_Subset_Annos$CHROM)),
  Count = as.numeric(table(Clustering_Subset_Annos$CHROM))
)

# Truth Set
Truth_Set_SNP_Distribution <- data.frame(
  Chromosome = names(table(Truth_Set_Annos$CHROM)),
  Count = as.numeric(table(Truth_Set_Annos$CHROM))
)

# Color scheme
vir_19 <- viridis(n = 19)
epsilon <- 1e-6

# Set different scales for the annotations
QD_breaks <- c(0,10,20,30,40)
QD_limits <- c(0,50)
logDP_breaks <- c(0,2,4,6,8,10,12)
logDP_limits <- c(0,13)
logFS_breaks <- c(-1,0,2,4)
logFS_limits <- c(-2,6)
MQ_breaks <- c(0,10,20,30,40,50,60)
MQ_limits <- c(0,70)
SOR_breaks <- c(0,2,4,6)
SOR_limits <- c(0,8)
MQRS_breaks <- c(-60,-40,-20,0,20,40,60)
MQRS_limits <- c(-80,80)
RPRS_breaks <- c(-30,-20,-10,0,10,20,30)
RPRS_limits <- c(-35,35)

## Combine data frames for joint plotting
Clustering_Subset_Annos$Dataset <- "Clustering Subset"
Truth_Set_Annos$Dataset        <- "Truth Set"

Combined_Annos <- rbind(
  Clustering_Subset_Annos,
  Truth_Set_Annos
)

# Function to plot annotation distributions in both sets
plot_global_density <- function(data, var, threshold, breaks, limits, label, show_legend = TRUE) {
	ggplot(data, aes(x = !!sym(var), colour = Dataset, fill = Dataset)) +
	geom_density(alpha = 0.6, linewidth = 1) +
	geom_vline(xintercept = threshold, color = "red", linewidth = 1, linetype = "dashed") +
	labs(title ="",
		x = paste(label), 
		y = "Density") +
	scale_y_continuous(expand=c(0,0))+
	scale_x_continuous(expand=c(0,0), breaks=breaks, limits=limits)+
	scale_color_manual(
		name = "Dataset",
		values = c("#00BFC4", "#7CAE00"),
		labels = c("1M Clustering Subset", "13k Truth Set")
		)+
	scale_fill_manual(
		name = "Dataset",
		values = c("#00BFC4", "#7CAE00"),
		labels = c("1M Clustering Subset", "13k Truth Set")
		)+
	theme(
		plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
		axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
		axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
		legend.position = if (show_legend) c(0.75, 0.85) else "none",
		legend.justification = "center",
		legend.background = element_blank(),
		legend.key = element_blank(),
		legend.key.width  = unit(1.2, "lines"),
		legend.key.height = unit(1.0, "lines"),
		legend.title = element_blank(),
		legend.text  = element_text(size = 14,face = "bold"),
		plot.title.position = "plot")
}

# Plot QD context annotation for all sets
QD_Plot_Global <- plot_global_density(Combined_Annos, "QD", 2, QD_breaks, QD_limits, "QualByDepth [QD]", show_legend = TRUE)
FS_Plot_Global <- plot_global_density(Combined_Annos, "log_FS", 4.1, logFS_breaks, logFS_limits, "log(FisherStrand) [FS]", show_legend = FALSE)
SOR_Plot_Global <- plot_global_density(Combined_Annos, "SOR", 3.0, SOR_breaks, SOR_limits, "StrandOddsRatio [SOR]", show_legend = FALSE)
MQ_Plot_Global <- plot_global_density(Combined_Annos, "MQ", 40.0, MQ_breaks, MQ_limits, "MappingQuality [MQ]", show_legend = FALSE)
MQRS_Plot_Global <- plot_global_density(Combined_Annos, "MQRankSum", -12.5, MQRS_breaks, MQRS_limits, "MQRankSum [MQRS]", show_legend = FALSE)
RPRS_Plot_Global <- plot_global_density(Combined_Annos, "ReadPosRankSum", -8, RPRS_breaks, RPRS_limits, "ReadPosRankSum [RPRS]", show_legend = FALSE)

## Combine in one plot
Global_Distribution_Plots <- cowplot::plot_grid(QD_Plot_Global, FS_Plot_Global,
												SOR_Plot_Global, MQ_Plot_Global,
												MQRS_Plot_Global, RPRS_Plot_Global,
												nrow=3, ncol=2, align="hv", scale = 0.9,
												labels = c("AUTO"), label_size = 24)
#png("./Figures/Global_Distribution_Plots.png", width=4000, height=4500, res=300)
Global_Distribution_Plots
#dev.off()

### Variant Context Annotation Clustering

# Label according to hard filtering
Clustering_Subset <- Clustering_Subset %>%
						mutate(PASS = !(
								(QD < 2.0 & !is.na(QD)) | 
								(FS > 60.0 & !is.na(FS)) |
								(SOR > 3.0 & !is.na(SOR)) | 
								(MQ < 40.0 & !is.na(MQ)) |
								(MQRankSum < -12.5 & !is.na(MQRankSum)) |
								(ReadPosRankSum < -8 & !is.na(ReadPosRankSum))
									)
								)

## Mark true sites in the subsets
# Complete Subset
Clustering_Subset <- Clustering_Subset %>%
	mutate(row_id = seq_len(nrow(Clustering_Subset)),
		is_truth = VariantID %in% Truth_Set$VariantID)

table(Clustering_Subset$is_truth)
# FALSE   TRUE 
# 999439    561

table(Clustering_Subset$PASS)
# FALSE   TRUE 
#273463 726537

# Hard filtering passing rate is ~ 72 %

table(Clustering_Subset$is_truth, Clustering_Subset$PASS)
#       
#         FALSE   TRUE
#  FALSE 273458 725981
#  TRUE       5    556

# 556 (99 %) of true SNPs passed hard filtering

## Calculate average annotation metrics

# Summary for the Truth Set
Truth_Summary <- Truth_Set %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary
#     QD    FS  SOR MQ MQRankSum ReadPosRankSum
#1 15.53 2.007 0.718 60         0          0.218

# Check the total number of true variants that passed hard filtering
Truth_Set_Filtered <- Truth_Set %>%
						mutate(PASS = !(
							(QD < 2.0 & !is.na(QD)) | 
							(FS > 60.0 & !is.na(FS)) |
							(SOR > 3.0 & !is.na(SOR)) | 
							(MQ < 40.0 & !is.na(MQ)) |
							(MQRankSum < -12.5 & !is.na(MQRankSum)) |
							(ReadPosRankSum < -8 & !is.na(ReadPosRankSum))
							)
						)

table(Truth_Set_Filtered$PASS)
#FALSE  TRUE 
#  186 13310

Truth_Summary_PASS <- Truth_Set_Filtered %>%
	filter(PASS == TRUE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary_PASS
#     QD    FS   SOR MQ MQRankSum ReadPosRankSum
#1 15.59 1.994 0.716 60         0         0.2175

Truth_Summary_FAIL <- Truth_Set_Filtered %>%
	filter(PASS == FALSE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary_FAIL
#     QD     FS SOR MQ MQRankSum ReadPosRankSum
#1 6.03 6.755 3.0095 59.63         0          0.301

# Summary for Subset
Subset_Summary <- Clustering_Subset %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Summary
#    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#1 13.6 2.136 0.887 55.92         0          0.198

# Summary for Subset PASS
Subset_Pass_Summary <- Clustering_Subset %>%
	filter(PASS == TRUE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Pass_Summary
#     QD    FS   SOR   MQ MQRankSum ReadPosRankSum
#1 14.82 2.119 0.805 58.94         0          0.188

# Summary for Subset FAIL
Subset_Fail_Summary <- Clustering_Subset %>%
	filter(PASS == FALSE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Fail_Summary
#    QD     FS   SOR    MQ MQRankSum ReadPosRankSum
#1 8.85 2.262 1.88 37.1    -0.778          0.239

## Calculate relationships among context annotations

# Pearson's Correlations
cor_matrix_subset <- Clustering_Subset %>%
  select(QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum) %>%
  cor(use="pairwise.complete.obs")

## Plot the metrics across the subsets

# Subset data frame with new category for plotting
plot_Subset <- Clustering_Subset %>%
			mutate(Category = case_when(
			is_truth == TRUE ~ "TRUTH",
			PASS == TRUE ~ "PASS",
			PASS == FALSE ~ "FAIL",
			TRUE ~ "NA"))
table(plot_Subset$Category)
#  FAIL   PASS  TRUTH 
#273458 725981    561
		
#####################################################################
############################ PART 2 #################################
#####################################################################

## Principal Component Analysis (PCA)
set.seed(123)

# Scale the features (Z-score normalization)
scaled_Subset <- scale(Clustering_Subset[, c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")])
# Convert to df
scaled_Subset <- as.data.frame(scaled_Subset)
# Add row ids
scaled_Subset$row_id <- Clustering_Subset$row_id
# Add VariantIDs
scaled_Subset$VariantID <- Clustering_Subset$VariantID

# Label hard filtered variants in the subsets
Clustering_Subset <- Clustering_Subset %>%
			  mutate(hardfiltered = case_when(
				PASS == TRUE ~ "Good",
				PASS == FALSE ~ "Bad",
				TRUE ~ "NA"
				))
				
# Merge scaled data to original subsets
Clustering_Subset <- Clustering_Subset %>%
	left_join(scaled_Subset, by="row_id") %>%
	rename(VariantID = VariantID.x,
		   QD_scaled = QD.y,
		   FS_scaled = FS.y,
		   MQ_scaled = MQ.y,
		   SOR_scaled = SOR.y,
		   MQRankSum_scaled = MQRankSum.y,
		   ReadPosRankSum_scaled = ReadPosRankSum.y,
		   QD = QD.x,
		   FS = FS.x,
		   MQ = MQ.x,
		   SOR = SOR.x,
		   MQRankSum = MQRankSum.x,
		   ReadPosRankSum = ReadPosRankSum.x) %>%
	select(row_id, VariantID, is_truth, hardfiltered,
		   QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, QD_scaled, FS_scaled, SOR_scaled, MQ_scaled,
		   MQRankSum_scaled, ReadPosRankSum_scaled)

# Create a subset of subsets for PCA - aim at ~100k variants with complete annotation profiles
Subset_idx <- sample(1:nrow(Clustering_Subset), 130000)

PCA_Subset <- na.omit(Clustering_Subset[Subset_idx, ])

# Select scaled annotation columns for PCA
Annotations <- c("QD_scaled", "FS_scaled", "SOR_scaled", 
                 "MQ_scaled", "MQRankSum_scaled", "ReadPosRankSum_scaled")
		 
# Check NAs
colSums(is.na(PCA_Subset[, Annotations]))	 
					 
# Run PCA on normalized data (six annotations)
PCA_res6_Subset <- prcomp(PCA_Subset[, Annotations], 
				  center = FALSE, scale. = FALSE)

# Scree Plot
Scree_Plot_Subset <- fviz_eig(PCA_res6_Subset,
	addlabels=TRUE, hjust = 0.5, ylim = c(0,40),
	main = "Clustering Subset (Sampled: n = 106,427)",
	ggtheme = theme(
    plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
  )) +
  scale_y_continuous(expand=c(0,0))

#png("./Figures/Scree_Plot_6Annotations_W03.png", width=2500, height=2000, res=300)
Scree_Plot_Subset
#dev.off()

# Create loading plots
Biplot_6Annos_Subset_PC1_PC2 <- fviz_pca_var(PCA_res6_Subset, axes = c(1,2),
					col.var = "contrib", # Color by contributions to the PC
					gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
					repel = TRUE,     # Avoid text overlapping
					title = " ",
					ggtheme = theme(
						plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
						axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
						legend.position = "right",
						legend.title = element_text(face = "bold", size = 12),
						legend.text = element_text(size = 12)
				))

Biplot_6Annos_Subset_PC1_PC3 <- fviz_pca_var(PCA_res6_Subset, axes = c(1,3),
					col.var = "contrib", # Color by contributions to the PC
					gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
					repel = TRUE,     # Avoid text overlapping
					title = " ",
					ggtheme = theme(
						plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
						axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
						legend.position = "right",
						legend.title = element_text(face = "bold", size = 12),
						legend.text = element_text(size = 12)
				))

Biplot_6Annos_Subset_PC2_PC3 <- fviz_pca_var(PCA_res6_Subset, axes = c(2,3),
					col.var = "contrib", # Color by contributions to the PC
					gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
					repel = TRUE,     # Avoid text overlapping
					title = " ",
					ggtheme = theme(
						plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
						axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
						axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
						legend.position = "right",
						legend.title = element_text(face = "bold", size = 12),
						legend.text = element_text(size = 12)
				))


Biplots_6Annotations_Subset <- cowplot::plot_grid(Biplot_6Annos_Subset_PC1_PC2, Biplot_6Annos_Subset_PC2_PC3,
												      Biplot_6Annos_Subset_PC1_PC3,
													  nrow=3, ncol=1, align="hv", labels = c("A","B","C"),
													  label_size = 18, scale = 1)
#png("./Figures/Biplots_6Annotations_Subset_W03.png", width=2000, height=3500, res=300)
Biplots_6Annotations_Subset
#dev.off()

# Check significance for FS and SOR
cor.test(Clustering_Subset$FS, Clustering_Subset$SOR, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation

#data:  Clustering_Subset$FS and Clustering_Subset$SOR
#t = 647.19, df = 999998, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.5421679 1.0000000
#sample estimates:
#      cor 
#0.5433282

# Create color palette and heatmap for clustering subset
coul <- circlize::colorRamp2(c(1,0.5,0,-0.5,-1),c("#172869FF","#088BBEFF","#F1F4EEFF","#E9A17CFF","#803233FF"))
cor_matrix_heatmap_subset <- Heatmap(cor_matrix_subset, col=coul,
								  show_row_names=TRUE,
								  row_dend_width=unit(1,"cm"),
								  column_dend_height=unit(.5,"cm"),
								  heatmap_legend_param=list(
								  title=expression(bold("R2")),
								  title_position="leftcenter-rot",
								  legend_height=unit(6,"cm"),
								  labels_gp=gpar(fontsize=10,fontface="bold"),
								  labels=c("-1","-0.5","0","0.5","1"),
								  border="black"),
								  border=TRUE,
								  column_title="",
								  row_title=NULL,
								  use_raster=TRUE,
								  raster_quality=10,						  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", cor_matrix_subset[i, j]),  # format numbers (2 decimals)
      x, y,
      gp = gpar(fontsize = 8, fontface = "bold", col = "black")
    )
  }
)
cor_matrix_heatmap_subset

#####################################################################
############################ PART 3 #################################
#####################################################################

### Gaussian-Mixture Model-based Clustering

## GATK hard filter thresholds
GATK_Thresholds_List <- list(
  QD  = 2,
  FS  = 60,
  MQ  = 40,
  SOR = 3,
  MQRankSum = -12.5,
  ReadPosRankSum = -8
)

# Function to calculate HF violations
count_hf_violations <- function(df, thresholds) {
  df %>%
    rowwise() %>%
    mutate(
      HF_Violations = sum(
        QD < thresholds$QD,
        FS > thresholds$FS,
        MQ < thresholds$MQ,
        SOR > thresholds$SOR,
        MQRankSum < thresholds$MQRankSum,
        ReadPosRankSum < thresholds$ReadPosRankSum,
        na.rm = TRUE
      )
    ) %>%
    ungroup()
}

# Summarise HF violations row-wise and add to dataframe
Clustering_Subset <- count_hf_violations(Clustering_Subset, GATK_Thresholds_List)

### Clustering analysis

# Prepare the raw set
Clustering_Subset_Complete <- na.omit(Clustering_Subset)
#[1] 819534     17

## GMM-based Clustering

# Prepare the data 
Opt_Cluster_Subset <- Clustering_Subset_Complete %>%
  select(QD_scaled, FS_scaled, SOR_scaled, MQ_scaled, MQRankSum_scaled, ReadPosRankSum_scaled)

# Identify opimal cluster number
Opt_GMM_Subset = Optimal_Clusters_GMM(Opt_Cluster_Subset, max_clusters = 10, criterion = "BIC", 
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               km_iter = 10, em_iter = 100, var_floor = 1e-10, 
                               plot_data = T)
			   
# Create data frames for plotting
GMM_DF_Subset <- data.frame(
  clusters = 1:length(Opt_GMM_Subset),
  BIC = Opt_GMM_Subset
)

GMM_deltaBIC_Subset <- bind_rows(
  data.frame(Clusters = 1:length(Opt_GMM_Subset), 
             BIC = Opt_GMM_Subset) %>%
  mutate(Delta_BIC = BIC - min(BIC)))

# Create the delta BIC plot
DeltaBIC_SubsetPlot <- ggplot(GMM_deltaBIC_Subset, aes(x = Clusters, y = Delta_BIC)) +
  geom_line(linewidth = 1.1, color="#1F77B4") +
  geom_point(size = 2, color="#1F77B4") +
  geom_vline(xintercept = 3, linetype = "dotted", size = 1) +
  geom_text_repel(aes(label = scales::comma(Delta_BIC)),
                  size = 4,
				  nudge_x = 0.75,
				  #nudge_y = -2,
                  show.legend = FALSE,
                  max.overlaps = Inf,   # allow many labels
                  box.padding = 0.5,
                  point.padding = 0.2,
                  segment.color = "black") +											 
  scale_x_continuous(breaks = 1:10, limits = c(1, 11)) +
  scale_y_continuous(labels = scales::label_number(big.mark=","), 
					breaks = seq(0,15000000, 2500000),
					limits=c(0,max(GMM_deltaBIC_Subset$Delta_BIC))) +
  labs(#title = "GMM Model Selection via ΔBIC",
	   title = "",
       x = "Number of Components (Gaussians)",
       y = expression(Delta*"BIC (rel. to last component)")) +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12),
	plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
DeltaBIC_SubsetPlot

## Set up the model 

# three cluster
GMM_K3_Subset <- GMM(
	data = Opt_Cluster_Subset,
	gaussian_comps = 3,
	dist_mode = "maha_dist",
	seed_mode = "random_subset",
	km_iter = 10,
	em_iter = 100,
	var_floor = 1e-10
	)

# --- Predict cluster memberships ---

# for three cluster
pred_GMM_K3_Subset <- predict_GMM(Opt_Cluster_Subset, 
								   GMM_K3_Subset$centroids, 
								   GMM_K3_Subset$covariance_matrices,
								   GMM_K3_Subset$weights)

# --- Add predictions to the main dataset ---
Clustering_Subset_Complete <- Clustering_Subset_Complete %>%
  mutate(
    GMM_K3 = as.factor(pred_GMM_K3_Subset$cluster_labels),
  )

# take the same random subset that was used for PCA
PCA_GMM_Labelled_Subset <- PCA_Subset %>%
  left_join(Clustering_Subset_Complete %>% select(VariantID, starts_with("GMM")),
            by = "VariantID")

# Check cluster size
table(Clustering_Subset_Complete$GMM_K3)
#     1      2      3 
#226925 406688 185921

# Check where truth variants fall
table(Clustering_Subset_Complete$GMM_K3, Clustering_Subset_Complete$is_truth)
#     FALSE   TRUE
#  1 226530    395
#  2 406608     80
#  3 185900     21

# Summarize metrics for K3
K3_GMM_Subset_Summary <- Clustering_Subset_Complete %>%
  group_by(GMM_K3) %>%
  summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))  
print(K3_GMM_Subset_Summary)

# A tibble: 3 × 7
#  GMM_K3    QD     FS   SOR    MQ MQRankSum ReadPosRankSum
#  <fct>  <dbl>  <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 1      13.4   2.37 0.717  60       0              0.143
#2 2      11.6   2.03 0.732  53.8    -0.954          0.208
#3 3       7.07 16.3  2.4    50.2    -2.24           0.266

table(Clustering_Subset_Complete$hardfiltered)
#   Bad   Good 
#204358 615176

table(Clustering_Subset_Complete$hardfiltered, Clustering_Subset_Complete$is_truth)
#        FALSE   TRUE
#  Bad  204354      4
#  Good 614684    492

# Summarize metrics for hard filtered
HF_Subset_Summary <- Clustering_Subset_Complete %>%
  group_by(hardfiltered) %>%
  summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))  
print(HF_Subset_Summary)
# A tibble: 2 × 7
#  hardfiltered    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <chr>        <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 Bad            5.2  6.01 1.58   38.2    -0.786          0.239
#2 Good          12.9  3.02 0.764  59.0     0              0.188

# Label cluster
Clustering_Subset_Complete$GMM_K3 <- factor(
  Clustering_Subset_Complete$GMM_K3,
  levels = c(1, 2, 3),
  labels = c("Good", "Medium", "Bad")
)

## Add probabilities of GMM K3 to dataframe
Clustering_Subset_Complete$GMM_K3_Good <- pred_GMM_K3_Subset$cluster_proba[, 1]
Clustering_Subset_Complete$GMM_K3_Med <- pred_GMM_K3_Subset$cluster_proba[, 2]
Clustering_Subset_Complete$GMM_K3_Bad <- pred_GMM_K3_Subset$cluster_proba[, 3]

## Illustrate variant flow with alluvial Plots

# Prepare dataframes
Subset_Alluvial_GMM <- Clustering_Subset_Complete %>%
  select(VariantID, hardfiltered, GMM_K3) %>%
  mutate(
	hardfiltered = factor(hardfiltered, levels = c("Good", "Bad")),
	GMM_K3 = factor(GMM_K3, levels = c("Good", "Medium", "Bad"))
	) %>%
	dplyr::count(hardfiltered, GMM_K3, name = "freq")

Subset_Alluvial_GMM <- Subset_Alluvial_GMM %>%
  mutate(across(
    c(hardfiltered, GMM_K3),
    ~ factor(recode(.,
      "Good" = "High",
      "Medium" = "Intermediate",
      "Bad" = "Low"
    ),
    levels = c("High", "Intermediate", "Low"))
  ))


# Color palette
my_colors <- c(
  "High"   = "#0E7175",
  "Intermediate" = "#808BC5",
  "Low"    = "#ED773C"
)	

# Alluvial plot
Alluvial_Plot_Subset_Two_Comparison <- ggplot(
  Subset_Alluvial_GMM,
  aes(
    axis1 = hardfiltered,
    axis2 = GMM_K3,
    y     = freq
  )
) +
  geom_alluvium(aes(fill = GMM_K3), width = 1/12) +
  geom_stratum(aes(fill = after_stat(stratum)), color = "black") +
  #geom_text(
  #  stat = "stratum",
  #  aes(label = after_stat(stratum)),
  #  size = 4
  #) +
  scale_fill_manual(values = my_colors, name = "Variant Quality Class") +
  scale_x_discrete(
    limits = c("hardfiltered", "GMM_K3"),
    labels = c("Hard Filters", "GMM (k = 3)")
  ) +
  scale_y_continuous(
    labels = scales::label_number(big.mark=","),
    #breaks = c(500000, 1000000, 1500000),
    limits = c(0, 950000)
  ) +
  labs(
    x = "",
    y = "Number of variants"
  ) +
  theme_gray(base_size = 14) +
  theme(
    axis.title.x  = element_text(face = "bold", color="black"),
    axis.title.y  = element_text(face = "bold", color="black"),
    axis.text.x   = element_text(face = "bold", color="black"),
    axis.text.y   = element_text(face = "bold", color="black"),
    legend.position = c(0.5, 1),
	legend.direction = "horizontal",
	legend.justification = c(0.5, 1),
    legend.title = element_text(face = "bold", color="black", size = 11),
    legend.text  = element_text(size = 10, color="black"),
	legend.background = element_blank(),
	legend.box.background = element_blank(),
	plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

### plot best model fit and alluvial plots together
Model_Output <- cowplot::plot_grid(DeltaBIC_SubsetPlot,
									  Alluvial_Plot_Subset_Two_Comparison,
									  nrow=2, ncol=1, align = "hv",
									  scale = 0.99, labels="AUTO", label_size=24)
#png("./Figures/Model_Output.png", width=2500, height=3500, res=300)
Model_Output
#dev.off()

#take the same random subset that was used for PCA
Clustering_Subset_PCA <- PCA_Subset %>%
  left_join(Clustering_Subset_Complete %>% select(VariantID, GMM_K3),
            by = "VariantID") %>%
  mutate(hardfiltered = factor(
			recode(hardfiltered,
			"Good" = "PASS",
			"Bad" = "FAIL"
		),
		levels = c("PASS", "FAIL"))
	) %>%
  mutate(GMM_K3 = factor(
			recode(GMM_K3,
			"Good" = "High",
			"Medium" = "Intermediate",
			"Bad" = "Low"
		),
		levels = c("High", "Intermediate", "Low"))
	)

# Color palette
my_colors <- c("#0E7175","#ED773C")	

# Biplots for hard filtering
Biplot_Ind_Var_HF_PC1_PC2 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(1,2),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$hardfiltered,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = "Hard Filtering [n = 106,427]",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
	labs(color = "Hard Filtering",
		 shape = "Hard Filtering",
		 fill = "Hard Filtering") +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = c(0.75,0.2),
	legend.background = element_blank(),
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))

Biplot_Ind_Var_HF_PC1_PC3 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(1,3),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$hardfiltered,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = " ",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = "none",
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))

Biplot_Ind_Var_HF_PC2_PC3 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(2,3),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$hardfiltered,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = " ",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = "none",
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))	

# Biplots for k = 3 clustering

my_colors <- c(
  "High"   = "#0E7175",
  "Intermediate" = "#808BC5",
  "Low"    = "#ED773C"
)	

### GMM Clustering
Biplot_Ind_Var_K3_GMM_PC1_PC2 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(1,2),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$GMM_K3,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = "GMM [k = 3; n = 106,427]",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
	labs(color = "Variant Quality Class",
		 shape = "Variant Quality Class",
		 fill = "Variant Quality Class") +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = c(0.72,0.2),
	legend.background = element_blank(),
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))

Biplot_Ind_Var_K3_GMM_PC1_PC3 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(1,3),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$GMM_K3,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = " ",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = "none",
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))

Biplot_Ind_Var_K3_GMM_PC2_PC3 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(2,3),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$GMM_K3,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = " ",
	col.var="black",
    palette = my_colors,
	repel = TRUE,        # <-- Important!
    labelsize = 3,       # reduce label size
    pointsize = 1.5) +
theme_gray(base_size = 14) +
theme(
	plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
	axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
	axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
	legend.position = "none",
	legend.title = element_text(face = "bold", size = 14),
	legend.text = element_text(size = 14))	

## Plot all together in one common plot (vertically):
Biplots_Var_Ind_GMM <- cowplot::plot_grid(Biplot_Ind_Var_HF_PC1_PC2, Biplot_Ind_Var_K3_GMM_PC1_PC2,
										 Biplot_Ind_Var_HF_PC2_PC3, Biplot_Ind_Var_K3_GMM_PC2_PC3,
										 Biplot_Ind_Var_HF_PC1_PC3, Biplot_Ind_Var_K3_GMM_PC1_PC3,
										 nrow=3, ncol=2, align="hv", labels = c("A","D","B","E","C","F"),
										 label_size = 20)
#png("./Figures/Biplots_Var_Ind_HF_vs_GMM.png",  width=4000, height=5000, res=350)
Biplots_Var_Ind_GMM
#dev.off()

# check the correlation among the annotations and the probs
cor_vars <- Clustering_Subset_Complete %>%
  select(QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum, GMM_K3_Good, GMM_K3_Med, GMM_K3_Bad) %>%
  cor(use="pairwise.complete.obs")

cor.test(Clustering_Subset_Complete$FS, Clustering_Subset_Complete$GMM_K3_Bad, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation
#
#data:  GMM_Subset$FS and GMM_Subset$GMM_K3_Bad
#t = 610.55, df = 819532, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.5578974 1.0000000
#sample estimates:
#      cor 
#0.5591476

cor.test(Clustering_Subset_Complete$SOR, Clustering_Subset_Complete$GMM_K3_Bad, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation
#
#data:  GMM_Subset$SOR and GMM_Subset$GMM_K3_Bad
#t = 924.76, df = 819532, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.7137039 1.0000000
#sample estimates:
#      cor 
#0.7145942

coul <- circlize::colorRamp2(c(1,0.5,0,-0.5,-1),c("#172869FF","#088BBEFF","#F1F4EEFF","#E9A17CFF","#803233FF"))

cor_vars_heatmap <- Heatmap(cor_vars, col=coul,
								  show_row_names=TRUE,
								  row_dend_width=unit(4,"cm"),
								  column_dend_height=unit(1,"cm"),
								  heatmap_legend_param=list(
								  title=expression(bold("R2")),
								  title_position="leftcenter-rot",
								  legend_height=unit(6,"cm"),
								  labels_gp=gpar(fontsize=10,fontface="bold"),
								  labels=c("-1","-0.5","0","0.5","1"),
								  border="black"),
								  border=TRUE,
								  column_title="",
								  row_title=NULL,
								  use_raster=TRUE,
								  raster_quality=10,						  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", cor_vars[i, j]),  # format numbers (2 decimals)
      x, y,
      gp = gpar(fontsize = 8, fontface = "bold", col = "black")
    )
  }
)
cor_vars_heatmap


Clustering_Subset_Complete$GMM_K3 <- factor(
  Clustering_Subset_Complete$GMM_K3,
  levels = c("Good", "Medium", "Bad")
)

## Run a multinomial/ordinal logistic regression to identify the features 
model <- multinom(GMM_K3 ~ QD + FS + MQ + SOR + MQRankSum + ReadPosRankSum, data = Clustering_Subset_Complete, Hess=TRUE)
# weights:  24 (14 variable)
#initial  value 900350.123394 
#iter  10 value 534637.259984
#iter  20 value 293095.731217
#iter  30 value 222515.594118
#iter  40 value 212684.886922
#iter  50 value 184162.454887
#final  value 184162.349585 
#converged

summary(model)
#Call:
#multinom(formula = GMM_K3 ~ QD + FS + MQ + SOR + MQRankSum + 
#    ReadPosRankSum, data = GMM_Subset, Hess = TRUE)

#Coefficients:
#       (Intercept)         QD         FS        MQ        SOR  MQRankSum ReadPosRankSum
#Medium    693.6063 0.03747538 -0.2034531 -11.58408 -0.9155304 -0.4054077     0.04248756
#Bad       687.8926 0.02808733  0.2518723 -11.63058  2.1111230 -0.3766939     0.01950559

#Std. Errors:
#       (Intercept)          QD          FS           MQ        SOR   MQRankSum ReadPosRankSum
#Medium  0.02079697 0.000725175 0.001585728 0.0004289557 0.01857107 0.008162318    0.005492842
#Bad     0.02080873 0.001036009 0.001767284 0.0006256066 0.01890414 0.008392811    0.007272866

#Residual Deviance: 368324.7 
#AIC: 368352.7 

exp(coef(model))
#         (Intercept)        QD       FS          MQ         SOR MQRankSum ReadPosRankSum
#Medium 1.695886e+301 1.038186 0.8159085 9.313166e-06 0.4003043 0.6667049       1.043403
#Bad    5.597118e+298 1.028485 1.2864317 8.890035e-06 8.2575094 0.6861260       1.019697



#####################################################################
############################ PART 4 #################################
#####################################################################

### Constructing resources for VQSR

## Stratify good variants
good_strat1_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.0 & GMM_K3_Good <= 0.25)
good_strat2_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.25 & GMM_K3_Good <= 0.5)
good_strat3_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.5 & GMM_K3_Good <= 0.6)
good_strat4_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.6 & GMM_K3_Good <= 0.7)
good_strat5_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.7 & GMM_K3_Good <= 0.8)
good_strat6_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.8 & GMM_K3_Good <= 0.9)
good_strat7_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Good" &
						  GMM_K3_Good >= 0.9)

## Check metrics
good_strat1_summary	<- good_strat1_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat1_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1     0    NA    NA    NA    NA        NA             NA
sum(good_strat1_prob$HF_Violations)
#[1] 0

good_strat2_summary	<- good_strat2_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat2_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1    191  13.0  241.  8.81  43.8     -6.14          0.309
sum(good_strat2_prob$HF_Violations)
#[1] 427

good_strat3_summary	<- good_strat3_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat3_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1   35  6.31  19.7  2.96    60         0          0.325
sum(good_strat3_prob$HF_Violations)
#[1] 17

good_strat4_summary	<- good_strat4_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat4_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1   40  8.40  23.0  2.75    60         0         -0.096
sum(good_strat4_prob$HF_Violations)
#[1] 18

good_strat5_summary	<- good_strat5_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat5_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1   40   7.7  21.9  2.78    60         0          0.143
sum(good_strat5_prob$HF_Violations)
#[1] 18

good_strat6_summary	<- good_strat6_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat6_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1   76  9.62  21.0  2.80    60         0          0.198
sum(good_strat6_prob$HF_Violations)
#[1] 31

good_strat7_summary	<- good_strat7_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(good_strat7_summary)
# A tibble: 1 × 7
#       n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#   <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 226543  13.4  2.36 0.716    60         0          0.143
sum(good_strat7_prob$HF_Violations)
#[1] 2990

# Prepare GATK thresholds
GATK_Thresholds <- data.frame(
  Metric = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"),
  Threshold = c(2.0, 60.0, 3.0, 40.0, -12.5, -8.0)
)

GATK_Thresholds$Metric <- factor(GATK_Thresholds$Metric, levels = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"))

## Analyse variant context annotations across the medium group
Good_Strats <- bind_rows(
  good_strat1_prob %>% mutate(Good_Prob = "0.0–0.25"),
  good_strat2_prob %>% mutate(Good_Prob = "0.25–0.5"),
  good_strat3_prob %>% mutate(Good_Prob = "0.5–0.6"),
  good_strat4_prob %>% mutate(Good_Prob = "0.6–0.7"),
  good_strat5_prob %>% mutate(Good_Prob = "0.7–0.8"),
  good_strat6_prob %>% mutate(Good_Prob = "0.8–0.9"),
  good_strat7_prob %>% mutate(Good_Prob = "0.9–1.0")
) %>%
pivot_longer(c(QD,FS,MQ,SOR,MQRankSum,ReadPosRankSum),
             names_to="Metric", values_to="Value")
Good_Strats$Metric <- factor(
  Good_Strats$Metric,
  levels = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")
)

Annotation_Metrics_Good_Probs_Stratification <- ggplot(Good_Strats, aes(Good_Prob, Value, fill=Good_Prob)) +
  geom_boxplot(outlier.size=0.3) +
  geom_hline(data = GATK_Thresholds,
             aes(yintercept = Threshold, group = Metric),
             color = "red",
             linetype = "dashed",
			 linewidth = 1) +
  facet_wrap(~Metric, scales="free") +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 9, face = "bold"),
	strip.text = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
	)
Annotation_Metrics_Good_Probs_Stratification


## Stratify medium variants
medium_strat1_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.0 & GMM_K3_Med <= 0.25)
medium_strat2_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.25 & GMM_K3_Med <= 0.5)
medium_strat3_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.5 & GMM_K3_Med <= 0.6)
medium_strat4_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.6 & GMM_K3_Med <= 0.7)
medium_strat5_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.7 & GMM_K3_Med <= 0.8)
medium_strat6_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.8 & GMM_K3_Med <= 0.9)
medium_strat7_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Medium" &
						  GMM_K3_Med >= 0.9)
## Check metrics
medium_strat1_summary	<- medium_strat1_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat1_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1     0    NA    NA    NA    NA        NA             NA
sum(medium_strat1_prob$HF_Violations)
#[1] 0

medium_strat2_summary	<- medium_strat2_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat2_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1   136  21.2  4.12 0.718  54.3     0.367          0.627
sum(medium_strat2_prob$HF_Violations)
#[1] 112

medium_strat3_summary	<- medium_strat3_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat3_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  5037  9.12  9.03  1.47  51.1     -1.39          0.218
sum(medium_strat3_prob$HF_Violations)
#[1] 1718

medium_strat4_summary	<- medium_strat4_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat4_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  5594  9.02  8.77  1.43  51.9     -1.65          0.234
sum(medium_strat4_prob$HF_Violations)
#[1] 1718

medium_strat5_summary	<- medium_strat5_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat5_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  7961  8.55  8.15  1.41  51.5     -1.82          0.219
sum(medium_strat5_prob$HF_Violations)
#[1]  2562

medium_strat6_summary	<- medium_strat6_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat6_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 16260  8.72  7.16  1.38  50.8     -1.39          0.212
sum(medium_strat6_prob$HF_Violations)
#[1] 6042

medium_strat7_summary	<- medium_strat7_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(medium_strat7_summary)
# A tibble: 1 × 7
#       n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#   <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 371700  11.9  1.78  0.71  54.0    -0.911          0.207
sum(medium_strat7_prob$HF_Violations)
#[1] 93886

## Analyse variant context annotations across the medium group
Medium_Strats <- bind_rows(
  medium_strat1_prob %>% mutate(Medium_Prob = "0.0–0.25"),
  medium_strat2_prob %>% mutate(Medium_Prob = "0.25–0.5"),
  medium_strat3_prob %>% mutate(Medium_Prob = "0.5–0.6"),
  medium_strat4_prob %>% mutate(Medium_Prob = "0.6–0.7"),
  medium_strat5_prob %>% mutate(Medium_Prob = "0.7–0.8"),
  medium_strat6_prob %>% mutate(Medium_Prob = "0.8–0.9"),
  medium_strat7_prob %>% mutate(Medium_Prob = "0.9–1.0")
) %>%
pivot_longer(c(QD,FS,MQ,SOR,MQRankSum,ReadPosRankSum),
             names_to="Metric", values_to="Value")
Medium_Strats$Metric <- factor(
  Medium_Strats$Metric,
  levels = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")
)

Annotation_Metrics_Medium_Probs_Stratification <- ggplot(Medium_Strats, aes(Medium_Prob, Value, fill=Medium_Prob)) +
  geom_boxplot(outlier.size=0.3) +
  geom_hline(data = GATK_Thresholds,
             aes(yintercept = Threshold, group = Metric),
             color = "red",
             linetype = "dashed",
			 linewidth = 1) +
  facet_wrap(~Metric, scales="free") +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
	strip.text = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
	)
Annotation_Metrics_Medium_Probs_Stratification

## Stratify bad variants
bad_strat1_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.0 & GMM_K3_Bad <= 0.25)
bad_strat2_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.25 & GMM_K3_Bad <= 0.5)
bad_strat3_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.5 & GMM_K3_Bad <= 0.6)
bad_strat4_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.6 & GMM_K3_Bad <= 0.7)
bad_strat5_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.7 & GMM_K3_Bad <= 0.8)
bad_strat6_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.8 & GMM_K3_Bad <= 0.9)
bad_strat7_prob <- Clustering_Subset_Complete %>%
				   filter(GMM_K3 == "Bad" &
						  GMM_K3_Bad >= 0.9)

## Check metrics
bad_strat1_summary	<- bad_strat1_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat1_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1     0    NA    NA    NA    NA        NA             NA
sum(bad_strat1_prob$HF_Violations)
#[1] 0

bad_strat2_summary	<- bad_strat2_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat2_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  3363  10.5  111.  6.28  47.6     -4.73          0.375
sum(bad_strat2_prob$HF_Violations)
#[1] 6666

bad_strat3_summary	<- bad_strat3_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat3_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  5111  9.23  9.86  1.55  50.9     -1.45          0.241
sum(bad_strat3_prob$HF_Violations)
#[1] 2459

bad_strat4_summary	<- bad_strat4_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat4_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  5415  9.01  10.2  1.59  51.1     -1.39          0.234
sum(bad_strat4_prob$HF_Violations)
#[1] 2553

bad_strat5_summary	<- bad_strat5_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat5_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  6369  9.03  10.2  1.70  50.4     -1.36          0.234
sum(bad_strat5_prob$HF_Violations)
#[1] 3290
						
bad_strat6_summary	<- bad_strat6_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)
print(bad_strat6_summary)
# A tibble: 1 × 7
#      n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1  8600  9.07  11.2  1.72  50.8     -1.60          0.220
sum(bad_strat6_prob$HF_Violations)
#[1] 4641

bad_strat7_summary	<- bad_strat7_prob %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)						
print(bad_strat7_summary)
# A tibble: 1 × 7
#       n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#   <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 157063  6.73  17.5  2.58  50.2     -2.33           0.27
sum(bad_strat7_prob$HF_Violations)
#[1] 124248

## Analyse variant context annotations across the medium group
Bad_Strats <- bind_rows(
  bad_strat1_prob %>% mutate(Bad_Prob = "0.0–0.25"),
  bad_strat2_prob %>% mutate(Bad_Prob = "0.25–0.5"),
  bad_strat3_prob %>% mutate(Bad_Prob = "0.5–0.6"),
  bad_strat4_prob %>% mutate(Bad_Prob = "0.6–0.7"),
  bad_strat5_prob %>% mutate(Bad_Prob = "0.7–0.8"),
  bad_strat6_prob %>% mutate(Bad_Prob = "0.8–0.9"),
  bad_strat7_prob %>% mutate(Bad_Prob = "0.9–1.0")
) %>%
pivot_longer(c(QD,FS,MQ,SOR,MQRankSum,ReadPosRankSum),
             names_to="Metric", values_to="Value")
Bad_Strats$Metric <- factor(
  Bad_Strats$Metric,
  levels = c("QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum")
)

Annotation_Metrics_Bad_Probs_Stratification <- ggplot(Bad_Strats, aes(Bad_Prob, Value, fill=Bad_Prob)) +
  geom_boxplot(outlier.size=0.3) +
  geom_hline(data = GATK_Thresholds,
             aes(yintercept = Threshold, group = Metric),
             color = "red",
             linetype = "dashed",
			 linewidth = 1) +
  facet_wrap(~Metric, scales="free") +
  theme_gray(base_size = 14) +
  theme(
    plot.title = element_text(colour = "black", size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 9, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
	strip.text = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 12)
	)
Annotation_Metrics_Bad_Probs_Stratification


### Prepare the training data set

## Filter good variants
VQSR_training_good <- Clustering_Subset_Complete %>%
	filter(GMM_K3 == "Good" & GMM_K3_Good >= 0.9)

## Filter medium variants
VQSR_training_medium <- Clustering_Subset_Complete %>%
	filter(GMM_K3 == "Medium" & GMM_K3_Med >= 0.9)

## Summary: good
VQSR_training_good_summary	<- VQSR_training_good %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)						
print(VQSR_training_good_summary)
# A tibble: 1 × 7
#       n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#   <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 226543  13.4  2.36 0.716    60         0          0.143

## Summary: medium
VQSR_training_medium_summary	<- VQSR_training_medium %>%
						summarise(
						n = n(),
						across(c(QD,FS,SOR,MQ,MQRankSum, ReadPosRankSum),
						\(x) median(x, na.rm = TRUE))
						)						
print(VQSR_training_medium_summary)
# A tibble: 1 × 7
#       n    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#   <int> <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 371700  11.9  1.78  0.71  54.0    -0.911          0.207

# Create a simple table for GATK to select good variants
VQSR_training_variants_good <- VQSR_training_good %>%
  select(VariantID)

# create BED file for good
VQSR_training_good_coords <- VQSR_training_variants_good %>%
  separate(VariantID, into=c("chr","pos_ref_alt"), sep=":") %>%
  separate(pos_ref_alt, into=c("pos","ref_alt"), sep="_") %>%
  select(chr,pos) %>%
  arrange(chr, as.numeric(pos))

# Save as BED
write.table(VQSR_training_good_coords, "./TRAINING/VQSR_training_good_W03.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Create a simple table for GATK to select medium variants
VQSR_training_variants_medium <- VQSR_training_medium %>%
  select(VariantID)

# create BED file for good
VQSR_training_medium_coords <- VQSR_training_variants_medium %>%
  separate(VariantID, into=c("chr","pos_ref_alt"), sep=":") %>%
  separate(pos_ref_alt, into=c("pos","ref_alt"), sep="_") %>%
  select(chr,pos) %>%
  arrange(chr, as.numeric(pos))
  

# Save as BED
write.table(VQSR_training_medium_coords, "./TRAINING/VQSR_training_medium_W03.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
