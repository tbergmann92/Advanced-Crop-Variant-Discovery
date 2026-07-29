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
setwd("PATH/TO/VCF/Round_1")

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
Clustering_Subset <- prep_VCF("./SUBSET/Samples_Bna_Minicore_1-135.subset.raw.SNPs.annotations.vcf.gz")
Truth_Set <- prep_VCF("./TRUTH_SET/Samples_Bna_Minicore_1-135.TRUTH_SET.annotations.vcf.gz")

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

# Set different scales for the annotations
QD_breaks <- c(0,10,20,30,40)
QD_limits <- c(0,45)
logFS_breaks <- c(-1,0,2,4)
logFS_limits <- c(-2,6)
MQ_breaks <- c(20,30,40,50,60)
MQ_limits <- c(18,65)
SOR_breaks <- c(0,2,4,6)
SOR_limits <- c(0,8)
MQRS_breaks <- c(-8,-4,0,4,8)
MQRS_limits <- c(-8.2,8.2)
RPRS_breaks <- c(-4,-2,0,2,4)
RPRS_limits <- c(-4.5,4.5)

## Combine data frames for joint plotting
Clustering_Subset_Annos$Dataset <- "Clustering Subset"
Truth_Set_Annos$Dataset        <- "Truth Set"

Combined_Annos <- rbind(
  Clustering_Subset_Annos,
  Truth_Set_Annos
)

# Function to plot global distribution in both sets
plot_global_density <- function(data, var, breaks, limits, label, show_legend = TRUE) {
	ggplot(data, aes(x = !!sym(var), colour = Dataset, fill = Dataset)) +
	geom_density(alpha = 0.6, linewidth = 1) +
	#labs(title = expression("Mini-Core Collection of " * italic("B. napus")), 
	labs(title ="",
		x = paste(label), 
		y = "Density") +
	scale_y_continuous(expand=c(0,0))+
	scale_x_continuous(expand=c(0,0), breaks=breaks, limits=limits)+
	scale_color_manual(
		name = "Dataset",
		values = c("#00BFC4", "#7CAE00"),
		labels = c("1M Clustering Subset", "33k Truth Set")
		)+
	scale_fill_manual(
		name = "Dataset",
		values = c("#00BFC4", "#7CAE00"),
		labels = c("1M Clustering Subset", "33k Truth Set")
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
		legend.text  = element_text(size = 13,face = "bold"),
		plot.title.position = "plot")
}

# Plot QD context annotation for all sets
QD_Plot_Global <- plot_global_density(Combined_Annos, "QD", QD_breaks, QD_limits, "QualByDepth [QD]", show_legend = FALSE)
FS_Plot_Global <- plot_global_density(Combined_Annos, "log_FS", logFS_breaks, logFS_limits, "log(FisherStrand) [FS]", show_legend = TRUE)
SOR_Plot_Global <- plot_global_density(Combined_Annos, "SOR", SOR_breaks, SOR_limits, "StrandOddsRatio [SOR]", show_legend = FALSE)
MQ_Plot_Global <- plot_global_density(Combined_Annos, "MQ", MQ_breaks, MQ_limits, "MappingQuality [MQ]", show_legend = FALSE)
MQRS_Plot_Global <- plot_global_density(Combined_Annos, "MQRankSum", MQRS_breaks, MQRS_limits, "MQRankSum [MQRS]", show_legend = FALSE)
RPRS_Plot_Global <- plot_global_density(Combined_Annos, "ReadPosRankSum", RPRS_breaks, RPRS_limits, "ReadPosRankSum [RPRS]", show_legend = FALSE)

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
# 998613   1387

table(Clustering_Subset$PASS)
# FALSE   TRUE 
#358147 641853

# Hard filtering passing rate is ~ 64.2 %

table(Clustering_Subset$is_truth, Clustering_Subset$PASS)
#       
#         FALSE   TRUE
#  FALSE 358132 640481
#  TRUE      15   1372

# 1372 (98.9 %) of true SNPs passed hard filtering


## Calculate average annotation metrics

# Summary for the Truth Set
Truth_Summary <- Truth_Set %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary
#     QD    FS  SOR MQ MQRankSum ReadPosRankSum
#1 23.41 0.604 0.67 60         0          0.101

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
#  305 33068

Truth_Summary_PASS <- Truth_Set_Filtered %>%
	filter(PASS == TRUE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary_PASS
#     QD    FS   SOR MQ MQRankSum ReadPosRankSum
#1 23.42 0.577 0.669 60         0          0.101

Truth_Summary_FAIL <- Truth_Set_Filtered %>%
	filter(PASS == FALSE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Truth_Summary_FAIL
#     QD     FS SOR MQ MQRankSum ReadPosRankSum
#1 19.66 51.258 3.535 55.72    -1.983          0.117

# Summary for Subset
Subset_Summary <- Clustering_Subset %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Summary
#    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#1 15.05 3.275 1.036 54.7    -1.741          0.112

# Summary for Subset PASS
Subset_Pass_Summary <- Clustering_Subset %>%
	filter(PASS == TRUE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Pass_Summary
#     QD    FS   SOR   MQ MQRankSum ReadPosRankSum
#1 19.7 2.005 0.811 56.62    -1.523          0.131

# Summary for Subset FAIL
Subset_Fail_Summary <- Clustering_Subset %>%
	filter(PASS == FALSE) %>%
	summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))
Subset_Fail_Summary
#    QD     FS   SOR    MQ MQRankSum ReadPosRankSum
#1 5.17 13.467 3.207 42.7    -2.102          0.047

## Calculate relationships among context annotations

# Pearson's Correlations
cor_matrix_subset <- Clustering_Subset %>%
  select(QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum) %>%
  cor(use="pairwise.complete.obs")
		
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
Subset_idx <- sample(1:nrow(Clustering_Subset), 115000)

PCA_Subset <- na.omit(Clustering_Subset[Subset_idx, ])

# Select scaled annotation columns for PCA
Annotations <- c("QD_scaled", "FS_scaled", "SOR_scaled", 
                 "MQ_scaled", "MQRankSum_scaled", "ReadPosRankSum_scaled")
	 
# Check NAs
colSums(is.na(PCA_Subset[, Annotations]))	 
					 
# Run PCA on normalized data (six annotations)
PCA_res6_Subset <- prcomp(PCA_Subset[, Annotations], 
				  center = FALSE, scale. = FALSE)

# Scree Plots
Scree_Plot_Subset <- fviz_eig(PCA_res6_Subset,
	addlabels=TRUE, hjust = 0.5, ylim = c(0,40),
	main = "Clustering Subset (Sampled: n = 103,185)",
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

#png("./Figures/Scree_Plot_6Annotations_Darmor.png", width=2500, height=2000, res=300)
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
#png("./Figures/Biplots_6Annotations_Subset.png", width=2000, height=3500, res=300)
Biplots_6Annotations_Subset
#dev.off()


# Check significance for FS and SOR
cor.test(Clustering_Subset$FS, Clustering_Subset$SOR, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation

#data:  Clustering_Subset$FS and Clustering_Subset$SOR
#t = 1071.1, df = 999998, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.7301895 1.0000000
#sample estimates:
#      cor 
#0.7309564

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
#[1] 897179     17

## GMM-based Clustering

# Prepare the data 
Opt_Cluster_Subset <- Clustering_Subset_Complete %>%
  select(QD_scaled, FS_scaled, SOR_scaled, MQ_scaled, MQRankSum_scaled, ReadPosRankSum_scaled)

# Identify opimal cluster number
Opt_Clustering_Subset_Complete = Optimal_Clusters_GMM(Opt_Cluster_Subset, max_clusters = 10, criterion = "BIC", 
                               dist_mode = "maha_dist", seed_mode = "random_subset",
                               km_iter = 10, em_iter = 100, var_floor = 1e-10, 
                               plot_data = T)
			   
# Create data frames for plotting
GMM_DF_Subset <- data.frame(
  clusters = 1:length(Opt_Clustering_Subset_Complete),
  BIC = Opt_Clustering_Subset_Complete
)

GMM_deltaBIC_Subset <- bind_rows(
  data.frame(Clusters = 1:length(Opt_Clustering_Subset_Complete), 
             BIC = Opt_Clustering_Subset_Complete) %>%
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
  scale_y_continuous(labels = scales::label_number(big.mark=","), breaks = c(0, 2500000, 5000000, 7500000, 10000000), limits=c(0,11500000)) +
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

## Set up the models and run them across the subsets 

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
    GMM_K3 = as.factor(pred_GMM_K3_Subset$cluster_labels)
  )

# take the same random subset that was used for PCA
PCA_GMM_Labelled_Subset <- PCA_Subset %>%
  left_join(Clustering_Subset_Complete %>% select(VariantID, starts_with("GMM")),
            by = "VariantID")

table(Clustering_Subset_Complete$GMM_K3)
#     1      2      3 
#338383 150403 408393

table(Clustering_Subset_Complete$GMM_K3, Clustering_Subset_Complete$is_truth)
#     FALSE   TRUE
#  1 338338     45
#  2 149251   1152
#  3 408215    178

# Summarize metrics for K3
K3_Clustering_Subset_Complete_Summary <- Clustering_Subset_Complete %>%
  group_by(GMM_K3) %>%
  summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))  
print(K3_Clustering_Subset_Complete_Summary)
# A tibble: 3 × 7
#  GMM_K3    QD     FS   SOR    MQ MQRankSum ReadPosRankSum
#  <fct>  <dbl>  <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 1       5.69 25.9   3.08   53       -2.98          0.054
#2 2      22.7   0.908 0.677  60        0             0.101
#3 3      14.2   1.83  0.766  53.6     -1.93          0.15

table(Clustering_Subset_Complete$hardfiltered)
#   Bad   Good 
#310557 586622

table(Clustering_Subset_Complete$hardfiltered, Clustering_Subset_Complete$is_truth)
#        FALSE   TRUE
#  Bad  310544     13
#  Good 585260   1362

# Summarize metrics for hard filtered
HF_Subset_Summary <- Clustering_Subset_Complete %>%
  group_by(hardfiltered) %>%
  summarise(across(QD:ReadPosRankSum, \(x) median(x, na.rm = TRUE)))  
print(HF_Subset_Summary)
# A tibble: 2 × 7
#  hardfiltered    QD    FS   SOR    MQ MQRankSum ReadPosRankSum
#  <chr>        <dbl> <dbl> <dbl> <dbl>     <dbl>          <dbl>
#1 Bad            4   19.6  3.26   46.6     -2.11          0.047
#2 Good          18.2  2.51 0.791  56.6     -1.52          0.131


# Label clusters
Clustering_Subset_Complete$GMM_K3 <- factor(
  Clustering_Subset_Complete$GMM_K3,
  levels = c(1, 2, 3),
  labels = c("Bad", "Good", "Medium")
)


## Add probabilities of GMM K3 to dataframe
Clustering_Subset_Complete$GMM_K3_Good <- pred_GMM_K3_Subset$cluster_proba[, 2]
Clustering_Subset_Complete$GMM_K3_Med <- pred_GMM_K3_Subset$cluster_proba[, 3]
Clustering_Subset_Complete$GMM_K3_Bad <- pred_GMM_K3_Subset$cluster_proba[, 1]

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
									  nrow=1, ncol=2, align = "hv",
									  scale = 0.99, labels="AUTO", label_size=24)
#png("./Figures/Model_Output.png", width=4500, height=2000, res=300)
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


# Biplots for k = 3 clustering

my_colors <- c(
  "High"   = "#0E7175",
  "Intermediate" = "#808BC5",
  "Low"    = "#ED773C"
)	

### ### GMM Clustering
Biplot_Ind_Var_K3_GMM_PC1_PC2 <- fviz_pca_biplot(PCA_res6_Subset, axes = c(1,2),
    geom.ind = "point",
    habillage = Clustering_Subset_PCA$GMM_K3,
    addEllipses = TRUE,
	ellipse.level = 0.95,
	title = "GMM [k = 3; n = 103,185]",
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

## Plot all together in one common plot (horizontally):
Biplots_Var_Ind_HF_vs_GMM_Horizontal <- cowplot::plot_grid(
										 Biplot_Ind_Var_K3_GMM_PC1_PC2,
										 Biplot_Ind_Var_K3_GMM_PC2_PC3, 
										 Biplot_Ind_Var_K3_GMM_PC1_PC3,
										 nrow=3, ncol=1, align="hv", labels = c("A","B","C"),
										 label_size = 20)
#png("./Figures/Biplots_Var_Ind_GMM.png",  width=2500, height=4500, res=350)
Biplots_Var_Ind_HF_vs_GMM_Horizontal
#dev.off()


# check the correlation among the annotations and the probs
cor_vars <- Clustering_Subset_Complete  %>%
  select(QD, FS, MQ, SOR, MQRankSum, ReadPosRankSum, GMM_K3_Good, GMM_K3_Med, GMM_K3_Bad) %>%
  cor(use="pairwise.complete.obs")

cor.test(Clustering_Subset_Complete$FS, Clustering_Subset_Complete$GMM_K3_Bad, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation
#
#data:  Clustering_Subset_Complete$FS and Clustering_Subset_Complete$GMM_K3_Bad
#t = 641.72, df = 897177, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.559698 1.000000
#sample estimates:
#      cor 
#0.5608894

cor.test(Clustering_Subset_Complete$SOR, Clustering_Subset_Complete$GMM_K3_Bad, alternative="greater", method="pearson", conf.level = 0.95)
#        Pearson's product-moment correlation
#
#data:  Clustering_Subset_Complete$SOR and Clustering_Subset_Complete$GMM_K3_Bad
#t = 935.83, df = 897177, p-value < 2.2e-16
#alternative hypothesis: true correlation is greater than 0
#95 percent confidence interval:
# 0.7019481 1.0000000
#sample estimates:
#      cor 
#0.7028279

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
#initial  value 985651.874551 
#iter  10 value 619473.521690
#iter  20 value 228717.998172
#iter  30 value 129167.948923
#iter  40 value 113198.762318
#iter  50 value 96726.321240
#final  value 96726.293666 
#converged

summary(model)
#Call:
#multinom(formula = GMM_K3 ~ QD + FS + MQ + SOR + MQRankSum + 
#    ReadPosRankSum, data = Clustering_Subset_Complete, Hess = TRUE)

#Coefficients:
#       (Intercept)          QD           FS        MQ        SOR  MQRankSum ReadPosRankSum
#Medium    300.4707 -0.04109829 -0.004160641 -5.039398 -0.3853911 -0.3312125     0.04978636
#Bad       288.3144 -0.08339560  0.744182891 -5.078179  5.0653570 -0.1945217    -0.20757528

#Std. Errors:
#       (Intercept)          QD          FS          MQ        SOR   MQRankSum ReadPosRankSum
#Medium  0.03778900 0.001013617 0.002669009 0.000803940 0.03042538 0.007359524     0.01190352
#Bad     0.03782281 0.001338984 0.004270833 0.001128099 0.03785423 0.008506358     0.01502286

#Residual Deviance: 193452.6 
#AIC: 193480.6

exp(coef(model))
#         (Intercept)        QD       FS          MQ         SOR MQRankSum ReadPosRankSum
#Medium 3.110139e+130 0.9597348 0.995848 0.006477648   0.6801846 0.7180525      1.0510465
#Bad    1.634382e+125 0.9199871 2.104721 0.006231244 158.4369983 0.8232283      0.8125521



#####################################################################
############################ PART 4 #################################
#####################################################################

### Constructing a training set

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
#1    31  9.38  5.98  1.90  59.7         0         -0.046
sum(good_strat2_prob$HF_Violations)
#[1] 4

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
#1   178  19.4  4.34 0.944  58.8         0           0.17
sum(good_strat3_prob$HF_Violations)
#[1] 3

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
#1   188  22.6  3.75 0.896  58.8         0              0
sum(good_strat4_prob$HF_Violations)
#[1] 8

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
#1   222  19.5  4.94 0.923  58.9         0          0.098
sum(good_strat5_prob$HF_Violations)
#[1] 3

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
#1   387  19.8  4.63 0.911  58.9         0          0.197
sum(good_strat6_prob$HF_Violations)
#[1] 15

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
#1 149397  22.8 0.897 0.676    60         0          0.101
sum(good_strat7_prob$HF_Violations)
#[1] 491

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
#1   148  10.7  5.39 0.574  50.1      1.73          0.779
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
#1  4852   7.8  9.74  1.50  52.7     -2.50          0.145
sum(medium_strat3_prob$HF_Violations)
#[1] 1545

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
#1  5362  7.90  9.48  1.47  53.0     -2.63          0.140
sum(medium_strat4_prob$HF_Violations)
#[1] 1663

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
#1  7298  8.04  9.09  1.43  53.0     -2.66          0.127
sum(medium_strat5_prob$HF_Violations)
#[1] 2277

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
#1 13193  7.05  8.19  1.39    53     -2.66          0.111
sum(medium_strat6_prob$HF_Violations)
#[1] 4360

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
#1 377540  14.9  1.62 0.746  53.6     -1.87          0.152
sum(medium_strat7_prob$HF_Violations)
#[1] 95536

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
#1  4873  15.1  193.  8.86  42.8      -2.1              0
sum(bad_strat2_prob$HF_Violations)
#[1] 10879

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
#1  4987  7.95  10.8   1.6  52.5     -2.56          0.116
sum(bad_strat3_prob$HF_Violations)
#[1] 2818

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
#1  5268  8.04  11.1  1.64  52.4     -2.52          0.126
sum(bad_strat4_prob$HF_Violations)
#[1] 2892

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
#1  6393  8.37  11.6  1.72  52.0     -2.48          0.133
sum(bad_strat5_prob$HF_Violations)
#[1] 3986
						
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
#1  9418  8.98  12.1  1.78  51.5     -2.35          0.126
sum(bad_strat6_prob$HF_Violations)
#[1] 6262

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
#1 307444  5.44  27.2  3.20  53.2     -3.03          0.049
sum(bad_strat7_prob$HF_Violations)
#[1] 320617

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
	filter(GMM_K3 == "Good" & GMM_K3_Good >= 0.8)

## Filter medium variants
VQSR_training_medium <- Clustering_Subset_Complete %>%
	filter(GMM_K3 == "Medium" & GMM_K3_Med >= 0.7)

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
#1 149784  22.8 0.902 0.676    60         0          0.101

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
#1 398031  14.4  1.76 0.759  53.6     -1.91          0.151

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
write.table(VQSR_training_good_coords, "./TRAINING/VQSR_training_SNP_good_DV10.txt",
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
write.table(VQSR_training_medium_coords, "./TRAINING/VQSR_training_SNP_medium_DV10.txt",
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
