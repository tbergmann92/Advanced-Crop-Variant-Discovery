### R Script 
### Variant Calling in a Dataset of Brassica napus (Mini-core collection)
### Plot Initial Base Quality Score Recalibration (BQSR)
### Author: Thomas Bergmann

### Prepare the environment
rm(list=ls())
set.seed(1)
setwd("C:/PATH/TO/BQSR/FILES")

### Load packages
library(ggrepel)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)
library(tidyr)
library(Metrics)

# Read your CSV file - BQSR 1
CRR379057_BQSR1 <- read.csv("./Round_1/CRR379057.recalibration.csv", stringsAsFactors = FALSE)
CRR379077_BQSR1 <- read.csv("./Round_1/CRR379077.recalibration.csv", stringsAsFactors = FALSE)
CRR379100_BQSR1 <- read.csv("./Round_1/CRR379100.recalibration.csv", stringsAsFactors = FALSE)
CRR379108_BQSR1 <- read.csv("./Round_1/CRR379108.recalibration.csv", stringsAsFactors = FALSE)
CRR379118_BQSR1 <- read.csv("./Round_1/CRR379118.recalibration.csv", stringsAsFactors = FALSE)
CRR379141_BQSR1 <- read.csv("./Round_1/CRR379141.recalibration.csv", stringsAsFactors = FALSE)
CRR379165_BQSR1 <- read.csv("./Round_1/CRR379165.recalibration.csv", stringsAsFactors = FALSE)
CRR379187_BQSR1 <- read.csv("./Round_1/CRR379187.recalibration.csv", stringsAsFactors = FALSE)
# Random sample
CRR379121_BQSR1 <- read.csv("./Round_1/CRR379121.recalibration.csv", stringsAsFactors = FALSE)

# Function to filter for Quality Score
filterQS <- function(recal_table) {
	filtered_table <- 	recal_table %>% 
						filter(CovariateName == "QualityScore") %>%
						mutate(ReadGroup = as.factor(ReadGroup),
							CovariateValue = as.factor(CovariateValue),
							CovariateName = as.factor(CovariateName),
							EventType = as.factor(EventType),
							Recalibration = as.factor(Recalibration)) %>%
						group_by(Recalibration, CovariateValue) %>%
						summarise(n = n(),
							EmpiricalQuality=mean(EmpiricalQuality),
							ReportedQuality=mean(AverageReportedQuality),
							Observations = sum(Observations, na.rm = TRUE))
	
	return(filtered_table)

}	

# Function to compute RMSE
comp_rmse <- function(table) {
  table %>%
    group_by(Recalibration) %>%
    summarise(RMSE = rmse(EmpiricalQuality, ReportedQuality), .groups = "drop")
}

# Filter the tables for QS
CRR379057_QS_BQSR1 <- filterQS(CRR379057_BQSR1)
CRR379077_QS_BQSR1 <- filterQS(CRR379077_BQSR1)
CRR379100_QS_BQSR1 <- filterQS(CRR379100_BQSR1)
CRR379108_QS_BQSR1 <- filterQS(CRR379108_BQSR1)
CRR379118_QS_BQSR1 <- filterQS(CRR379118_BQSR1)
CRR379141_QS_BQSR1 <- filterQS(CRR379141_BQSR1)
CRR379165_QS_BQSR1 <- filterQS(CRR379165_BQSR1)
CRR379187_QS_BQSR1 <- filterQS(CRR379187_BQSR1)

# Compute RMSEs
comp_rmse(CRR379057_QS_BQSR1)
comp_rmse(CRR379077_QS_BQSR1)
comp_rmse(CRR379100_QS_BQSR1)
comp_rmse(CRR379108_QS_BQSR1)
comp_rmse(CRR379118_QS_BQSR1)
comp_rmse(CRR379141_QS_BQSR1)
comp_rmse(CRR379165_QS_BQSR1)
comp_rmse(CRR379187_QS_BQSR1)

# Add RMSE labels (too lazy to paste into a function here)
# Feel free to automize this
label_df_CRR379057_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 2.36)
  )
)
label_df_CRR379077_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.6),
    paste0("BQSR 1, RMSE = ", 2.61)
  )
)
label_df_CRR379100_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.4),
    paste0("BQSR 1, RMSE = ", 2.19)
  )
)
label_df_CRR379108_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 2.13)
  )
)
label_df_CRR379118_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.4),
    paste0("BQSR 1, RMSE = ", 2.37)
  )
)
label_df_CRR379141_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 1.76)
  )
)
label_df_CRR379165_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 2.11)
  )
)
label_df_CRR379187_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 1.87)
  )
)

# Function to plot the overall effect of bias for BQSR1
plotBias <- function(table, sample, sample_label) {
	ggplot(table, aes(x=ReportedQuality, 
				      y=EmpiricalQuality, 
					  size = Observations, 
					  color=Recalibration, 
					  fill=Recalibration
					  ))+
		geom_point(shape = 21, alpha = 0.8, stroke = 0.4) +
		geom_abline(slope = 1, intercept=0, linetype="dashed", linewidth=1)+
		scale_color_manual(
				values = c("Before" = "#C51B7D", "After" = "#2166AC"), 
				aesthetics = c("color", "fill"),
				labels = c("Before" = "Original", "After" = "BQSR 1")
				) +
		scale_size_continuous(
				name = "Observations",
				range = c(2, 6),      # adjust this to your dataset
				trans = "sqrt"        # keeps large counts readable
				) +
		labs(title = paste0("Overall effect of bias for ", sample),
			x = "Reported Quality",
			y = "Empirical Quality",
			color = NULL, fill = NULL,
			size = "Observations")+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3, hjust = 0) +
		scale_y_continuous(expand=c(0,0), limits=c(5,45), breaks=c(0,10,20,30,40))+
		scale_x_continuous(expand=c(0,0), limits=c(5,45), breaks=c(0,10,20,30,40))+
		theme(
			legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 8, face = "bold"),
			legend.title = element_text(colour = "black", size = 9, face = "bold"),
			legend.background = element_blank(),
			legend.spacing.y = unit(0.05, "cm"), # Increase horizontal space
			plot.title = element_text(colour = "black", size = 12, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 10, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 10, face = "bold"))
}

# Prepare single plots
CRR379057_QS_BQSR1_Scatter <- plotBias(CRR379057_QS_BQSR1, "CRR379057", label_df_CRR379057_BQSR1)
CRR379077_QS_BQSR1_Scatter <- plotBias(CRR379077_QS_BQSR1, "CRR379077", label_df_CRR379077_BQSR1)
CRR379100_QS_BQSR1_Scatter <- plotBias(CRR379100_QS_BQSR1, "CRR379100", label_df_CRR379100_BQSR1)
CRR379108_QS_BQSR1_Scatter <- plotBias(CRR379108_QS_BQSR1, "CRR379108", label_df_CRR379108_BQSR1)
CRR379118_QS_BQSR1_Scatter <- plotBias(CRR379118_QS_BQSR1, "CRR379118", label_df_CRR379118_BQSR1)
CRR379141_QS_BQSR1_Scatter <- plotBias(CRR379141_QS_BQSR1, "CRR379141", label_df_CRR379141_BQSR1)
CRR379165_QS_BQSR1_Scatter <- plotBias(CRR379165_QS_BQSR1, "CRR379165", label_df_CRR379165_BQSR1)
CRR379187_QS_BQSR1_Scatter <- plotBias(CRR379187_QS_BQSR1, "CRR379187", label_df_CRR379187_BQSR1)

# Plot all together in one common plot:
QS_BQSR1_Scatter <- cowplot::plot_grid(CRR379057_QS_BQSR1_Scatter, CRR379077_QS_BQSR1_Scatter, 
									   CRR379100_QS_BQSR1_Scatter, CRR379108_QS_BQSR1_Scatter,
									   CRR379118_QS_BQSR1_Scatter,CRR379141_QS_BQSR1_Scatter, 
									   CRR379165_QS_BQSR1_Scatter, CRR379187_QS_BQSR1_Scatter,
									   nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_Overall_Effect_DV10_BQSR1.png", width=3000, height=3500, res=300)
QS_BQSR1_Scatter
#dev.off()

# Function to filter for Dinucleotide Context
filterContext <- function(recal_table) {
	filtered_table <- 	recal_table %>% 
						filter(CovariateName == "Context") %>%
						mutate(ReadGroup = as.factor(ReadGroup),
							CovariateValue = as.factor(CovariateValue),
							CovariateName = as.factor(CovariateName),
							EventType = as.factor(EventType),
							Recalibration = as.factor(Recalibration)) %>%
						group_by(Recalibration, CovariateValue) %>%
						summarise(n = n(), 
							EmpiricalQuality=mean(EmpiricalQuality),
							ReportedQuality=mean(AverageReportedQuality),
							meanAccuracy=mean(Accuracy),
							AccuracySD=sd(Accuracy),
							.groups = "drop")	
	return(filtered_table)

}	

# Filter the tables for Context
CRR379057_Nucl_BQSR1 <- filterContext(CRR379057_BQSR1)
CRR379077_Nucl_BQSR1 <- filterContext(CRR379077_BQSR1)
CRR379100_Nucl_BQSR1 <- filterContext(CRR379100_BQSR1)
CRR379108_Nucl_BQSR1 <- filterContext(CRR379108_BQSR1)
CRR379118_Nucl_BQSR1 <- filterContext(CRR379118_BQSR1)
CRR379141_Nucl_BQSR1 <- filterContext(CRR379141_BQSR1)
CRR379165_Nucl_BQSR1 <- filterContext(CRR379165_BQSR1)
CRR379187_Nucl_BQSR1 <- filterContext(CRR379187_BQSR1)

# Compute RMSEs
comp_rmse(CRR379057_Nucl_BQSR1)
comp_rmse(CRR379077_Nucl_BQSR1)
comp_rmse(CRR379100_Nucl_BQSR1)
comp_rmse(CRR379108_Nucl_BQSR1)
comp_rmse(CRR379118_Nucl_BQSR1)
comp_rmse(CRR379141_Nucl_BQSR1)
comp_rmse(CRR379165_Nucl_BQSR1)
comp_rmse(CRR379187_Nucl_BQSR1)

# Add RMSE labels
label_df_CRR379057_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 6.52),
    paste0("BQSR 1, RMSE = ", 0.30)
  )
)
label_df_CRR379077_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 6.12),
    paste0("BQSR 1, RMSE = ", 0.41)
  )
)
label_df_CRR379100_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 6.87),
    paste0("BQSR 1, RMSE = ", 0.49)
  )
)
label_df_CRR379108_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 7.36),
    paste0("BQSR 1, RMSE = ", 0.42)
  )
)
label_df_CRR379118_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 7.04),
    paste0("BQSR 1, RMSE = ", 0.34)
  )
)
label_df_CRR379141_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 7.87),
    paste0("BQSR 1, RMSE = ", 0.51)
  )
)
label_df_CRR379165_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 7.27),
    paste0("BQSR 1, RMSE = ", 0.32)
  )
)
label_df_CRR379187_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 7.86),
    paste0("BQSR 1, RMSE = ", 0.45)
  )
)

# Plot the dinucleotide bias
plotContext <- function(table, sample, sample_label) {
	ggplot(table, aes(x=CovariateValue, y=meanAccuracy, group=Recalibration, color=Recalibration))+
		geom_abline(slope=0, intercept=0, linetype="dashed", linewidth=1)+
		geom_point(size=15, shape="-")+
		labs(
			title=paste0("Bias on nucleotide context for ", sample),
			x="Dinucleotide",
			y="Error (Empirical - Reported Quality)"
			)+
		scale_color_manual(
			values = c("blue", "deeppink3"),
			labels = c("After" = "BQSR1", "Before" = "Original")) +
		scale_y_continuous(expand=c(0,0), limits=c(-14,6), breaks=seq(-14,4,2))+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 5, hjust = 0)+
		theme(
			legend.position = "none",
			#legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 10, face = "bold"),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.spacing.y = unit(8, "pt"),            # Adjust vertical spacing
			legend.key.height = unit(10, "pt"),          # Adjust symbol height
			plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 13, face = "bold"))+
		guides(color = guide_legend(
			keywidth = 11,
			keyheight = 11,
			default.unit = "pt"))

}

# Create the single plots
CRR379057_Dinucl_BQSR1 <- plotContext(CRR379057_Nucl_BQSR1, "CRR379057", label_df_CRR379057_BQSR1)
CRR379077_Dinucl_BQSR1 <- plotContext(CRR379077_Nucl_BQSR1, "CRR379077", label_df_CRR379077_BQSR1)
CRR379100_Dinucl_BQSR1 <- plotContext(CRR379100_Nucl_BQSR1, "CRR379100", label_df_CRR379100_BQSR1)
CRR379108_Dinucl_BQSR1 <- plotContext(CRR379108_Nucl_BQSR1, "CRR379108", label_df_CRR379108_BQSR1)
CRR379118_Dinucl_BQSR1 <- plotContext(CRR379118_Nucl_BQSR1, "CRR379118", label_df_CRR379118_BQSR1)
CRR379141_Dinucl_BQSR1 <- plotContext(CRR379141_Nucl_BQSR1, "CRR379141", label_df_CRR379141_BQSR1)
CRR379165_Dinucl_BQSR1 <- plotContext(CRR379165_Nucl_BQSR1, "CRR379165", label_df_CRR379165_BQSR1)
CRR379187_Dinucl_BQSR1 <- plotContext(CRR379187_Nucl_BQSR1, "CRR379187", label_df_CRR379187_BQSR1)

# Plot all together in one common plot:
Dinucl_Plots_BQSR1 <- cowplot::plot_grid(CRR379057_Dinucl_BQSR1, CRR379077_Dinucl_BQSR1,
										 CRR379100_Dinucl_BQSR1, CRR379108_Dinucl_BQSR1,
										 CRR379118_Dinucl_BQSR1, CRR379141_Dinucl_BQSR1,
										 CRR379165_Dinucl_BQSR1, CRR379187_Dinucl_BQSR1,
										 nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_Dinucleotide_Effect_DV10_BQSR1.png", width=4500, height=5000, res=300)
Dinucl_Plots_BQSR1
#dev.off()

# Function to filter for Dinucleotide Context
filterCycle <- function(recal_table) {
	filtered_table <- 	recal_table %>% 
						filter(CovariateName == "Cycle") %>%
						mutate(ReadGroup = as.factor(ReadGroup),
							CovariateValue = as.integer(CovariateValue),
							CovariateName = as.factor(CovariateName),
							EventType = as.factor(EventType),
							Recalibration = as.factor(Recalibration)) %>%
						group_by(Recalibration, CovariateValue) %>%
						summarise(n = n(), 
							EmpiricalQuality=mean(EmpiricalQuality),
							ReportedQuality=mean(AverageReportedQuality),
							meanAccuracy=mean(Accuracy),
							AccuracySD=sd(Accuracy),
							.groups = "drop")	
	return(filtered_table)

}	

# Filter the tables for Cycle
CRR379057_Cycle_BQSR1 <- filterCycle(CRR379057_BQSR1)
CRR379077_Cycle_BQSR1 <- filterCycle(CRR379077_BQSR1)
CRR379100_Cycle_BQSR1 <- filterCycle(CRR379100_BQSR1)
CRR379108_Cycle_BQSR1 <- filterCycle(CRR379108_BQSR1)
CRR379118_Cycle_BQSR1 <- filterCycle(CRR379118_BQSR1)
CRR379141_Cycle_BQSR1 <- filterCycle(CRR379141_BQSR1)
CRR379165_Cycle_BQSR1 <- filterCycle(CRR379165_BQSR1)
CRR379187_Cycle_BQSR1 <- filterCycle(CRR379187_BQSR1)

# Compute RMSEs
comp_rmse(CRR379057_Cycle_BQSR1)
comp_rmse(CRR379077_Cycle_BQSR1)
comp_rmse(CRR379100_Cycle_BQSR1)
comp_rmse(CRR379108_Cycle_BQSR1)
comp_rmse(CRR379118_Cycle_BQSR1)
comp_rmse(CRR379141_Cycle_BQSR1)
comp_rmse(CRR379165_Cycle_BQSR1)
comp_rmse(CRR379187_Cycle_BQSR1)

# Add RMSE labels
label_df_CRR379057_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.57),
    paste0("BQSR 1, RMSE = ", 0.32)
  )
)
label_df_CRR379077_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.11),
    paste0("BQSR 1, RMSE = ", 0.39)
  )
)
label_df_CRR379100_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.26),
    paste0("BQSR 1, RMSE = ", 0.43)
  )
)
label_df_CRR379108_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ",  7.96),
    paste0("BQSR 1, RMSE = ", 0.43)
  )
)
label_df_CRR379118_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.94),
    paste0("BQSR 1, RMSE = ", 0.44)
  )
)
label_df_CRR379141_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 8.40),
    paste0("BQSR 1, RMSE = ", 0.54)
  )
)
label_df_CRR379165_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 8.23),
    paste0("BQSR 1, RMSE = ", 0.39)
  )
)
label_df_CRR379187_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 8.82),
    paste0("BQSR 1, RMSE = ", 0.44)
  )
)

# Plot the cycle covariate
plotCycle <- function(table, sample, sample_label) {
	ggplot(table, aes(x=CovariateValue, 
					  y=meanAccuracy, 
					  group=Recalibration, 
					  color=Recalibration,
					  fill=Recalibration))+
		geom_point(shape = 21, alpha = 0.8, stroke = 0.4) +
		geom_abline(slope = 0, intercept=0, linetype="dashed", linewidth=1)+
		scale_color_manual(
				values = c("Before" = "#C51B7D", "After" = "#2166AC"), 
				aesthetics = c("color", "fill"),
				labels = c("Before" = "Original", "After" = "BQSR 1")
				) +
		labs(
			title=paste0("Bias on cycle covariate for ", sample),
			x="Cycle Covariate",
			y="Error (Empirical - Reported Quality)"
			)+
		scale_y_continuous(expand=c(0,0), limits=c(-14,6), breaks=seq(-14,6,2))+
		scale_x_continuous(expand=c(0,0), limits=c(-160,160), breaks=seq(-150,150,50))+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 5, hjust = 0)+
		theme(
			legend.position = "none",
			#legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 10, face = "bold"),
			legend.title = element_blank(),
			legend.background = element_blank(),
			legend.key = element_blank(),
			legend.spacing.y = unit(8, "pt"),            # Adjust vertical spacing
			legend.key.height = unit(10, "pt"),          # Adjust symbol height
			plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 13, face = "bold"))+
		guides(color = guide_legend(
			keywidth = 10,
			keyheight = 10,
			default.unit = "pt"))

}

# Create the single plots
CRR379057_CycleCov_BQSR1 <- plotCycle(CRR379057_Cycle_BQSR1, "CRR379057", label_df_CRR379057_BQSR1)
CRR379077_CycleCov_BQSR1 <- plotCycle(CRR379077_Cycle_BQSR1, "CRR379077", label_df_CRR379077_BQSR1)
CRR379100_CycleCov_BQSR1 <- plotCycle(CRR379100_Cycle_BQSR1, "CRR379100", label_df_CRR379100_BQSR1)
CRR379108_CycleCov_BQSR1 <- plotCycle(CRR379108_Cycle_BQSR1, "CRR379108", label_df_CRR379108_BQSR1)
CRR379118_CycleCov_BQSR1 <- plotCycle(CRR379118_Cycle_BQSR1, "CRR379118", label_df_CRR379118_BQSR1)
CRR379141_CycleCov_BQSR1 <- plotCycle(CRR379141_Cycle_BQSR1, "CRR379141", label_df_CRR379141_BQSR1)
CRR379165_CycleCov_BQSR1 <- plotCycle(CRR379165_Cycle_BQSR1, "CRR379165", label_df_CRR379165_BQSR1)
CRR379187_CycleCov_BQSR1 <- plotCycle(CRR379187_Cycle_BQSR1, "CRR379187", label_df_CRR379187_BQSR1)

# Plot all together in one common plot:
CycleCov_Plots_BQSR1 <- cowplot::plot_grid(CRR379057_CycleCov_BQSR1, CRR379077_CycleCov_BQSR1, 
										   CRR379100_CycleCov_BQSR1, CRR379108_CycleCov_BQSR1, 
										   CRR379118_CycleCov_BQSR1, CRR379141_CycleCov_BQSR1,
										   CRR379165_CycleCov_BQSR1, CRR379187_CycleCov_BQSR1,
										   nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_CycleCov_Effect_DV10_BQSR1.png", width=4500, height=5000, res=300)
CycleCov_Plots_BQSR1
#dev.off()


####################################################################################################
############################# Publication Plot for CRR379121 #######################################
####################################################################################################

# prepare data for overall bias
CRR379121_QS_BQSR1 <- filterQS(CRR379121_BQSR1)
# calculate RMSE
comp_rmse(CRR379121_QS_BQSR1)
# label the RMSE
label_df_CRR379121_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 10.5),
    paste0("BQSR 1, RMSE = ", 1.65)
  )
)

# Function to plot the overall effect of bias for BQSR1
plotBias <- function(table, sample, sample_label) {
	
	ggplot(table, aes(x=ReportedQuality, 
				      y=EmpiricalQuality, 
					  size = Observations, 
					  color=Recalibration, 
					  fill=Recalibration
					  ))+
		geom_point(shape = 21, alpha = 0.8, stroke = 0.4) +
		geom_abline(slope = 1, intercept=0, linetype="dashed", linewidth=1)+
		scale_color_manual(
				values = c("Before" = "#C51B7D", "After" = "#2166AC"), 
				aesthetics = c("color", "fill"),
				labels = c("Before" = "Original", "After" = "BQSR 1")
				) +
		scale_size_continuous(
				name = "Observations",
				range = c(2, 6),      # adjust this to your dataset
				trans = "sqrt"        # keeps large counts readable
				) +
		labs(title = paste0("Overall effect of bias for ", sample),
			x = "Reported Quality",
			y = "Empirical Quality",
			color = NULL, fill = NULL,
			size = "Observations")+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4, hjust = 0) +
		scale_y_continuous(expand=c(0,0), limits=c(5,45), breaks=c(0,10,20,30,40))+
		scale_x_continuous(expand=c(0,0), limits=c(5,45), breaks=c(0,10,20,30,40))+
		theme(
			legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 12, face = "bold"),
			legend.title = element_text(colour = "black", size = 12, face = "bold"),
			legend.background = element_blank(),
			legend.spacing.y = unit(0.05, "cm"), # Increase horizontal space
			plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 13, face = "bold"))
}
# generate plot
CRR379121_QS_BQSR1_Scatter <- plotBias(CRR379121_QS_BQSR1, "CRR379121", label_df_CRR379121_BQSR1)

# prepare data for dinucleotide context
CRR379121_Nucl_BQSR1 <- filterContext(CRR379121_BQSR1)
# calculate RMSE
comp_rmse(CRR379121_Nucl_BQSR1)
# label the RMSE
label_df_CRR379121_BQSR1 <- tibble(
  x = 2, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.77),
    paste0("BQSR 1, RMSE = ", 0.43)
  )
)

# Plot the dinucleotide bias
plotContext <- function(table, sample, sample_label) {
	ggplot(table, aes(x=CovariateValue, y=meanAccuracy, group=Recalibration, color=Recalibration))+
		geom_abline(slope=0, intercept=0, linetype="dashed", linewidth=1)+
		geom_point(size=15, shape="-")+
		labs(
			title=paste0("Bias on nucleotide context for ", sample),
			x="Dinucleotide",
			y="Error (Empirical - Reported Quality)"
			)+
		scale_color_manual(
			values = c("#2166AC", "#C51B7D"),
			labels = c("After" = "BQSR1", "Before" = "Original")) +
		scale_y_continuous(expand=c(0,0), limits=c(-14,6), breaks=seq(-14,4,2))+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4, hjust = 0)+
		theme(
			legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 12, face = "bold"),
			legend.title = element_text(colour = "black", size = 12, face = "bold"),
			legend.background = element_blank(),
			legend.spacing.y = unit(0.01, "cm"), # Increase horizontal space
			plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 12, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 13, face = "bold"))
}

# generate the plot
CRR379121_Dinucl_BQSR1 <- plotContext(CRR379121_Nucl_BQSR1, "CRR379121", label_df_CRR379121_BQSR1)

# prepare data for machine cycle context
CRR379121_Cycle_BQSR1 <- filterCycle(CRR379121_BQSR1)
# calculate RMSE
comp_rmse(CRR379121_Cycle_BQSR1)
# label the RMSE
label_df_CRR379121_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 8.46),
    paste0("BQSR 1, RMSE = ", 0.48)
  )
)

# Plot the cycle covariate
plotCycle <- function(table, sample, sample_label) {
	ggplot(table, aes(x=CovariateValue, 
					  y=meanAccuracy, 
					  group=Recalibration, 
					  color=Recalibration,
					  fill=Recalibration))+
		geom_point(shape = 21, alpha = 0.8, stroke = 0.4) +
		geom_abline(slope = 0, intercept=0, linetype="dashed", linewidth=1)+
		scale_color_manual(
				values = c("Before" = "#C51B7D", "After" = "#2166AC"), 
				aesthetics = c("color", "fill"),
				labels = c("Before" = "Original", "After" = "BQSR 1")
				) +
		labs(
			title=paste0("Bias on cycle covariate for ", sample),
			x="Cycle Covariate",
			y="Error (Empirical - Reported Quality)"
			)+
		scale_y_continuous(expand=c(0,0), limits=c(-14,6), breaks=seq(-14,6,2))+
		scale_x_continuous(expand=c(0,0), limits=c(-160,160), breaks=seq(-150,150,50))+
		geom_text(data = sample_label, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 4, hjust = 0)+
		theme(
			legend.position.inside = c(0.15,0.60),
			legend.text = element_text(colour = "black", size = 12, face = "bold"),
			legend.title = element_text(colour = "black", size = 12, face = "bold"),
			legend.background = element_blank(),
			legend.spacing.y = unit(0.05, "cm"), # Increase horizontal space
			plot.title = element_text(colour = "black", size = 15, face = "bold", hjust = 0.5),
			axis.title.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.title.y = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.x = element_text(colour = "black", size = 13, face = "bold"),
			axis.text.y = element_text(colour = "black", size = 13, face = "bold"))
}

# generate the plot
CRR379121_CycleCov_BQSR1 <- plotCycle(CRR379121_Cycle_BQSR1, "CRR379121", label_df_CRR379121_BQSR1)


### Plot altogether
BQSR_Plots_CRR379121 <- cowplot::plot_grid(CRR379121_QS_BQSR1_Scatter,
										   CRR379121_Dinucl_BQSR1,
										   CRR379121_CycleCov_BQSR1,
										   nrow=3, ncol=1, align ="hv", scale=0.9,
										   labels="AUTO", label_size= 24)

#png("./Round_1/BQSR_Covariates_CRR379121_DV10_BQSR.png", width=2500, height=3500, res=300)
BQSR_Plots_CRR379121
#dev.off()

