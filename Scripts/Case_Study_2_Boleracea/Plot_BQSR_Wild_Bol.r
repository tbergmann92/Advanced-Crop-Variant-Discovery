### R Script 
### Variant Calling in a Dataset of Brassica oleracea
### Plot Initial Base Quality Score Recalibration (BQSR)
### Author: Thomas Bergmann

### Prepare the environment
rm(list=ls())
set.seed(1)
setwd("PATH/TO/BQSR")

### Load packages
library(ggrepel)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)
library(tidyr)
library(Metrics)

# Read your CSV file - BQSR 1
SRR24208602_BQSR1 <- read.csv("./Round_1/SRR24208602.recalibration.csv", stringsAsFactors = FALSE)
SRR24208603_BQSR1 <- read.csv("./Round_1/SRR24208603.recalibration.csv", stringsAsFactors = FALSE)
SRR24208604_BQSR1 <- read.csv("./Round_1/SRR24208604.recalibration.csv", stringsAsFactors = FALSE)
SRR24208605_BQSR1 <- read.csv("./Round_1/SRR24208605.recalibration.csv", stringsAsFactors = FALSE)
SRR24208607_BQSR1 <- read.csv("./Round_1/SRR24208607.recalibration.csv", stringsAsFactors = FALSE)
SRR24208608_BQSR1 <- read.csv("./Round_1/SRR24208608.recalibration.csv", stringsAsFactors = FALSE)
SRR24208616_BQSR1 <- read.csv("./Round_1/SRR24208616.recalibration.csv", stringsAsFactors = FALSE)

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
SRR24208602_QS_BQSR1 <- filterQS(SRR24208602_BQSR1)
SRR24208603_QS_BQSR1 <- filterQS(SRR24208603_BQSR1)
SRR24208604_QS_BQSR1 <- filterQS(SRR24208604_BQSR1)
SRR24208605_QS_BQSR1 <- filterQS(SRR24208605_BQSR1)
SRR24208607_QS_BQSR1 <- filterQS(SRR24208607_BQSR1)
SRR24208608_QS_BQSR1 <- filterQS(SRR24208608_BQSR1)
SRR24208616_QS_BQSR1 <- filterQS(SRR24208616_BQSR1)

# Compute RMSEs
comp_rmse(SRR24208602_QS_BQSR1)
comp_rmse(SRR24208603_QS_BQSR1)
comp_rmse(SRR24208604_QS_BQSR1)
comp_rmse(SRR24208605_QS_BQSR1)
comp_rmse(SRR24208607_QS_BQSR1)
comp_rmse(SRR24208608_QS_BQSR1)
comp_rmse(SRR24208616_QS_BQSR1)

# Add RMSE labels (too lazy to paste into a function here)
label_df_SRR24208602_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 12.4),
    paste0("BQSR 1, RMSE = ", 0.679)
  )
)
label_df_SRR24208603_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 12.4),
    paste0("BQSR 1, RMSE = ", 0.655)
  )
)
label_df_SRR24208604_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 11.9),
    paste0("BQSR 1, RMSE = ", 0.433)
  )
)
label_df_SRR24208605_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 11.2),
    paste0("BQSR 1, RMSE = ", 0.577)
  )
)
label_df_SRR24208607_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 11.7),
    paste0("BQSR 1, RMSE = ", 0.683)
  )
)
label_df_SRR24208608_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 11.2),
    paste0("BQSR 1, RMSE = ", 0.661)
  )
)
label_df_SRR24208616_BQSR1 <- tibble(
  x = 8, y = c(42, 39),
  label = c(
    paste0("Original, RMSE = ", 7.70),
    paste0("BQSR 1, RMSE = ", 0.707)
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
SRR24208602_QS_BQSR1_Scatter <- plotBias(SRR24208602_QS_BQSR1, "SRR24208602", label_df_SRR24208602_BQSR1)
SRR24208603_QS_BQSR1_Scatter <- plotBias(SRR24208603_QS_BQSR1, "SRR24208603", label_df_SRR24208603_BQSR1)
SRR24208604_QS_BQSR1_Scatter <- plotBias(SRR24208604_QS_BQSR1, "SRR24208604", label_df_SRR24208604_BQSR1)
SRR24208605_QS_BQSR1_Scatter <- plotBias(SRR24208605_QS_BQSR1, "SRR24208605", label_df_SRR24208605_BQSR1)
SRR24208607_QS_BQSR1_Scatter <- plotBias(SRR24208607_QS_BQSR1, "SRR24208607", label_df_SRR24208607_BQSR1)
SRR24208608_QS_BQSR1_Scatter <- plotBias(SRR24208608_QS_BQSR1, "SRR24208608", label_df_SRR24208608_BQSR1)
SRR24208616_QS_BQSR1_Scatter <- plotBias(SRR24208616_QS_BQSR1, "SRR24208616", label_df_SRR24208616_BQSR1)

# Plot all together in one common plot:
QS_BQSR1_Scatter <- cowplot::plot_grid(SRR24208602_QS_BQSR1_Scatter, SRR24208603_QS_BQSR1_Scatter, 
									   SRR24208604_QS_BQSR1_Scatter, SRR24208605_QS_BQSR1_Scatter,
									   SRR24208607_QS_BQSR1_Scatter,SRR24208608_QS_BQSR1_Scatter, 
									   SRR24208616_QS_BQSR1_Scatter,
									   nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_Overall_Effect_W03_BQSR1.png", width=3000, height=3500, res=300)
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
	
# Filter the tables for QS
SRR24208602_Nucl_BQSR1 <- filterContext(SRR24208602_BQSR1)
SRR24208603_Nucl_BQSR1 <- filterContext(SRR24208603_BQSR1)
SRR24208604_Nucl_BQSR1 <- filterContext(SRR24208604_BQSR1)
SRR24208605_Nucl_BQSR1 <- filterContext(SRR24208605_BQSR1)
SRR24208607_Nucl_BQSR1 <- filterContext(SRR24208607_BQSR1)
SRR24208608_Nucl_BQSR1 <- filterContext(SRR24208608_BQSR1)
SRR24208616_Nucl_BQSR1 <- filterContext(SRR24208616_BQSR1)

# Compute RMSEs
comp_rmse(SRR24208602_Nucl_BQSR1)
comp_rmse(SRR24208603_Nucl_BQSR1)
comp_rmse(SRR24208604_Nucl_BQSR1)
comp_rmse(SRR24208605_Nucl_BQSR1)
comp_rmse(SRR24208607_Nucl_BQSR1)
comp_rmse(SRR24208608_Nucl_BQSR1)
comp_rmse(SRR24208616_Nucl_BQSR1)

# Add RMSE labels
label_df_SRR24208602_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 9.95),
    paste0("BQSR 1, RMSE = ", 0.405)
  )
)
label_df_SRR24208603_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 9.54),
    paste0("BQSR 1, RMSE = ", 0.400)
  )
)
label_df_SRR24208604_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 8.30),
    paste0("BQSR 1, RMSE = ", 0.425)
  )
)
label_df_SRR24208605_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 8.30),
    paste0("BQSR 1, RMSE = ", 0.422)
  )
)
label_df_SRR24208607_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 8.48),
    paste0("BQSR 1, RMSE = ", 0.449)
  )
)
label_df_SRR24208608_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 8.04),
    paste0("BQSR 1, RMSE = ", 0.442)
  )
)
label_df_SRR24208616_BQSR1 <- tibble(
  x = 2, y = c(4, 2),
  label = c(
    paste0("Original, RMSE = ", 4.18),
    paste0("BQSR 1, RMSE = ", 0.337)
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
SRR24208602_Dinucl_BQSR1 <- plotContext(SRR24208602_Nucl_BQSR1, "SRR24208602", label_df_SRR24208602_BQSR1)
SRR24208603_Dinucl_BQSR1 <- plotContext(SRR24208603_Nucl_BQSR1, "SRR24208603", label_df_SRR24208603_BQSR1)
SRR24208604_Dinucl_BQSR1 <- plotContext(SRR24208604_Nucl_BQSR1, "SRR24208604", label_df_SRR24208604_BQSR1)
SRR24208605_Dinucl_BQSR1 <- plotContext(SRR24208605_Nucl_BQSR1, "SRR24208605", label_df_SRR24208605_BQSR1)
SRR24208607_Dinucl_BQSR1 <- plotContext(SRR24208607_Nucl_BQSR1, "SRR24208607", label_df_SRR24208607_BQSR1)
SRR24208608_Dinucl_BQSR1 <- plotContext(SRR24208608_Nucl_BQSR1, "SRR24208608", label_df_SRR24208608_BQSR1)
SRR24208616_Dinucl_BQSR1 <- plotContext(SRR24208616_Nucl_BQSR1, "SRR24208616", label_df_SRR24208616_BQSR1)

# Plot all together in one common plot:
Dinucl_Plots_BQSR1 <- cowplot::plot_grid(SRR24208602_Dinucl_BQSR1, SRR24208603_Dinucl_BQSR1,
										 SRR24208604_Dinucl_BQSR1, SRR24208605_Dinucl_BQSR1,
										 SRR24208607_Dinucl_BQSR1, SRR24208608_Dinucl_BQSR1,
										 SRR24208616_Dinucl_BQSR1,
										 nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_Dinucleotide_Effect_W03_BQSR1.png", width=4500, height=5000, res=300)
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
SRR24208602_Cycle_BQSR1 <- filterCycle(SRR24208602_BQSR1)
SRR24208603_Cycle_BQSR1 <- filterCycle(SRR24208603_BQSR1)
SRR24208604_Cycle_BQSR1 <- filterCycle(SRR24208604_BQSR1)
SRR24208605_Cycle_BQSR1 <- filterCycle(SRR24208605_BQSR1)
SRR24208607_Cycle_BQSR1 <- filterCycle(SRR24208607_BQSR1)
SRR24208608_Cycle_BQSR1 <- filterCycle(SRR24208608_BQSR1)
SRR24208616_Cycle_BQSR1 <- filterCycle(SRR24208616_BQSR1)

# Compute RMSEs
comp_rmse(SRR24208602_Cycle_BQSR1)
comp_rmse(SRR24208603_Cycle_BQSR1)
comp_rmse(SRR24208604_Cycle_BQSR1)
comp_rmse(SRR24208605_Cycle_BQSR1)
comp_rmse(SRR24208607_Cycle_BQSR1)
comp_rmse(SRR24208608_Cycle_BQSR1)
comp_rmse(SRR24208616_Cycle_BQSR1)

# Add RMSE labels
label_df_SRR24208602_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 9.72),
    paste0("BQSR 1, RMSE = ", 0.431)
  )
)
label_df_SRR24208603_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 9.44),
    paste0("BQSR 1, RMSE = ", 0.389)
  )
)
label_df_SRR24208604_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.87),
    paste0("BQSR 1, RMSE = ", 0.380)
  )
)
label_df_SRR24208605_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ",  8.19),
    paste0("BQSR 1, RMSE = ", 0.399)
  )
)
label_df_SRR24208607_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 8.38),
    paste0("BQSR 1, RMSE = ", 0.412)
  )
)
label_df_SRR24208608_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 7.78),
    paste0("BQSR 1, RMSE = ", 0.401)
  )
)
label_df_SRR24208616_BQSR1 <- tibble(
  x = -130, y = c(4, 2.5),
  label = c(
    paste0("Original, RMSE = ", 4.16),
    paste0("BQSR 1, RMSE = ", 0.381)
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
SRR24208602_CycleCov_BQSR1 <- plotCycle(SRR24208602_Cycle_BQSR1, "SRR24208602", label_df_SRR24208602_BQSR1)
SRR24208603_CycleCov_BQSR1 <- plotCycle(SRR24208603_Cycle_BQSR1, "SRR24208603", label_df_SRR24208603_BQSR1)
SRR24208604_CycleCov_BQSR1 <- plotCycle(SRR24208604_Cycle_BQSR1, "SRR24208604", label_df_SRR24208604_BQSR1)
SRR24208605_CycleCov_BQSR1 <- plotCycle(SRR24208605_Cycle_BQSR1, "SRR24208605", label_df_SRR24208605_BQSR1)
SRR24208607_CycleCov_BQSR1 <- plotCycle(SRR24208607_Cycle_BQSR1, "SRR24208607", label_df_SRR24208607_BQSR1)
SRR24208608_CycleCov_BQSR1 <- plotCycle(SRR24208608_Cycle_BQSR1, "SRR24208608", label_df_SRR24208608_BQSR1)
SRR24208616_CycleCov_BQSR1 <- plotCycle(SRR24208616_Cycle_BQSR1, "SRR24208616", label_df_SRR24208616_BQSR1)

# Plot all together in one common plot:
CycleCov_Plots_BQSR1 <- cowplot::plot_grid(SRR24208602_CycleCov_BQSR1, SRR24208603_CycleCov_BQSR1, 
										   SRR24208604_CycleCov_BQSR1, SRR24208605_CycleCov_BQSR1, 
										   SRR24208607_CycleCov_BQSR1, SRR24208608_CycleCov_BQSR1,
										   SRR24208616_CycleCov_BQSR1,
										   nrow=4, ncol=2, align="hv")
#png("./Round_1/BQSR_CycleCov_Effect_W03_BQSR1.png", width=4500, height=5000, res=300)
CycleCov_Plots_BQSR1
#dev.off()

