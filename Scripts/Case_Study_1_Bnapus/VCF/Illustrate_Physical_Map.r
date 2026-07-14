# Load required libraries
library(ggplot2)
library(cowplot)

# Load your data
known_set <- read.table("PATH/TO/KNOWN/SITES/TXT")
truth_set <- read.table("PATH/TO/TRUE/SITES/TXT")
subset <- read.table("PATH/TO/SUBSET/SITES/TXT")
colnames(known_set) <- c("CHROM", "POS")
colnames(truth_set) <- c("CHROM", "POS")
colnames(subset) <- c("CHROM", "POS")

setwd("PATH/TO/WORK/DIR/")

# convert columns to appropriate formats
known_set$Chromosome <- as.factor(known_set$CHROM)
known_set$POS <- as.numeric(known_set$POS)
truth_set$CHROM <- as.factor(truth_set$CHROM)
truth_set$POS <- as.numeric(truth_set$POS)
subset$CHROM <- as.factor(subset$CHROM)
subset$POS <- as.numeric(subset$POS)

# change the factor order of the chromosomes for plotting
rev_order <- c("C09", "C08", "C07", "C06", "C05", "C04", "C03", "C02", "C01",
		   "A10", "A09", "A08", "A07", "A06", "A05", "A04", "A03", "A02", "A01")
known_set$CHROM <- factor(known_set$CHROM, levels = rev_order)
truth_set$CHROM <- factor(truth_set$CHROM, levels = rev_order)
subset$CHROM <- factor(subset$CHROM, levels = rev_order)

# check
print(levels(known_set$CHROM))
print(levels(truth_set$CHROM))
print(levels(subset$CHROM))

### Known Sites
Known_Sites_Map <- ggplot(known_set, aes(x = POS, y = CHROM)) +
geom_point(color = "dodgerblue4", size = 1) +
scale_x_continuous(breaks=c(0,10000000,20000000,30000000,40000000,50000000,60000000,70000000),
	labels=c("0", "10","20","30","40","50","60","70")) +
annotate("segment", x=0, y=1, xend=66465249, yend=1, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=2, xend=41681856, yend=2, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=3, xend=55656957, yend=3, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=4, xend=50218839, yend=4, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=5, xend=56382805, yend=5, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=6, xend=65837619, yend=6, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=7, xend=73669886, yend=7, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=8, xend=62297340, yend=8, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=9, xend=48239358, yend=9, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=10, xend=20778245, yend=10, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=11, xend=53549826, yend=11, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=12, xend=26309499, yend=12, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=13, xend=29390523, yend=13, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=14, xend=45146386, yend=14, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=15, xend=42112164, yend=15, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=16, xend=23101715, yend=16, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=17, xend=39685748, yend=17, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=18, xend=33432960, yend=18, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=19, xend=32958928, yend=19, linewidth=0.5, color="black", linetype=1)+
labs(
title = " ",
x = "Length (Mbp)",
y = "Chromosome"
) +
theme_gray() +
theme(
	panel.grid.major = element_line(color = "gray80"),
	panel.grid.minor = element_blank(),
	legend.position="none",
	plot.title=element_text(hjust=0.5, colour="black", size=12, face="bold"),
	legend.text = element_text(colour="black", size=12, face="bold"),
	legend.title = element_blank(),
	legend.direction = "horizontal",
	axis.title.x = element_text(colour="black", size=11, face="bold"),
	axis.title.y = element_text(colour="black", size=11, face="bold", margin = margin(r = 10)),
	axis.text.x = element_text(colour="black", size=11, face="bold", angle = 45, hjust = 1),
	axis.text.y = element_text(colour="black", size=11, face="bold"),
	plot.margin=unit(c(0.25,0.25,0.25,0.25), "lines"))
#png("./Brassica_60k_SNP_Array/Known_Sites_Physical_Map.png", width=3000, height=1500, type="cairo", res=300)
Known_Sites_Map  
#dev.off()

### Truth Set Sites
Truth_Set_Map <- ggplot(truth_set, aes(x = POS, y = CHROM)) +
geom_point(color = "dodgerblue4", size = 1) +
scale_x_continuous(breaks=c(0,10000000,20000000,30000000,40000000,50000000,60000000,70000000),
	labels=c("0", "10","20","30","40","50","60","70")) +
annotate("segment", x=0, y=1, xend=66465249, yend=1, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=2, xend=41681856, yend=2, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=3, xend=55656957, yend=3, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=4, xend=50218839, yend=4, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=5, xend=56382805, yend=5, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=6, xend=65837619, yend=6, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=7, xend=73669886, yend=7, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=8, xend=62297340, yend=8, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=9, xend=48239358, yend=9, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=10, xend=20778245, yend=10, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=11, xend=53549826, yend=11, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=12, xend=26309499, yend=12, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=13, xend=29390523, yend=13, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=14, xend=45146386, yend=14, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=15, xend=42112164, yend=15, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=16, xend=23101715, yend=16, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=17, xend=39685748, yend=17, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=18, xend=33432960, yend=18, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=19, xend=32958928, yend=19, linewidth=0.5, color="black", linetype=1)+
labs(
title = "Physical Map of the 33k Truth Set",
x = "Physical Position (Mbp)",
y = "Chromosome"
) +
theme_gray() +
theme(
	panel.grid.major = element_line(color = "gray80"),
	panel.grid.minor = element_blank(),
	legend.position="none",
	plot.title=element_text(hjust=0.5, colour="black", size=12, face="bold"),
	legend.text = element_text(colour="black", size=12, face="bold"),
	legend.title = element_blank(),
	legend.direction = "horizontal",
	axis.title.x = element_text(colour="black", size=11, face="bold"),
	axis.title.y = element_text(colour="black", size=11, face="bold", margin = margin(r = 10)),
	axis.text.x = element_text(colour="black", size=11, face="bold", angle = 45, hjust = 1),
	axis.text.y = element_text(colour="black", size=11, face="bold"),
	plot.margin=unit(c(0.25,0.25,0.25,0.25), "lines"))
#png("./VCF/Round_1/Figures/Truth_Set_Physical_Map.png", width=3000, height=1500, type="cairo", res=300)
Truth_Set_Map  
#dev.off()

### Subset Sites
Subset_Map <- ggplot(subset, aes(x = POS, y = CHROM)) +
geom_point(color = "dodgerblue4", size = .5) +
scale_x_continuous(breaks=c(0,10000000,20000000,30000000,40000000,50000000,60000000,70000000),
	labels=c("0", "10","20","30","40","50","60","70")) +
annotate("segment", x=0, y=1, xend=66465249, yend=1, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=2, xend=41681856, yend=2, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=3, xend=55656957, yend=3, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=4, xend=50218839, yend=4, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=5, xend=56382805, yend=5, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=6, xend=65837619, yend=6, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=7, xend=73669886, yend=7, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=8, xend=62297340, yend=8, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=9, xend=48239358, yend=9, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=10, xend=20778245, yend=10, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=11, xend=53549826, yend=11, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=12, xend=26309499, yend=12, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=13, xend=29390523, yend=13, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=14, xend=45146386, yend=14, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=15, xend=42112164, yend=15, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=16, xend=23101715, yend=16, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=17, xend=39685748, yend=17, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=18, xend=33432960, yend=18, linewidth=0.5, color="black", linetype=1)+
annotate("segment", x=0, y=19, xend=32958928, yend=19, linewidth=0.5, color="black", linetype=1)+
labs(
title = " ",
x = "Length (Mbp)",
y = "Chromosome"
) +
theme_gray() +
theme(
	panel.grid.major = element_line(color = "gray80"),
	panel.grid.minor = element_blank(),
	legend.position="none",
	plot.title=element_text(hjust=0.5, colour="black", size=12, face="bold"),
	legend.text = element_text(colour="black", size=12, face="bold"),
	legend.title = element_blank(),
	legend.direction = "horizontal",
	axis.title.x = element_text(colour="black", size=11, face="bold"),
	axis.title.y = element_text(colour="black", size=11, face="bold", margin = margin(r = 10)),
	axis.text.x = element_text(colour="black", size=11, face="bold", angle = 45, hjust = 1),
	axis.text.y = element_text(colour="black", size=11, face="bold"),
	plot.margin=unit(c(0.25,0.25,0.25,0.25), "lines"))
#png("./VCF/Round_1/Figures/Clustering_Subset_Physical_Map.png", width=2500, height=1500, type="cairo", res=300)
Subset_Map  
#dev.off()





