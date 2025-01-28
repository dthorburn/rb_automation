#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~				  Resurrect Bio Cell Culture Ouptut Plotting              ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <miles@resurrect.bio>
## Last Update: 27/01/25

## Loading libraries
options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk/"))
## Loading libraries
if (suppressWarnings(!require("dplyr", quietly = TRUE, warn.conflicts = FALSE))) {
  install.packages("dplyr")
}
if (!require("ggplot2", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("ggplot2")
}
if (!require("optparse", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("optparse")
}
if (!require("patchwork", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("patchwork")
}
if (!require("data.table", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("data.table")
}

## Handling data passed to R
option_list <- list(
  make_option(c("--input"), type="character", default=NULL,
    help="Input file in csv format. Must be the provided template file."),
  make_option(c("--log_trans"), action="store_true", default=FALSE,
    help="Log transform Y axis [default: %default]"),
  make_option(c("--add_mean"), action="store_true", default=FALSE,
    help="Adds mean line to plot [default: %default]"),
  make_option(c("--reorder_x"), action="store_true", default=FALSE,
    help="Readjusts X axis by row [default: %default]"),
  make_option(c("--draw_control_means"), action="store_true", default=FALSE,
    help="Adds horizontal mean lines for controls [default: %default]"),
  make_option(c("--pos_ctl_well"), type="character", default=31,
    help="Which number well is 1st positive control [default: %default]."),
  make_option(c("--neg_ctl_well"), type="character", default=29,
    help="Which number well is 1st negative control [default: %default]."),
  make_option(c("--p1"), action="store_true", default=FALSE,
    help="Just plot plate 1 [default: %default]"),
  make_option(c("--p2"), action="store_true", default=FALSE,
    help="Just plot plate 2 [default: %default]"),
  make_option(c("--p3"), action="store_true", default=FALSE,
    help="Just plot plate 3 [default: %default]")
)

parser <- OptionParser(option_list = option_list,
               description = "\nThis script processes platereader data from in a template file.\n
                              \nNote: The output file will be placed in the same directory as the input file.")

## For testing
#opt <- ""
#opt$input <- "C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_od_calculator/templates/plotting_input_Plates4_5_6.csv"
#opt$input <- "C:/Users/miles/Documents/Resurrect_Bio/Projects/06_MegaScreen/Results/csvs_250120/25.01.20_N32_EffectorPlate1.csv"
opt <- parse_args(parser)
newname <- opt$input %>% gsub(pattern = "csv", replacement = "jpg")

## Reading data
plate_names <- fread(opt$input, header = TRUE) %>% head(., n = 8)
plate_1 <- fread(opt$input, header = TRUE, skip = 9) %>% head(., n = 8)
plate_2 <- fread(opt$input, header = TRUE, skip = 18) %>% head(., n = 8)
plate_3 <- fread(opt$input, header = TRUE, skip = 27) %>% head(., n = 8)
names(plate_names)[1] <- "Row"
names(plate_1)[1] <- "Row"
names(plate_2)[1] <- "Row"
names(plate_3)[1] <- "Row"

pn_melted <- plate_names %>% melt(., id.vars = "Row")
p1_melted <- plate_1 %>% melt(., id.vars = "Row")
p2_melted <- plate_2 %>% melt(., id.vars = "Row")
p3_melted <- plate_3 %>% melt(., id.vars = "Row")

p1_melted[, "Effector" := pn_melted$value[ 1:32] %>% rep(., times = 3)]
p2_melted[, "Effector" := pn_melted$value[33:64] %>% rep(., times = 3)]
p3_melted[, "Effector" := pn_melted$value[65:96] %>% rep(., times = 3)]

## Preparing for plot
p1_melted <- p1_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup()
p2_melted <- p2_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup()
p3_melted <- p3_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup()

p1_melted <- p1_melted %>% mutate(Effector = factor(Effector, levels = unique(Effector)))
p2_melted <- p2_melted %>% mutate(Effector = factor(Effector, levels = unique(Effector)))
p3_melted <- p3_melted %>% mutate(Effector = factor(Effector, levels = unique(Effector)))

if(opt$reorder_x == TRUE){
	new_order <- c()
	for (i in 1:8) {
	  new_order <- c(new_order, seq(i, 28, by = 8))
	}
	new_order <- c(new_order, 29:32)
	p1_melted$Effector <- factor(p1_melted$Effector, levels = c(p1_melted$Effector[new_order]))
	p2_melted$Effector <- factor(p2_melted$Effector, levels = c(p2_melted$Effector[new_order]))
	p3_melted$Effector <- factor(p3_melted$Effector, levels = c(p3_melted$Effector[new_order]))
}

## Plotting
p1 <- ggplot(data = p1_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.85) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + 
		coord_cartesian(clip = "off") 

p2 <- ggplot(data = p2_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.85) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + 
		coord_cartesian(clip = "off") 

p3 <- ggplot(data = p3_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.75) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") + 
		coord_cartesian(clip = "off") 

## Adjustable parameters
if(opt$reorder_x == TRUE){
#	p1.1 <- p1 + geom_vline(xintercept = 1:4, linetype = "dashed", color = "grey30") +
#		 geom_vline(xintercept = 5:8, linetype = "dashed", color = "grey60") 
#	p2.1 <- p2 + geom_vline(xintercept = 1:4, linetype = "dashed", color = "grey30") +
#		 geom_vline(xintercept = 5:8, linetype = "dashed", color = "grey60") 
#	p3.1 <- p3 + geom_vline(xintercept = 1:4, linetype = "dashed", color = "grey30") +
#		 geom_vline(xintercept = 5:8, linetype = "dashed", color = "grey60") 
	p1.1 <- p1 + annotate("text", x = seq_along(unique(p1_melted$Effector)), y = max(as.numeric(p1_melted$value)) + 1,
	           label = c(rep(seq(1, 4), each = 4), rep(seq(5,8), each = 3), 5:8), vjust = -1, hjust = 0.5, size = 2.5)
	p2.1 <- p2 + annotate("text", x = seq_along(unique(p2_melted$Effector)), y = max(as.numeric(p2_melted$value)) + 1,
	           label = c(rep(seq(1, 4), each = 4), rep(seq(5,8), each = 3), 5:8), vjust = -1, hjust = 0.5, size = 2.5)
	p3.1 <- p3 + annotate("text", x = seq_along(unique(p3_melted$Effector)), y = max(as.numeric(p3_melted$value)) + 1,
	           label = c(rep(seq(1, 4), each = 4), rep(seq(5,8), each = 3), 5:8), vjust = -1, hjust = 0.5, size = 2.5)
	newname <- gsub(newname, pattern = ".jpg", replacement = "_rx.jpg")
} else {
	p1.1 <- p1 + geom_vline(xintercept = c(1, 9, 17, 25), linetype = "dashed", color = "grey30") +
		 geom_vline(xintercept = c(8, 16, 24, 32), linetype = "dashed", color = "grey60") 
	p2.1 <- p2 + geom_vline(xintercept = c(1, 9, 17, 25), linetype = "dashed", color = "grey30") +
		 geom_vline(xintercept = c(8, 16, 24, 32), linetype = "dashed", color = "grey60") 
	p3.1 <- p3 + geom_vline(xintercept = c(1, 9, 17, 25), linetype = "dashed", color = "grey30") +
		 geom_vline(xintercept = c(8, 16, 24, 32), linetype = "dashed", color = "grey60") 
	p1.1 <- p1.1 + annotate("text", x = seq_along(unique(p1_melted$Effector)), y = max(as.numeric(p1_melted$value)) + 1,
	           label = rep(seq(1,8),4), vjust = -1, hjust = 0.5, size = 2.5)
	p2.1 <- p2.1 + annotate("text", x = seq_along(unique(p2_melted$Effector)), y = max(as.numeric(p2_melted$value)) + 1,
	           label = rep(seq(1,8),4), vjust = -1, hjust = 0.5, size = 2.5)
	p3.1 <- p3.1 + annotate("text", x = seq_along(unique(p3_melted$Effector)), y = max(as.numeric(p3_melted$value)) + 1,
		           label = rep(seq(1,8),4), vjust = -1, hjust = 0.5, size = 2.5)
}

if(opt$draw_control_means == TRUE){
	newname <- gsub(newname, pattern = ".jpg", replacement = "_cm.jpg")
	p1_pos_ctl <- subset(p1_melted, Effector == p1_melted$Effector[opt$pos_ctl_well %>% as.numeric])
	p1_neg_ctl <- subset(p1_melted, Effector == p1_melted$Effector[opt$neg_ctl_well %>% as.numeric])
	p1.1 <- p1.1 + geom_hline(yintercept = mean(p1_pos_ctl$value %>% as.numeric), linetype = "dashed", color = "chartreuse3") +
		geom_hline(yintercept = mean(p1_neg_ctl$value %>% as.numeric), linetype = "dashed", color = "tomato1") 
	p2_pos_ctl <- subset(p2_melted, Effector == p2_melted$Effector[opt$pos_ctl_well %>% as.numeric])
	p2_neg_ctl <- subset(p2_melted, Effector == p2_melted$Effector[opt$neg_ctl_well %>% as.numeric])
	p2.1 <- p2.1 + geom_hline(yintercept = mean(p2_pos_ctl$value %>% as.numeric), linetype = "dashed", color = "chartreuse3") +
		geom_hline(yintercept = mean(p2_neg_ctl$value %>% as.numeric), linetype = "dashed", color = "tomato1") 
	p3_pos_ctl <- subset(p3_melted, Effector == p3_melted$Effector[opt$pos_ctl_well %>% as.numeric])
	p3_neg_ctl <- subset(p3_melted, Effector == p3_melted$Effector[opt$neg_ctl_well %>% as.numeric])
	p3.1 <- p3.1 + geom_hline(yintercept = mean(p3_pos_ctl$value %>% as.numeric), linetype = "dashed", color = "chartreuse3") +
		geom_hline(yintercept = mean(p3_neg_ctl$value %>% as.numeric), linetype = "dashed", color = "tomato1") 
}

if(opt$log_trans == TRUE){
	p1.1 <- p1.1 + scale_y_log10(expand = c(0,0))
	p2.1 <- p2.1 + scale_y_log10(expand = c(0,0))
	p3.1 <- p3.1 + scale_y_log10(expand = c(0,0))
	newname <- gsub(newname, pattern = ".jpg", replacement = "_lt.jpg")
} else {
	p1.1 <- p1.1 + scale_y_continuous(expand = c(0,0))
	p2.1 <- p2.1 + scale_y_continuous(expand = c(0,0))
	p3.1 <- p3.1 + scale_y_continuous(expand = c(0,0))
}
if(opt$add_mean == TRUE){
	p1.1 <- p1.1 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	p2.1 <- p2.1 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	p3.1 <- p3.1 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	newname <- gsub(newname, pattern = ".jpg", replacement = "_ml.jpg")
}

if(opt$p1 == TRUE){
	newname_p1 <- gsub(newname, pattern = ".jpg", replacement = "_p1.jpg")
	ggsave(newname_p1, p1.1, width = 200, height = 112, dpi = 300, units = "mm")
}
if(opt$p2 == TRUE){
	newname_p2 <- gsub(newname, pattern = ".jpg", replacement = "_p2.jpg")
	ggsave(newname_p2, p2.1, width = 200, height = 112, dpi = 300, units = "mm")	
}
if(opt$p3 == TRUE){
	newname_p3 <- gsub(newname, pattern = ".jpg", replacement = "_p3.jpg")
	ggsave(newname_p3, p3.1, width = 200, height = 112, dpi = 300, units = "mm")	
}
if(opt$p1 == FALSE & opt$p2 == FALSE & opt$p3 == FALSE){
  combined_plot <- p1.1 / p2.1 / p3.1
	ggsave(newname, combined_plot, width = 200, height = 350, dpi = 300, units = "mm")
}