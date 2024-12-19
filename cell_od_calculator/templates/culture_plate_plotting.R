#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~				  Resurrect Bio Cell Culture Ouptut Plotting              ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <miles@resurrect.bio>
## Last Update: 20/11/24

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
  make_option(c("--p1"), action="store_true", default=FALSE,
    help="Just plot plate 1 [default: %default]"),
  make_option(c("--p2"), action="store_true", default=FALSE,
    help="Just plot plate 2 [default: %default]"),
  make_option(c("--p3"), action="store_true", default=FALSE,
    help="Just plot plate 3 [default: %default]")
)

parser <- OptionParser(option_list = option_list,
               description = "This script processes platereader data from in a template file.\n
                              \nNote: The output file will be placed in the same directory as the input file.")

## For testing
#opt <- ""
#opt$input <- "C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_od_calculator/templates/example1.csv"
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

## Plotting
p1 <- ggplot(data = p1_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.5) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")
p2 <- ggplot(data = p2_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.5) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")
p3 <- ggplot(data = p3_melted, aes(x = Effector, y = as.numeric(value), colour = effector_num)) +
		geom_point(size = 3, alpha = 0.5) +
		theme_minimal(base_size = 20) +
		labs(y = "Luminosity", x = NULL) +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none")

## Adjustable parameters
if(opt$log_trans == TRUE){
	p1 <- p1 + scale_y_log10(expand = c(0,0))
	p2 <- p2 + scale_y_log10(expand = c(0,0))
	p3 <- p3 + scale_y_log10(expand = c(0,0))
	newname <- gsub(newname, pattern = ".jpg", replacement = "_lt.jpg")
} else {
	p1 <- p1 + scale_y_continuous(expand = c(0,0))
	p2 <- p2 + scale_y_continuous(expand = c(0,0))
	p3 <- p3 + scale_y_continuous(expand = c(0,0))
}
if(opt$add_mean == TRUE){
	p1 <- p1 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	p2 <- p2 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	p3 <- p3 + stat_summary(fun = mean, aes(group = 1),geom = "line", colour = "black")
	newname <- gsub(newname, pattern = ".jpg", replacement = "_ml.jpg")
}

if(opt$p1 == TRUE){
	newname_p1 <- gsub(newname, pattern = ".jpg", replacement = "_p1.jpg")
	ggsave(newname_p1, p1, width = 200, height = 112, dpi = 300, units = "mm")
}
if(opt$p2 == TRUE){
	newname_p2 <- gsub(newname, pattern = ".jpg", replacement = "_p2.jpg")
	ggsave(newname_p2, p2, width = 200, height = 112, dpi = 300, units = "mm")	
}
if(opt$p3 == TRUE){
	newname_p3 <- gsub(newname, pattern = ".jpg", replacement = "_p3.jpg")
	ggsave(newname_p3, p3, width = 200, height = 112, dpi = 300, units = "mm")	
}
if(opt$p1 == FALSE & opt$p2 == FALSE & opt$p3 == FALSE){
	combined_plot <- p1 / p2 / p3
	ggsave(newname, combined_plot, width = 200, height = 350, dpi = 300, units = "mm")
}