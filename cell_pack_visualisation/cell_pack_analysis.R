#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~	    			    Resurrect Bio Cell Pack Ouptut Analysis                 ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <miles@resurrect.bio>
## Last Update: 05/02/25

## Workflow
## 1. Read in data and stop if no files found ✓
## 2. Process each file
##   i.    Parse and attach well and sample names ✓
##   ii.   Calcualte control wells ✓
##   iii.  Calculate every other well mean ✓ 
##   iv.   Calculate normalisation ✓
##   v.    Calculate z-score ✓
##   vi.   Calculate skew index
##   vii.  Calculate mean-of-means and median-of-means
##   viii. Highlight wells with deviations from expecations

## Loading libraries
options(repos = c(CRAN = "https://cran.ma.imperial.ac.uk/"))
## Loading libraries
if (suppressWarnings(!require("dplyr", quietly = TRUE, warn.conflicts = FALSE))) {
  install.packages("dplyr")
}
if (!require("optparse", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("optparse")
}
if (!require("data.table", quietly = TRUE, warn.conflicts = FALSE)) {
  install.packages("data.table")
}

## Handling data passed to R
option_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL,
    help="Input directory of csv files. Must be the provided template file."),
  make_option(c("--output_name"), type="character", default="rb_cell_pack_results",
    help="Output file name [default: %default]."),
  make_option(c("--pos_ctl_well"), type="character", default=31,
    help="Which number well is 1st positive control [default: %default]."),
  make_option(c("--neg_ctl_well"), type="character", default=29,
    help="Which number well is 1st negative control [default: %default]."),
)

parser <- OptionParser(option_list = option_list,
               description = "\nThis script processes platereader data in a template file format.\n
                              \nNote: The output file will be placed in the same directory as the input file.")
opt <- parse_args(parser)

## Functions
plate_mean <- function(temp_plate){
  ## removing controls - this is so hideous
  td_expt <- subset(temp_plate, !(well_num >= opt$neg_ctl_well & well_num <=32))
  t_out <- data.table(plate_mean = mean(td_expt$value %>% as.numeric),
                      plate_sd = sd(td_expt$value %>% as.numeric),
                      pos_ctl_mean = temp_plate[well_num == opt$pos_ctl_well, value] %>% as.numeric %>% mean,
                      pos_ctl_sd = temp_plate[well_num == opt$pos_ctl_well, value] %>% as.numeric %>% sd,
                      neg_ctl_mean = temp_plate[well_num == opt$neg_ctl_well, value] %>% as.numeric %>% mean,
                      neg_ctl_sd = temp_plate[well_num == opt$neg_ctl_well, value] %>% as.numeric %>% sd)
  return(t_out)
}

well_means <- function(temp_eff, temp_plate){
  temp_dat <- subset(temp_plate, Effector == temp_eff)
  well_nums <- paste0(unique(temp_dat$Row), gsub(pattern = "C", "", x = temp_dat$variable))
  temp_out <- data.table(
                file = gsub(pattern = "^.*/(.*).csv", "\\1", x = temp_file),
                effector = temp_eff,
                plate = ifelse(sum(temp_plate[1,] == p1_withz[1,]) == ncol(p1_withz), "plate1", 
                          ifelse(sum(temp_plate[1,] == p2_withz[1,]) == ncol(p1_withz), "plate2", 
                            ifelse(sum(temp_plate[1,] == p3_withz[1,]) == ncol(p1_withz), "plate3", "error"))),
                wells = paste0(well_nums, collapse = ","),
                well_type = ifelse(temp_dat$well_num[1] == opt$pos_ctl_well, "pos_ctl", 
                              ifelse(temp_dat$well_num[1] == opt$neg_ctl_well, "neg_ctl", 
                                ifelse(temp_dat$well_num[1] > opt$neg_ctl_well, "ctl", "expt"))),
                raw_mean = mean(temp_dat$value %>% as.numeric),
                raw_sd = sd(temp_dat$value %>% as.numeric),
                z_mean = mean(temp_dat$z_score %>% as.numeric),
                z_sd = sd(temp_dat$z_score %>% as.numeric),
                z_neg_mean = mean(temp_dat$z_score_negctl %>% as.numeric),
                z_neg_sd = sd(temp_dat$z_score_negctl %>% as.numeric),
                z_both_mean = mean(temp_dat$z_score_bothctls %>% as.numeric),
                z_both_sd = sd(temp_dat$z_score_bothctls %>% as.numeric))
  return(temp_out)
}

## For testing
#opt <- ""
#opt$input <- "C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_pack_visualisation/example_data"
#opt$pos_ctl_well <- 31
#opt$neg_ctl_well <- 29

## Checkpoint - File Numbers
file_list <- grep(pattern = "*.csv", x = dir(opt$input), value = TRUE)
if(length(file_list) > 0){
  message(paste0("Found ", length(file_list), " csvs in ", opt$input, "."))
} else {
  stop(paste0("Found 0 csv files in ", opt$input, ".\nExiting."))
}

temp_file <- paste0(opt$input, "/", file_list[1])
plate_names <- fread(temp_file, header = TRUE) %>% head(., n = 8)
plate_1 <- fread(temp_file, header = TRUE, skip = 9) %>% head(., n = 8)
plate_2 <- fread(temp_file, header = TRUE, skip = 18) %>% head(., n = 8)
plate_3 <- fread(temp_file, header = TRUE, skip = 27) %>% head(., n = 8)
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

p1_melted[, "well_num" := rep(1:32,3)]
p2_melted[, "well_num" := rep(1:32,3)]
p3_melted[, "well_num" := rep(1:32,3)]

## Only needed if plotting to check
p1_melted <- p1_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup() %>% as.data.table
p2_melted <- p2_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup() %>% as.data.table
p3_melted <- p3_melted %>%
  group_by(Effector) %>%
  mutate(effector_num = paste0(row_number())) %>%
  ungroup() %>% as.data.table


## Getting plate means and sd
p1_means <- plate_mean(p1_melted)
p2_means <- plate_mean(p2_melted)
p3_means <- plate_mean(p3_melted)

## First z-score
p1_withz <- p1_melted[,"z_score" := (as.numeric(value) - p1_means$plate_mean) / p1_means$plate_sd]
p2_withz <- p2_melted[,"z_score" := (as.numeric(value) - p2_means$plate_mean) / p2_means$plate_sd]
p3_withz <- p3_melted[,"z_score" := (as.numeric(value) - p3_means$plate_mean) / p3_means$plate_sd]

## z-score with negative controls
p1_withz[,"z_score_negctl" := (as.numeric(value) - p1_means$neg_ctl_mean) / p1_means$neg_ctl_sd]
p2_withz[,"z_score_negctl" := (as.numeric(value) - p2_means$neg_ctl_mean) / p2_means$neg_ctl_sd]
p3_withz[,"z_score_negctl" := (as.numeric(value) - p3_means$neg_ctl_mean) / p3_means$neg_ctl_sd]

## z-score with both controls
p1_withz[,"z_score_bothctls" := (as.numeric(value) - p1_means$neg_ctl_mean) / (p1_means$pos_ctl_mean - p1_means$neg_ctl_mean)]
p2_withz[,"z_score_bothctls" := (as.numeric(value) - p2_means$neg_ctl_mean) / (p2_means$pos_ctl_mean - p2_means$neg_ctl_mean)]
p3_withz[,"z_score_bothctls" := (as.numeric(value) - p3_means$neg_ctl_mean) / (p3_means$pos_ctl_mean - p3_means$neg_ctl_mean)]

## Getting all well means
p1_results <- lapply(FUN = well_means, temp_plate = p1_withz, X = p1_withz$Effector[1:32]) %>% rbindlist
p2_results <- lapply(FUN = well_means, temp_plate = p2_withz, X = p2_withz$Effector[1:32]) %>% rbindlist
p3_results <- lapply(FUN = well_means, temp_plate = p3_withz, X = p3_withz$Effector[1:32]) %>% rbindlist

