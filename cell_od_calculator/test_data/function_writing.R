suppressMessages(library("data.table"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("shiny"))
suppressMessages(library("DT"))

setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_od_calculator/test_data")
dir()
dat <- fread("example_felix_input_28012025.txt")
dat[, "col" := gsub(stock_well, pattern = "^([0-9][0-9]?)([A-H])", replacement = "\\1") %>% as.numeric]
dat[, "row" := gsub(stock_well, pattern = "^([0-9][0-9]?)([A-H])", replacement = "\\2")]
setorderv(dat, c("col", "row"), c(1,1))

br_num <- "3"

beckman_stock_out <- data.table(
	"Source" = paste0("source1.", br_num),
	"Source well" = dat$stock_well,
	"Destination" = paste0("destination1.", br_num),
	"Destination well" = dat$stock_well,
	"volume" = dat$stock_vol %>% round(.,0))

beckman_ai_out <- data.table(
	"Source" = ifelse(br_num < 3, "source2", "source3"),
	"Source well" = 1,
	"Destination" = paste0("destination1.", br_num),
	"Destination well" = dat$stock_well,
	"volume" = dat$ai_vol %>% round(.,0))