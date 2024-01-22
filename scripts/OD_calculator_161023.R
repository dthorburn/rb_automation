suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))

setwd("C:/Users/miles/Documents/Resurrect_Bio/Projects/03_robot/01_od_protocol")
dir()

## Paramaters 
run_name  <- "solo_agroprep_run2"
blank_col <- 10
og_df     <- 5               ## original diluation factor. ## Change to 1:5
target_ul <- 1000
drift     <- 3.95            ## Plate reads :: spectrophotomoter drift factor.

## Reading in files
ods <- fread("example_ods.csv")
#prr <- read.xlsx("Endpoint Abs @ 600 2023-10-18 14-38-37.xlsx",
prr <- read.xlsx("Endpoint Abs @ 600 2023-12-13 12-18-11.xlsx",
					startRow = 8, 
					colNames = FALSE) %>% as.data.table

## Removing redundant rows
dat <- prr[!c(9,10),!1]
#dat <- prr[!c(9,10),!c(1,11)]

## Removing blank rows and getting original concentration
no_bs <- dat[,!..blank_col] - mean(dat[,..blank_col] %>% unlist)
w_drift <- no_bs * drift
#og_co <- w_drift * og_df  ## Not taking from this original stock well anymore.
og_co <- w_drift  
fwrite(file = paste0(run_name, "_adjustedODs.csv"),  round(og_co, digits = 3), col.names = FALSE)
## ul of stock
ul_stock <- (ods / og_co) * target_ul
## ul of water
ul_water <- target_ul - ul_stock

## Writing results
fwrite(file = paste0(run_name, "_stock_ul_fromStock.csv"), round(ul_stock, digits = 3), col.names = FALSE)
fwrite(file = paste0(run_name, "_water_ul_fromStock.csv"), round(ul_water, digits = 3), col.names = FALSE)


