suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))

setwd("C:/Users/miles/Documents/Resurrect_Bio/Projects/03_robot/01_od_protocol")
dir()

calc_stock <- function(row_num){
	temp_val <- (dat$target_od[row_num] / dat$diluted_plate_od_dc[row_num]) * target_ul
	round(temp_val, digits = 3) %>% return
}

## Paramaters 
run_name  <- "FeliX_test2_pl1"
blank_col <- 12
og_df     <- 5               ## original diluation factor. ## Change to 1:5
target_ul <- 1000
drift     <- 3.95            ## Plate reads :: spectrophotomoter drift factor. 5.70, new drift 23.01.24

## Reading in files
ods <- fread("example_ods.csv")
prr <- read.xlsx("C:/Users/miles/Documents/Resurrect_Bio/Projects/03_robot/01_od_protocol/Endpoint Abs PL1 @ 600 2024-01-23 15-47-50.xlsx",
					startRow = 8, 
					colNames = FALSE) %>% as.data.table
## Changing the names for easier manipulation later
names(prr) <- c("row_name", paste0("V", 1:(ncol(prr)-1)))

## Removing redundant rows
dat <- prr[!c(9,10),!1]
#dat <- prr[!c(9,10),!c(1,11)]



## Removing blank rows and getting original concentration
no_bs <- dat[,!..blank_col] - mean(dat[,..blank_col] %>% unlist)
no_bs <- no_bs[,1]
## Adjust for drift
w_drift <- no_bs * drift
#og_co <- w_drift * og_df  ## Not taking from this original stock well anymore.
#og_co <- w_drift

## Manipulating data for FeliX
w_drift[,"row_name" := letters[1:8] %>% toupper()]
dat_melted <- melt(w_drift, id.vars=c("row_name")) 
names(dat_melted) <- c("row_name", "col_name", "diluted_plate_od_dc")
dat_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]


ods[,"row_name" := letters[1:8] %>% toupper()]
ods_melted <- melt(ods, id.vars=c("row_name")) 
names(ods_melted) <- c("row_name", "col_name", "target_od")
ods_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]

dat <- merge(dat_melted, ods_melted, by = c("row_name", "col_name", "stock_well"))
## Ordering the data for single pipette work. 
setorderv(dat, c("col_name", "row_name"), c(1, -1))

## Calculating the volumes
dat[,"stock_vol" := calc_stock(1:nrow(dat))]
dat[,"ai_vol" := target_ul - stock_vol]

## TO ADD:
## 1. Annotated problematic columns
## 2. Add a way to skip problematic cells

## Writing results
#fwrite(file = paste0(run_name, "_stock_ul_fromStock.csv"), round(ul_stock, digits = 3), col.names = FALSE)
#fwrite(file = paste0(run_name, "_water_ul_fromStock.csv"), round(ul_water, digits = 3), col.names = FALSE)

fwrite(file = paste0(run_name, "_volumes_", format(Sys.Date(), "%d%m%Y"), ".txt"), dat[,c("stock_well", "stock_vol", "ai_vol")])