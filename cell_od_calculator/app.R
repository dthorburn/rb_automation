#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~				  Resurrect Bio Cell Culture AgroPrep OD Calculator for FeliX automation              ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <miles@resurrect.bio>
## Last Update: 20/11/24

## Loading in required libraries
if (!require("DT", quietly = TRUE))
	install.packages("DT")
if (!require("shiny", quietly = TRUE))
	install.packages("shiny")
if (!require("dplyr", quietly = TRUE))
	install.packages("dplyr")
if (!require("ggplot2", quietly = TRUE))
	install.packages("ggplot2")
if (!require("openxlsx", quietly = TRUE))
	install.packages("openxlsx")
if (!require("data.table", quietly = TRUE))
	install.packages("data.table")

suppressMessages(library("data.table"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("shiny"))
suppressMessages(library("DT"))

## shinyapps.io through rsconnect -- This is the one for https://dthorburn.shinyapps.io/cell_od_calculator/
#rsconnect::deployApp('C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_od_calculator')
## For debugging 
##setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/cell_od_calculator")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~ Define UI
ui <- fluidPage(
	titlePanel("RB Cell Culture OD Calculator: A tool to calculate aliquot volumes for the FeliX."),
	sidebarLayout(
		sidebarPanel(
			width = 3,
			h3("Mandatory input data:"),
			fileInput("xslx_file", "Platereader XSLX File", multiple = FALSE),
			
			h3("Optional input data:"),
			fileInput("od_csv", "Target OD CSV upload", multiple = FALSE),
			radioButtons("in_od_table_source", "Source of target OD table:",
				choices = list("Homogeneous table" = 1, "Upload CSV" = 0), selected = 1, inline = TRUE),

			h3("Calculation parameters:"),
			textInput("run_name", "Name of Run:", value = "", placeholder = "psojae_batch1"),
 			selectInput("num_out_plates", "Number of output plates?", choices = c("1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5"), selected = "1"),
 			selectInput("to_platereader", "Add volume for output platereader test?", choices = c("Yes" = "Yes", "No" = "No"), selected = "Yes"),
			textInput("last_well", "Last well in use:", value = "A12"),
			textInput("target_ul", "Target volume per replicate (ul):", value = "200"),
			textInput("min_vol", "Minimum volume of stock to dilute (ul):", value = "25"),
			textInput("blank_mean", "Blank mean:", value = "0.04"),
			textInput("drift", "Platereader drift:", value = "3.95"),
			textInput("max_vol", "Max output volume per well or channel (ul):", value = "1000"),
 			selectInput("remove_control_wells", "Remove +/- control wells from calculation?", choices = c("Yes" = "Yes", "No" = "No"), selected = "Yes"),
			radioButtons("default_od", "Default target OD for homogeneous table:",
				choices = list("0.05" = "0.05", "0.1" = "0.1", "0.2" = "0.2", "0.3" = "0.3"), selected = "0.1", inline = TRUE),

			hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
			helpText("To report bugs or request functions to be added, please contact Miles Thorburn <miles@resurrect.bio>")
		),
		
		mainPanel(
			tabsetPanel(
				tabPanel("Target ODs",
					hr(),
					h3("Target OD table:"),
					h5("Note: Ignore values in unused wells."),
					textOutput('num_plate_output'),
					tags$head(
						tags$style("#num_plate_output{color: red;
							font-size: 20px;
						}")
					),
					hr(),
					DTOutput("all_od_tabs"),
				),

				tabPanel("Platereader table",
					hr(),
					h3("Input XLSX table from platereader:"),
					h5("Estimated ODs of platereader data based on current drift level."),
					DTOutput("fin_od_tab"),
					hr(),
					h5("The raw platereader plate input file. Mostly for debugging purposes."),
					DTOutput("pr_tab")

				),

				tabPanel("Output",
					hr(),
					h3("Summary Output Table:"),
					downloadButton("downloadsumData", "Download FeliX Input Table"),
					textOutput('warns'),
					tags$head(
						tags$style("#warns{color: red;
							font-size: 20px;
							font-style: bold;
						}"),
					),
					hr(),
					DTOutput("out_tab_summ"),
					hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
					h3("Full Output Table:"),
					downloadButton("downloadfullData", "Download Full Output Table"),
					hr(),
					DTOutput("out_tab_full")
				),
				
				tabPanel("Instructions",
					hr(),
					h3("How to use the OD calculator:"),
					h5("1. Select local XLSX file emitted from the platereader in the shared lab."),
					h5("2. Add a run name, and if necessary, change any of the default calculation parameters."),
					h5("3. Alter editable target OD table, upload target OD table, or use homogeneous target OD table."),
					h5("4. Click on Output tab and inspect data."),
					h5("5. Download summary table which is the input for FeliX."),
					hr(),
					h3("Important Considerations:"),
					h5("1. Wells cannot be skipped. These will have to be manually removed after FeliX is complete."),
					h5("2. Due to how FeliX operated, wells need to be filled from H -> A for each column."),
					h5("3. If you forget to adjust order or need to skip wells, you can ignore those wells, only 40ul of AI will be dispensed into them."),
				),
				
				tabPanel("Debug",
					DTOutput("debug_out_tab"),
					textOutput("debug_out_print")
				),
			),
		)
	)
)

#shinyApp(ui = ui, server = server)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define server logic
server <- function(input, output) {
	## ~~ Setting up callable reactive variables, in same order as selection above. 
	run_name         <- reactive({ as.character(trimws(input$run_name)) })
	num_out_plates   <- reactive({ as.numeric(input$num_out_plates) })
	to_platereader   <- reactive({ as.character(input$to_platereader) })
	last_well        <- reactive({ as.character(trimws(input$last_well)) })
	target_ul        <- reactive({ as.numeric(trimws(input$target_ul)) })
	min_vol          <- reactive({ as.numeric(trimws(input$min_vol)) })
	blank_mean       <- reactive({ as.numeric(input$blank_mean) })
	drift            <- reactive({ as.numeric(trimws(input$drift)) })
	max_vol          <- reactive({ as.numeric(trimws(input$max_vol)) })
	dod              <- reactive({ as.numeric(input$default_od) })
	rcw              <- reactive({ as.character(input$remove_control_wells) })
	in_od_table_val  <- reactive({ as.numeric(input$in_od_table_source) })

	## ~~ Step 1: Error Checking for number of output plates
	total_volume <- reactive({
		output_vol <- num_out_plates() * target_ul()
		if(to_platereader() == "Yes"){
			output_vol <- output_vol + 200
		}
		as.numeric(output_vol)
	})
	output$num_plate_output <- renderText({
		if(total_volume() > max_vol()){
			paste0("WARN! Total volumes per well outside of acceptable range: ", total_volume(), "ul. Reduce number of output plates (including platereader plate).")
		} else {
			paste0("")
			#paste0("Total Output volume per well acceptable: ", total_volume())
		}
	})

	## ~~ Step 2: Reading in platereader output table
	xlsx_in <- reactive({
		req(input$xslx_file)
		df <- read.xlsx(input$xslx_file$datapath,
					startRow = 8, 
					colNames = FALSE) %>% as.data.table
		names(df) <- c("Row", paste0("V", 1:12))
		df
	})
	output$pr_tab <- renderDT({
		datatable(xlsx_in(), rownames= FALSE)
	})

	## ~~ Step 3: Target OD tables
	read_input_ods <- reactive({
		req(input$od_csv)
		df <- fread(input$od_csv$datapath)
		num_rows <- nrow(df)
		num_cols <- ncol(df)
		if(num_rows != 8){
			message(paste0("WARN: Incorrect number of rows identified - ", num_rows, "identified. There should be 8."))
		}
		if(num_cols != 12){
			message(paste0("WARN: Incorrect number of columns identified - ", num_cols, "identified. There should be 12."))
		}
		df
	})
	set_default_tab <- reactive({
		dod1 <- dod()
		od_table <- data.table( Row = letters[1:8] %>% toupper(),
					V1 = rep(dod1, 8), V2 = rep(dod1, 8),
					V3 = rep(dod1, 8), V4 = rep(dod1, 8),
					V5 = rep(dod1, 8), V6 = rep(dod1, 8),
					V7 = rep(dod1, 8), V8 = rep(dod1, 8),
					V9 = rep(dod1, 8), V10 = rep(dod1, 8),
					V11 = rep(dod1, 8), V12 = rep(dod1, 8))
		od_table
	})
	output$all_od_tabs <- renderDT({
		if(input$in_od_table_source == 0){
			read_input_ods() %>% datatable(., rownames = FALSE)
		} else if(input$in_od_table_source == 1){
			set_default_tab() %>% datatable(., rownames = FALSE)
		}
	})


	## ~~ Step 4: Main OD calculation function
	drift_cor <- reactive({		
		## Removing redundant data
		dat <- xlsx_in()[!c(9,10),!1] %>% as.data.frame
		
		## Removing daily blank value from wells - Handled differently from previous version
		b_mean <- blank_mean() %>% as.numeric
		no_bs <- dat - b_mean

		## Adjust for drift
		temp_drift <- drift() %>% as.numeric
		w_drift <- no_bs * temp_drift
		w_drift %>% as.data.table
	})
	calc_dat <- reactive({
		w_drift <- drift_cor()
		## Manipulating data 
		w_drift[,"row_name" := letters[1:8] %>% toupper()]
		dat_melted <- melt(w_drift, id.vars=c("row_name")) 
		names(dat_melted) <- c("row_name", "col_name", "diluted_plate_od_dc")
		dat_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]

 		# Applying df to OD value
		iotv <- in_od_table_val() %>% as.numeric
		if(iotv == 0){
			ods <- read_input_ods()
		} else if(iotv == 1){
			ods <- set_default_tab()
		}
		ods[,"Row"] <- NULL
		ods[,"row_name" := letters[1:8] %>% toupper()]
		ods_melted <- melt(ods, id.vars=c("row_name")) 
		names(ods_melted) <- c("row_name", "col_name", "target_od")
		ods_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]
		
		## Merging melted datasets
		dat <- merge(dat_melted, ods_melted, by = c("row_name", "col_name", "stock_well"))
		
		## Ordering the data for single pipette work. 
		setorderv(dat, c("col_name", "row_name"), c(1, -1))

		## Calculating the volumes
		calc_stock <- function(row_num){
			temp_val <- (as.numeric(dat$target_od[row_num]) / as.numeric(dat$diluted_plate_od_dc[row_num])) * total_volume()
			round(temp_val, digits = 1) %>% return
		}
		dat[,"stock_vol" := calc_stock(1:nrow(dat))]
		dat[,"ai_vol" := total_volume() - stock_vol %>% round(., digits = 1)]
		dat
	})
	output$out_tab_full <- renderDT({
		datatable(calc_dat(), rownames= FALSE) %>% formatRound(c(4), 3) %>% formatRound(c(6,7), 3)
	})
	output$downloadfullData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_full_table_", format(Sys.Date(), "%d%m%Y"), ".csv")
		},
		content = function(file) {
			write.csv(calc_dat(), file, row.names = FALSE, quote = FALSE)
		}
	)
	output$fin_od_tab <- renderDT({
		w_drift <- drift_cor()
		w_drift[,"Row" := letters[1:8] %>% toupper()]
		datatable(w_drift[,c(1:12)], rownames= FALSE) %>% formatRound(c(2:12), 3)
	})

	## ~~ Step 6: Creating FeliX input table by adjusting values based on last 
	adjust_dat <- reactive({
		dat <- calc_dat() %>% as.data.table
		lw <- last_well() %>% as.character
		lw_char <- paste0(gsub(pattern = "([A-z]?)([0-9]+)([A-z]?)", replacement = "\\2", lw),
			gsub(pattern = "([A-z]?)([0-9]+)([A-z]?)", replacement = "\\1\\3", lw)
		)
		lw_num <- which(dat$stock_well == lw_char)
		dat[1:lw_num,c("stock_well", "stock_vol", "ai_vol")]		
	})

	## ~~ Step 7: Adjusting wells with out of range values and reporting.
	fix_anoms <- reactive({
		dat <- adjust_dat()
		adjust_vols <- function(row_num){
			temp_row <- dat[row_num,]
			## Adjusting if values higher than limits
			temp_row[,"stock_vol" := ifelse((stock_vol > total_volume()), 20, stock_vol)]
			temp_row[,"ai_vol" := ifelse((ai_vol > total_volume()), 20, ai_vol)]
			## Adjusting if values higher than limits
			temp_row[,"stock_vol" := ifelse((stock_vol < min_vol()), 20, stock_vol)]
			temp_row[,"ai_vol" := ifelse((ai_vol > max_vol()), 20, ai_vol)]
			## Adjusting if values are negative
			temp_row[,"stock_vol" := ifelse((stock_vol < 0), 20, stock_vol)]
			temp_row[,"ai_vol" := ifelse((ai_vol < 0), 20, ai_vol)]
			# Adjusting second volume if first was changed
			if(temp_row$stock_vol == 20){
				temp_row[,"ai_vol" := 20]
			} else if(temp_row$ai_vol == 20){
				temp_row[,"stock_vol" := 20]
			}
			return(temp_row)
		}
		out <- lapply(FUN = adjust_vols, 1:nrow(dat)) %>% rbindlist
		if(rcw() == "Yes"){
			target_wells <- c("4G", "4H", "8G", "8H", "12G", "12H")
			out[stock_well %in% target_wells, `:=`(stock_vol = 0, ai_vol = 0)] ## Can this be NULL instead? Does FeliX skip these?
		}
		out
	})
	anoms_ids <- reactive({
		fin_dat <- fix_anoms() %>% as.data.table
		errs <- subset(fin_dat, stock_vol == 20 & ai_vol == 20)
		errs
	})
	output$warns <- renderText({
		errs_out <- anoms_ids() %>% as.data.table
		paste0("WARN! Volumes outside of acceptable range for wells: ", paste(errs_out$stock_well, collapse = " "))
	})
	output$out_tab_summ <- renderDT({
		datatable(fix_anoms(), rownames= FALSE, options = list(
			paging =TRUE,
			pageLength = 8
		)) %>% formatRound(c(2:3), 1)
	})
	output$downloadsumData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_", format(Sys.Date(), "%d%m%Y"), ".txt")
		},
		content = function(file) {
			write.csv(fix_anoms(), file, row.names = FALSE, quote = FALSE)
		}
	)
}

shinyApp(ui = ui, server = server)
