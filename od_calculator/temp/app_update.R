#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~				  Resurrect Bio AgroPrep OD Calculator for FeliX automation						  ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <dthorburn@imperial.ac.uk>
## Last Update: 30/01/24

## Loading in required libraries
if (!require("DT", quietly = TRUE))
	install.packages("DT")
if (!require("shiny", quietly = TRUE))
	install.packages("shiny")
if (!require("dplyr", quietly = TRUE))
	install.packages("dplyr")
#if (!require("httpuv", quietly = TRUE))
#	install.packages("httpuv")
if (!require("ggplot2", quietly = TRUE))
	install.packages("ggplot2")
if (!require("openxlsx", quietly = TRUE))
	install.packages("openxlsx")
#if (!require("rsconnect", quietly = TRUE))
#	install.packages('rsconnect')
#if (!require("shinylive", quietly = TRUE))
#	install.packages("shinylive")
if (!require("data.table", quietly = TRUE))
	install.packages("data.table")

suppressMessages(library("data.table"))
#suppressMessages(library("shinylive"))
#suppressMessages(library("rsconnect"))
suppressMessages(library("openxlsx"))
suppressMessages(library("ggplot2"))
#suppressMessages(library("httpuv"))
suppressMessages(library("dplyr"))
suppressMessages(library("shiny"))
suppressMessages(library("DT"))


## shinyapps.io through rsconnect
#rsconnect::deployApp('C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/od_calculator')

## Just for setting the working dir when adjusting things. 
##setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/od_calculator")
##logo_src <- "https://resurrect.bio/assets/images/logo/colour.png"

## Setting up the shinylive server
#setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation")
#shinylive::export(appdir = "od_calculator", destdir = "od_calculator_shinylive")
## Checking it works locally
#httpuv::runStaticServer("od_calculator_shinylive", port=8008)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~ Define UI
#shinyUI(fluidPage(
ui <- fluidPage(
	titlePanel("RB OD Calculator: A tool to calculate aliquot volumes for the FeliX."),
	sidebarLayout(
		sidebarPanel(
			width = 3,
			h3("Mandatory input data:"),
			fileInput("xslx_file", "Platereader XSLX File", multiple = FALSE),
			
			h3("Optional input data:"),
			fileInput("od_csv", "Target OD CSV upload", multiple = FALSE),
			
			h3("Calculation parameters:"),
			textInput("run_name", "Name of Run:", value = "", placeholder = "N77_AA13"),
			textInput("blank_col", "Blank column:", value = "12"),
			textInput("target_ul", "Target volume (ul):", value = "1000"),
			textInput("drift", "Platereader drift:", value = "3.95"),
			textInput("min_vol", "Max stock volume to dilute (ul):", value = "975"),
			textInput("max_vol", "Max output volume per well or channel (ul):", value = "1000"),
 			selectInput("in_type", "Infection assay?", choices = c(Yes = "Yes", No = "No")),
			#radioButtons("in_type", h4("Infection assay?"),
			#			  choices = list("Yes" = 0, "No" = 1), selected = character(0), inline = TRUE),
			## For the life of me I cannot figure out why this is not working...
			textInput("last_well", "Last well in use:", value = "A11"),
 			selectInput("dilf", "Select Dilution Factor", choices = c("None", "1:5", "1:10", "1:20"), selected = "None"),
#			conditionalPanel(
#				condition = "input.in_type == 'Y'",
#				textInput("last_well", "Last well in use:", value = "A11"),
#				conditionalPanel(
#					condition = "input.in_type == 'N'",
#					selectInput("dilf", "Select Dilution Factor", choices = c("None", "1:5", "1:10", "1:20"), selected = "None"),
#				)
#			),

			hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
			helpText("To report bugs or request functions to be added, please contact Miles Thorburn <miles@resurrect.bio>")
		),
		
		mainPanel(
			tabsetPanel(
				tabPanel("Target ODs",
					hr(),
					h3("Target OD table:"),
					radioButtons("in_od_table_source", h4("Source of target OD table:"),
						  choices = list("Homogeneous table" = 2, "Editable table" = 1, "Upload CSV" = 0), selected = 1, inline = TRUE),
					radioButtons("default_od", h4("Default target OD for homogeneous table:"),
						  choices = list("0.05" = "0.05", "0.1" = "0.1", "0.2" = "0.2", "0.3" = "0.3"), selected = "0.1", inline = TRUE),
					hr(),
					h5("Please adjust values to reflect target OD FeliX will aim to achieve."),
					h5("Note: Ignore values in unused wells."),
					#DTOutput("homogeneous_od_table"),
					#DTOutput("input_csv_table"),
					DTOutput("all_od_tabs"),
					DTOutput("reactive_od_table")
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

				tabPanel("Single Effector Output",
					hr(),
					h3("Summary Output Table:"),
					#actionButton("calc_go",label = "Calculate"),
					#hr(),
					downloadButton("downloadsumData", "Download FeliX Input Table"),
					textOutput('warns'),
					tags$head(
						tags$style("#warns{color: red;
							font-size: 20px;
							font-style: bold;
						}")
					),
					#uiOutput("warns"),
					hr(),
					DTOutput("out_tab_summ"),
					hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
					h3("Full Output Table:"),
					downloadButton("downloadfullData", "Download Full Output Table"),
					hr(),
					DTOutput("out_tab_full")
				),

				## New panel for infection 
				tabPanel("Infection Assay Output",
					hr(),
					h3("Summary Output Table:"),
					#actionButton("calc_go",label = "Calculate"),
					#hr(),
					downloadButton("downloadsumData", "Download FeliX Input Table"),
					textOutput('ai_floor_val'),
					textOutput('stock_floor_val'),
					textOutput('warns'),
					tags$head(
						tags$style("#warns{color: red;
							font-size: 20px;
							font-style: bold;
						}")
					),
					#uiOutput("warns"),
					hr(),
					DTOutput("out_tab_summInf"),
					hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
					h3("Full Output Table:"),
					downloadButton("downloadfullInfData", "Download Full Output Table"),
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
#)
#shinyApp(ui = ui, server = server)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define server logic
#shinyServer(function(input, output) {
server <- function(input, output) {
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Redundancy around adjustable input values
	## ~~ Setting up callable reactive variables, might be overkill. 
	run_name <- reactive({
		run_name <- input$run_name %>% trimws
	})
	blank_col <- reactive({
		blank_col <- input$blank_col %>% trimws
	})
	target_ul <- reactive({
		target_ul <- input$target_ul %>% trimws
	})
	drift <- reactive({
		drift <- input$drift %>% trimws
	})
	last_well <- reactive({
		last_well <- input$last_well %>% trimws
	})
	in_od_table_val <- reactive({
		in_od_table_val <- input$in_od_table_source %>% as.numeric
	})
	dod <- reactive({
		dod <- input$default_od %>% as.numeric
	})
	min_vol <- reactive({
		min_vol <- input$min_vol %>% as.numeric
	})
	max_vol <- reactive({
		max_vol <- input$max_vol %>% as.numeric
	})
	dilf <- reactive({
		dilf <- input$dilf %>% as.character
	})
	in_type <- reactive({
		in_type <- input$in_type %>% as.character
	})

	## ~~ Reading data in platereader data
	xlsx_in <- reactive({
		req(input$xslx_file)
		df <- read.xlsx(input$xslx_file$datapath,
					startRow = 8, 
					colNames = FALSE) %>% as.data.table
		names(df) <- c("Row", paste0("V", 1:12))
		df
	})

	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Calculations
	## ~~ Main OD calculation function
	calc_dat <- reactive({		
		## Removing redundant data
		dat <- xlsx_in()[!c(9,10),!1] %>% as.data.frame
		
		## Removing blank rows and getting original concentration - how is this working? Surely removing column 1 means the black_col number is wrong... 
		bc <- blank_col() %>% as.numeric
		b_mean <- dat[,bc] %>% mean
		dat[,bc] <- NULL
		no_bs <- dat - b_mean

		## Adjust for drift
		temp_drift <- drift() %>% as.numeric
		w_drift <- no_bs * temp_drift
		w_drift <- w_drift %>% as.data.table

		## Adjusting for dilution factor
		temp_dilf <- dilf() %>% as.character
		if(temp_dilf == "None"){
			w_dilf <- w_drift
		} else if(temp_dilf == "1:5"){
			w_dilf <- w_drift * 5
		} else if(temp_dilf == "1:10"){
			w_dilf <- w_drift * 10
		} else if(temp_dilf == "1:20"){
			w_dilf <- w_drift * 20
		}

		## Manipulating data 
		w_dilf[,"row_name" := letters[1:8] %>% toupper()]
		dat_melted <- melt(w_dilf, id.vars=c("row_name")) 
		names(dat_melted) <- c("row_name", "col_name", "diluted_plate_od_dc")
		dat_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]

		#choices = list("Homogeneous table" = 2, 
 		#"Editable table" = 1, "Upload CSV" = 0)
		iotv <- in_od_table_val() %>% as.numeric
		if(iotv == 0){
			ods <- read_input_ods()
		} else if(iotv == 1){
			ods <- isolate(target_od_dat$data) %>% as.data.table
		} else if(iotv == 2){
			ods <- set_default_tab()
		}
		#ods <- ods %>% as.data.table
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
		temp_target_vol <- target_ul() %>% as.numeric
		calc_stock <- function(row_num){
			temp_val <- (dat$target_od[row_num] / dat$diluted_plate_od_dc[row_num]) * temp_target_vol
			round(temp_val, digits = 1) %>% return
		}
		dat[,"stock_vol" := calc_stock(1:nrow(dat))]
		dat[,"ai_vol" := temp_target_vol - stock_vol %>% round(., digits = 1)]
		dat
		#ods
	})

	## ~~ Calculating final OD from platereader -- A lot repeated in above function. Could separate when I have 30 minutes.
	calc_fin_od <- reactive({
		dat <- xlsx_in()[!c(9,10),!1] %>% as.data.frame
		bc <- blank_col() %>% as.numeric
		b_mean <- dat[,bc] %>% mean
		dat[,bc] <- NULL
		no_bs <- dat - b_mean

		## Adjust for drift
		temp_drift <- drift() %>% as.numeric
		w_drift <- no_bs * temp_drift
		w_drift <- w_drift %>% as.data.table
		w_drift[,"Row" := letters[1:8] %>% toupper()]
		w_drift[,c(12,1:11)]
	})

	## ~~ Creating FeliX input table
	adjust_dat <- reactive({
		#req(input$calc_go)
		dat <- calc_dat() %>% as.data.table
		lw <- last_well() %>% as.character
		it <- in_type() %>% as.character
		## Implemented the last well filter
		#print(lw)
		#lw <- "B1"
		lw_char <- paste0(gsub(pattern = "([A-z]?)([0-9]+)([A-z]?)", replacement = "\\2", lw),
			gsub(pattern = "([A-z]?)([0-9]+)([A-z]?)", replacement = "\\1\\3", lw)
		)
		lw_num <- which(dat$stock_well == lw_char)
		#print(lw_char); print(lw_num)
		if(it == "Yes"){
			dat[c(21:24,9:16,1:8),c("stock_well", "stock_vol", "ai_vol")]
		} else if(it == "No"){
			dat[1:lw_num,c("stock_well", "stock_vol", "ai_vol")]
		}
	})

	## Processing data for infection assays
	process_infdat <- reactive({
		it <- in_type() %>% as.character
		if(it == "Yes"){
			dat_inf <- adjust_dat()
			ai_floor <- floor((dat_inf$ai_vol/1000)) %>% min
			stock_floor <- floor((dat_inf$stock_vol/1000)) %>% min
			if(ai_floor < 0){
				ai_floor <- 0
			}
			if(stock_floor < 0){
				stock_floor <- 0
			}
			dat_inf[,"res_well" := c("A12", "A11", "A10", "A9", "A8", "A7", "A6", "A5", "A4", "A3", "A2", "A1", "H1", "G1", "F1", "E1", "D1", "C1", "B1", "A1")]
			dat_inf[,"original_ai_vol" := ai_vol]
			dat_inf[,"original_stock_vol" := stock_vol]
			dat_inf[,"ai_vol" := ai_vol - (ai_floor * 1000)]
			dat_inf[,"stock_vol" := stock_vol - (stock_floor * 1000)]
			dat_inf
		}
	}) 

	## ~~ Making a text output for the floor volumes to add to the FeliX parameters
	ai_floor_proc <- reactive({
		it <- in_type() %>% as.character
		if(it == "Yes"){
			dat_inf <- adjust_dat()
			ai_floor <- floor((dat$ai_vol/1000)) %>% min
			if(ai_floor < 0){
				ai_floor <- 0
			}
			ai_floor
		}
	})
	stock_floor_proc <- reactive({
		it <- in_type() %>% as.character
		if(it == "Yes"){
			dat_inf <- adjust_dat()
			stock_floor <- floor((dat$stock_vol/1000)) %>% min
			if(stock_floor < 0){
				stock_floor <- 0
			}
			stock_floor
		}
	})
	output$ai_floor_val <- renderText({
		ai_floor_out <- ai_floor_proc() %>% as.character
		paste0("AI Floor Volume (ml): ", ai_floor_out)
	})
	output$stock_floor_val <- renderText({
		stock_floor_out <- stock_floor_proc() %>% as.character
		paste0("Stock Floor Volume (ml): ", stock_floor_out)
	})

	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fixing erroneous values based on input limit
	## Making new datatable
	fix_anoms <- reactive({
		it <- in_type() %>% as.character
		if(it == "Yes"){
			dat <- process_infdat() %>% as.data.table
		} else if(it == "No"){
			dat <- adjust_dat() %>% as.data.table
		}
		min_vol_num <- min_vol() %>% as.numeric
		max_vol_num <- max_vol() %>% as.numeric
		adjust_vols <- function(row_num){
			temp_row <- dat[row_num,]
			## Adjusting if values higher than limits
			temp_row[,"stock_vol" := ifelse((stock_vol > min_vol_num), 20, stock_vol)]
			temp_row[,"ai_vol" := ifelse((ai_vol > max_vol_num), 20, ai_vol)]
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
		if(it == "Yes"){
			out[,c("stock_well", "stock_vol", "ai_vol", "res_well")]
		} else if(it == "No"){
			out
		}
	})

	## Identifying problems and rendering output
	anoms_ids <- reactive({
		fin_dat <- fix_anoms() %>% as.data.table
		errs <- subset(fin_dat, stock_vol == 20 & ai_vol == 20)
		errs
	})
	output$warns <- renderText({
		errs_out <- anoms_ids() %>% as.data.table
		paste0("WARN! Volumes outside of acceptable range for wells: ", paste(errs_out$stock_well, collapse = " "))
	})
	#output$warns_html <- renderUI({
	#	errs_out <- anoms_ids() %>% as.data.table
	#	HTML(as.character(div(style="color: red;", 
	#		paste0("WARN! Volumes outside of acceptable range for wells: ", paste(errs_out$stock_well, collapse = " ")))
	#})

	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Printing tables and download handling
	## ~~ Rendering data tables for inspection
	output$pr_tab <- renderDT({
		datatable(xlsx_in(), rownames= FALSE)
	})
	output$fin_od_tab <- renderDT({
		datatable(calc_fin_od(), rownames= FALSE) %>% formatRound(c(2:12), 3)
	})
	output$out_tab_full <- renderDT({
		datatable(calc_dat(), rownames= FALSE) %>% formatRound(c(4), 3) %>% formatRound(c(6,7), 3)
	})
	output$out_tab_summ <- renderDT({
		#datatable(adjust_dat(), rownames= FALSE)
		datatable(fix_anoms(), rownames= FALSE, options = list(
			paging =TRUE,
			pageLength = 8
		)) %>% formatRound(c(2:3), 1)
	})
	output$out_tab_summInf <- renderDT({
		#datatable(adjust_dat(), rownames= FALSE)
		datatable(fix_anoms(), rownames= FALSE, options = list(
			paging =TRUE,
			pageLength = 8
		)) %>% formatRound(c(2:3), 1)
	})


	## ~~ Handing download for datasets
	output$downloadsumData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_", format(Sys.Date(), "%d%m%Y"), ".txt")
		},
		content = function(file) {
			#write.csv(adjust_dat(), file, row.names = FALSE, quote = FALSE)
			write.csv(fix_anoms(), file, row.names = FALSE, quote = FALSE)
		}
	)
	output$downloadfullData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_full_table_", format(Sys.Date(), "%d%m%Y"), ".csv")
		},
		content = function(file) {
			write.csv(calc_dat(), file, row.names = FALSE, quote = FALSE)
		}
	)
	output$downloadfullInfData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_full_table_", format(Sys.Date(), "%d%m%Y"), ".csv")
		},
		content = function(file) {
			write.csv(process_infdat(), file, row.names = FALSE, quote = FALSE)
		}
	)
	## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Target OD tables
	## ~~ Define default depending on input source data
	#od_table_reactive <- 
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

	## ~~ Outputting the target OD tables from different sources
 	#choices = list("Homogeneous table" = 2, 
 	#"Editable table" = 1, "Upload CSV" = 0), selected = 1, inline = TRUE),
	#output$input_csv_table <- renderDT({
	#	req(input$in_od_table_source == 0)
	#	read_input_ods() %>% datatable(., rownames = FALSE)
	#})
	output$reactive_od_table <- renderDT({
		req(input$in_od_table_source == 1)
		datatable(target_od_dat$data, editable = TRUE, rownames= FALSE)
	})
	#output$homogeneous_od_table <- renderDT({
	#	req(input$in_od_table_source == 2)
	#	set_default_tab() %>% datatable(., rownames = FALSE)
	#})
	output$all_od_tabs <- renderDT({
		if(input$in_od_table_source == 0){
			read_input_ods() %>% datatable(., rownames = FALSE)
		} else if(input$in_od_table_source == 2){
			set_default_tab() %>% datatable(., rownames = FALSE)
		}
	})

	## ~~ Reactive OD table elements
	target_od_dat <- reactiveValues(data = {
		dat <- data.table( Row = letters[1:8] %>% toupper(),
			V1 = rep(0.1, 8), V2 = rep(0.1, 8),
			V3 = rep(0.1, 8), V4 = rep(0.1, 8),
			V5 = rep(0.1, 8), V6 = rep(0.1, 8),
			V7 = rep(0.1, 8), V8 = rep(0.1, 8),
			V9 = rep(0.1, 8), V10 = rep(0.1, 8),
			V11 = rep(0.1, 8), V12 = rep(0.1, 8))
		dat
	})
	## Note, if the input table name is n, the input here is n_cell_edit
	observeEvent(input$reactive_od_table_cell_edit, {
		## Getting values
		info = input$reactive_od_table_cell_edit
		i = as.numeric(info$row)
		j = as.numeric(info$col)
		k = as.numeric(info$value)

		## Writing values to reactive object
		target_od_dat$data[i,j] <- k
	})

	## ~~ Debugging lines. Was used for reactive table inputs. 
	debug_table <- reactive({
		iotv <- in_od_table_val() %>% as.numeric
		if(iotv == 0){
			ods <- read_input_ods() %>% as.data.table
		} else if(iotv == 1){
			ods <- isolate(target_od_dat$data) %>% as.data.table
			ods[,"Row"] <- NULL
		} else if(iotv == 2){
			ods <- set_default_tab() %>% as.data.table
			ods[,"Row"] <- NULL
		}
		ods[,"row_name" := letters[1:8] %>% toupper()]
		ods_melted <- melt(ods, id.vars=c("row_name")) 
		names(ods_melted) <- c("row_name", "col_name", "target_od")
		ods_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]
		ods
	})
	output$debug_out_tab <- renderDT({ 
		debug_table()
  	})
	#output$debug_out_print <- renderText({ 
	#	input$reactive_od_table_cell_edit %>% str %>% print
	#})
}
#)
shinyApp(ui = ui, server = server)
#runApp()