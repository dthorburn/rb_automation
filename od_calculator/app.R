#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~				  Resurrect Bio AgroPrep OD Calculator for FeliX automation						  ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <dthorburn@imperial.ac.uk>
## Last Update: 23/01/24

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

## Just for setting the working dir when adjusting things. 
setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/od_calculator")

## Define functions


## Define default data
od_table <- data.table(	Row = letters[1:8] %>% toupper(),
						V1 = rep(0.1, 8), V2 = rep(0.1, 8),
						V3 = rep(0.1, 8), V4 = rep(0.1, 8),
						V5 = rep(0.1, 8), V6 = rep(0.1, 8),
						V7 = rep(0.1, 8), V8 = rep(0.1, 8),
						V9 = rep(0.1, 8), V10 = rep(0.1, 8),
						V11 = rep(0.1, 8), V12 = rep(0.1, 8))

# Define UI ----
##img(src='rb_logo.png', align = "top", height = 30, width = 100),
ui <- fluidPage(
	titlePanel(strong("RB OD Calculator: A tool to calculate aliquot volumes for the FeliX.")),
    sidebarLayout(
    	sidebarPanel(
	    	h2("Input data:"),
		    # Input: Select a file ----
		    fileInput("xslx_file", "Platereader XSLX File",
        	    multiple = FALSE),

	    	h4("Calculation parameters:"),
		    textInput("run_name", "Name of Run:", value = "", placeholder = "N77_AA13"),
		    textInput("blank_col", "Blank column:", value = "12"),
		    textInput("target_ul", "Target volume (ul):", value = "1000"),
		    textInput("drift", "Platereader drift:", value = "3.95"),
		    textInput("last_well", "Last well in use:", value = "A11"),
 
    		hr(style="border-color: #b2b2b3; margin-bottom: 0px"),
    		helpText("To report bugs or request functions to be added, please contact Miles Thorburn <miles@resurrect.bio>")
    	),
	    
	    mainPanel(
		    tabsetPanel(
	    		tabPanel("Input",
			    	hr(),
		    		h2("Target OD table:"),
		    		h5("Please adjust values to reflect target OD FeliX will aim to achieve."),
		    		h5("Note: Leave unused wells at default 0.1."),
		    		DTOutput("reactive_od_table")
	    		),

	    		tabPanel("Platereader table",
	    			hr(),
		    		h2("Input XLSX table from platereader:"),
		    		h5("An uneditable version of the platereader plate input file. Mostly for debugging purposes."),
		    		h5("If no results are showing, please press Submit again once file path has been selected."),
		    		DTOutput("pr_tab"),
		    		#textOutput("debug_out")
				),

	    		tabPanel("Output",
	    			hr(),
		    		h2("Summary Output Table:"),
		    		downloadButton("downloadsumData", "Download FeliX Input Table"),
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
	    			h2("How to use the OD calculator:"),
	    			h5("1. Select local XLSX file emitted from the platereader in the shared lab."),
	    			h5("2. If necessary, change any of the default calculation parameters."),
	    			h5("3. Adjust the target OD table."),
	    			h5("4. Press Submit button."),
	    			h5("4. Download FeliX input table as a .txt file by pressing Download Data button."),
	    			hr(),
	    			h2("Important Considerations:"),
	    			h5("1. Due to how FeliX operated, wells need to be filled from H -> A for each column."),
	    			h5("2. Wells cannot be skipped. These will have to be manually removed after FeliX is complete."),
				),
	    	),
		)
	)
)
#shinyApp(ui = ui, server = server)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Define server logic
server <- function(input, output) {
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Redundancy around adjustable input values
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

    ## Reading data in platereader data
    xlsx_in <- reactive({
    	req(input$xslx_file)
    	df <- read.xlsx(input$xslx_file$datapath,
					startRow = 8, 
					colNames = FALSE) %>% as.data.table
    	names(df) <- c("Row", paste0("V", 1:12))
    	df
    })

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start of calculations
    calc_dat <- reactive({
    	## Removing redundant data
    	dat <- xlsx_in()[!c(9,10),!1] %>% as.data.frame
    	
    	## Removing blank rows and getting original concentration
    	bc <- blank_col() %>% as.numeric
    	b_mean <- dat[,bc] %>% mean
    	dat[,bc] <- NULL
		no_bs <- dat - b_mean

		## Adjust for drift
		temp_drift <- drift() %>% as.numeric
		w_drift <- no_bs * temp_drift
		w_drift <- w_drift %>% as.data.table

		## Manipulating data 
		w_drift[,"row_name" := letters[1:8] %>% toupper()]
		dat_melted <- melt(w_drift, id.vars=c("row_name")) 
		names(dat_melted) <- c("row_name", "col_name", "diluted_plate_od_dc")
		dat_melted[,"stock_well" := paste0(gsub(pattern = "V", replacement = "", col_name), row_name)]

		ods <- target_od_dat$data %>% as.data.table
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
    })

    adjust_dat <- reactive({
    	dat <- calc_dat() %>% as.data.table
    	dat[,c("stock_well", "stock_vol", "ai_vol")]
    })

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Printing tables and download handling
    ## Just to be able to inspect them manually if results are weird
    output$pr_tab <- renderDT({
    	#datatable(calc_dat(), rownames= FALSE)
    	datatable(xlsx_in(), rownames= FALSE)
    })
    
    output$out_tab_full <- renderDT({
    	datatable(calc_dat(), rownames= FALSE)
    })

	output$out_tab_summ <- renderDT({
    	datatable(adjust_dat(), rownames= FALSE)
    })

	output$downloadsumData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_", format(Sys.Date(), "%d%m%Y"), ".txt")
		},
		content = function(file) {
			write.csv(adjust_dat(), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
		}
	)
	output$downloadfullData <- downloadHandler(
		filename = function() {
			paste0(run_name(), "_full_table_", format(Sys.Date(), "%d%m%Y"), ".csv")
		},
		content = function(file) {
			write.csv(calc_dat(), file, row.names = FALSE, col.names = FALSE, quote = FALSE)
		}
	)
    #output$debug_out <- renderText({ 
    #	calc_dat()
  	#})


    ## Reactive OD table
    ## Needs to be fixed!!!! - can't see why the values are not updating. 
    ## Something about inputs is clearly being lost
	target_od_dat <- reactiveValues(data = od_table)
    output$reactive_od_table <- renderDT({
        datatable(target_od_dat$data, editable = TRUE, rownames= FALSE)
    })
    observeEvent(input$reactive_od_table_edit, {
        #get values
        info = input$reactive_od_table_edit
        i = as.numeric(info$row)
        j = as.numeric(info$col)
        k = as.numeric(info$value)

        #write values to reactive object
        target_od_dat$data[i,j] <- k
    })
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Run the app
shinyApp(ui = ui, server = server)