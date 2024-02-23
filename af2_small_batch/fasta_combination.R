if (!require("seqinr", quietly = TRUE))
	install.packages("seqinr")
if (!require("dplyr", quietly = TRUE))
	install.packages("dplyr")
suppressMessages(library(seqinr))
suppressMessages(library(dplyr))

## args order should be: nlr fasta, effector fasta, AA length cutoff, and output file name
args <- commandArgs(trailingOnly = TRUE)
#print(args)
nlr_fasta 		<- args[1] %>% as.character
eff_fasta 		<- args[2] %>% as.character
length_cutoff 	<- args[3] %>% as.numeric
output_name 	<- args[4] %>% as.character

## reading data
effs <- read.fasta(file = nlr_fasta, as.string = TRUE, forceDNAtolower = FALSE)
nlrs <- read.fasta(file = eff_fasta, as.string = TRUE, forceDNAtolower = FALSE)

## output name handling
if(grepl(pattern = "(\\.fa$)|(\\.faa$)|(\\.fasta$)", output_name)){
		long_out_name <- gsub(
			pattern = "(^.*)(\\.fa(a)?(sta)?)", 
			replacement = "\\1_too_long\\2", 
			x = output_name)
} else {
	long_out_name <- paste0(output_name, "too_long.fasta")
	out_name <- paste0(output_name, ".fasta")
}

## Loop to mix sequences - could vectorise, but seems unecessary.
cor_length <- 0
incor_length <- 0
for(temp_nlr in names(nlrs)){
	nlr_seq <- nlrs[[temp_nlr]][1]
	for(temp_eff in names(effs)){
		eff_seq <- effs[[temp_eff]][1]
		tot_len <- sum(nchar(nlr_seq), nchar(eff_seq))
		temp_name <- paste0(">", temp_nlr, "_", temp_eff, "_", tot_len,  "\n")
		temp_seq <- paste0(nlr_seq, ":", eff_seq, "\n")

		## Length checkpoint
		if(tot_len >= length_cutoff){
			incor_length <- incor_length + 1
			if(incor_length == 1){
				cat(file = paste0("./", long_out_name), temp_name, append = FALSE)
			} else {
				cat(file = paste0("./", long_out_name), temp_name, append = TRUE)
			}
			cat(file = paste0("./", long_out_name), temp_seq, append = TRUE)
		} else if(tot_len < length_cutoff){
			cor_length <- cor_length + 1
			if(cor_length == 1){
				cat(file = paste0("./", out_name), temp_name, append = FALSE)
			} else {
				cat(file = paste0("./", out_name), temp_name, append = TRUE)
			}
			cat(file = paste0("./", out_name), temp_seq, append = TRUE)
		}
	}
}
wd <- getwd()
message(paste0("Correct length sequences: ", cor_length,"\nSequences too long: ", incor_length, "\nFiles located: ", wd))