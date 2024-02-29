if (!require("data.table", quietly = TRUE))
	install.packages("data.table")
if (!require("seqinr", quietly = TRUE))
	install.packages("seqinr")
if (!require("dplyr", quietly = TRUE))
	install.packages("dplyr")
suppressMessages(library(data.table))
suppressMessages(library(seqinr))
suppressMessages(library(dplyr))

get_seq <- function(seqname){
    data.table(seqname = seqname, 
               sequence = fas[[seqname]][1]) %>% return
}

## args order should be: PATHOGEN pipeline output scored csv (no seqs), PATHOGEN output fasta files, include short seqs (boolean), output path
args <- commandArgs(trailingOnly = TRUE)
#print(args)
path_in 			<- args[1] %>% as.character
fasta_in 			<- args[2] %>% as.character
keep_short_bool 	<- args[3] %>% as.character
path_out			<- args[4] %>% as.character

#path_in <- "~/Resurrect_Bio/Projects/01_Candidate_Search/02_Effectors/11_fgraminearum/01_PATHOGEN/PATHOGEN-OUT/Fgraminearum_PH1v3_peptides_processed_effs_scored.csv"
#fasta_in <- "~/Resurrect_Bio/Projects/01_Candidate_Search/02_Effectors/11_fgraminearum/01_PATHOGEN/SEQUENCE-OUT/Fgraminearum_PH1v3_peptides.passed.fasta"
#keep_short_bool <- "TRUE"
#path_out <- "~/Resurrect_Bio/Projects/01_Candidate_Search/02_Effectors/11_fgraminearum/01_PATHOGEN/SEQUENCE-OUT/Fgraminearum_PH1v3_peptides.scored.NoSP.fasta"

dat <- fread(path_in)
fas <- read.fasta(file = fasta_in, as.string = TRUE, forceDNAtolower = FALSE)

dat <- subset(dat, Score != "Failed_TP" & Score != "Failed_TM" & Score != "Failed_SP")
if(!grepl(keep_short_bool, pattern = "true", ignore.case = TRUE)){
	dat <- subset(dat, Score != "Failed_Short")
}
temp_seqs <- lapply(FUN = get_seq, dat$seqname %>% unique) %>% rbindlist
dat_out <- merge(dat, temp_seqs, by = "seqname")
dat_out[,"trunc_seq" := paste0("M", substr(sequence, start = (as.numeric(SP_cleavage_site)+1), stop = as.numeric(AA_length)))]
dat_out[,"trunc_length" := nchar(trunc_seq)]

## Using a loop over a vector to ensure writing is in correct order. 
for(i in 1:nrow(dat)){
	temp_row <- dat_out[i,]
	temp_name <- paste0(">", temp_row$seqname, "_NoSP_", temp_row$OrthoGroup, "_", temp_row$trunc_length, "\n")
	temp_seq  <- paste0(temp_row$trunc_seq, "\n")
	if(i == 1){
		cat(file = path_out, temp_name, append = FALSE)
	} else {
		cat(file = path_out, temp_name, append = TRUE)
	}
	cat(file = path_out, temp_seq, append = TRUE)
}

message(paste0("Emitted NoSP sequences: ", i,"\nOutput files: ", path_out))
