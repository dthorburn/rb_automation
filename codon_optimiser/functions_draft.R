#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~             Resurrect Bio Codon Optimisation and Processing App             ~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Author: Doko-Miles Thorburn <miles@resurrect.bio>
## 


                                                                    ################################################
                                                                    ##                Libraries                   ##
                                                                    ################################################
## Loading in required libraries
if (!require("DT", quietly = TRUE))
  install.packages("DT")
if (!require("shiny", quietly = TRUE))
  install.packages("shiny")
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("seqinr", quietly = TRUE))
  install.packages("seqinr")
if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!require("data.table", quietly = TRUE))
  install.packages("data.table")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")

suppressMessages(library("data.table"))
suppressMessages(library("Biostrings"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
#suppressMessages(library("seqinr"))
suppressMessages(library("dplyr"))
suppressMessages(library("shiny"))
suppressMessages(library("DT"))

                                                                    ################################################
                                                                    ##                Parameters                  ##
                                                                    ################################################
seed           <- 2021
min_length     <- 300
#input_seq      <- "AQDLARQPNCGQNCLLTAFTTLSNGCTQTDFACLCKSQKFTGTSNLRYSSSCTSTEAATTKSWGAKTCASVGVNVTATLGNTTSVGAGGLVPPLLEQADSDVVNVDEQVEHELAKISDSLQSNGFTSEE"
x5_flank       <- "CACTCTGTGGTCTCAA"
x3_flank       <- "GGTTCGTGAGACCACGAAGTG"
cut_site_avoid <- "GGTCTC" ## Default: BSAI
cut_iterations <- 20
filler_seq     <- "aatcagttcaacgccgtcttacaacaggcagctaagtgccatattaacgagaatgaggtgcttcgggaactggtcaagaaaatcagaactgtagtgaacagcgcggaagacgccatagacaaatttgtaagagaatgaggtgcttcgggaactggtcaagaa"
filler_loc     <- "3_prime" ## Where to add filler sequence, 5' or 3'
remove_X       <- TRUE
add_M          <- TRUE
lower_flanks   <- TRUE

#setwd("C:/Users/miles/Documents/Resurrect_Bio/Scripts/rb_automation/codon_optimiser/test_data")

## Maybe implement
# best_codon vs sampling based on distribution.
# Checklist to click through to ensure data is processed successfully.
#   - Removal of signal peptide

                                                                    ################################################
                                                                    ##               Codon table                  ##
                                                                    ################################################

## Codon usage table taken from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100&aa=1&style=N
#("UAA", "UGA", "UAG", "GCU", "GCC", "GCA", "GCG", "UGU", "UGC", "GAU", "GAC", "GAA", "GAG", "UUU", 
#      "UUC", "GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "AAA", "AAG", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", 
#      "AUG", "AAU", "AAC", "CCU", "CCC", "CCA", "CCG", "CAA", "CAG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "UCU", "UCC", "UCA", 
#      "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "GUU", "GUC", "GUA", "GUG", "UGG", "UAU", "UAC")
codon_table <- structure(
  list(
    codon = c("TAA", "TGA", "TAG", "GCT", "GCC", "GCA", "GCG", "TGT", "TGC", "GAT", "GAC", "GAA", "GAG", "TTT", "TTC", "GGT", "GGC", 
      "GGA", "GGG", "CAT", "CAC", "ATT", "ATC", "ATA", "AAA", "AAG", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATG", "AAT", "AAC", 
      "CCT", "CCC", "CCA", "CCG", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
      "ACT", "ACC", "ACA", "ACG", "GTT", "GTC", "GTA", "GTG", "TGG", "TAT", "TAC"), 
    aa = c("X", "X", "X", "A", "A", "A", "A", "C", "C", "D", "D", "E", "E", "F", "F", "G", "G", "G", "G", "H", "H", "I", "I", "I", "K", "K", 
      "L", "L", "L", "L", "L", "L", "M", "N", "N", "P", "P", "P", "P", "Q", "Q", "R", "R", "R", "R", "R", "R", "S", "S", "S", "S", "S", "S", 
      "T", "T", "T", "T", "V", "V", "V", "V", "W", "Y", "Y"), 
    frequency = c(0.32, 0.39, 0.29, 0.44, 0.17, 0.31, 0.09, 0.56, 0.44, 0.7, 0.3, 0.53, 0.47, 0.57, 0.43, 0.35, 0.16, 
      0.32, 0.16, 0.61, 0.39, 0.51, 0.26, 0.23, 0.43, 0.57, 0.13, 0.25, 0.26, 0.13, 0.1, 0.12, 1, 0.63, 0.37, 0.39, 0.13, 0.35, 0.13, 
      0.51, 0.49, 0.15, 0.08, 0.11, 0.1, 0.3, 0.25, 0.28, 0.13, 0.21, 0.07, 0.18, 0.13, 0.36, 0.21, 0.32, 0.11, 0.42, 0.17, 0.16, 0.25, 
      1, 0.55, 0.45)), 
    row.names = c(NA, -67L)) %>% as.data.table()
#    class = c("data.table", "data.frame"))

                                                                    ################################################
                                                                    ##                 Functions                  ##
                                                                    ################################################

## Convert stop codons from * to X - Just for ease with donwstream work.
convert_stops <- function(data_row){
  new_input <- data.table(seqname = data_row$seqname, 
                          og_aa_seq = gsub(pattern = "\\*", replacement = "X", data_row$seqs))
  return(new_input)
}



## Function to optimize codon usage based on frequency
optimise_codon_seq <- function(data_row, rand_seed, remove_X, add_M, codon_table){
  optimised_seq <- ""
  set.seed(rand_seed)
  tokenised_seq <- strsplit(data_row$og_aa_seq, split = "")[[1]]
  if(remove_X == TRUE & tokenised_seq[length(tokenised_seq)] == "X"){
    tokenised_seq <- tokenised_seq[1:(length(tokenised_seq) - 1)]
  }
  if(add_M == TRUE & tokenised_seq[1] != "M"){
    tokenised_seq <- append("M", tokenised_seq)
  }
  for(temp_aa in tokenised_seq){
    aa_codons <- subset(codon_table, aa == temp_aa)
    codon_out <- sample(aa_codons$codon, size = 1, replace = FALSE, prob = aa_codons$frequency)
    #best_codon <- aa_codons[which.max(aa_codons$Frequency), "Codon"]
    optimised_seq <- paste0(optimised_seq, codon_out)
  }
  output <- data.table(seqname = data_row$seqname, opt_dna_seqs = optimised_seq, og_aa_seq = data_row$og_aa_seq) 
  return(output)
}

## Function that checks if site to avoid is present in intial optimised sequence. If so, it increments seed by 1 and tries again. 
## If after `cut_iterations` no solution without restriction site is found, it fails, but still emits the sequence. 
#data_row <- seq2[2,]
restriction_site_check <- function(data_row, seed, cut_iterations, cut_site_avoid, remove_X, add_M){
  contains_site <- grepl(pattern = cut_site_avoid, data_row$opt_dna_seqs, ignore.case = TRUE)
  if(contains_site == FALSE){
    return(data_row)
  } else if(contains_site == TRUE){
    counter  <- 0
    while(contains_site == TRUE){
      counter   <- counter + 1
      temp_seed <- seed + 1
      temp_seq  <- optimise_codon_seq(data_row[,c(1,3)], temp_seed, remove_X, add_M, codon_table)
      contains_site <- grepl(pattern = cut_site_avoid, temp_seq$opt_dna_seqs)
      if(counter == cut_iterations){
        break
        warn_messages <- paste0("ERROR: Completed ", counter , " iterations of random codon optimisation in ", data_row$seqname, ". Cannot remove cut site ", cut_site_avoid, "\nIncrease iterations or change seed and try again. If problem persists, manually inspect.")
        ## If this becomes an issue, add debug details with coordinates and output all iterations to log file.  
      }
    }
    return(temp_seq)
  }
}

## Function to add flanks and add filler if needed. Filler will be lowercase just to easily identify. 
data_row <- seq3[9,]
add_flanks <- function(data_row, x5_flank, x3_flank, filler_seq, min_length, filler_loc, lower_flanks){
  if(lower_flanks == TRUE){
    x5_flank <- tolower(x5_flank)
    x3_flank <- tolower(x3_flank)
  } else if(lower_flanks == FALSE){
    x5_flank <- toupper(x5_flank)
    x3_flank <- toupper(x3_flank)    
  }
  with_flanks <- paste0(x5_flank, data_row$opt_dna_seqs, x3_flank)

  ## Second length check and filler addition if required.
  len_diff <- nchar(with_flanks) - as.numeric(min_length)
  if(len_diff < 0){
    to_fill <- substring(filler_seq, first = 1, last = abs(len_diff))
    if(filler_loc == "5_prime"){
      out_seq <- paste0(to_fill, with_flanks)
    } else if(filler_loc == "3_prime"){
      out_seq <- paste0(with_flanks, to_fill)
    }
  } else {
    out_seq <- with_flanks
  }
  output <- data.table(seqname = data_row$seqname, opt_dna_seqs_wflanks = out_seq, 
                      opt_dna_seqs = data_row$opt_dna_seqs, og_aa_seq = data_row$og_aa_seq)
  return(output)
}
data_row <- seq4[5,]
#output_check <- function(data_row, seq_w_flanks_dt, original_seq_dt, cut_site_avoid){
sanity_check <- function(data_row, cut_site_avoid, min_length, x5_flank, x3_flank){
  #warn_messages <- ""

  ## Length
  min_length_wflanks <- min_length - nchar(x5_flank) - nchar(x3_flank)
  leftover_length <- (nchar(data_row$og_aa_seq) * 3) - min_length_wflanks
  data_row[,"DNA_length" := nchar(data_row$opt_dna_seqs_wflanks)]
  if(leftover_length < 0){
    data_row[,"DNA_filler" := abs(leftover_length)]
    data_row[,"filler_warn" := 1]
    #warn_messages <- paste0(warn_messages, "WARN: ", data_row$seqname, " is ", abs(leftover_length), " short of minimum length. Making up length with filler.\n")
  } else {
    data_row[,"DNA_filler" := 0]
    data_row[,"filler_warn" := 0]
  }

  ## Starting with M 
  starting_aa <- strsplit(data_row$og_aa_seq, split = "")[[1]][1]
  if(starting_aa != "M"){
    data_row[,"start_m_warn" := 1]
    #warn_messages <- paste0(warn_messages, "WARN: ", data_row$seqname, " starts with ", starting_aa, ", not M. Please ensure this is correct.\n")
  } else {
    data_row[,"start_m_warn" := 0]
  }

  ## Overlap with second check. Moved relevant aspects. 
  #tokenised_seq <- strsplit(data_row$og_aa_seq, split = "")[[1]]
  #if(grepl(pattern = "X", data_row$og_aa_seq)){
  #  stops <- which(tokenised_seq == "X")
  #  data_row[,"internal_x_warn" := 1]
  #  #warn_messages <- paste0(warn_messages, "WARN: Stop codon identified in ", data_row$seqname, " position(s) ", paste(stops, collapse = " + "), " of ", nchar(data_row$seqs), ". Please ensure this is correct.\n")
  #} else {
  #  data_row[,"internal_x_warn" := 0]
  #}

  ## Second rescriction site check
  #seq_w_flanks <- seq_w_flanks_dt[which(data_row$seqname == seq_w_flanks_dt$seqname), 2]
  num_restrictions <- str_count(pattern = toupper(cut_site_avoid), toupper(data_row$opt_dna_seqs_wflanks))
  data_row[,"num_restriction_sites" := num_restrictions]
  if(num_restrictions > 1){
    data_row[,"restriction_sites_warn" := 1]
    #warn_messages <- paste0(warn_messages, "ERROR: ", num_restrictions, " sites to avoid detected with sequence ", cut_site_avoid, ". There should be 1.\n")
  } else if(num_restrictions == 0){
    data_row[,"restriction_sites_warn" := 1]
    #warn_messages <- paste0(warn_messages, "ERROR: ", num_restrictions, " sites to avoid detected with sequence ", cut_site_avoid, ". There should be 1.\n")
  } else {
    data_row[,"restriction_sites_warn" := 0]
  }

  ## Internal stop codon check
  translated <- DNAString(data_row$opt_dna_seqs) %>% Biostrings::translate(.) %>% as.character %>% gsub(pattern = "\\*", replacement = "X", .)
  tokenised_seq <- strsplit(translated, split = "")[[1]]
  num_int_stops <- str_count(pattern = "X", translated)
  if(num_int_stops > 0){
    data_row[,"internal_x_warn" := 1]
    stops <- which(tokenised_seq == "X")
    data_row[,"internal_x_positions" := paste(stops, collapse = ", ")]
    #warn_messages <- paste0(warn_messages, "ERROR: ", num_int_stops, " stop codons detected in ", data_row$seqname,".\n")
  } else {
    data_row[,"internal_x_warn" := 0]
    data_row[,"internal_x_positions" := 0]
  }

  ## Sanity check against input
  #og_seq <- original_seq_dt[which(data_row$seqname == original_seq_dt$seqname), 2]
  match_input <- grepl(pattern = translated, data_row$og_aa_seq, ignore.case = TRUE)
  if(match_input == FALSE){
    ## Removing start and end to see if pattern can be found since these can be added/removed during optimisation. If not, a problem might be present.
    match_attempt2 <- grepl(pattern = str_sub(translated, start = 2, end = -2), data_row$og_aa_seq)
    if(match_attempt2 == FALSE){
      data_row[,"sequence_match_warn" := 1]
      #warn_messages <- paste0(warn_messages, "ERROR: Input ", data_row$seqname," sequence and translated optimised sequence do not match.\n")
    } else {
      data_row[,"sequence_match_warn" := 0]
    }
  } else {
    data_row[,"sequence_match_warn" := 0]
  }

  return(data_row)
}


## Testing funcitons - works well enough.
# Main functionality
seq1 <- convert_stops(input_seq) 
#check_input_seq(seq1, min_length, x5_flank, x3_flank)
seq2 <- optimise_codon_seq(seq1, seed, remove_X, add_M)
seq3 <- restriction_site_check(seq2, seed, cut_iterations, cut_site_avoid)
seq4 <- add_flanks(seq3, x5_flank, x3_flank, filler_seq, min_length, filler_loc, lower_flanks)
output_check(seq3, seq4, seq1, cut_site_avoid)



## Updating to batch process fastas
fas <- seqinr::read.fasta(file = "effectors.fasta", as.string = TRUE, forceDNAtolower = FALSE)
fas_dt <- data.table(seqname = names(fas), seqs = fas %>% unlist %>% as.data.table)
names(fas_dt)[2] <- "seqs" ## This is odd behaviour. 

seq1 <- convert_stops(fas_dt) 
#check_input_seq(seq1[1,], min_length, x5_flank, x3_flank)
#optimise_codon_seq(seq1[2,], seed, remove_X, add_M, codon_table)[,2] %>% as.character %>% DNAString(.) %>% Biostrings::translate(.)

for(i in 1:nrow(seq1)){
  temp_out <- optimise_codon_seq(seq1[eval(i),], seed, remove_X, add_M, codon_table)
  if(i == 1){
    seq2 <- temp_out
  } else {
    seq2 <- rbind(seq2, temp_out)
  }
}

restriction_site_check(seq2[2,], seed, cut_iterations, cut_site_avoid, remove_X, add_M)
for(i in 1:nrow(seq2)){
  temp_out <- restriction_site_check(seq2[eval(i),], seed, cut_iterations, cut_site_avoid, remove_X, add_M)
  if(i == 1){
    seq3 <- temp_out
  } else {
    seq3 <- rbind(seq3, temp_out)
  }
}

for(i in 1:nrow(seq3)){
  temp_out <- add_flanks(seq3[eval(i),], x5_flank, x3_flank, filler_seq, min_length, filler_loc, lower_flanks)
  if(i == 1){
    seq4 <- temp_out
  } else {
    seq4 <- rbind(seq4, temp_out)
  }
}


for(i in 1:nrow(seq4)){
  temp_out <- sanity_check(seq4[eval(i),], cut_site_avoid, min_length, x5_flank, x3_flank)
  if(i == 1){
    seq5 <- temp_out
  } else {
    seq5 <- rbind(seq5, temp_out)
  }
}

## Melting relevant columns of seq5
to_plot <- melt(seq5[,c(1, 7, 8, 10, 11, 13)], id.vars = c("seqname"))

ggplot(to_plot, aes(y = seqname, x = variable, fill = value)) +
  geom_tile(color = "black") +
  labs(x = "Warnings", y = "Input Sequence Name") +
  scale_fill_gradient(low = "#09b81b", high = "#cf0a0a") +
  theme_classic(base_size = 20) +
  scale_x_discrete(labels = c("Length Filler", "Input 5' M", "Sites to Avoid", "Internal Stops", "Input :: Output\nSequence Match"), expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "None")




## Testing sample distributions
aa_codon <- subset(codon_table, aa == "L")
for(i in seq(1,1000,1)){
  message(i)
  temp_row <- data.table(iteration = i, codon = sample(aa_codons$codon, size = 1, replace = FALSE, prob = aa_codons$frequency))
  if(i == 1){
    output <- temp_row
  } else {
    output <- rbind(output, temp_row)
  }
} 
output$codon %>% table/1000 %>% sort
aa_codons

