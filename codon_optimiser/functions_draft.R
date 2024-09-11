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
if (!require("stringr", quietly = TRUE))
  install.packages("stringr")
if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!require("data.table", quietly = TRUE))
  install.packages("data.table")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("Biostrings", quietly = TRUE))
    install.packages("Biostrings")

suppressMessages(library("data.table"))
suppressMessages(library("Biostrings"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("shiny"))
suppressMessages(library("DT"))

                                                                    ################################################
                                                                    ##                Parameters                  ##
                                                                    ################################################
seed           <- 2021
min_length     <- 300
input_seq      <- "AQDLARQPNCGQNCLLTAFTTLSNGCTQTDFACLCKSQKFTGTSNLRYSSSCTSTEAATTKSWGAKTCASVGVNVTATLGNTTSVGAGGLVPPLLEQADSDVVNVDEQVEHELAKISDSLQSNGFTSEE"
x5_flank       <- "CACTCTGTGGTCTCAA"
x3_flank       <- "GGTTCGTGAGACCACGAAGTG"
cut_site_avoid <- "GGTCTC" ## Default: BSAI
cut_iterations <- 20
filler_seq     <- "aatcagttcaacgccgtcttacaacaggcagctaagtgccatattaacgagaatgaggtgcttcgggaactggtcaagaaaatcagaactgtagtgaacagcgcggaagacgccatagacaaatttgtaagagaatgaggtgcttcgggaactggtcaagaa"
filler_loc     <- "3_prime" ## Where to add filler sequence, 5' or 3'
remove_X       <- TRUE
add_M          <- TRUE
lower_flanks   <- TRUE
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
convert_stops <- function(input_sequence){
  new_input <- gsub(pattern = "\\*", replacement = "X", input_sequence)
  return(new_input)
}

## Input sequence health check
check_input_seq <- function(input_sequence, min_len, flank_5, flank_3){
  warnings <- ""
  ## Length
  min_length_wflanks <- min_len - nchar(flank_5) - nchar(flank_3)
  leftover_length <- (nchar(input_sequence) * 3) - min_length_wflanks
  if(leftover_length < 0){
    warnings <- paste0(warnings, "WARN: Input sequence is ", abs(leftover_length), " short of minimum length. Making up length with filler.\n")
  }

  ## Starting with M
  starting_aa <- strsplit(input_sequence, split = "")[[1]][1]
  if(starting_aa != "M"){
    warnings <- paste0(warnings, "WARN: Input sequence starts with ", starting_aa, ", not M. Please ensure this is correct.\n")
  }

  ## Internal stops
  #grepl(pattern = "(X\\*)", "MSXAV")
  tokenised_seq <- strsplit(input_sequence, split = "")[[1]]
  if(grepl(pattern = "X", input_sequence)){
    stops <- which(tokenised_seq == "X")
    warnings <- paste0(warnings, "WARN: Stop codon identified in input sequence position(s) ", stops, " of ", nchar(input_sequence), ". Please ensure this is correct.\n")
  }
  if(warnings != ""){
    message(warnings)
  }
}

## Function to optimize codon usage based on frequency
optimise_codon_seq <- function(input_sequence, rand_seed, remove_X, add_M){
  optimised_seq <- ""
  set.seed(rand_seed)
  tokenised_seq <- strsplit(input_sequence, split = "")[[1]]
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
  return(optimised_seq)
}

## Function that checks if site to avoid is present in intial optimised sequence. If so, it increments seed by 1 and tries again. If after `cut_iterations` no solution without restriction site is found, it fails, but still emits the sequence. 
restriction_site_check <- function(input_sequence, rand_seed, stop_counter, cut_site_avoid){
  contains_site <- grepl(pattern = cut_site_avoid, input_sequence)
  if(contains_site == FALSE){
    return(input_sequence)
  } else if(contains_site == TRUE){
    counter  <- 0
    while(contains_site == TRUE){
      counter   <- counter + 1
      temp_seed <- rand_seed + 1
      temp_seq  <- optimise_codon_seq(input_sequence, rand_seed)
      contains_site <- grepl(pattern = cut_site_avoid, temp_seq)
      if(counter == stop_counter){
        break
        warnings <- paste0("ERROR: Completed ", counter , " iterations of random codon optimisation. Cannot remove cut site ", cut_site_avoid, "\nIncrease iterations or change seed and try again. If problem persists, manually inspect.")
        ## If this becomes an issue, add debug details with coordinates and output all iterations to log file.  
      }
    }
    return(temp_seq)
  }
}

## Function to add flanks and add filler if needed. Filler will be lowercase just to easily identify. 
add_flanks <- function(input_sequence, x5_flank, x3_flank, filler_seq, min_length, filler_loc, lower_flanks){
  if(lower_flanks == TRUE){
    x5_flank <- tolower(x5_flank)
    x3_flank <- tolower(x3_flank)
  } else if(lower_flanks == FALSE){
    x5_flank <- toupper(x5_flank)
    x3_flank <- toupper(x3_flank)    
  }
  with_flanks <- paste0(x5_flank, input_sequence, x3_flank)

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
  return(out_seq)
}

output_check <- function(input_sequence, cut_site_avoid){
  warnings <- ""

  ## Second rescriction site check
  num_restrictions <- str_count(pattern = cut_site_avoid, input_sequence)
  if(num_restrictions > 1){
    warnings <- paste0(warnings, "ERROR: ", num_restrictions, " sites to avoid detected with sequence ", cut_site_avoid, "\n")
  }

  ## Internal stop codon check
  translated <- DNAString(input_sequence) %>% translate
  num_int_stops <- str_count(pattern = "X", as.character(translated))
  if(num_int_stops > 0){
    warnings <- paste0(warnings, "ERROR: ", num_int_stops, " stop codons detected with sequence.\n")
  }

  ## Sanity check against input
  match_input <- grepl(pattern = translated, input_sequence)
  if(match_input == FALSE){
    ## Removing start and end to see if pattern can be found since these can be added/removed during optimisation. If not, a problem might be present.
     match_attempt2 <- grepl(pattern = str_sub(translated, start = 2, end = -2), input_sequence)
     if(match_attempt2 == FALSE){
      warnings <- paste0(warnings, "ERROR: Input sequence and translated optimised sequence do not match.")
     }
  }

  if(warnings != ""){
    message(warnings)
  }
}





## Testing funcitons - works well enough.
# Main functionality
seq1 <- convert_stops(input_seq) 
check_input_seq(seq1, min_length, x5_flank, x3_flank)
seq2 <- optimise_codon_seq(seq1, seed, remove_X, add_M)
seq3 <- restriction_site_check(seq2, seed, cut_iterations, cut_site_avoid)
seq4 <- add_flanks(seq3, x5_flank, x3_flank, filler_seq, min_length, filler_loc, lower_flanks)
output_check(seq4, cut_site_avoid)




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