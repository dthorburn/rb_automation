# Create the codon table with single-letter amino acids and frequencies
codon_table <- data.frame(
  Amino_Acid = c("A", "A", "A", "A", "R", "R", "R", "R", "N", "N",
                 "D", "D", "C", "C", "Q", "Q", "E", "E", "G", "G",
                 "G", "G", "H", "H", "I", "I", "I", "L", "L", "L",
                 "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
                 "P", "P", "S", "S", "S", "S", "S", "S", "T", "T",
                 "T", "T", "W", "Y", "Y", "V", "V", "V", "V"),
  Codon = c("GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG", "AAU", "AAC",
            "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG", "GGU", "GGC",
            "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA", "UUG", "CUU",
            "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC", "CCU", "CCC",
            "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "ACU", "ACC",
            "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA", "GUG"),
  Frequency = c(0.18, 0.36, 0.22, 0.24, 0.07, 0.29, 0.12, 0.52, 0.47, 0.53,
                0.44, 0.56, 0.43, 0.57, 0.28, 0.72, 0.26, 0.74, 0.16, 0.34,
                0.26, 0.24, 0.42, 0.58, 0.47, 0.39, 0.14, 0.13, 0.14, 0.13,
                0.07, 0.53, 0.53, 0.28, 0.72, 1.00, 0.46, 0.54, 0.23, 0.35,
                0.20, 0.22, 0.17, 0.31, 0.11, 0.12, 0.15, 0.14, 0.29, 0.24,
                0.31, 0.16, 1.00, 0.48, 0.52, 0.17, 0.30, 0.22, 0.31)
)

# Function to optimize codon usage based on frequency
optimize_codon_sequence <- function(aa_sequence) {
  # Initialize an empty string for the codon-optimized sequence
  codon_optimized_sequence <- ""
  
  # Loop through each amino acid in the input sequence
  for (aa in strsplit(aa_sequence, "")[[1]]) {
    # Subset the codon table for the current amino acid
    aa_codons <- subset(codon_table, Amino_Acid == aa)
    
    # Get the codon with the highest frequency for the current amino acid
    best_codon <- aa_codons[which.max(aa_codons$Frequency), "Codon"]
    
    # Append the best codon to the codon-optimized sequence
    codon_optimized_sequence <- paste0(codon_optimized_sequence, best_codon)
  }
  
  # Return the codon-optimized sequence
  return(codon_optimized_sequence)
}

# Example usage
aa_sequence <- "MVKTA"
optimized_sequence <- optimize_codon_sequence(aa_sequence)
print(optimized_sequence)
