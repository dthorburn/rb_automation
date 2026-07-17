# Codon Optimiser for *Nicotiana benthamiana*

Codon Optimiser is a bioinformatics tool designed to translate protein sequences into DNA sequences optimised specifically for expression in *Nicotiana benthamiana*. The goal of this project is to provide sequences that balance high expression potential with biological stability, ensuring they are ready for synthesis and *in planta* testing.

---

## Table of Contents
1. [Usage](#usage)
2. [Features & Methodology](#features--methodology)
3. [Parameters](#parameters)
4. [Results Interpretation](#results-interpretation)

---

## Usage

To use the optimiser, run the script via Python 3. You must provide a FASTA file with protein sequences as the primary input. Below is an example command using several custom constraints:

```bash
python3 CodonOptimiser.py Effectors.fasta --window 50 --gc_floor 28 --homopolymer 10 --save_f --use_ime
```



## Features & Methodology

### 1. Codon Optimisation via Stochastic Sampling
Instead of a rigid "one-codon-one-amino-acid" replacement (which can cause metabolic stress and repetitive sequences), this script utilises stochastic sampling:

* **Natural Frequency:** Codons are sampled based on their actual usage frequency in the *N. benthamiana* genome.
*Obtained from: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4100*
* **Iterative Refinement:** The algorithm performs up to 2,000 iterations per protein to find a candidate that satisfies all constraints.
* **Balanced Diversity:** This approach ensures the sequence is diverse enough for DNA synthesis while remaining highly readable by the host.

### 2. Codon Adaptation Index (CAI)
The script calculates the CAI to estimate translation efficiency based on the work of Sharp and Li (1987).


It uses the geometric mean of relative adaptiveness values ($w$):

$$CAI = \left( \prod_{i=1}^{L} w_i \right)^{1/L}$$

* **The Score:** CAI ranges from 0 to 1.
* **Why 0.7–0.8?** A score of 1.0 is mathematically "perfect" but biologically unrealistic, as it would use only the single most frequent codon for every amino acid, potentially depleting specific tRNA pools. A score between 0.7 and 0.8 is considered optimal for *Benthi*.

### 3. GC Content Regulation & Spliceosome Avoidance
Plants identify introns primarily by high AT content and specific AU-rich repeats. To prevent your mRNA from being incorrectly processed and degraded, the script "hides" it from the spliceosome:

* **GC Floor:** Maintains a minimum GC percentage (Target: ~44% for *Benthi*) to prevent the sequence from appearing "intron-like."
* **Sliding Window:** Users can adjust the `--window` size to increase the stringency of the GC check.
* **Homopolymer Removal:** Long repeats (e.g., AAAAAA or TTTTTT) are automatically removed to ensure stability.

### 4. Algorithmic Sequence Repair
The script ensures sequences are ready for **Golden Gate Assembly**:

* **Automated Removal:** Scans and removes BsaI sites (`GGTCTC`/`GAGACC`) by default.
* **Customisation:** Additional sites (e.g., NotI, XbaI) can be excluded via the `--exclude` parameter.
* **Synonymous Swaps:** The script breaks restriction sites using synonymous mutations, ensuring the protein sequence remains 100% identical.

### 5. IV2 Intron Insertion (IME)
For difficult-to-express sequences, the script can implement **Intron-Mediated Enhancement (IME)** by inserting the IV2 Intron from the *ST-LS1* gene (Vancanneyt et al., 1990).

* **Optimal Targeting:** It searches the first 60bp for the optimal `AG|G` motif.
* **Fallback:** If no motif is found, the intron is inserted immediately following the Start Codon (`ATG`).

### 6. Audit Logs & Verification
* **Failure Breakdown:** If a sequence fails within the 2,000-iteration limit, the script logs exactly why (e.g., `Effector_581 | {'Low GC (26.0%)': 1700}`).
* **Sanity Checks:** Every generated sequence is back-translated to confirm 100% protein identity. If IME is active, the script "virtually splices" the DNA before translation to verify correctness.

---

## Parameters

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--window` | `int` | `50` | Size of the sliding window for GC content checks. |
| `--gc_floor` | `float` | `35.0` | Minimum GC percentage required within the window. |
| `--exclude` | `str` | `GGTCTC,GAGACC` | Comma-separated restriction sites to avoid. |
| `--homopolymer` | `int` | `6` | Maximum allowed length of a single-nucleotide repeat. |
| `--use_ime` | `flag` | `False` | Enables IV2 intron insertion (Intron-Mediated Enhancement). |
| `--save_f` | `flag` | `False` | Saves sequences that fail constraints but pass translation check. |

---

## Results Interpretation

* **Success:** The sequence passed all biological (GC, Homopolymer) and synthesis (Restriction Sites) constraints.
* **Best Effort:** If using `--save_f`, these sequences maintain the correct protein identity but may contain a minor constraint violation that could not be resolved.
* **CAI Target:** Aim for results in the **0.75** range for the best balance of expression and reliability in *N. benthamiana*.