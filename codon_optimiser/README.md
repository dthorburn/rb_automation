## RB Codon Optimiser Project
### Outline
The Twist BioScience website codon optimisation tool isn’t great. You can only process 1 sample at a time, and it chooses the sequence randomly which makes continuity difficult to achieve. 
The idea here is to solve that by launching a hosted RShiny App that provides repeatability measures like a seed.

### Updates
- 11/09/2024: Project initialised and draft R script to write and test functions begins

### App Outline
**Input:** Amino acid sequences
**Output:** *Nicotiana benthamiana* codon optimised DNA sequences with appropriate flanks. 

**Parameters:** 

1. Seed (default: 2021)
2. 5’ and 3’ flanks (default: current use)
3. Sequences to avoid (default: BSAI cut site)
4. Button for batch addition of flanks (is this necessary?)

**Complications:**

1. Unsure how seed will interact with hosted app. Can set seed at start of *every* calculation. Is the randomness necessary?
2. Identify if BSAI site is located in region with enough codon variability to be solved. 
3. How to solve when BSAI is identified: (i) use `while` loop to iterate through each codon in site incrementally changing until site disappears, or (ii) increment seed when randomly selecting sites and keep iterating until solution is found (limit to *n* iterations before throwing error).
4. Visualisation will not be as good as Twist’s website (or are there libraries for this?). 
5. Warnings for common mistakes like missing stop codon, missing `ATG`, length, etc…
6. If sequences are too short (<300bp), add consistent filler sequence at either 5’ or 3’ end.
7. How much compute does hosted service have? Can i use a vectorised function (i.e., `lapply`), or should I play it safe with a `for` loop? Can I implement both and use `tryCatch`?

