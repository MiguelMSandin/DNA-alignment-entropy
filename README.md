# DNA-alignment-entropy
Calculate Shannon entropy for every position of an aligned fasta file.

# Requirements
1. R: libraries: "vegan" & "seqinr"
   - optional library: "ggplot2"
2. An aligned fasta file.

# Usage
Basically, this R script calculates the Shannon entropy of every position in an aligned fasta file. 
For the entropy calculation it uses the *diversity* function from the library *vegan* (H = -sum p_i log(b) p_i).

For **a simple usage**, copy the R script [entropy.R](https://github.com/MiguelMSandin/DNA-alignment-entropy/blob/master/entropy.R) file to the preferred directory and run: 
```
Rscript path/to/entropy.R path/to/FILE_NAME.fasta 
```

This will take the aligned fasta file (```FILE_NAME.fasta```) and create a TSV file (```FILE_NAME_positionInfo.tsv```) with the following columns for the given position:
1. posi:  Position
2. shan:  Shannon entropy
3. shanc: Shannon entropy removing gaps ("-")
4. rich:  Richness
5. richc: Richness removing gaps ("-")
6. uniq:  Unique bases
7. repe:  Repetitions of the unique bases

For **a more elaborated usage** check the script [entropy_plotting.R](https://github.com/MiguelMSandin/DNA-alignment-entropy/blob/master/entropy_plotting.R). **Recommended to open it in *R-studio***. 
With this script you can plot different entropy curves for different group of sequences within the same aligned fasta file.
The script has comments for main functions and sentences to ease the understanding. Otherwise don't hesitate to contact me.

**Lastly**, if you want to upgrade/improve or efficiently rewrite the code (in any other language) I'll be very happy to collaborate.

