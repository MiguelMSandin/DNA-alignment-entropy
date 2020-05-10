# DNA-alignment-entropy
Calculate Shannon entropy for every position of an aligned fasta file.

# Requirements
1. R (or R-studio):
   - libraries: "vegan" & "seqinr"
   - optional library: "ggplot2"
2. An aligned fasta file.

# Usage
Basically, this R script calculates the Shannon entropy of every position in an aligned fasta file. 
For the entropy calculation it takes the * *diversity* * function from the library * *vegan* * (H = -sum p_i log(b) p_i).

For **a simple usage**, copy the Rscript [entropy.R](https://github.com/MiguelMSandin/DNA-alignment-entropy/blob/master/entropy.R) file to the working directory where you have your aligned fasta file (```FILE_NAME.fasta```) and run: 
```
Rscript entropy.R FILE_NAME.fasta 
```

This will create a TSV file (```FILE_NAME_positionInfo.tsv```) with the following columns:
1. posi:  Position
2. shan:  Shannon entropy of given position
3. shanc: Shannon entropy of given position removing gaps ("-")
4. rich:  Richness in given position
5. richc: Richness in given position removing gaps ("-")
6. uniq:  Unique bases in given position
7. repe:  Repetitions of the unique bases in given position

For **a more elaborated usage** check the script [entropy_plotting.R](https://github.com/MiguelMSandin/DNA-alignment-entropy/blob/master/entropy_plotting.R). Recommended to open it in R-studio. 
With this script you can plot different entropy curves for different group of sequences within the same aligned fasta file.
The script has comments for main functions and sentences to ease the understanding. Otherwise don't hesitate to contact me.


**Lastly**, I'm not a (bio)informatitian, so if you want to upgrade/improve or efficiently rewrite the code (in any other language) I'll be very happy to collaborate.

