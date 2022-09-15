# Genome_comparison_tool

 With Genome_Comparison_tool.R microbial protein annotation files of a referece is compared to as many other microbial files. 
 The tool is used to determine orthology via reciprocal blasting  (Tools uses BLAST version 2.10.0+ with the default settings)[1] . 
 The tool works by creating separate databases(protein) for the reference sequence and queries. Then it performs reciprocal blast(blastp), 
 by blasting the reference sequence against the database of the query. The highest scoring hit is then taken and blasted against the reference database. 
 If the result is the same reference sequence, then that sequence is considered a putative orthologue.
 
 Example of tool output 
 ![Example long](https://github.com/nvecchini/Genome_comparison_tool/blob/main/data/example_images/long_formate.png)
 
 
 Wide output that can be generated 
 ![Example wide](https://github.com/nvecchini/Genome_comparison_tool/blob/main/data/example_images/Wide%20format.png)
 
# Before Use
1. Download the Genome comparison project and it should run in linux or windows
2. [Blast+] must be installed and setup to run globally
3. The packages below should be installed in R.
``` r
install.packages("tidyverse")
install.packages("seqinr")
```

# How to Use

Their are three available ways to run the tool: `Run_start_markdown.Rmd` , `running_script.R`, and `full_directory_run.R`

1. `Run_start_markdown.Rmd` allow you to run the tool directly from the Markdown
 - Instructions on how to use are in the markdown
 - Recommend to test in markdown before running the `.R` scripts 
 - More memory intensive 
 
 2. `running_script.R` Run the tool from cmd prompt or bash
 3. `full_directory_run.R` Run the tool from cmd prompt or bash
 - goes one by one using each file as a reference and running it against all other files 
 
 
 
 
# Packages used
 

 1.tidyverse[2]
 
 2.seqinr[3]


# Reference
1. Camacho, C. et al. (2009) BLAST+: architecture and applications. BMC Bioinformatics 10, 421. 10.1186/1471-2105-10-421
2. Wickham, H. et al. (2019) Welcome to the tidyverse. Journal of Open Source Software 4, 1686. 10.21105/joss.01686 
3. Charif D, Lobry J (2007). “SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.” In Bastolla U, Porto M, Roman H, Vendruscolo M (eds.), Structural approaches to sequence evolution: Molecules, networks, populations, series Biological and Medical Physics, Biomedical Engineering, 207-232. Springer Verlag, New York. ISBN : 978-3-540-35305-8.
