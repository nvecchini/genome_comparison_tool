


############ script will compare all files with each other in genome folder ####

# Extremly resource intensive.

# source code to run 
source(file = "scripts/Genome_comparison_tool.R")

full_comparison<-full_folder_genome_comparison("genome_files")

saveRDS(full_comparison,"rds/full_comparison.rds")