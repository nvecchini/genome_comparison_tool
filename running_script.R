####### Runs the genome comparison tool directly from script ###########


# run the code below to get all the necessary tools 
source(file = "scripts/Genome_comparison_tool.R")

# make a list of all genomes 
genomes <- list.files("genome_files", full.names = T)

# To run the analysis you are required to provide ref genome and the query genomes,
genome_ref <- genomes[str_detect("genome_file/reference.faa", genomes)]
genomes_query <-genomes[!str_detect(genome_ref, genomes)]
#create list of locus tags, can be whole genome or subsections can be selected 
locus_tags<-locus_tag_creation_vector(genome_ref)

# how to select a subset: locus_tags[1:?]

#run genome comparison 
mult_genome_comp<-Genome_comparison_tool(genomes_query,genome_ref,locus_tags)
saveRDS(mult_genome_comp,"rds/mult_genome_comp.rds")