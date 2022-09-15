library(tidyverse)

### Extra functions to extract unmatched genes in query genomes and reference
# genomes 


unmatch_gene_retrival<-function(gene_comp_db,genome_db,length_largest_genome=5000){
  # this function needs the function fasta_converters to operate
  # This function obtains the unmatch genes of a genome comparison 
  
  #First create a vector of unique genomes used
  
  genome_files<- as_vector(unique(gene_comp_db$genome_query))
  # determine a number to create a dataframe of that size 
  n_rows<-length_largest_genome
  # creates the final DB to place the iterations 
  not_available_db<-data.frame(a=1:n_rows)

  
  # length of unmatch genes to cut out extra rows that are not important 
  lengths_col<-NULL
  # iterates through each genome to determine which genes did not match with any gene in our reference genome 
  #i<- genome_files[2]
  for (i in genome_files){
    #print(i)
    #file path to the genomes
    file<- paste(i,".faa",sep = "")
    file_path<- file.path(genome_db,file)
    
   # converts the fasta files to be used to determine gen nemes 
    genome_anno<-fasta_converster(file_path) %>% 
      distinct(id,.keep_all = TRUE)
    
    # extract all matched genes 
    inquiry<- gene_comp_db%>% 
      filter(str_detect(genome_query, i)) %>%
      select(query_sequence) %>% 
      separate(query_sequence,c("id","name"),sep = " ")
    #not_match<-data.frame(i=intersect(genome_file$id,inquiry$id)) 
    # creates a vector of genes that did not match  
    not_match<-intersect(genome_anno$id,inquiry$id)
   
    # Extract names of unmatched genes 
    unmatch_genes<-genome_anno %>% 
      filter(!id%in%not_match) %>%
      ungroup() %>% 
      select(sequence) %>% 
      as_vector()
    
    # store length of non matching genes 
    len_col<-length(unmatch_genes)
    # creates a vector of all length of unmatch genes 
    lengths_col<-append(lengths_col,len_col)
    
    #create the data frame with names unmateched genes 
    max_len<-n_rows
    not_match.1<-data.frame(unmatch_genes[1:max_len])
    colnames(not_match.1)[1]<-i
    not_available_db<-cbind(not_available_db,not_match.1) 
  }
  # removes the na and places a blank space 
  not_available_db[is.na(not_available_db)]<-""
  # removes the dummy first row 
  not_available_db<-not_available_db %>%
    select(-1) 
  #Determine which genome had the most unmatched genes 
  max_len<-max(lengths_col)
  # creates a final dataframe with all unmatched genes 
  unmatched_db_final<-not_available_db[1:max_len,]
  
  return(unmatched_db_final)
}
#unmatch_genes_amp_comp<-unmatch_gene_retrival(gene_comp_db,"genome_files")


#ref_genome<-"M_521"
unmatch_referece_gene_retrival<-function(gene_comp_db,genome_db,ref_genome,length_largest_genome=5000){
  # this function needs the function fasta_converters to operate
  # This function obtains the unmatch genes of a genome comparison 
  
  #First create a vector of unique genomes used
  
  genome_files<- as_vector(unique(gene_comp_db$genome_query))
  # determine a number to create a dataframe of that size 
  n_rows<-length_largest_genome
  # creates the final DB to place the iterations 
  not_available_db<-data.frame(a=1:n_rows)
  
  
  # length of unmatch genes to cut out extra rows that are not important 
  lengths_col<-NULL
  # iterates through each genome to determine which genes did not match with any gene in our reference genome 
  #i<- genome_files[1]
  for (i in genome_files){
    #print(i)
    #file path to the genomes
    file<- paste(ref_genome,".faa",sep = "")
    file_path<- file.path(genome_db,file)
    
    # converstes the fasta files to be used to determine gen nemes 
    genome_anno<-fasta_converster(file_path) %>% 
      distinct(id,.keep_all = TRUE)
    
    # extract all matched genes 
    inquiry<- gene_comp_db%>% 
      filter(str_detect(genome_query, i)) %>%
      na_if("") %>% 
      filter(is.na(query_sequence)) %>% 
      select(ref_sequence) %>% 
      as_vector()
    
    # store length of non matching genes 
    len_col<-length(inquiry)
    # creates a vector of all length of unmatch genes 
    lengths_col<-append(lengths_col,len_col)
    
    #create the data frame with names unmateched genes 
    max_len<-n_rows
    not_match.1<-data.frame(inquiry[1:max_len])
    colnames(not_match.1)[1]<-i
    not_available_db<-cbind(not_available_db,not_match.1) 
  }
  # removes the na and places a blank space 
  not_available_db[is.na(not_available_db)]<-""
  # removes the dummy first row 
  not_available_db<-not_available_db %>%
    select(-1) 
  #Determine which genome had the most unmatched genes 
  max_len<-max(lengths_col)
  # creates a final dataframe with all unmatched genes 
  unmatched_db_final<-not_available_db[1:max_len,]
  
  return(unmatched_db_final)
}



