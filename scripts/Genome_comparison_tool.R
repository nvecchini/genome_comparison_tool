#packages for the tools 

#Loading Libraries 
library(tidyverse)
library(seqinr)

###################### Functions used by the Genome comparison tool###################


# Converts fasta files protein sequences to dataframes

fasta_converster<-function(fasta_file){
    db_rt<- read.delim(fasta_file,header = FALSE,quote="") %>%
    mutate(id=str_extract(V1,"[>][^ ]+"))%>% fill(id) %>% 
    mutate(id2=ifelse(str_detect(V1,">"),1,2)) %>% 
    mutate(id=str_replace(id,">","")) %>%
    group_by(id,id2) %>% 
    summarise(sequence=paste(V1,collapse = "")) %>% 
    select(id,sequence) %>%
    mutate(sequence=str_replace(sequence,">",""))
  return(db_rt)
}

# Extract the locus tags of a fasta file for query. 

locus_tag_creation_vector<- function(ref_genome_file){
  locus<-fasta_converster(ref_genome_file)
  locus_tag<-locus %>% 
    distinct(id) %>%  
    pull(id)
  return(locus_tag)
}

#Create blast database 

#path<- file.path("./fasta_files","C_auto.faa")
create_blast_db<-function(path){
  file.name <- str_split(path, "/") %>% map(., last)
  sample.name <- str_split(file.name, "\\.") %>% map(., dplyr::first) %>% unlist()
  out.path<- paste("blast_db/",sample.name,sep = "")
  if(file.exists(paste(out.path,".pdb",sep = ""))=="TRUE"){
    print("db is available")
  } else { 
    system2('makeblastdb', args = c('-in',path, '-out', out.path, '-dbtype', "prot"))
  }
  return(sample.name)
}

# Converts the blast+ output into a dataframe

out_file_converter<-function(file_path){
  x<- read.delim(file_path,skip = 23)
  length<- str_extract(x[1,1],"[0-9]+")
  names(x)[1]<- "col1"
  x2<- x %>% 
    filter(str_detect(col1,"[^ ]+_[0-9]{5}")) %>%
    filter(!str_detect(col1,"Length")) %>% 
    separate(col1,c("Sequence","Score(Bits)","Evalue"),"\\s\\s+") %>%
    mutate(id=str_extract(Sequence,"[^ ]+")) %>% 
    select(4,1:3)
  x3<- x %>%
    filter(str_detect(col1,"[>][^ ]+")|str_detect(col1,"Iden")) %>% 
    mutate(id=str_extract(col1,"[>][^ ]+")) %>%
    fill(id) %>% 
    select(id,col1) %>% 
    filter(str_detect(col1,"Iden")) %>%
    mutate(id=str_replace(id,">","")) %>% 
    separate(col1,c("Identities","Positives","Gaps"),",") %>% 
    mutate(Similarity = str_extract(Identities,"[0-9]?[0-9]?[0-9]?%")) %>% 
    mutate(Similarity = str_replace(Similarity,"%",""),Identities = str_replace(Identities,"Identities = ",""),Positives = str_replace(Positives,"Positives = ",""),Gaps = str_replace(Gaps,"Gaps = ","")) %>% 
    left_join(x2,by="id") %>% 
    select(id,Sequence,7,8,2:5)
  return(x3)
}

######################### Genome Comparison tool ########################


# Function to iterate through each genome 


Genome_comparison_tool <- function(genome_entry, genome_ref, locus_tag){
  print(Sys.time())
  result=NULL
  print(genome_entry)
  print(genome_ref)
  print(locus_tag)
  genome.name <- str_split(genome_ref, "/") %>% map(., last)
  genome_name.1 <- str_split(genome.name, "/") %>% map(., dplyr::last) %>% unlist()
  genome_name<- str_remove(genome_name.1,"[.][A-z]{1,}")
  logfile_name<-paste(genome_name,"_",Sys.Date(),".Rout",sep = "")
  # log code
  log_out <- file(paste("outputs/",logfile_name,sep = ""), open = "wt")
  sink(log_out)
  sink(log_out, type = "message")
  try(log("a"))
  
  for(genome in genome_entry){
    result_1 <- Locus_comparison_tool(genome,genome_ref,locus_tag)
    result<-rbind(result,result_1)  

  }
  sink()
  print(Sys.time())
   return(result)
}

# Function to iterates through each locus_tag that will be compared 


Locus_comparison_tool<-function(genome_2, genome_1, locus_tags){
  
  #creates the database if necessary of genomes and then get the db name
  db_name_gen1<-create_blast_db(genome_1) 
  db_name_gen2<-create_blast_db(genome_2)
  db_path2<-paste("blast_db/", db_name_gen2, sep = "")
  db_path1<-paste("blast_db/", db_name_gen1, sep = "")
  
  
  # convert genome_1 to fasta
  db_gen_1<- fasta_converster(genome_1)%>% 
    rename(ref_sequence=sequence)
  db_gen_2<- fasta_converster(genome_2)
  # creates a vector 
  final_table_4=NULL
  locus_tag<-"Ckluy_00001"
for (locus_tag in locus_tags){ 
  print("start")
  #add the genome name
  genome.name <- str_split(genome_2, "/") %>% map(., last)
  genome_name.1 <- str_split(genome.name, "/") %>% map(., dplyr::last) %>% unlist()
  #genome_name <- str_split(genome.name.1, "\\.",fixed=TRUE) %>% map(., dplyr::first) %>% unlist()
  genome_name<- str_remove(genome_name.1,"[.][A-z]{1,}")
  
  db_path<-paste("blast_db/db_name_gen2")
  
  print(locus_tag)
  sequence_1<- db_gen_1 %>% 
    filter(id %in% locus_tag)
  #write a fasta output to be used for the out file 
  write.fasta(sequence_1[2,2], names =sequence_1[2,1] , file.out = "./fastas/text.fasta")
  #runs blast on the locus tag
  system2("blastp", args = c('-query',"./fastas/text.fasta", '-out',locus_tag, '-db', db_path2, "-num_alignments","100") ,stderr = "stderr.txt")
  #move blast file into a temporary file in unix
  system2("mv",args = c(locus_tag,"./temp_blast/"),wait = F)
  #Moves in windows 
  cmd_1<-paste("MOVE",as.character(locus_tag),"./temp_blast/")
  shell(cmd_1)
  
  
  locus_path<-file.path(".","temp_blast",locus_tag)
  #read the out file 
  out_file_1<-out_file_converter(locus_path)
  # final table with first hit  
  final_table<-sequence_1[1,] %>% 
    append(out_file_1[1,]) %>% 
    as.data.frame() %>% 
    rename(query_sequence=Sequence) %>% 
    mutate(Evalue=as.numeric(Evalue),Score.Bits.=as.numeric(Score.Bits.))
  
  system2("rm",args = c("-f","./temp_blast/*"))
  #paste("del","/q",)
  shell("del /q .\\temp_blast ")
  
  #Reciprocal check 
  #extract locus tag
  locus_tag_2.1<-final_table[1,3]
  #filter sequences
  sequence_2<- db_gen_2 %>% 
    filter(id %in% locus_tag_2.1)
  write.fasta(sequence_2[2,2], names =sequence_2[2,1] , file.out = "./fastas/text.fasta")
  #runs blast on the locus tag
  system2("blastp", args = c('-query',"./fastas/text.fasta", '-out',locus_tag_2.1, '-db', db_path1, "-num_alignments","100") ,stderr = "stderr.txt")
  #move blast file into a temporary file 
  system2("mv",args = c(locus_tag_2.1,"./temp_blast/"),wait = F)
  # MOve files in windows 
  cmd_2<-paste("MOVE",as.character(locus_tag_2.1),"./temp_blast/")
  shell(cmd_2)
  
  
  locus_path_2<-file.path(".","temp_blast",locus_tag_2.1)
  #read the out file 
  out_file_2<-out_file_converter(locus_path_2)
  a<-as_vector(final_table[1,2])
  b<-as_vector(final_table[1,4])
  c<-as_vector(out_file_2[1,2])
  print(a)
  print(b)
  print(c)
  
  system2("rm",args = c("-f","./temp_blast/*"))
  #paste("del","/q",)
  shell("del /q .\\temp_blast ")
  
  
  if(final_table[1,1]==out_file_2[1,1]|is.na(out_file_2[1,1])==TRUE){
    final_table_2<-final_table %>% 
      select(ref_sequence,query_sequence,Score.Bits.,Evalue,Identities,Positives,Gaps,Similarity) 
   
  }else{
    final_table[1,3]<-""
    final_table_2<-final_table%>% 
      select(ref_sequence,query_sequence,Score.Bits.,Evalue,Identities,Positives,Gaps,Similarity)
  }
  if(final_table_2[1,3]<50|is.na(out_file_2[1,3])==TRUE){
    final_table_2[1,2]<-""
    final_table_3<-final_table_2
  }else if(final_table_2[1,4]>0.005|is.na(out_file_2[1,4])==TRUE){
    final_table_2[1,2]<-""
    final_table_3<-final_table_2
    
  }else{
    final_table_3<-final_table_2
  }
  final_table_3<-final_table_3 %>% 
    mutate(genome_query = genome_name)
 
  #remove locus file 
  #system2("rm",args = c("-f","./temp_blast/*"))
  final_table_4<-rbind(final_table_4,final_table_3)
   print("end")
}

  return(final_table_4)
}

# New function full genome comparison compares a whole directory to itself
# every genome is treated as a new query

full_folder_genome_comparison<-function(direc){
  # creates a directory with all genomes 
  genomes <- list.files(direc , full.names = T)
  result=NULL
  
  log_out <- file(paste("outputs/",logfile_name,sep = ""), open = "wt")
  sink(log_out)
  sink(log_out, type = "message")
  try(log("a"))
  
  
  for(genome_1 in genomes){
    
    genome_ref <- genomes[str_detect(genome_1, genomes)]
    genomes_query <-genomes[!str_detect(genome_ref, genomes)]
    #create locus tags 
    locus_tags<-locus_tag_creation_vector(genome_ref)
    
    # run genome comparison 
    result_1<-Genome_comparison_tool(genomes_query,genome_ref,locus_tags)
    result<-rbind(result,result_1)
  }
  sink()
  return(result)
}
