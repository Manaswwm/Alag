# callpoint for getting the longest transcript for a given gene from its gff file
#this function requires that the gff file (subset) is loaded by default

get_largest_transcript_cds_info = function(geneID_list){
  
  #converting data type
  geneID = as.character(geneID_list)
  
  #have to change the datatype
  ingroup_gff_subset$attributes = as.character(ingroup_gff_subset$attributes)
  
  ## since I have done most of the heavy lifting in identifying the gene specific transcript - I will not just use the gene conversion table to subset the gff
  #keeping the original structure as this helps me in getting the cds_info as list of dataframes
  get_cds_intervals = function(geneID){
    
    #accessing the peptide IDs per gene - this will be used to subset the gff
    transcriptIDs = gene_conversion_table$ensembl_transcript_id[gene_conversion_table$ensembl_gene_id %in% geneID]
    
    #subsetting the gff file now
    df = ingroup_gff_subset[ingroup_gff_subset$attributes %in% transcriptIDs,]
    
    #subsetting columns to keep only those of interest
    df = subset(df, select = -c(source, type, score, phase))
    
    #giving back the geneIDs and transcriptIDs per peptideIDs
    df = merge(df, gene_conversion_table, by.x = "attributes", by.y = "ensembl_transcript_id")
    
    #changing column names and orientation
    df = df[,c("seqid", "start", "end", "strand", "ensembl_gene_id", "attributes")]
    
    #changing column name of the last column to peptide id
    colnames(df)[6] = "ensembl_transcript_id"
    
    #returning
    return(df)
  }
  
  #giving call to every gene sequentially and accessing its cds info
  cds_info = mclapply(geneID_list, function(x) get_cds_intervals(x), mc.cores = 1)
  
  #### in some rare cases annotated genes do not have a CDS - in this case removing them ####
  cds_info = cds_info[sapply(cds_info, function(x){!dim(x)[1] == 0})]
  
  #retunring longest transcript and cds_info as a single dataframe
  return(cds_info)
}

## code snippet 1 begins ##

#### what I would have to run if I did not have the access to ensembl MANE transcript for humans -- instead would have to choose the longest transcript ####
# #going over all ids sequentially and storing their info on longest transcript
# #and their coding coordinates
# get_cds_intervals = function(geneID){
# 
#   #checking for all the available transcripts - using stringr dirty trick
#   
#   ## changed from the original code;  
#   ## grep(paste(geneID,"_.[[:digit:]]",sep=""), ingroup_gff_subset$attributes)])
#   ## because hirsuta's transcript ID format= Chir01Ox_b00050.1_T001 - with a '_'
#   
#   transcriptIDs = unique(ingroup_gff_subset$attributes[
#     grep(paste(geneID,".[[:digit:]]",sep=""), ingroup_gff_subset$attributes)])
#   
#   
#   #getting the largest transcript - condensed function
#   longest_transcript = transcriptIDs[which.max(lapply(transcriptIDs, function(x){
#     sum(ingroup_gff_subset$end[ingroup_gff_subset$attributes == x] - ingroup_gff_subset$start[ingroup_gff_subset$attributes == x])}))]
#   
#   #getting the subset of the thaliana gff for the longest transcript
#   df = ingroup_gff_subset[ingroup_gff_subset$attributes == longest_transcript,]
#   
#   #subsetting the dataframe to contin required info only
#   df = subset(df, select = -c(source, type, score))
#   
#   #giving back the geneID
#   df$geneID = rep(geneID)
#   
#   #removing rownames
#   rownames(df) = NULL
#   
#   #chaning column names for consistency
#   colnames(df) = c("seqid", "start", "end", "strand", "phase", "transcriptID", "geneID")
#   
#   #collecting information on the longest transcript here
#   return(df)
# 
#   }#closing the for-loop - gone over all geneIDs - collected longest transcript and cds coords
# 
# #calling the function here
# cds_info = mclapply(geneID_list, function(x) get_cds_intervals(x), mc.cores = 1)


##code snippet 1 ends ##


  
  # #accounting for phasing
  # cds_info$start = cds_info$start + cds_info$phase
  # 
  # #processing the dataframe before processing
  # fa1 = readDNAStringSet("dummy_input_files/Athaliana_seq/Athaliana_1.fa")
  # 
  # #using forloop for going over all the rows - may be room for optimization
  # cds_seq = ""
  # 
  # for(item in 1:nrow(cds_info)){
  #   
  #   row = cds_info[item,]
  # 
  #   seq = toString(subseq(fa1, start = row$start, end = row$end))
  #   
  #   cds_seq = paste(cds_seq, seq, sep = "")
  #   
  # }
  # 
  # #getting protein sequence
  # cds_seq = s2c(cds_seq)
  # 
  # cds_seq = paste(translate(cds_seq), collapse = "")  