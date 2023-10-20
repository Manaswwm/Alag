#call point for getting the alignment of each row of the cds coordinates
get_alignment_cds = function(local_region){
  
  ### dummy input -- delete before launching the batch! ###
  cds_df = cds_info_sub[1,]
  
  ## part 1 - I need to get the sequence info of the local region in consideration ##
  
  #declaring the fasta sequence file from which I want to extract the sequence
  fa_file = readDNAStringSet(paste("dummy_input_files/Athaliana_seq/Athaliana_",as.character(cds_df$seqid),".fa", sep =""))
  
  #subsetting the sequence of interest
  fa_seq = toString(subseq(fa_file, start = cds_df$start, end = cds_df$end))
  
}