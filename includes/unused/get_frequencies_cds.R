# call point for getting frequencies_cds dataframe
get_frequencies_cds = function(path_to_vcf_folder, transcriptID){
  
  # #declaring an empty dataframe to fill in frequencies_cds data
  # frequencies_cds = data.frame()
  
  #print check
  print(paste("Making frequencies_cds df for transcript ID : ",transcriptID, sep = ""))
  
  #loading transcript information now
  load(paste(backup_file,"transcript_alignments",sep = ""))
  
  ## calculating frequency of synonymous variants ##
  #checking if the positions file exists
  if(file.exists(paste(path_to_vcf_folder,transcriptID,"_synonymous_positions",sep=""))){
    
    # giving call to vcftools for calculating the frequency of variants
    calculate_freq_command_syn = paste("vcftools --vcf ",paste(path_to_vcf_folder,transcriptID,".vcf",sep = "")," --positions ",
                                       paste(path_to_vcf_folder,transcriptID,"_synonymous_positions",sep = "")," --max-alleles 2 --max-missing 0.8 --freq2 --stdout", sep="")
    
    #print check
    print(calculate_freq_command_syn)
    
    #giving system call
    command_return = system(command = calculate_freq_command_syn, intern = T)
    
    #if there are variants within this region
    if(length(command_return)>1){
      
      #reading the output as table
      poly_syn_freq = read.table(text = command_return, skip = 1)
      
      #adding column names to the table
      colnames(poly_syn_freq) = c("chr", "start", "alleles", "nchr", "freq_ref", "freq_alt")
      
      #getting maf
      poly_syn_freq$maf = apply(poly_syn_freq[,c(5,6)],1 ,min)
      
      #loading vcf file to give the reference and alternate allele
      syn_position = read.table(paste(path_to_vcf_folder,transcriptID,"_synonymous_positions",sep = ""),header = TRUE)
      
      #adding back information on positions into the freq dataframe to get ref and alt allele
      poly_syn_freq = merge(poly_syn_freq, syn_position, by.x = "start", by.y = "pos")
      
      #cleaning the output dataframe
      poly_syn_freq = subset(poly_syn_freq, select = -c(chr.y))
      colnames(poly_syn_freq)[2] = "chr"
      
      ###### --- removing variants which are already fixed or lost within the population --- ######
      ## question - why was this not removed by default ##
      poly_syn_freq = poly_syn_freq[!poly_syn_freq$freq_alt %in% c(0,1),]
      
      #checking if the dataframe is empty
      if(!dim(poly_syn_freq)[1] == 0){
      
        #creating empty dataframe to store the ingroup and outgroup information
        ancestral_state_df = data.frame()
        
        #getting the ingroup and outgroup of variants -- using the transcript alignments here
        for(item in 1:nrow(poly_syn_freq)){
          
          #accessing individual row here
          row = poly_syn_freq[item,]
          
          #accessing the transcript alignment belonging to this transcript
          transcript_info = transcript_alignments[sapply(transcript_alignments, function(x) x$id == transcriptID)]
          
          #outgroup information -- divergence species
          outgroup = transcript_info[[1]]$alignment[,which(colnames(transcript_info[[1]]$alignment) == row$start)][2]
          
          #ingroup information - species of interest
          ingroup = transcript_info[[1]]$alignment[,which(colnames(transcript_info[[1]]$alignment) == row$start)][1]
          
          #storing start position -- to merge later
          start = row$start
          
          #storing df
          df = data.frame(start = start, ingroup = ingroup, outgroup = outgroup)
          
          #merging with df
          ancestral_state_df = rbind(ancestral_state_df, df)
          
        }
        
        #merging ancestral state with original dataframe
        poly_syn_freq = merge(poly_syn_freq, ancestral_state_df, by = "start")
        
        #subsetting this dataframe to remove positions for which there is a gap in either in or outgroup
        poly_syn_freq = poly_syn_freq[!(poly_syn_freq$ingroup == "-" | poly_syn_freq$outgroup == "-"),]
        
        #cheking if the dataframe is now empty
        if(dim(poly_syn_freq)[1] > 0){
          
          #polarizing the data
          poly_syn_freq = ldply(apply(poly_syn_freq, 1,function(x) polarize_data(x)), rbind)
          
          #checking for dimensions
          if(dim(poly_syn_freq)[1] > 0){
            
            #cleaning up a bit
            poly_syn_freq = subset(poly_syn_freq, select = -c(.id))
            
            #changing date type here
            poly_syn_freq$freq_der = as.numeric(as.character(poly_syn_freq$freq_der))
            poly_syn_freq$freq_anc = as.numeric(as.character(poly_syn_freq$freq_anc))
            
            #adding annotation here
            poly_syn_freq$type = rep("syn")
            
          }else{poly_syn_freq = data.frame()}
        }else{poly_syn_freq = data.frame()}
      }else{poly_syn_freq = data.frame()}
    }else{poly_syn_freq = data.frame()}
    
  }else{poly_syn_freq = data.frame()}#gone over synonymous positions for a transcript
  
  ## calculating frequency of nonsynonymous variants ##
  #checking if the positions file exists
  if(file.exists(paste(path_to_vcf_folder,transcriptID,"_nonsynonymous_positions",sep=""))){
    
    # giving call to vcftools for calculating the frequency of variants
    calculate_freq_command_nonsyn = paste("vcftools --vcf ",paste(path_to_vcf_folder,transcriptID,".vcf",sep = "")," --positions ",
                                          paste(path_to_vcf_folder,transcriptID,"_nonsynonymous_positions",sep = "")," --max-alleles 2 --max-missing 0.8 --freq2 --stdout", sep="")
    
    #giving system call
    command_return = system(command = calculate_freq_command_nonsyn, intern = T)
    
    #if there are variants within this region
    if(length(command_return)>1){
      
      #reading the output as table
      poly_nonsyn_freq = read.table(text = command_return, skip = 1)
      
      #adding column names to the table
      colnames(poly_nonsyn_freq) = c("chr", "start", "alleles", "nchr", "freq_ref", "freq_alt")
      
      #getting maf
      poly_nonsyn_freq$maf = apply(poly_nonsyn_freq[,c(5,6)],1 ,min)
      
      #loading vcf file to give the reference and alternate allele
      nonsyn_position = read.table(paste(path_to_vcf_folder,transcriptID,"_nonsynonymous_positions",sep = ""),header = TRUE)
      
      #adding back information on positions into the freq dataframe to get ref and alt allele
      poly_nonsyn_freq = merge(poly_nonsyn_freq, nonsyn_position, by.x = "start", by.y = "pos")
      
      #cleaning the output dataframe
      poly_nonsyn_freq = subset(poly_nonsyn_freq, select = -c(chr.y))
      colnames(poly_nonsyn_freq)[2] = "chr"
      
      ###### --- removing variants which are already fixed or lost within the population --- ######
      ## question - why was this not removed by default ##
      poly_nonsyn_freq = poly_nonsyn_freq[!poly_nonsyn_freq$freq_alt %in% c(0,1),]
      
      #checking if the dataframe is empty
      if(!dim(poly_nonsyn_freq)[1] == 0){
      
        #creating empty dataframe to store the ingroup and outgroup information
        ancestral_state_df = data.frame()
        
        #getting the ingroup and outgroup of variants -- using the transcript alignments here
        for(item in 1:nrow(poly_nonsyn_freq)){
          
          #accessing individual row here
          row = poly_nonsyn_freq[item,]
          
          #accessing the transcript alignment belonging to this transcript
          transcript_info = transcript_alignments[sapply(transcript_alignments, function(x) x$id == transcriptID)]
          
          #outgroup information -- divergence species
          outgroup = transcript_info[[1]]$alignment[,which(colnames(transcript_info[[1]]$alignment) == row$start)][2]
          
          #ingroup information - species of interest
          ingroup = transcript_info[[1]]$alignment[,which(colnames(transcript_info[[1]]$alignment) == row$start)][1]
          
          #storing start position -- to merge later
          start = row$start
          
          #storing df
          df = data.frame(start = start, ingroup = ingroup, outgroup = outgroup)
          
          #merging with df
          ancestral_state_df = rbind(ancestral_state_df, df)
          
        }
        
        #merging ancestral state with original dataframe
        poly_nonsyn_freq = merge(poly_nonsyn_freq, ancestral_state_df, by = "start")
        
        #subsetting this dataframe to remove positions for which there is a gap in either in or outgroup
        poly_nonsyn_freq = poly_nonsyn_freq[!(poly_nonsyn_freq$ingroup == "-" | poly_nonsyn_freq$outgroup == "-"),]
        
        #cheking if the dataframe is now empty
        if(dim(poly_nonsyn_freq)[1] > 0){
        
          #polarizing the data
          poly_nonsyn_freq = ldply(apply(poly_nonsyn_freq, 1,function(x) polarize_data(x)), rbind)
        
          #checking for dimensions
          if(dim(poly_nonsyn_freq)[1] > 0){
        
            #cleaning up a bit
            poly_nonsyn_freq = subset(poly_nonsyn_freq, select = -c(.id))
          
            #changing date type here
            poly_nonsyn_freq$freq_der = as.numeric(as.character(poly_nonsyn_freq$freq_der))
            poly_nonsyn_freq$freq_anc = as.numeric(as.character(poly_nonsyn_freq$freq_anc))
          
            #adding annotation here
            poly_nonsyn_freq$type = rep("nonsyn")
        
          }else{poly_nonsyn_freq = data.frame()}
        }else{poly_nonsyn_freq = data.frame()}
      }else(poly_nonsyn_freq = data.frame())
    }else{poly_nonsyn_freq = data.frame()}
    
  }else{poly_nonsyn_freq = data.frame()}#gone over nonsynonymous positions for a transcript
  
  
  ## done computing syn and nonsyn stats -- now combining them into a dataframe ##
  frequencies_cds = rbind(poly_syn_freq, poly_nonsyn_freq)
  
  #checking if the dataframe is empty before assiging the transcriptID
  if(!dim(frequencies_cds)[1]==0){
  
    #giving back the transcript ID here
    frequencies_cds$transcriptID = rep(transcriptID)
    
    #print check
    print(frequencies_cds)
    
    #returning frequencies cds 
    return(frequencies_cds)
  }else{return(frequencies_cds = data.frame())}
  
  
}