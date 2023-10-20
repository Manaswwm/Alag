# call point for getting the synonymous and nonsynonymous positions for every transcript
get_syn_nonsyn_positions_poly = function(path_to_vcf_folder, transcriptID, transcript_alignments){
  
  #loading cds_info to keep only those positions which are a part of this transcript
  load(paste(backup_file,"cds_info",sep=""))
  
  #loading transcript specific cds_info
  transcript_alignment_pos = as.integer(colnames(transcript_alignments[sapply(transcript_alignments, function(x) x$id == transcriptID)][[1]]$alignment))
  
  ## subset the syn and nonsyn positions to contain only those in cds_info
  
  #checking if there is a vcf file corresponding to the transcript ID
  if(file.exists(paste(path_to_vcf_folder,transcriptID,"_annotated.vcf",sep=""))){
  
    #print check
    print(paste("Making syn and nonsyn positions for transcript ID - ",transcriptID,sep=""))
    
    #reading the annotated vcf file from snpEff
    file = readLines(paste(path_to_vcf_folder,transcriptID,"_annotated.vcf",sep=""))
  
    #checking if read file is empty
    if(length(file)>0){
      
      #given that read file is non-empty -- re-reading this file with other packages
      vcf_pos = as.data.frame(info(readVcf(paste(path_to_vcf_folder,transcriptID,"_annotated.vcf",sep=""))))
      
      #accessing meta-information which contains position and chr information
      vcf_pos$meta_pos = row.names(vcf_pos)
      
      #declaring dataframes to store the syn and nonsyn positions
      syn_pos_df = data.frame()
      nonsyn_pos_df = data.frame()
      
      #going over all rows
      for(item in 1:nrow(vcf_pos)){
        
        #going over every row of the annotated vcf
        row = vcf_pos[item,]
        
        #extracting information on chromosome and position
        chr = str_extract_all(pattern = "[[:digit:]]+", string = row$meta_pos)[[1]][1]
        pos = str_extract_all(pattern = "[[:digit:]]+", string = row$meta_pos)[[1]][2]
        
        #extracting information on reference and alternate allele
        ref = unlist(str_split(string = sub(".*_","",row$meta_pos), pattern = "/"))[1]
        alt = unlist(str_split(string = sub(".*_","",row$meta_pos), pattern = "/"))[2]
        
        #extracting rows which are syn annotated
        if(length(grep("synonymous_variant", row$ANN)) > 0){
          
          #preparing syn variant position information for a transcript
          df = data.frame(chr = chr, pos = pos, ref = ref, alt = alt)
          syn_pos_df = rbind(syn_pos_df, df)
          
        }
        
        #extracting rows which are nonsyn annotated
        else if(length(grep("missense_variant", row$ANN)) > 0){
          
          #preparing nonsyn variant position information for transcript
          df = data.frame(chr = chr, pos = pos, ref = ref, alt = alt)
          nonsyn_pos_df = rbind(nonsyn_pos_df, df)
          
        }
        
      } #went over all rows of the annotated vcf -- now writing the position information to files
      
      #writing position information to a file now - only if information exists -- syn
      if(dim(syn_pos_df)[1] > 0){
        
        ############## -- the gff problem -- ###############
        #since the gff contains all transcripts - at times the variants reported here are from other transcripts
        #I correct for this by removing variants that do not fall within the longest transcript
        
        #changing class of the position in syn_pos_df
        syn_pos_df$pos = as.integer(as.character(syn_pos_df$pos))
        
        #keeping only those positions which are there within the transcript alignment with the outgroup
        syn_pos_df = syn_pos_df[syn_pos_df$pos %in% transcript_alignment_pos,]
        
        ####################################################
        
        
        #writing syn positions dataframe to a table
        write.table(file = paste(path_to_vcf_folder,transcriptID,"_synonymous_positions",sep=""), x = syn_pos_df,
                    quote = FALSE, row.names = FALSE)
        
      }
      
      #writing position information to a file now - only if information exists -- nonsyn
      if(dim(nonsyn_pos_df)[1] > 0){
        
        ############## -- the gff problem -- ###############
        #since the gff contains all transcripts - at times the variants reported here are from other transcripts
        #I correct for this by removing variants that do not fall within the longest transcript
        
        #changing class of the position in nonsyn_pos_df
        nonsyn_pos_df$pos = as.integer(as.character(nonsyn_pos_df$pos))
        
        #keeping only those positions which are there within the transcript alignment with the outgroup
        nonsyn_pos_df = nonsyn_pos_df[nonsyn_pos_df$pos %in% transcript_alignment_pos,]
        
        ####################################################
        
        #writing nonsyn positions dataframe to a table
        write.table(file = paste(path_to_vcf_folder,transcriptID,"_nonsynonymous_positions",sep=""), x = nonsyn_pos_df,
                    quote = FALSE, row.names = FALSE)
        
      }
      
    }#gone over the annotated file and done writing syn and nonsyn positions -- if there are any
    
    }#check if there is indeed an annotation file present
}