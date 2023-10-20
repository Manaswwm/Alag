#callpoint for calculating divergence statistics 
get_divergence_stats = function(divergent_sites, df_syn_nonsyn_lengths){
  
  #declaring list of codons which are pre-built within BioStrings
  gen_code = as.data.frame(GENETIC_CODE)
  gen_code$codon = rownames(gen_code)
  colnames(gen_code) = c("aa", "codon")
  row.names(gen_code) = NULL
  
  #empty dataframe to fill in all the divergent stats
  div_stats = data.frame()

  ## sequentially going over every geneID
  for(id in unique(divergent_sites$id)){
    
    #@print check
    print(paste("Calculating divergence stats for transcript id : ",id, sep = ""))
    
    #first accessing the transcript sequence
    transcript_seq = as.character(as.character(
      readDNAStringSet(file = paste(backup_file,"/transcript_sequences/",id,"_sequence.fa",sep=""))))
    
    #taking subset of the divergent sites dataframe for this transcript
    divergent_sites_subset = divergent_sites[divergent_sites$id == id,]
    
    #chopping the transcript sequence into codons
    codons_list = regmatches(transcript_seq, gregexpr(".{3}", transcript_seq))[[1]]
    
    ## counts first - taking count of syn and nonsyn variants ##
    
    #dataframe for storing syn and nonsyn information for every variant
    syn_nonsyn_info = data.frame()
    
    #going over all variants instead of all codons
    for(item in divergent_sites_subset$cds_pos){
      
      #if variant occurs on a 3rd position then consider the nth codon
      #if variant does not occur on 3rd position then consider the n+1th codon
      if(item%%3 == 0){cod_pos = (item%/%3)}else{cod_pos = (item%/%3) + 1}
      
      #accessing the codon under question
      var_cod = codons_list[cod_pos]
      
      #wild type amino acid
      aa_wild = as.character(gen_code$aa[gen_code$codon == var_cod])
      
      #the position of the variant within the codon
      var_pos = item%%3
      
      ### bug - if the codon is at the third position then it gives position as 0 on modulus - changing this to 3 - MJ11122022 ###
      ### usually third position substitutions are synonymous - but some could be nonsyn - our test transcripts did not seem to have such mutations before ###
      ### if the variant falls in position 1 or 2 then I will be able to catch them any ways ###
      if(var_pos == 0){var_pos = 3}
      
      #substituting the wild type nucleotide with the variant nucleotide within the codon
      substr(var_cod, var_pos, var_pos) = divergent_sites_subset$list_of_divergent_alleles_alt[divergent_sites_subset$cds_pos == item]
      
      #mutated amino acid
      aa_mut = gen_code$aa[gen_code$codon == var_cod]
      
      # #checkin if the wild type and mutated amino acid are the same - if yes then syn, else nonsyn
      # if(aa_wild == aa_mut){
      #   df = data.frame(pos = item, type = "syn")
      #   }else if(!aa_wild == aa_mut){df = data.frame(pos = item, type = "nonsyn")}
      
      #checking if the wild type and mutated amino acid are the same - if yes then syn, else nonsyn
      #additional layer of information - checking for premature stop codons - more on these here - https://en.wikipedia.org/wiki/Nonsense_mutation
      #in terms of technicality - I need to identify if there are any mutations on any of the codons barring from the last codon
      #that change the codon to either TAA (UAA), TGA (UGA) or TAG (UAG)
      #checking if the mutated amino acid is * (termination codon) --- MJ 07062022
      #nonsyn mutations are further classifed into nonsyn and nonsense
      if(aa_wild == aa_mut){#a plain syn mutation
        df = data.frame(pos = item, type = "syn")
      }else if(!aa_wild == aa_mut){
        #if the mutated amino acid is nonsyn but not * (a termination codon)
        if(!aa_mut == "*"){ # a plain nonsyn mutation
          df = data.frame(pos = item, type = "nonsyn")
        }else if(item < (nchar(transcript_seq) - 2) & aa_mut == "*"){ # a nonsyn mutation which is a nonsense mutation
          df = data.frame(pos = item, type = "nonsense")
        }      
      }
      
      
      #binding this information into the big dataframe
      syn_nonsyn_info = rbind(syn_nonsyn_info, df)  
      
    }
    
    ## count ends ##
    
    ## taking ratios from the counts of variants and length of the variants ##
    
    #making a small jump here - taking the geneID to work with the lengths dataframe
    geneID = as.character(unique(divergent_sites_subset$geneID))
    
    #accessing lengths here
    length_df_sub = df_syn_nonsyn_lengths[df_syn_nonsyn_lengths$geneID == geneID,]
    
    #constructing the final dataframe for div stats per gene
    final = data.frame(div_syn_count = (dim(syn_nonsyn_info[syn_nonsyn_info$type == "syn",])[1]) ,
                       div_nonsyn_count = (dim(syn_nonsyn_info[syn_nonsyn_info$type == "nonsyn",])[1]),
                       div_nonsense_count = (dim(syn_nonsyn_info[syn_nonsyn_info$type == "nonsense",])[[1]]))
    
    #making ratios
    final$div_syn = final$div_syn_count/length_df_sub$synonymous_length
    final$div_nonsyn = final$div_nonsyn_count/length_df_sub$nonsynonymous_length
    final$div_nonsense = final$div_nonsense_count/length_df_sub$nonsynonymous_length
    
    #making statistics
    final$KnKs = final$div_nonsyn/final$div_syn
    final$KnonsenseKs = final$div_nonsense/final$div_syn
    
    #giving back information on geneID and transcriptID
    final$geneID = geneID
    final$transcriptID = id
    
    #binding information with the big dataframe
    div_stats = rbind(div_stats, final)
    
  }
  
  #returning the dataframe containing the divergence stats
  return(div_stats)
  
}
