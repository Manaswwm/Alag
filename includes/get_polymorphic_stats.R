#callpoint for calculating polymorphism statistics
#local version of differentiating syn and nonsyn mutations from locally available sequence information
#inspired from get_divergence_stats.R for C. hirsuta

get_polymorphic_stats = function(polymorphic_sites, list_of_regions){
  
  
  #declaring list of codons which are pre-built within BioStrings
  gen_code = as.data.frame(GENETIC_CODE)
  gen_code$codon = rownames(gen_code)
  colnames(gen_code) = c("aa", "codon")
  row.names(gen_code) = NULL
  
  ## Part 1 - first creating transcript level positions from cds positions ##
  cds_pos = c()
  
  #going over all rows sequentially
  for(item in 1:nrow(list_of_regions)){
    
    #listing a single row here
    row = list_of_regions[item,]
    
    #all pos in one row
    pos = unlist(row$start):unlist(row$end)
    
    #joining
    cds_pos = c(cds_pos,pos)
    
  }
  
  #determination of the polymorphism information is strand specific - this is different as compared to the divergence information
  #in divergence we directly compare sequence to sequence and identify the impact of variants
  #in case of polymorphism - I dont have to luxury to work with transcripts directly but rather the ref and alt alleles are reverse complemented
  #determination of the impact of variants is strand specific
  
  ### -- for genes that are on negative strand -- ###
  if(unique(list_of_regions$strand == "-")){

    #first determining the position of the allele within the alignment - as if the strand is positive
    polymorphic_sites$transcript_pos_neg_strand = match(polymorphic_sites$mapped.start, cds_pos)
    
    #next determining the position of the allele when the strand will be reverse complemented
    polymorphic_sites$transcript_pos_neg_strand = (length(cds_pos) - polymorphic_sites$transcript_pos_neg_strand) + 1
    
    #reverse complementing the reference
    polymorphic_sites$ref_allele_reverse_complement = sapply(polymorphic_sites$list_of_divergent_alleles_ref, function(x) as.character(reverseComplement(DNAString(x))))
    
    #reverse complementing the alternate
    polymorphic_sites$alt_allele_reverse_complement = sapply(polymorphic_sites$list_of_divergent_alleles_alt, function(x) as.character(reverseComplement(DNAString(x))))

    #subsetting dataframe to keep only the required information for syn and nonsyn splitting
    syn_nonsyn_split_df = polymorphic_sites[,c("transcript_pos_neg_strand", "ref_allele_reverse_complement", "alt_allele_reverse_complement")]
    colnames(syn_nonsyn_split_df) = c("cds_pos", "ref", "alt")
    
    ### discretizing variants to be either syn or nonsyn on the basis of the reference transcript ###
    #first accessing the transcript sequence
    transcript_seq = as.character(as.character(
      readDNAStringSet(file = paste(backup_file,"/transcript_sequences/",unique(list_of_regions$ensembl_transcript_id),"_sequence.fa",sep=""))))  
    
    #chopping the transcript sequence into codons
    codons_list = regmatches(transcript_seq, gregexpr(".{3}", transcript_seq))[[1]]
    
    ## determining syn and nonsyn locally ##
    
    #dataframe for storing syn and nonsyn information for every variant
    syn_nonsyn_info = data.frame()
    
    #going over all variants instead of all codons
    for(item in syn_nonsyn_split_df$cds_pos){
      
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
      substr(var_cod, var_pos, var_pos) = syn_nonsyn_split_df$alt[syn_nonsyn_split_df$cds_pos == item]
      
      #mutated amino acid
      aa_mut = gen_code$aa[gen_code$codon == var_cod]
      
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
      # 
      # 
      # if(item < (nchar(transcript_seq) - 2) & aa_mut == "*"){
      #   
      #   #attempt 1 - writing this to a seperate file
      #   log_file_name = paste("nonsense_mutation_file.txt")
      #   cat(paste("identified a nonsense mutation for transcript :",unique(list_of_regions$transcriptID), " for position : ",item, "\n", sep = ""),
      #       file = log_file_name, append = TRUE)
      #   
      # }
      
      #binding this information into the big dataframe
      syn_nonsyn_info = rbind(syn_nonsyn_info, df)  
      
    }
    
    #merging with the original dataframe
    syn_nonsyn_split_df = merge(syn_nonsyn_split_df, syn_nonsyn_info, by.x = "cds_pos", "pos")
    
    #merging with the other polymorphic information
    polymorphic_sites = merge(polymorphic_sites, syn_nonsyn_split_df[,c("cds_pos", "type")], by.x = "transcript_pos_neg_strand", by.y = "cds_pos")
    
    #renaming the columns to keep the "true" ref and "reverse complement" ref
    colnames(polymorphic_sites) = c("mapped.start", "mapped.seq_region_name", "mapped.start_true", "snpID_vector", "ref_allele_true", "alt_allele_true",
                                    "geneID", "ref_allele", "alt_allele", "type")
    
    #rearranging the columns
    polymorphic_sites = polymorphic_sites[,c(2,3,5,6,4,1,8,9,7,10)]
    
    #splitting syn and nonsyn
    polymorphic_sites_syn = polymorphic_sites[polymorphic_sites$type == "syn",]
    polymorphic_sites_nonsyn = polymorphic_sites[polymorphic_sites$type == "nonsyn",]
    polymorphic_sites_nonsense = polymorphic_sites[polymorphic_sites$type == "nonsense",]
    
    #returning information collected so far
    return(list(syn = polymorphic_sites_syn, nonsyn = polymorphic_sites_nonsyn, nonsense = polymorphic_sites_nonsense))
    
    
    }
  
  ### -- for genes that are on positive strand -- ###
  
  if(unique(list_of_regions$strand == "+")){
    
    #first determining the position of the allele within the alignment - as if the strand is positive
    polymorphic_sites$transcript_pos_pos_strand = match(polymorphic_sites$mapped.start, cds_pos)
    
    #next determining the position of the allele when the strand will be reverse complemented
    #polymorphic_sites$transcript_pos_neg_strand = (length(cds_pos) - polymorphic_sites$transcript_pos_neg_strand) + 1
    
    #no reverse complementing as this is the positive strand
    polymorphic_sites$ref_allele_reverse_complement = polymorphic_sites$list_of_divergent_alleles_ref
    
    #no reverse complementing as this is the positive strand
    polymorphic_sites$alt_allele_reverse_complement = polymorphic_sites$list_of_divergent_alleles_alt
    
    #subsetting dataframe to keep only the required information for syn and nonsyn splitting
    syn_nonsyn_split_df = polymorphic_sites[,c("transcript_pos_pos_strand", "list_of_divergent_alleles_ref", "list_of_divergent_alleles_alt")]
    colnames(syn_nonsyn_split_df) = c("cds_pos", "ref", "alt")
    
    ### discretizing variants to be either syn or nonsyn on the basis of the reference transcript ###
    #first accessing the transcript sequence
    transcript_seq = as.character(as.character(
      readDNAStringSet(file = paste(backup_file,"/transcript_sequences/",unique(list_of_regions$ensembl_transcript_id),"_sequence.fa",sep=""))))  
    
    #chopping the transcript sequence into codons
    codons_list = regmatches(transcript_seq, gregexpr(".{3}", transcript_seq))[[1]]
    
    ## determining syn and nonsyn locally ##
    
    #dataframe for storing syn and nonsyn information for every variant
    syn_nonsyn_info = data.frame()
    
    #going over all variants instead of all codons
    for(item in syn_nonsyn_split_df$cds_pos){
      
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
      substr(var_cod, var_pos, var_pos) = syn_nonsyn_split_df$alt[syn_nonsyn_split_df$cds_pos == item]
      
      #mutated amino acid
      aa_mut = gen_code$aa[gen_code$codon == var_cod]
      
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
    
    #merging with the original dataframe
    syn_nonsyn_split_df = merge(syn_nonsyn_split_df, syn_nonsyn_info, by.x = "cds_pos", "pos")
    
    #merging with the other polymorphic information
    polymorphic_sites = merge(polymorphic_sites, syn_nonsyn_split_df[,c("cds_pos", "type")], by.x = "transcript_pos_pos_strand", by.y = "cds_pos")
    
    #renaming the columns to keep the "true" ref and "reverse complement" ref
    colnames(polymorphic_sites) = c("mapped.start", "mapped.seq_region_name", "mapped.start_true", "snpID_vector", "ref_allele_true", "alt_allele_true",
                                    "geneID", "ref_allele", "alt_allele", "type")
    
    #rearranging the columns
    polymorphic_sites = polymorphic_sites[,c(2,3,5,6,4,1,8,9,7,10)]
    
    #splitting syn and nonsyn
    polymorphic_sites_syn = polymorphic_sites[polymorphic_sites$type == "syn",]
    polymorphic_sites_nonsyn = polymorphic_sites[polymorphic_sites$type == "nonsyn",]
    polymorphic_sites_nonsense = polymorphic_sites[polymorphic_sites$type == "nonsense",]
    
    
    #returning information collected so far
    return(list(syn = polymorphic_sites_syn, nonsyn = polymorphic_sites_nonsyn, nonsense = polymorphic_sites_nonsense))
    
  }
  
  else{
    
    print(paste("Invalid strand - check ID : ",unique(list_of_regions$transcriptID),sep = ""))
    return(NULL)
    
  }
}