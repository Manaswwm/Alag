#center call-point from key.R, contains all the main functions relating to local alag
main = function(geneID_list, backup_file, counter){

  ## Part 1 - make list of transcripts and construct the information on the cds co-ordinates of the transcripts
  # due to limited annotation - our selection of transcript is limited to the largest transcript
  
  #print check
  print("Part 1 - extracting coordinates information on the representative transcript per gene")
  
  #checking if information on the largest transcript and cds info is present
  if(!file.exists(paste(backup_file,"/cds_info", sep=""))){
    
    #callpoint for function to get the largest transcript and the cds info
    cds_info = get_largest_transcript_cds_info(geneID_list)
    
    #saving cds_information in backup directory
    save(cds_info, file = paste(backup_file,"/cds_info", sep = ""))
    
  }
  
  #loading cds_info here
  load(paste(backup_file,"/cds_info",sep = ""))
  
  
  ## Part 2 - accessing the transcript sequences from the reference sequence file
  # and performing alignments to the outgroup - we use local blast for performing the alignment for now
  
  #print check
  print("Part 2 - extracting transcript alignments with the outgroup species")
  
  #callpoint for function to get alignments of the largest transcripts
  #first creating a directory to store all the transcript sequences
  if(!dir.exists(paste(backup_file,"/transcript_sequences", sep=""))){
    dir.create(paste(backup_file,"/transcript_sequences", sep=""), showWarnings = TRUE)
  }
  #second creating a directory to store transcript alignments
  if(!dir.exists(paste(backup_file,"/transcript_alignments", sep=""))){
    dir.create(paste(backup_file,"/transcript_alignments", sep=""), showWarnings = TRUE)
  }
  
  #checking if the information on transcript alignments is present
  if(!file.exists(paste(backup_file,"/transcript_alignments_all", sep=""))){
    
    ## mclapply - give one transcript (one list of df) at a time
    transcript_alignments = mclapply(cds_info, function(x) get_transcript_reciprocal_alignments(x), mc.cores = 1)
    
    #checking if any of the alignments are null
    transcript_alignments = transcript_alignments[!sapply(transcript_alignments, function(x) is.null(x$alignment))]  
    
    #saving transcript alignments in backupdirectory
    save(transcript_alignments, file = paste(backup_file,"/transcript_alignments_all",sep="")) 
    
  }
  
  #loading transcript alignments here
  load(paste(backup_file,"/transcript_alignments_all", sep = ""))
  
  #removing transcript_seq object here
  #rm(transcript_seq)
  
  ###### --- Divergence stats calculation part begins --- ######
  
  
  ## Part 3 - getting information on divergence sites and syn and nonsyn lengths from the transcript alignments
  
  #print check
  print("Part 3 - extracting divergent sites")
  
  #checking if the information on divergent sites and lengths is present
  if(!file.exists(paste(backup_file,"/divergent_sites",sep = ""))){
    
    #callpoint for getting divergent sites
    divergent_sites=mclapply(transcript_alignments,function(x) get_divergent_sites2(x), mc.cores = 1)
    
    #saving divergent sites in backup directory
    save(divergent_sites, file = paste(backup_file,"/divergent_sites", sep = ""))
  }
  
  #loading divergent sites here
  load(paste(backup_file,"/divergent_sites",sep = ""))
  
  #extract df with syn and nonsyn length
  df_syn_nonsyn_lengths = lapply(divergent_sites, function(x) x$lengths)
  
  #convert to a dataframe
  df_syn_nonsyn_lengths = ldply(df_syn_nonsyn_lengths, rbind)
  
  #extract df with list of divergent sites
  divergent_sites = lapply(divergent_sites, function(x) x$variant_list)
  
  #convert to a dataframe
  divergent_sites = ldply(divergent_sites, rbind)
  
  
  ## Part 4 - getting divergence stats from the divergence information collected so far
  
  #print check
  print("Part 4 - calculating divergence stats")
  
  #checking if information on divergence stats is present
  if(!file.exists(paste(backup_file,"/div_per_gene", sep = ""))){
    
    #call point for getting divergence stats per gene
    div_per_gene = get_divergence_stats(divergent_sites, df_syn_nonsyn_lengths)
    
    #saving the divergence stats per gene in backup directory
    save(div_per_gene, file = paste(backup_file,"/div_per_gene", sep = ""))
  }
  
  #loading divergence per gene stats
  load(paste(backup_file,"/div_per_gene", sep=""))
  
  
  # ##### -- Polymorphism stats calculation part begins -- #####
  
  ## Part 5 - getting population-specific frequencies of variants falling within the cds
  
  #print check
  print("Part 5 - extracting population-specific variants")
  
  #creating folder to write the vcf information if the file does not exist already
  if(!dir.exists(paste(backup_file,"/vcf_cds",sep=""))){
         dir.create(paste(backup_file,"/vcf_cds",sep=""))
    }

  #re-loading cds_info
  load(paste(backup_file,"/","cds_info", sep=""))
  
  #keep track of all geneIDs present at this stage
  #geneIDs_with_data=lapply(transcript_alignments, function(x) return(x$id))
  
  #### to account for the presence of characters other than A,G,C and T - I have changed the get_divergent_sites2.R code to drop those IDs
  #### so from this point on - only those IDs that are present in the div_per_gene dataframe will be kept
  geneIDs_with_data = div_per_gene$transcriptID
  
  #removing cds_info of transcripts which do not have an alignment
  cds_info = lapply(cds_info, function(x) {subset(x, unique(unlist(x$ensembl_transcript_id)) %in% geneIDs_with_data)})
  cds_info = cds_info[sapply(cds_info, function(x) dim(x)[1] > 0)]
  
  #getting vcf locally
  if(!file.exists(paste(backup_file,"/","frequencies_cds", sep=""))){
    
    #### the vcf file from Danijela is not split for populations, supplying different population samples here
    #### pop_info = read.delim("/netscratch/dep_tsiantis/grp_laurent/joshi/alag_local_files/vcf_for_Cardamine_Danijela_Azores/imputed_version/pop_files/az_pure_0.78", header = FALSE)
    #### I read this file in key.R; line 111
    
    #formatting this to contain samples in the syntax that is required by bcftools
    pop_info = paste(pop_info$V1, collapse = ",")
    
    #giving a call for calculating the frequencies cds
    frequencies_cds = ldply(mclapply(cds_info, function(x){
      get_frequencies_cds_for_hirsuta(x, paste(backup_file,"/","vcf_cds",sep=""), transcript_alignments, pop_info)
    }, mc.cores = 1), rbind)

    save(frequencies_cds, file =paste(backup_file,"/","frequencies_cds", sep=""))
  }
  
  #loading divergence per gene stats
  load(paste(backup_file,"/frequencies_cds", sep=""))
  
  #giving back the transcript ID
  frequencies_cds = merge(frequencies_cds, gene_conversion_table[,c("ensembl_gene_id", "ensembl_transcript_id")], by.x = "geneID", by.y = "ensembl_gene_id")
  
  ## Part 6 - getting polymorphism stats per gene 
  
  #@print check
  print("Part 6 - calculating polymorphism stats per gene")
  
  #extracting per gene polymorphism information
  if(!file.exists(paste(backup_file,"/","polymorphism_per_gene", sep=""))){
    
    #calculate diversity indices for synonymous variants
    polymorphism_syn_per_gene=ldply(lapply(unlist(geneIDs_with_data), function(x) get_polymorphism_stats_per_gene2(transcriptID = x, allele_frequencies = frequencies_cds, min_freq, type = "syn")), rbind)
    
    #calculate diversity indices for non-synonymous variants
    polymorphism_nonsyn_per_gene=ldply(lapply(unlist(geneIDs_with_data), function(x) get_polymorphism_stats_per_gene2(transcriptID = x, allele_frequencies = frequencies_cds, min_freq, type = "nonsyn")), rbind)
    
    #calculate diversity indices for nonsense variants
    polymorphism_nonsense_per_gene=ldply(lapply(unlist(geneIDs_with_data), function(x) get_polymorphism_stats_per_gene2(transcriptID = x, allele_frequencies = frequencies_cds, min_freq, type = "nonsense")), rbind)
    
    #merge syn and nonsyn results
    polymorphism_per_gene=merge(polymorphism_syn_per_gene, polymorphism_nonsyn_per_gene, by="id")
    
    #merge syn and nonsyn and nonsense results
    polymorphism_per_gene = merge(polymorphism_per_gene, polymorphism_nonsense_per_gene, by = "id")
    
    #reformat geneID to allow merging
    #polymorphism_per_gene$id=gsub(".[[:digit:]]*$","",polymorphism_per_gene$id)
    
    #calculating the pi_nonsyn/pi_syn
    polymorphism_per_gene$pin_pis = polymorphism_per_gene$pi_nonsyn/polymorphism_per_gene$pi_syn
    
    #calculating the pi_nonsense/pi_syn
    polymorphism_per_gene$pinonsense_pis = polymorphism_per_gene$pi_nonsense/polymorphism_per_gene$pi_syn 
    
    #save to backup folder 
    save(polymorphism_per_gene, file =paste(backup_file,"/","polymorphism_per_gene", sep=""))
    
    }
  load(paste(backup_file,"/","polymorphism_per_gene", sep=""))
  

  ## Part 7 - the climax
  #merging the divergence and polymorphism information
  
  #@print check
  print("Part 7 - final step of merging polymorphism and divergence stats")
  
  #merging on transcript ID
  poly_div_gene_merge = merge(polymorphism_per_gene, div_per_gene, by.x = "id", by.y = "transcriptID")
  
  #writing this to a file
  write.table(poly_div_gene_merge, file = paste("poly_div_gene_stats_batch_",counter,".txt", sep=""), row.names = FALSE,
              col.names = TRUE, sep = "\t", quote = FALSE)
}
  
  # ##### -- Divergence stats calculation part concludes -- #####
  # 
  # 
  # ##### -- Polymorphism stats calculation part begins -- #####
  # 
  # 
  # ## Part 5 - subsetting the vcf file according to set of transcript coordinates and population accessions
  # 
  # #print check
  # print("Part 5")
  # 
  # #this part will be inserted once the vcf files are available
  # 
  # 
  # ## Part 6 - annotating the polymorphic variants and subsetting syn and nonsyn variants on their positions
  # 
  # #print check
  # print("Part 6")
  # 
  # #sequential call to annotate_polymorphic_variants_function
  # #no need to use the next variable as all the information is already stored  
  # annotate_variants = mclapply(div_per_gene$transcriptID, function(x) {
  #   annotate_polymorphic_variants(path_to_vcf_folder = paste(backup_file,"vcf_cds/",x,"/",sep=""), transcriptID = x)}, mc.cores = 1)
  # 
  # 
  # #sequential call for getting the synonymous and nonsynonymous positions
  # #no need to use the next variable as all the information is already stored
  # syn_nonsyn_pos_list = mclapply(div_per_gene$transcriptID, function(x) {
  #   get_syn_nonsyn_positions_poly(path_to_vcf_folder = paste(backup_file,"vcf_cds/",x,"/",sep=""), transcriptID = x, transcript_alignments)}, mc.cores = 1)
  # 
  # 
  # ## Part 7 - calculating the frequencies of the polymorphic variants - constructing frequencies_cds dataframe
  # 
  # #print check
  # print("Part 7")
  # 
  # #sequentially giving call to the calculation of frequencies_cds calculation
  # if(!file.exists(paste(backup_file,"/frequencies_cds",sep = ""))){
  #   
  #   #callpoint for getting the frequencies_cds dataframe
  #   frequencies_cds = ldply(mclapply(div_per_gene$transcriptID, function(x) {
  #     get_frequencies_cds(path_to_vcf_folder = paste(backup_file,"vcf_cds/",x,"/",sep=""), transcriptID = x)}, mc.cores = 1), rbind)
  # 
  #   #saving frequencies_cds in backup directory
  #   save(frequencies_cds, file = paste(backup_file,"/frequencies_cds",sep = ""))
  # }
  # 
  # #loading frequencies_cds dataframe here
  # load(paste(backup_file,"/frequencies_cds",sep = ""))
  # 
  # 
  # ## Part 8 - calculating the polymorphism stats per gene from frequencies_cds dataframe
  # 
  # #print check
  # print("Part 8")
  # 
  # #sequentially giving call to the polymorphism stats calculation
  # if(!file.exists(paste(backup_file,"/polymorphism_per_gene",sep = ""))){
  # 
  #   #callpoint for getting the poly syn stats per gene
  #   polymorphism_syn_per_gene = ldply(lapply(unique(frequencies_cds$transcriptID), function(x){
  #     get_polymorphism_stats_per_gene2(geneID = x, allele_frequencies = frequencies_cds, min_freq = min_freq, variant_type = "syn")}))
  #   
  #   #callpoint for getting the poly nonsyn stats per gene
  #   polymorphism_nonsyn_per_gene = ldply(lapply(unique(frequencies_cds$transcriptID), function(x){
  #     get_polymorphism_stats_per_gene2(geneID = x, allele_frequencies = frequencies_cds, min_freq = min_freq, variant_type = "nonsyn")}))
  #   
  #   #merging the syn and nonsyn stats
  #   polymorphism_per_gene = merge(polymorphism_syn_per_gene, polymorphism_nonsyn_per_gene, by = "id")
  #   
  #   #saving poly syn and nonsyn stats per gene to backup directory
  #   save(polymorphism_per_gene, file = paste(backup_file,"/polymorphism_stats_per_gene",sep=""))
  #   
  # }
  # 
  # #loading frequencies_cds dataframe here
  # load(paste(backup_file,"/polymorphism_stats_per_gene",sep = "")) 
  # 
  # ##### -- Polymorphism stats calculation part concludes -- #####
  # 
  # 
  # ## Part 9 - merging the divergence and polymorphism stats now
  # 
  # #print check
  # print("Part 9")
  # 
  # poly_div_per_gene = merge(polymorphism_per_gene, div_per_gene, by.x = "id", by.y = "transcriptID")
  # 
  # #writing the output to a file now
  # write.table(poly_div_per_gene, file = paste("poly_div_per_gene_stats_",counter,".txt",sep = ""), row.names = FALSE,
  #             col.names = TRUE, quote = FALSE, sep = "\t")


