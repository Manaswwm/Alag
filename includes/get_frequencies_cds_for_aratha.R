#callpoint for getting the population specific subset of the vcf file

#sequence of events - inputs - maf, accessions path, path where stuff should be stored, cds_info (one transcript at a time)
# 1. preparing the tabix/bcftools command for system
# 2. getting the vcf file and storing it 
# 3. giving this output to vcf tools to calculate the frequency
# 4. formatting the output of vcftools - removing row 1, adding the column names and calculating the maf
# 5. getting the ancestral state

## storing the vcf files locally - the requests to ftp using tabix/bcftools gives random errors
## this is likely caused due to unstable connectivity on my side or the ftp side
get_frequencies_cds_for_hirsuta = function(list_of_regions, vcf_folder, transcript_alignments, pop_info){
  
  #### -- pre-processing section -- ####
  # 
  # list_of_regions = list_of_regions[[1]]
  
  #declaring transcriptID here
  transcriptID = as.character(unique(unlist(list_of_regions$ensembl_transcript_id)))
  
  #### for thaliana - since there is a single vcf file - I will hardcode the path here directly ####
  bcftools_vcf_file = paste("bcftools view ", vcf_file_path, " ", sep = "")
  
  #have to add the command for subsetting on individuals
  bcftools_polymorphic = paste(" -c 1 ")

  #the imputed version of the file is single for all populations, I have to specify populations for every pop-run
  bcftools_pop_info = paste("-s ", pop_info, " ",sep = "")
  
  #only for C. hirsuta -- I have to rename chromosomes from Chr1 to 1 ..
  #list_of_regions$seqid = as.character(gsub("Chr", "", list_of_regions$seqid))
  
  #common sections of commands for every chromosome - remaining part of the command will be 
  #chr dependent
  
  #-- hardcoded vcf file
  #-- bcftools_main = paste("/opt/share/software/bin/bcftools view /netscratch/dep_tsiantis/grp_laurent/joshi/alag_local_files/vcf_for_Cardamine_Danijela_Azores/imputed_version/chirsuta_imputed_version_all_pops_reheader.vcf.gz ")
  #-- CHANGED TO:
  #-- where the vcf file is set in key.R; line 108
  # bcftools_main = paste("/opt/share/software/bin/bcftools view ", vcf_file_path)
  # 
  # 
  # #have to add the command for subsetting on individuals
  # bcftools_polymorphic = paste(" -c 1 ")
  # 
  # #the imputed version of the file is single for all populations, I have to specify populations for every pop-run
  # bcftools_pop_info = paste("-s ", pop_info, " ",sep = "")
  # 
  #print check - which id am I currently working on
  print(paste("Calculating vcf for id : ",as.character(unique(list_of_regions$ensembl_gene_id))),sep="")
  
  #creating an empty dataframe to fill in all the information
  polymorphism_info = data.frame()
  vcf_data = data.frame()
  
  #subsetting and changing the datatype of the list of regions here
  list_of_regions=list_of_regions[,c("seqid","strand", "start","end","ensembl_gene_id", "ensembl_transcript_id")]
  list_of_regions$seqid = as.character(list_of_regions$seqid)
  
  #### -- pre-processing section concludes -- ####
  
  #### -- going over all cds intervals -- ####
  
  #going over every row in the list of regions - coordinates
  for (item in 1:nrow(list_of_regions)){
    
    #going over every row sequentially
    row = list_of_regions[item,]
    
    #defining the region
    #for humans - adding chr to the chromosome - different for different species
    bcftools_region = paste("-r ",row$seqid, ":", row$start, "-", row$end, sep = "")
    
    #print check
    print(bcftools_region)
    
    #print check
    print(paste(bcftools_vcf_file, bcftools_polymorphic, bcftools_pop_info, bcftools_region, sep=""))
    
    #merging the commands
    bcftools_command = system(paste(bcftools_vcf_file, bcftools_pop_info, bcftools_region, sep=""), intern = TRUE)
    
    #checking if there is any output to process - if yes then merging it with dataframe
    if(length(bcftools_command)>18){
      bcftools_command = bcftools_command[-grep(pattern = "^##", x = bcftools_command)]
      bcftools_command[1] = str_remove(string = bcftools_command[1], pattern = "#")
      bcftools_command = read.table(text = bcftools_command, sep="\t", header = TRUE, colClasses = "character")
      vcf_data = rbind(vcf_data, bcftools_command)
    }
    #}
  } # gone over all rows of the list of regions
  
  #### -- going over all cds intervals concludes -- ####
  
  #### -- vcf processing -- ####
  
  #checking if I have at all collected data in the vcf-file
  if(length(vcf_data)>0){
    path_to_vcf = paste(vcf_folder,"/",transcriptID,".vcf",sep="")
    colnames(vcf_data)[1] = "#CHROM"
    write.table(x = vcf_data, file = path_to_vcf,quote = FALSE, row.names = FALSE, sep ="\t")
    vcf_command = system(paste("vcftools --vcf ",path_to_vcf, " --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf ",maf," --recode --stdout | cut -f 1-5", sep=""), intern = TRUE)
    
    #in here because I have indeed collected some data
    #processing the vcf command now segregating into syn and nonsyn - if the output exists
    if(length(grep("#",vcf_command, invert=T))>0){
      
      #reading the vcf command output
      polymorphic_sites = read.table(text = vcf_command, skip = 1, sep = "\t", colClasses = c("character","numeric","character","character","character"))
      colnames(polymorphic_sites) = c("mapped.seq_region_name","mapped.start","snpID_vector","list_of_divergent_alleles_ref","list_of_divergent_alleles_alt")
      
      #bringing back the geneIDs
      geneID = unique(list_of_regions$ensembl_gene_id)
      polymorphic_sites$geneID=rep(geneID,nrow(polymorphic_sites))
      
      #vep for all the variants in the vcf file
      polymorphic_sites_syn_and_nonsyn=get_polymorphic_stats(polymorphic_sites, list_of_regions)
      
      #### -- break point into syn and nonsyn and nonsense -- ####
      
      ##### first for synonymous sites #####
      #@print check
      print("Calculating polymorphism syn stats ..")
      #extracting syn data
      polymorphic_synonymous_sites=polymorphic_sites_syn_and_nonsyn$syn
      
      #writing the syn data to file
      write.table(t(do.call(rbind, polymorphic_synonymous_sites[,c(1,2)])), quote=F, row.names = F, col.names = F, file=paste(vcf_folder,"/",geneID,"_positions_synonymous_polymorphims",sep=""))
      
      #calculating frequency here
      calculate_freq_command=paste("vcftools --vcf ",path_to_vcf, " --min-alleles 2 --max-alleles 2 --positions ", vcf_folder,"/",geneID, "_positions_synonymous_polymorphims --freq2 --stdout", sep="")
      command_return = system(command = calculate_freq_command, intern = T)
      
      #processing the syn information - frequency, ancestral, ingroup/outgroup
      if(length(command_return)>1){
        
        #format vcftools output
        poly_syn_freq = read.table(text=command_return, skip = 1)
        
        #add names to columns
        colnames(poly_syn_freq)=c("chr","mapped.start_true","alleles","nchr","freq_ref","freq_alt")
        
        #get maf
        poly_syn_freq$maf=apply(poly_syn_freq[,c(5,6)],1, min)
        
        #extract information about positions and allelic states
        #vcf_alleles=polymorphic_sites[,c(2,4,5)]
        
        # #add column names
        # colnames(vcf_alleles)=c("start","ref","alt")
        
        #merge with frequency information
        #poly_syn_freq=merge(poly_syn_freq, vcf_alleles, by.x = "mapped.start_true", by.y = "mapped.start")
        
        ## my headache with the reverse allele ##
        #adding the "true" ref and alt along with the "reverse complement" ref and alt
        poly_syn_freq = merge(poly_syn_freq[,c("mapped.start_true", "alleles", "nchr", "freq_ref", "freq_alt", "maf")], polymorphic_synonymous_sites, by = "mapped.start_true")
        
        ## since C. hirsuta and C. oligospermaalignments are not available on Ensembl -- I polarize the data here from local alignments ##
        #checking if the dataframe is empty and if the transcript alignment from divergence has coverage for polymorphic variants
        if(!dim(poly_syn_freq)[1] == 0){
          
          #creating empty dataframe to store the ingroup and outgroup information
          ancestral_state_df = data.frame()
          
          #getting the ingroup and outgroup of variants -- using the transcript alignments here
          for(item in 1:nrow(poly_syn_freq)){
            
            #accessing individual row here
            row = poly_syn_freq[item,]
            
            #accessing the transcript alignment belonging to this transcript
            transcript_info = transcript_alignments[sapply(transcript_alignments, function(x) x$id == transcriptID)]
            
            ## -- important -- checking if the alignment from divergence is in range with the polymorphic sites
            ## -- specifically - in some cases the alignments from divergence covers only subset of transcript
            ## -- however polymorphic variants can fall outside of this alignments but within the cds of gene
            #### -- checking the absolute position of the variant within the alignment - not the transcript-specific with genomic-specific position ####
            if(row$mapped.start_true %in% colnames(transcript_info[[1]]$alignment)){
              
              #outgroup information -- divergence species - the second entry is always the outgroup
              outgroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][2]
              
              #ingroup information - species of interest - the first entry is always the ingroup
              ingroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][1]
            }else{
              
              outgroup = "-"
              ingroup = "-"
              
            }
            
            ##$$$ old polarization code $$##
            # if(length(colnames(transcript_info[[1]]$alignment)) >= row$mapped.start){ 
            # 
            #   #outgroup information -- divergence species
            #   outgroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][2]
            #   
            #   #ingroup information - species of interest
            #   ingroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][1]
            # }else{
            #   
            #   outgroup = "-"
            #   ingroup = "-"
            #   
            # }
            
            #storing start position -- to merge later
            start = row$mapped.start_true
            
            #storing df
            df = data.frame(start = start, ingroup = ingroup, outgroup = outgroup)
            
            #merging with df
            ancestral_state_df = rbind(ancestral_state_df, df)
            
          }
          
          #merging ancestral state with original dataframe
          poly_syn_freq = merge(poly_syn_freq, ancestral_state_df, by.x = "mapped.start_true", by.y = "start")
          
          #subsetting this dataframe to remove positions for which there is a gap in either in or outgroup
          poly_syn_freq = poly_syn_freq[!(poly_syn_freq$ingroup == "-" | poly_syn_freq$outgroup == "-"),]
          
          ## IMPORTANT -- not all the pop-specific variants have ancestral state -- maybe because this is done on local alignments ##
          #removing positions for which I am not able to obtain an ancestral state
          poly_syn_freq = poly_syn_freq[!is.na(poly_syn_freq$ingroup),]
        }
        
        
        #here i obtain outgroup information from the alignment using chr/pos information
        #prepare necessary information
        #variant_coordinates=poly_syn_freq[,c(2,1,1)]
        
        #get outgroup infromation
        #outgroup_info=ldply(apply(variant_coordinates,1, function(x)get_ancestral_state(x)), rbind)
        
        #merge outgroup information with ingroup frequency information
        #poly_syn_freq=merge(poly_syn_freq,outgroup_info)
        
        if(dim(poly_syn_freq)[1]>0)
        {
          #polarize data
          #how much data am i loosing because it not polarizable?
          poly_syn_freq=ldply(apply(poly_syn_freq, 1, function(x) polarize_data(x)), rbind)
          
          #checking if I lose all the data
          if(dim(poly_syn_freq)[2] > 0){
            
            # #removing the first column
            poly_syn_freq = subset(poly_syn_freq, select = -c(.id))
            
            #setting the formats right again
            poly_syn_freq$freq_der=as.numeric(as.character(poly_syn_freq$freq_der))
            poly_syn_freq$freq_anc=as.numeric(as.character(poly_syn_freq$freq_anc))
            
            #adding new variable indicating whether SNP is syn or nonsyn
            poly_syn_freq$type=rep("syn",dim(poly_syn_freq)[1])
            
          }else{poly_syn_freq = data.frame()}
        }else{poly_syn_freq=data.frame()}
      }else{poly_syn_freq=data.frame()}
      
      ##### second for non-synonymous sites #####
      print("Calculating polymorphism non-syn stats ..")
      #extracting nonsyn data
      polymorphic_nonsynonymous_sites=polymorphic_sites_syn_and_nonsyn$nonsyn
      
      #writing nonsyn info to a file
      write.table(t(do.call(rbind, polymorphic_nonsynonymous_sites[,c(1,2)])), quote=F, row.names = F, col.names = F, file=paste(vcf_folder,"/",geneID,"_positions_nonsynonymous_polymorphims",sep=""))
      
      #format vcftools command
      calculate_freq_command=paste("vcftools --vcf ",path_to_vcf, " --min-alleles 2 --max-alleles 2 --positions ", vcf_folder,"/",geneID, "_positions_nonsynonymous_polymorphims --max-missing 0.8 --freq2 --stdout", sep="")
      
      #call vcftools
      command_return = system(command = calculate_freq_command, intern = T)
      
      #if variants are found in this region
      if(length(command_return)>1){
        
        #format vcftools output
        poly_nonsyn_freq = read.table(text=command_return, skip = 1)
        
        #add names to columns
        colnames(poly_nonsyn_freq)=c("chr","mapped.start_true","alleles","nchr","freq_ref","freq_alt")
        
        #get maf
        poly_nonsyn_freq$maf=apply(poly_nonsyn_freq[,c(5,6)],1, min)
        
        #extract information about positions and allelic states
        #vcf_alleles=polymorphic_sites[,c(2,3,4)]
        
        # #add column names
        # colnames(vcf_alleles)=c("start","ref","alt")
        
        #merge with frequency information
        #poly_nonsyn_freq=merge(poly_nonsyn_freq, vcf_alleles)
        
        ## my headache with the reverse allele ##
        #adding the "true" ref and alt along with the "reverse complement" ref and alt
        poly_nonsyn_freq = merge(poly_nonsyn_freq[,c("mapped.start_true", "alleles", "nchr", "freq_ref", "freq_alt", "maf")], polymorphic_nonsynonymous_sites, by = "mapped.start_true")
        
        #here i obtain outgroup information from the alignment using chr/pos information
        #prepare necessary information
        ## since D. simulans to D. melanogaster alignments are not available on Ensembl -- I polarize the data here from local alignments ##
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
            
            ## -- important -- checking if the alignment from divergence is in range with the polymorphic sites
            ## -- specifically - in some cases the alignments from divergence covers only subset of transcript
            ## -- however polymorphic variants can fall outside of this alignments but within the cds of gene
            #### -- checking the absolute position of the variant within the alignment - not the transcript-specific with genomic-specific position ####
            if(row$mapped.start_true %in% colnames(transcript_info[[1]]$alignment)){
              
              #outgroup information -- divergence species - the second entry is always the outgroup
              outgroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][2]
              
              #ingroup information - species of interest - the first entry is always the ingroup
              ingroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][1]
            }else{
              
              outgroup = "-"
              ingroup = "-"
              
            }
            
            ##$$$ old polarization code $$$##
            # if(length(colnames(transcript_info[[1]]$alignment)) >= row$mapped.start){ 
            #   
            #   # #outgroup information -- divergence species
            #   # outgroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][2]
            #   # 
            #   # #ingroup information - species of interest
            #   # ingroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][1]
            # }else{
            #   
            #   outgroup = "-"
            #   ingroup = "-"
            #   
            # }
            
            #storing start position -- to merge later
            start = row$mapped.start_true
            
            #storing df
            df = data.frame(start = start, ingroup = ingroup, outgroup = outgroup)
            
            #merging with df
            ancestral_state_df = rbind(ancestral_state_df, df)
            
          }
          
          #merging ancestral state with original dataframe
          poly_nonsyn_freq = merge(poly_nonsyn_freq, ancestral_state_df,by.x = "mapped.start_true", by.y = "start")
          
          #subsetting this dataframe to remove positions for which there is a gap in either in or outgroup
          poly_nonsyn_freq = poly_nonsyn_freq[!(poly_nonsyn_freq$ingroup == "-" | poly_nonsyn_freq$outgroup == "-"),]
          
          ## IMPORTANT -- not all the pop-specific variants have ancestral state -- maybe because this is done on local alignments ##
          #removing positions for which I am not able to obtain an ancestral state
          poly_nonsyn_freq = poly_nonsyn_freq[!is.na(poly_nonsyn_freq$ingroup),]
        }
        
        if(dim(poly_nonsyn_freq)[1]>0)
        {
          #polarize data
          #how much data am i loosing because it not polarizable?
          poly_nonsyn_freq=ldply(apply(poly_nonsyn_freq, 1, function(x) polarize_data(x)), rbind)
          
          #checking if I lose all the data in polarization
          if(dim(poly_nonsyn_freq)[2] > 0){
            
            #removing the first column
            poly_nonsyn_freq = subset(poly_nonsyn_freq, select = -c(.id))
            
            #setting the formats right again
            poly_nonsyn_freq$freq_der=as.numeric(as.character(poly_nonsyn_freq$freq_der))
            poly_nonsyn_freq$freq_anc=as.numeric(as.character(poly_nonsyn_freq$freq_anc))
            
            #adding new variable indicating whether SNP is syn or nonsyn
            poly_nonsyn_freq$type=rep("nonsyn",dim(poly_nonsyn_freq)[1])
            
          }else{poly_nonsyn_freq = data.frame()}
        }else{poly_nonsyn_freq=data.frame()}
      }else{poly_nonsyn_freq=data.frame()} 
      
      
      ##### third for nonsense mutations (if there exist any per gene - they are few) #####
      ###first for synonymous sites####
      #@print check
      print("Calculating polymorphism nonsense stats ..")
      #extracting syn data
      polymorphic_nonsense_sites=polymorphic_sites_syn_and_nonsyn$nonsense
      
      #writing the syn data to file
      write.table(t(do.call(rbind, polymorphic_nonsense_sites[,c(1,2)])), quote=F, row.names = F, col.names = F, file=paste(vcf_folder,"/",geneID,"_positions_nonsense_polymorphims",sep=""))
      
      #calculating frequency here
      calculate_freq_command=paste("vcftools --vcf ",path_to_vcf, " --min-alleles 2 --max-alleles 2 --positions ", vcf_folder,"/",geneID, "_positions_nonsense_polymorphims --freq2 --stdout", sep="")
      command_return = system(command = calculate_freq_command, intern = T)
      
      #processing the syn information - frequency, ancestral, ingroup/outgroup
      if(length(command_return)>1){
        
        #format vcftools output
        poly_nonsense_freq = read.table(text=command_return, skip = 1)
        
        #add names to columns
        colnames(poly_nonsense_freq)=c("chr","mapped.start_true","alleles","nchr","freq_ref","freq_alt")
        
        #get maf
        poly_nonsense_freq$maf=apply(poly_nonsense_freq[,c(5,6)],1, min)
        
        #extract information about positions and allelic states
        #vcf_alleles=polymorphic_sites[,c(2,4,5)]
        
        # #add column names
        # colnames(vcf_alleles)=c("start","ref","alt")
        
        #merge with frequency information
        #poly_syn_freq=merge(poly_syn_freq, vcf_alleles, by.x = "mapped.start_true", by.y = "mapped.start")
        
        ## my headache with the reverse allele ##
        #adding the "true" ref and alt along with the "reverse complement" ref and alt
        poly_nonsense_freq = merge(poly_nonsense_freq[,c("mapped.start_true", "alleles", "nchr", "freq_ref", "freq_alt", "maf")], polymorphic_nonsense_sites, by = "mapped.start_true")
        
        ## since C. hirsuta and C. oligospermaalignments are not available on Ensembl -- I polarize the data here from local alignments ##
        #checking if the dataframe is empty and if the transcript alignment from divergence has coverage for polymorphic variants
        if(!dim(poly_nonsense_freq)[1] == 0){
          
          #creating empty dataframe to store the ingroup and outgroup information
          ancestral_state_df = data.frame()
          
          #getting the ingroup and outgroup of variants -- using the transcript alignments here
          for(item in 1:nrow(poly_nonsense_freq)){
            
            #accessing individual row here
            row = poly_nonsense_freq[item,]
            
            #accessing the transcript alignment belonging to this transcript
            transcript_info = transcript_alignments[sapply(transcript_alignments, function(x) x$id == transcriptID)]
            
            ## -- important -- checking if the alignment from divergence is in range with the polymorphic sites
            ## -- specifically - in some cases the alignments from divergence covers only subset of transcript
            ## -- however polymorphic variants can fall outside of this alignments but within the cds of gene
            #### -- checking the absolute position of the variant within the alignment - not the transcript-specific with genomic-specific position ####
            if(row$mapped.start_true %in% colnames(transcript_info[[1]]$alignment)){
              
              #outgroup information -- divergence species - the second entry is always the outgroup
              outgroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][2]
              
              #ingroup information - species of interest - the first entry is always the ingroup
              ingroup = transcript_info[[1]]$alignment[,as.character(row$mapped.start_true)][1]
            }else{
              
              outgroup = "-"
              ingroup = "-"
              
            }
            
            ##$$$ old polarization code $$##
            # if(length(colnames(transcript_info[[1]]$alignment)) >= row$mapped.start){ 
            #   
            #   #outgroup information -- divergence species
            #   outgroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][2]
            #   
            #   #ingroup information - species of interest
            #   ingroup = transcript_info[[1]]$alignment[,colnames(transcript_info[[1]]$alignment)[row$mapped.start]][1]
            # }else{
            #   
            #   outgroup = "-"
            #   ingroup = "-"
            #   
            # }
            
            
            #storing start position -- to merge later
            start = row$mapped.start_true
            
            #storing df
            df = data.frame(start = start, ingroup = ingroup, outgroup = outgroup)
            
            #merging with df
            ancestral_state_df = rbind(ancestral_state_df, df)
            
          }
          
          #merging ancestral state with original dataframe
          poly_nonsense_freq = merge(poly_nonsense_freq, ancestral_state_df, by.x = "mapped.start_true", by.y = "start")
          
          #subsetting this dataframe to remove positions for which there is a gap in either in or outgroup
          poly_nonsense_freq = poly_nonsense_freq[!(poly_nonsense_freq$ingroup == "-" | poly_nonsense_freq$outgroup == "-"),]
          
          ## IMPORTANT -- not all the pop-specific variants have ancestral state -- maybe because this is done on local alignments ##
          #removing positions for which I am not able to obtain an ancestral state
          poly_nonsense_freq = poly_nonsense_freq[!is.na(poly_nonsense_freq$ingroup),]
        }
        
        
        #here i obtain outgroup information from the alignment using chr/pos information
        #prepare necessary information
        #variant_coordinates=poly_syn_freq[,c(2,1,1)]
        
        #get outgroup infromation
        #outgroup_info=ldply(apply(variant_coordinates,1, function(x)get_ancestral_state(x)), rbind)
        
        #merge outgroup information with ingroup frequency information
        #poly_syn_freq=merge(poly_syn_freq,outgroup_info)
        
        if(dim(poly_nonsense_freq)[1]>0)
        {
          #polarize data
          #how much data am i loosing because it not polarizable?
          poly_nonsense_freq=ldply(apply(poly_nonsense_freq, 1, function(x) polarize_data(x)), rbind)
          
          #checking if I lose all the data
          if(dim(poly_nonsense_freq)[2] > 0){
            
            # #removing the first column
            poly_nonsense_freq = subset(poly_nonsense_freq, select = -c(.id))
            
            #setting the formats right again
            poly_nonsense_freq$freq_der=as.numeric(as.character(poly_nonsense_freq$freq_der))
            poly_nonsense_freq$freq_anc=as.numeric(as.character(poly_nonsense_freq$freq_anc))
            
            #adding new variable indicating whether SNP is syn or nonsyn
            poly_nonsense_freq$type=rep("nonsense",dim(poly_nonsense_freq)[1])
            
          }else{poly_nonsense_freq = data.frame()}
        }else{poly_nonsense_freq=data.frame()}
      }else{poly_nonsense_freq=data.frame()}
      
      
      ##### compilation #####
      
      #sanity check if any information on the nonsyn and syn variants in present at all
      if(dim(poly_syn_freq)[1]+dim(poly_nonsyn_freq)[1]+dim(poly_nonsense_freq)[1] > 0){
        
        #concatenating syn and nonsyn dataframes
        poly_syn_nonsyn_freq=rbind(poly_syn_freq, poly_nonsyn_freq, poly_nonsense_freq)
        
        #adding back geneID
        #i have to add the geneID here 
        poly_syn_nonsyn_freq$geneID=rep(rep(geneID),dim(poly_syn_nonsyn_freq)[1])
        
        ##formatting the vcf data table
        #changing datatypes
        #colnames(vcf_data)[1] ="CHROM"
        #vcf_data$CHROM = as.character(vcf_data$CHROM)
        #vcf_data$POS = as.character(vcf_data$POS)
        
        #changing the first column name for clarity
        colnames(poly_syn_nonsyn_freq)[1] = "start"
        
        #giving strand for reference
        poly_syn_nonsyn_freq$strand = unique(list_of_regions$strand)
        
        #giving transcriptID
        poly_syn_nonsyn_freq$transcriptID = unique(list_of_regions$transcriptID)
        
        #@ print check 
        print(poly_syn_nonsyn_freq)
        
      }else{poly_syn_nonsyn_freq = data.frame()}
    }else{poly_syn_nonsyn_freq = data.frame()}
  }else{poly_syn_nonsyn_freq = data.frame()}
  
  #return statement
  return(poly_syn_nonsyn_freq)
  
}