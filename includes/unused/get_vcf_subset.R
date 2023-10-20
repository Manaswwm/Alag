get_vcf_subset=function(list_of_regions, path_to_strains, vcf_folder, maf){
  
  #this is very important. other wise this function seems to be requesting to much
  Sys.sleep(1)
  
  #store all accessions' names
  accessions=read.table(file=path_to_strains)
  
  #store gene id
  transcriptID=unlist(list_of_regions[1,"transcriptID"])
  #convert to vector of formatted factors
  list_of_regions=list_of_regions[,c("seqid","start","end","transcriptID")]
  list_of_regions_formatted=as.factor(gsub(" ","",apply(list_of_regions,1, format_gene_coordinate)))

  #preparing the URL (for GET)
  #server="http://tools.1001genomes.org"
  #ext="/api/v1/vcfsubset"
  #formatted_partial_url_strains=paste("/strains/",toString(accessions$V1), sep="")
  #formatted_partial_url_regions=paste("/regions/",gsub("-","..",toString(list_of_regions)), sep="")
  #additional_options="/type/fullgenome/format/vcf"
  #final_url=gsub(" ","",paste(server,ext,formatted_partial_url_strains, formatted_partial_url_regions, additional_options, sep=""))
  
  #preparing the URL and body (for POST)
  server="http://tools.1001genomes.org/api/v1/vcfsubset/"
  my_body=list(strains=paste(accessions$V1, collapse=","), regions=paste(list_of_regions_formatted, collapse=","), type="fullgenome", format="vcf")

  
  #POSTING
  r=POST(url = server, body=my_body, encode = "multipart", verbose(), config = httr_config)

  if(r$status==200){
    

    #get vcf from content
    vcf=content(r)
    
    #path to vcf file (need to get the gene id back
    path_to_vcf=paste(vcf_folder,"/",geneID,".vcf", sep="")
    #write vcf to file for further analysis with vcftools (extraction of allele frequencies). could be done in R (improve!)
    cat(vcf, file=path_to_vcf)
    #Exclude indels and non-biallelic SNPs and get list and nucleotide state of remaining SNPs
    #print("filter")
    filtering_command=paste("vcftools --vcf ",path_to_vcf, " --remove-indels --min-alleles 2 --max-alleles 2 --max-missing 0.8 --maf ",maf," --recode --stdout | cut -f 1-5", sep="")
    command_return = system(command = filtering_command, intern = T)
    #print("command return")
    #if number of lines (excluding comments which always start with "#") is larger than 0
    if(length(grep("#",command_return, invert=T))>0){
      
      #print("inside if")
      polymorphic_sites = read.table(text = command_return, colClasses = c("numeric","numeric","character","character","character"))
      rm(command_return)
      
      #Rename columns. give the same name as used for the divergence data (see get_divergent_sites2.R)
      #print("changing col names")
      colnames(polymorphic_sites)<-c("mapped.seq_region_name","mapped.start","snpID_vector","list_of_divergent_alleles_ref","list_of_divergent_alleles_alt")
      
      #add back geneID
      #print("adding bacl gene IDs")
      polymorphic_sites$geneID=rep(geneID,nrow(polymorphic_sites))

      #splitting and storing syn and nonsyn variants
      #print("splitting syn and nonsyn")
      polymorphic_sites_syn_and_nonsyn=split_syn_nonsyn_sites_poly(polymorphic_sites)
      #print("done splitting syn and nonsyn")
      
      #UPDATING SYNONYMOUS INFORMATION IF VARIANTS ARE FOUND
      #extract synonymous variants
      polymorphic_synonymous_sites=polymorphic_sites_syn_and_nonsyn$syn
      #print("extracted syn variants")
      #.write to external file for vcftools
      write.table(t(do.call(rbind, polymorphic_synonymous_sites[,c(1,2)])), quote=F, row.names = F, col.names = F, file=paste(vcf_folder,"/",geneID,"_positions_synonymous_polymorphims",sep=""))
      #format vcftools command
      #print("wrote table")
      calculate_freq_command=paste("vcftools --vcf ",path_to_vcf, " --positions ", vcf_folder,"/",geneID, "_positions_synonymous_polymorphims --freq2 --stdout", sep="")
      #call vcftools
      #print("giving call to vcftools")
      command_return = system(command = calculate_freq_command, intern = T)
      #if variants are found in this region
      if(length(command_return)>1){
        
        #format vcftools output
        poly_syn_freq = read.table(text=command_return, skip = 1)
        #add names to columns
        colnames(poly_syn_freq)=c("chr","start","alleles","nchr","freq_ref","freq_alt")
        
        #get maf
        poly_syn_freq$maf=apply(poly_syn_freq[,c(5,6)],1, min)
        
        #extract information about positions and allelic states
        vcf_alleles=polymorphic_sites[,c(2,4,5)]
        #add column names
        colnames(vcf_alleles)=c("start","ref","alt")
        
        #merge with frequency information
        poly_syn_freq=merge(poly_syn_freq, vcf_alleles)
        
        #here i obtain outgroup information from the alignment using chr/pos information
        #prepare necessary information
        variant_coordinates=poly_syn_freq[,c(2,1,1)]
        #get outgroup infromation
        outgroup_info=ldply(apply(variant_coordinates,1, function(x)get_ancestral_state(x)), rbind)
        
        #merge outgroup information with ingroup frequency information
        poly_syn_freq=merge(poly_syn_freq,outgroup_info)
        
        if(dim(poly_syn_freq)[1]>0)
        {
          #polarize data
          #how much data am i loosing because it not polarizable?
          poly_syn_freq=ldply(apply(poly_syn_freq, 1, function(x) polarize_data(x)), rbind)
          
          #setting the formats right again
          poly_syn_freq$freq_der=as.numeric(as.character(poly_syn_freq$freq_der))
          poly_syn_freq$freq_anc=as.numeric(as.character(poly_syn_freq$freq_anc))
          
          #adding new variable indicating whether SNP is syn or nonsyn
          poly_syn_freq$type=rep("syn",dim(poly_syn_freq)[1])} else {poly_syn_freq=data.frame()}
        
        } else {poly_syn_freq=data.frame()}
        


      #UPDATING NONSYNONYMOUS INFORMATION IF VARIANTS ARE FOUND
      #extract nonsynonymous sites
      polymorphic_nonsynonymous_sites=polymorphic_sites_syn_and_nonsyn$nonsyn
      #.write to external file for vcftools
      write.table(t(do.call(rbind, polymorphic_nonsynonymous_sites[,c(1,2)])), quote=F, row.names = F, col.names = F, file=paste(vcf_folder,"/",geneID,"_positions_nonsynonymous_polymorphims",sep=""))
      #format vcftools command
      calculate_freq_command=paste("vcftools --vcf ",path_to_vcf, " --positions ", vcf_folder,"/",geneID, "_positions_nonsynonymous_polymorphims --max-missing 0.8 --freq2 --stdout", sep="")
      #call vcftools
      command_return = system(command = calculate_freq_command, intern = T)
      #if variants are found in this region
      if(length(command_return)>1){
        
        #format vcftools output
        poly_nonsyn_freq = read.table(text=command_return, skip = 1)
        #add names to columns
        colnames(poly_nonsyn_freq)=c("chr","start","alleles","nchr","freq_ref","freq_alt")
        
        #get maf
        poly_nonsyn_freq$maf=apply(poly_nonsyn_freq[,c(5,6)],1, min)
        
        #extract information about positions and allelic states
        vcf_alleles=polymorphic_sites[,c(2,4,5)]
        #add column names
        colnames(vcf_alleles)=c("start","ref","alt")
        #merge with frequency information
        poly_nonsyn_freq=merge(poly_nonsyn_freq, vcf_alleles)
        
        
        #here i obtain outgroup information from the alignment using chr/pos information
        #prepare necessary information
        variant_coordinates=poly_nonsyn_freq[,c(2,1,1)]
        #get outgroup infromation
        outgroup_info=ldply(apply(variant_coordinates,1, function(x)get_ancestral_state(x)), rbind)
        
        #merge outgroup information with ingroup frequency information
        poly_nonsyn_freq=merge(poly_nonsyn_freq,outgroup_info)
        
        if(dim(poly_nonsyn_freq)[1]>0)
        {
          #polarize data
          #how much data am i loosing because it not polarizable?
          poly_nonsyn_freq=ldply(apply(poly_nonsyn_freq, 1, function(x) polarize_data(x)), rbind)
          
          #setting the formats right again
          poly_nonsyn_freq$freq_der=as.numeric(as.character(poly_nonsyn_freq$freq_der))
          poly_nonsyn_freq$freq_anc=as.numeric(as.character(poly_nonsyn_freq$freq_anc))
          
          #adding new variable indicating whether SNP is syn or nonsyn
          poly_nonsyn_freq$type=rep("nonsyn",dim(poly_nonsyn_freq)[1])} else { poly_nonsyn_freq=data.frame() }
        
        } else {
        poly_nonsyn_freq=data.frame()}
          
       
        if(dim(poly_syn_freq)[1]+dim(poly_nonsyn_freq)[1])
        {
          #concatenating syn and nonsyn dataframes
          poly_syn_nonsyn_freq=rbind(poly_syn_freq, poly_nonsyn_freq)
          
          #adding back geneID
          #i have to add the geneID here 
          poly_syn_nonsyn_freq$geneID=rep(rep(geneID),dim(poly_syn_nonsyn_freq)[1])
          return(poly_syn_nonsyn_freq)}
        
    
      }
    } else {
          log_file_name=paste(analysisID,"_get_vcf_subset_file.txt",sep="")
          cat(paste("Problem in get_vcf_subset line 54", "status_code:", r$status_code, "geneID:", geneID, "\n"),file = log_file_name, append = T)
        }
  }

