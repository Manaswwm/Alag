#I have a serious bug here tha leads to rejection of genes for which diversity is 0,
#also a bug when div_syn or div_non_syn is 0. 
#I will have to keep track of geneIDs before this function, then check for genes for which div is 0.
#Also I will have to take care of SNPs first (vep-syn/non_syn) and later assign number of SNPs with total syn and non_syn total length (in bp)

#alignment matrix is a  list returned by the function get_transcript_alignment
get_divergent_sites2=function(alignment_matrix)
{
  
  #@print check
  print(paste("Estimating divergent sites of ",alignment_matrix$id))
  print("Sub-part 1 - extracting the variant sites and calculating lengths")
  
  #get alignment between reference sequence and outgroup sequence
  m=alignment_matrix$alignment
  #keep track of the transcript ID
  transcript_id=alignment_matrix$id
  #keep track of geneID
  geneID=gene_conversion_table$ensembl_gene_id[gene_conversion_table$ensembl_transcript_id == transcript_id]
  #keep tracl of chromosome
  chr=alignment_matrix$chr

  #PART 1: Get the lappropriate length for syn/non-syn divergence and polymorphism
  #this bloc is to calculate total syn and non-syn length and will go to another function
  #delete all gaps
  m_no_gaps=del.colgapsonly(m, threshold = 0.0000001)
  #re-format the ungapped reference sequence for the next function
  seq_ref_no_gaps=paste(m_no_gaps[1,], collapse="")#this brakes the reference coordinate system (because nucleotides in reference for which no outgroup info is available are removed)!!
  
  #### important point that needs to be looked at for alag local - if there are other characters present in the sequence of transcript besides A,C,G and T then such transcripts are to be dropped ####
  #### MJ-11122022
  #checking if the total length of the sequence now is equal to the sum of the occurances of A,G,C and T in the sequence
  if(sum(str_count(string = toupper(seq_ref_no_gaps), c("A", "G", "C", "T"))) == nchar(seq_ref_no_gaps)){
  
    #get length of ungapped reference sequence
    ungapped_length=nchar(seq_ref_no_gaps)
    #here I should check whether I still have an ORF or at least whether the length is multiple of 3
    #calculate total possible number of synonymous/non-synonymous sites (Nei and Gojobori)
    synonymous_length=get_number_of_synonymous_sites(toupper(seq_ref_no_gaps))
    nonsynonymous_length=ungapped_length-synonymous_length
    
  
    #delete all sites with a gap in thaliana (to get positions back into TAIR10 reference system)
    #identify which columns have a gap "-" in thaliana
    if(length(which(m[1,]=="-"))>0)
    {
      m_no_gaps_in_hirsuta=m[,-which(m[1,]=="-")]
    } else if(length(which(m[1,]=="-"))==0){
      m_no_gaps_in_hirsuta=m
    } else(stop("line 27 in get_divergent_sites: problem with filtering gapped sites in A. thaliana"))
    
    #extract sequences
    seq1=m_no_gaps_in_hirsuta[1,]
    seq2=m_no_gaps_in_hirsuta[2,]
    
    #re-format the ungapped reference sequence for the next function
    seq_ref_cds=paste(seq1, collapse="")
    #get length of ungapped reference sequence
    cds_length=nchar(seq_ref_cds)
    #calculate total possible number of synonymous/non-synonymous sites on total cds (for polymorphism diversity)
    synonymous_length_cds=get_number_of_synonymous_sites(toupper(seq_ref_cds))
    nonsynonymous_length_cds=cds_length - synonymous_length_cds
    
    #prepare first table for output: all 6 lengths and geneID
    df_syn_nonsyn_lengths_for_div_poly=data.frame(geneID, ungapped_length,synonymous_length,nonsynonymous_length)
    #PART 1 is done
    
    #Now PART 2: Get lis of divergent sites to be further classified into synonymous and non-synonymous
    print("Sub-part 2 - binding the collected information on variants and lengths together")
    
    #identify divergent sites
    list_of_divergent_sites=which(seq1!=seq2)
    
    #MJ-29052020 - In some cases of DNABD, there seems to be noe divergent sites in between thaliana and lyrata
    #I try to catch for this by modifying the code - I put this now in an if..else loop
    if(length(list_of_divergent_sites) == 0){
      divergent_sites = NULL
    }
    
    else{
    
    #identify thaliana (ref) alleles at divergent sites
    list_of_divergent_alleles_ref=as.array(seq1)[which(seq1!=seq2)]
    
    #identify lyrata (alt) alleles at divergent sites
    list_of_divergent_alleles_alt=as.array(seq2)[which(seq1!=seq2)]
    
    #PREPARING COORDINATE CONVERSION PROCESS
    #creating vector with replicated chromosome ID and transcript ID
    id_vector=rep(transcript_id, length(list_of_divergent_sites))
    #put in dataframe together with cds positions
    cds_positions=data.frame(id=id_vector, cds_pos=list_of_divergent_sites)
    #now use endpoint map/cds/:id/:region (GET) to convert cds positions back into genomic coordinates
    #cds_position_genomic_coord=apply(cds_positions,1, function(x) convert_cds_to_genomic_coordinates(x))
    
    #### writing code for local alag here #### 
    cds_positions$original.start = rownames(cds_positions)
    
    #convert list of df to df
    #cds_positions=ldply(cds_position_genomic_coord,rbind)
    
    #creating vector of replicated (useless) dots (to conform to format requirements for json body in VEP POST)
    snpID_vector=rep(".",length(list_of_divergent_sites))
    
    #creating vector with gene IDs 
    geneID=rep(geneID,length(list_of_divergent_sites))
    
    #pack everything together
    div_sites=cbind.data.frame(list_of_divergent_sites, snpID_vector, list_of_divergent_alleles_ref, list_of_divergent_alleles_alt, geneID)
    
    #giving the genomic coordinates
    div_sites$list_of_divergent_sites = rownames(div_sites)
    
    #merge with table containing genomic positions
    divergent_sites=merge(div_sites, cds_positions, by.x="list_of_divergent_sites", by.y="original.start")
    
    #exclude sites that are different because of a deletion in A. lyrata
    divergent_sites=subset(divergent_sites,divergent_sites$list_of_divergent_alleles_alt!="-")
    
    #checking if there are any divergent sites after removing the gapped divergent sites - 20082021
    if(!dim(divergent_sites)[1] == 0){
    
      #giving back the chromosome information
      divergent_sites$mapped.seq_region_name = rep(alignment_matrix$chr)
      
      #removing cds_pos
      #divergent_sites = subset(divergent_sites, select = -c(cds_pos))
      
      #changing column names
      colnames(divergent_sites) = c("mapped.start", "snpID_vector", "list_of_divergent_alleles_ref",
                                    "list_of_divergent_alleles_alt", "geneID", "id", "cds_pos", "mapped.seq_region_name")
      
      divergent_sites = divergent_sites[,c(8,7,1,2,3,4,5,6)]
      
      #giving back the transcript positions of the variants here
      #adjusting the position here to -1 -- MJ12122022
      divergent_sites$transcript_pos = divergent_sites$cds_pos + alignment_matrix$transcript_start_pos - 1
      
      #reduce dimension
      #divergent_sites=divergent_sites[,c("mapped.seq_region_name","mapped.start","snpID_vector","list_of_divergent_alleles_ref","list_of_divergent_alleles_alt","geneID")]
      
      #print check
      #######################print(divergent_sites)
    }else{(divergent_sites = NULL)}
    
    return(list(lengths=df_syn_nonsyn_lengths_for_div_poly, variant_list=divergent_sites))
    
    }
  }
}
