summarize_frequencies=function(x=NULL, min_freq=NULL, variant_type=NULL, transcriptID=NULL){
  
  #calculating total heterozygosity
  pi=sum(2*x$freq_der*(1-x$freq_der))
  
  #filter out frequencies lower than user-specified value (min_freq)
  derived_frequencies_filtered=x$freq_der[x$freq_der>=min_freq]
  
  #calculating S (or P in the MK framework) : the number of variant sites
  P=length(derived_frequencies_filtered)
  
  
  #if UPSTREAM DATA
  if(variant_type=="upstream"){
    
    
    #heterozygosity per base pair
    pi_upstream=pi/window_size #CAUTION: this is a global variable
    
    #formatting upstream-specific results
    results=data.frame(id=transcriptID, P_upstream=P, pi_upstream=pi_upstream, window_size=window_size)
    
    #if SYNONYMOUS DATA  
    
  } else if(variant_type=="syn"){  
    
    #get size (denominator) for this region
    #print(get_transcript_cds(geneID)$seq)
    
    #accessing fasta file containing the sequence of the transcript
    transcript_seq = toupper(read.fasta(file = paste(backup_file,"/transcript_sequences/",transcriptID,"_sequence.fa",sep=""), as.string = TRUE))[1]
    
    P_syn_size=get_number_of_synonymous_sites(transcript_seq)
    #P_syn_size = syn_length
    #heterozygosity per base pair
    pi_syn=pi/P_syn_size
    
    #formatting syn-specific results
    results=data.frame(id=transcriptID, P_syn=P, pi_syn=pi_syn, syn_length=P_syn_size)
    
    #if NON-SYNONYMOUS DATA
  } else if(variant_type=="nonsyn") {
    
    #accessing fasta file containing the sequence of the transcript
    transcript_seq = toupper(read.fasta(file = paste(backup_file,"/transcript_sequences/",transcriptID,"_sequence.fa",sep=""), as.string = TRUE))[1]
    
    P_nonsyn_size= nchar(transcript_seq) - get_number_of_synonymous_sites(transcript_seq)
    #P_non_syn_size = nonsyn_length
    #heterozygosity per base pair
    pi_nonsyn=pi/P_nonsyn_size
    
    #formatting syn-specific results
    results=data.frame(id=transcriptID, P_nonsyn=P, pi_nonsyn=pi_nonsyn, nonsyn_length=P_nonsyn_size)
  
  }else if(variant_type=="nonsense"){
    
    #accessing fasta file containing the sequence of the transcript
    transcript_seq = toupper(read.fasta(file = paste(backup_file,"/transcript_sequences/",transcriptID,"_sequence.fa",sep=""), as.string = TRUE))[1]
    
    P_nonsense_size= nchar(transcript_seq) - get_number_of_synonymous_sites(transcript_seq)
    #P_non_syn_size = nonsyn_length
    #heterozygosity per base pair
    pi_nonsense=pi/P_nonsense_size
    
    #formatting syn-specific results
    results=data.frame(id=transcriptID, P_nonsense=P, pi_nonsense=pi_nonsense, nonsense_length=P_nonsense_size)
  }
  return(results)
}

get_polymorphism_stats_per_gene=function(allele_frequencies=NULL, min_freq, variant_type=NULL, window_size=window_size){
  
  #extracting SNPs of the specified type (upstream ,syn, or nonsyn)
  allele_frequencies=subset(allele_frequencies, allele_frequencies$type==variant_type)
  by(data = allele_frequencies,INDICES = allele_frequencies$geneID,FUN = function(x) summarize_frequencies(x, min_freq = min_freq, variant_type, window_size))
  
}

get_polymorphism_stats_per_gene2=function(transcriptID=NULL, allele_frequencies=NULL, min_freq, type=NULL){
  
  #extracting SNPs of the specified type (upstream ,syn, or nonsyn)
  allele_frequencies=allele_frequencies[allele_frequencies$type==type & allele_frequencies$ensembl_transcript_id==as.character(transcriptID),]
  #by(data = allele_frequencies,INDICES = allele_frequencies$geneID,FUN = function(x) summarize_frequencies(x, min_freq = min_freq, variant_type, window_size))
  summarize_frequencies(allele_frequencies, min_freq = min_freq, type, transcriptID)
  
}
