polarize_data=function(one_variant){
  
  one_variant=as.data.frame(t(one_variant))

  
  #print(as.character(one_variant$outgroup) %in% c(as.character(one_variant$ref), as.character(one_variant$alt)))
  
  if(as.character(one_variant$outgroup) %in% c(as.character(one_variant$ref_allele), as.character(one_variant$alt_allele)))
  {
    
    if(one_variant$ref_allele==one_variant$outgroup){
      
      one_variant$freq_anc=one_variant$freq_ref
      one_variant$freq_der=one_variant$freq_alt
      
    } else if(one_variant$alt_allele==one_variant$outgroup) {
      
      one_variant$freq_anc=one_variant$freq_alt
      one_variant$freq_der=one_variant$freq_ref
      
    }
    
    return(one_variant)
    
  } else {
    log_file_name=paste("polarize_data_log_file.txt",sep="")
    cat(paste("identified multiple hit i.e more than 2 states at mapped.start: ", one_variant$mapped.start," ; for gene: ",one_variant$geneID,"\n"),file = log_file_name, append = T)
  }
  

  
  
}