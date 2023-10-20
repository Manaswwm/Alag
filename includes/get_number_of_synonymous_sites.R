get_number_of_synonymous_sites=function(seq1_orf){
  
  nei_gojobori_table=get_nei_gojobori_table()
  
  #split input sequence into array of codon
  all_codons=substring(seq1_orf,seq(1,nchar(seq1_orf)-2,3),seq(3,nchar(seq1_orf),3))

  #calculate number of synonymous sites
  codon_specific_distances=sapply(all_codons, function(x) return(nei_gojobori_table[[x]]))
  #print(class(codon_specific_distances))
  return(sum(codon_specific_distances))
    
  
}