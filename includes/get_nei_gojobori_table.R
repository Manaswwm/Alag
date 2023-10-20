  get_nei_gojobori_table=function()
{
  #add code to check if table is present!!!!!!!

  #if table is not present
  #declare hash table that will be filled with codon specific synonymous distances
  codon_specific_number_of_syn_sites=new.env(hash = T)
  
  #GENETIC code is a hash table defined in Biostrings  
  #loop on each codon in the genetic code
  lapply(names(GENETIC_CODE),function(x) {
    
    #initialize current codon in hash table
    codon_specific_number_of_syn_sites[[x]]=0
    
    #loop over each nucleotide
    lapply(c(1:3), function(y) {
      
      #identify current codon
      current_codon=unlist(strsplit(x, split=""))
      #identify current nucleotide
      current_nucleotide=current_codon[y]
      #identify current alternative nucleotides (and change format to string)
      alternative_nucleotides=c("A","T","C","G")[-which(c("A","T","C","G")==current_nucleotide)]
      #create possible alternative codons
      #loop over possible alternative nucleotides
      lapply(c(1:3), function(z) {
        #create alternative codon for codon x at position y with alternative allele z (! it is x that gets modified here)
        substr(x,y,y+0)<-alternative_nucleotides[z]
        #(format: current codon is formatted back to a single string)
        current_codon=paste(current_codon, collapse = "")
        #test if current change is synonymous, if yes increment 
        if(GENETIC_CODE[current_codon]==GENETIC_CODE[x]){
          codon_specific_number_of_syn_sites[[current_codon]]=codon_specific_number_of_syn_sites[[current_codon]]+(1/3)
#          if(current_codon=="TTA"){
#            print(paste("debug", current_codon, y, current_nucleotide))
#          }
        }
      })
    })
  })
    
    return(codon_specific_number_of_syn_sites)
  
}