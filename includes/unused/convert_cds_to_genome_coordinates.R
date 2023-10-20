# call point for converting cds to genomic co-ordinates

convert_cds_to_genome_coordinates = function(){
  
  #dummy input -- please delete later!!! #
  item = transcript_alignments[[1]]
  
  #getting a subset of the gff file for the give transcript first
  transcript_gff_subset = thaliana_gff_subset[thaliana_gff_subset$attributes == item$id,]
  
  #collecting all the cds positions here
  cds_positions = c()
  
  #going over all rows in gff subset
  for(item in 1:nrow(transcript_gff_subset)){
    
    row = transcript_gff_subset[item,]
    
    cds_positions = append(cds_positions, seq(row$start, row$end, by = 1))
    
  }

  
    
}
