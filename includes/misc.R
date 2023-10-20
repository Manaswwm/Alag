
#simple function to convert the three component of genomic coordinates (chr,start,end) in to a string that can be embeded in url for GET or POST

format_gene_coordinate=function(coordinate_vector)
{
  coordinate_string=paste(coordinate_vector[1],":",coordinate_vector[2],"-",coordinate_vector[3], sep="")
  return(coordinate_string)
  
}

#this is to format the key query strings for the VEP endpoint POST vep/:species/region
#the argument should be a table with four columns chr, pos, ref, alt
format_allele_coordinate=function(allele_info_table){
  
  allele_info_table=apply(allele_info_table, 1, toString)
  allele_info_table=gsub(",","",allele_info_table)

  return(allele_info_table)
}
